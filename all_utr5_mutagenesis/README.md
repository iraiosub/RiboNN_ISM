# ORF-start-focused 5′UTR mutagenesis

This is a separate, additive workflow for saturating selected 5′UTR positions
with:

- all three possible single-nucleotide substitutions; and
- a one-nucleotide deletion.

By default it no longer mutates every 5′UTR base. It uses the ORF prediction
table to keep only start codons fully inside the 5′UTR, then mutates exactly
`orf_start`, `orf_start + 1`, and `orf_start + 2`. The catalog stage checks that
each retained ORF start codon is `ATG`.

It predicts RiboNN TE changes without producing per-transcript heatmaps. The
main outputs are:

- `output/<species>_orf_starts_<te_label>/final/all_utr5_position_scores.tsv.gz`: one row
  per targeted 5′UTR base
- `output/<species>_orf_starts_<te_label>/final/orf_start_codon_scores.tsv.gz`: one row
  per retained ORF start codon, averaging the three start-codon positions

Important columns are:

- `transcript_id`, `utr5_position_1based`, `offset_from_cds_start`, `ref_base`
- `orf_start_1based`, `orf_position_in_start_codon`, `target_codon`
- `wt_mean_predicted_TE`
- `substitution_te_change_mean`: mean of the three alternative-base scores
  minus the transcript WT score
- `substitution_direction_vs_wt`: `increase`, `decrease`, or `unchanged`
- `substitution_effect_pattern`: whether all three alternatives agree or are
  mixed
- `deletion_te_change`: single-base deletion score minus WT
- `deletion_direction_vs_wt`

The ORF-level table adds:

- `orf_substitution_te_change_mean_3nt`
- `orf_deletion_te_change_mean_3nt`
- direction columns for those three-position averages

Here, `mean_predicted_TE` matches the repository's prediction wrapper: it is
computed after native multitask prediction as the mean across all
`predicted_TE_*` outputs, then across folds. Human defaults to
`predicted_TE_normal_brain_tissue` averaged across folds, while mouse defaults
to `mean_predicted_TE`. Pass `--te-column mean_predicted_TE` to use this
all-output aggregate for any species.

## Input choice

The default is the same species-specific genome FASTA plus
`longest_cds_transcripts.gtf.gz` reference used by the existing SCN2A RiboNN
analyses. Human ORF-start mode uses the human ORF prediction file under
`/camp/lab/ulej/home/shared/oscar_ira_riboloco/ref/human/orfs`. Mouse ORF-start
mode uses the unified cross-tissue master table:
`/camp/lab/ulej/home/shared/oscar_ira_riboloco/analysis_results/cross_tissue.unmixing.master_table.with_below_tpm_threshold.tsv.gz`.

```bash
bash all_utr5_mutagenesis/submit_all_utr5_mutagenesis.sh \
  --species human
```

For the mouse cross-tissue run:

```bash
bash all_utr5_mutagenesis/submit_all_utr5_mutagenesis.sh \
  --species mouse
```

Mouse defaults to averaging all mouse RiboNN predicted TE outputs across folds,
so the default output directory is
`output/mouse_orf_starts_mean_predicted_TE`.

Run the launcher from the checkout with `bash`; it submits the prep, GPU-array,
and merge jobs itself. If a site wrapper or accidental `sbatch` call executes a
copy from `/tmp/slurmd`, pass the checkout explicitly:

```bash
bash all_utr5_mutagenesis/submit_all_utr5_mutagenesis.sh \
  --repo-root /path/to/RiboNN_ISM \
  --species human
```

To provide a different ORF table:

```bash
bash all_utr5_mutagenesis/submit_all_utr5_mutagenesis.sh \
  --species human \
  --orf-predictions /path/to/orf_predictions.csv.gz
```

The ORF table may be a CSV ORF-prediction file with
`transcript_id,orf_start,orf_stop,orf_frame,annotated`, or a TSV master table
with `orf_start_1based`/`orf_stop_1based` aliases or an `orf_id` formatted as
`<transcript_id>_<orf_start>_<orf_stop>_<orf_frame>`.

To use a full transcript FASTA directly, provide the matching GTF so that the
workflow can locate the 5′UTR/CDS boundaries:

```bash
bash all_utr5_mutagenesis/submit_all_utr5_mutagenesis.sh \
  --species human \
  --transcript-fasta /path/to/gencode.transcripts.fa.gz \
  --gtf /path/to/matching.annotation.gtf.gz
```

An existing RiboNN-format transcript table is also accepted:

```bash
bash all_utr5_mutagenesis/submit_all_utr5_mutagenesis.sh \
  --species human \
  --input-table /path/to/transcripts.tsv
```

The old full-5′UTR screen is still available, but it is intentionally opt-in:

```bash
bash all_utr5_mutagenesis/submit_all_utr5_mutagenesis.sh \
  --species human \
  --all-utr5
```

The RiboScanner file `master_table.context_m40_p40.fa.gz` used elsewhere in
this repository contains fixed 81-nt start-codon contexts, not full
transcripts. It cannot support whole-5′UTR RiboNN mutagenesis because it lacks
complete 5′UTR, CDS, and 3′UTR sequences.

## HPC scaling

The workflow uses a weighted SLURM array. Transcripts remain intact within a
shard, while shard boundaries are chosen by targeted variant count, so GPU
tasks have similar loads and each transcript WT is predicted only once.

```bash
bash all_utr5_mutagenesis/submit_all_utr5_mutagenesis.sh \
  --species human \
  --num-shards 128 \
  --max-concurrent 8 \
  --batch-size 256
```

Each array task:

1. materializes only its own compressed targeted-variant TSV;
2. runs RiboNN;
3. immediately reduces the model outputs to one scalar per variant;
4. writes the per-position summary; and
5. removes the reproducible variant input/scores after successful summary.

Use `--keep-intermediates` to retain those shard files. Completed shards have a
`work/shard_NNN/complete` marker and are skipped on a rerun.

Useful tuning:

- increase `--num-shards` if a task is still too slow or uses too much RAM;
- reduce `--batch-size` if GPU memory is limiting;
- set `--max-concurrent` to the number of GPUs you want to occupy;
- use `--top-k 5` for parity with the existing analyses.

## QC and reproducibility

The catalog stage excludes transcripts that RiboNN cannot encode, including
non-ACGT sequence, empty/over-limit 5′UTRs, noncanonical CDS starts, invalid
CDS lengths, or missing terminal stop codons. In default ORF mode it further
excludes transcripts without any ATG ORF start codon fully inside the 5′UTR.
Long 3′UTRs are truncated by default to preserve the full CDS and satisfy the
model limit.

Inspect:

- `catalog/catalog_audit.tsv.gz` for inclusion/exclusion reasons
- `catalog/orf_target_audit.tsv.gz` for ORF-start ATG and 5′UTR checks
- `catalog/shard_manifest.tsv` for array load balance
- `workflow_manifest.json` for total transcript, position, and variant counts
- `final/run_summary.json` for merged-output validation
