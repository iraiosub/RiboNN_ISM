# Transcriptome-wide 5′UTR mutagenesis

This is a separate, additive workflow for saturating every retained 5′UTR
position with:

- all three possible single-nucleotide substitutions; and
- a one-nucleotide deletion.

It predicts RiboNN TE changes without producing per-transcript heatmaps. The
main output is one compressed table row per 5′UTR position:

`output/<species>/final/all_utr5_position_scores.tsv.gz`

Important columns are:

- `transcript_id`, `utr5_position_1based`, `offset_from_cds_start`, `ref_base`
- `wt_mean_predicted_TE`
- `substitution_te_change_mean`: mean of the three alternative-base scores
  minus the transcript WT score
- `substitution_direction_vs_wt`: `increase`, `decrease`, or `unchanged`
- `substitution_effect_pattern`: whether all three alternatives agree or are
  mixed
- `deletion_te_change`: single-base deletion score minus WT
- `deletion_direction_vs_wt`

Here, mean predicted TE is averaged across the RiboNN output cell types and
test folds, matching the scalar `mean_predicted_TE` idea used by the existing
SCN2A workflow.

## Input choice

The default is the same species-specific genome FASTA plus
`longest_cds_transcripts.gtf.gz` reference used by the existing SCN2A RiboNN
analyses:

```bash
bash all_utr5_mutagenesis/submit_all_utr5_mutagenesis.sh \
  --species human
```

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

The RiboScanner file `master_table.context_m40_p40.fa.gz` used elsewhere in
this repository contains fixed 81-nt start-codon contexts, not full
transcripts. It cannot support whole-5′UTR RiboNN mutagenesis because it lacks
complete 5′UTR, CDS, and 3′UTR sequences.

## HPC scaling

The workflow uses a weighted SLURM array. Transcripts remain intact within a
shard, while shard boundaries are chosen by predicted variant count, so GPU
tasks have similar loads and each transcript WT is predicted only once.

```bash
bash all_utr5_mutagenesis/submit_all_utr5_mutagenesis.sh \
  --species human \
  --num-shards 128 \
  --max-concurrent 8 \
  --batch-size 256
```

Each array task:

1. materializes only its own compressed variant TSV;
2. runs RiboNN;
3. immediately reduces the model outputs to one scalar per variant;
4. writes the per-position summary; and
5. removes the reproducible variant input/scores after successful summary.

Use `--keep-intermediates` to retain those shard files. Completed shards have a
`work/shard_NNN/complete` marker and are skipped on a rerun.

Useful tuning:

- increase `--num-shards` if a task uses too much RAM;
- reduce `--batch-size` if GPU memory is limiting;
- set `--max-concurrent` to the number of GPUs you want to occupy;
- use `--top-k 5` for parity with the existing analyses.

## QC and reproducibility

The catalog stage excludes transcripts that RiboNN cannot encode, including
non-ACGT sequence, empty/over-limit 5′UTRs, noncanonical CDS starts, invalid
CDS lengths, or missing terminal stop codons. Long 3′UTRs are truncated by
default to preserve the full CDS and satisfy the model limit.

Inspect:

- `catalog/catalog_audit.tsv.gz` for inclusion/exclusion reasons
- `catalog/shard_manifest.tsv` for array load balance
- `workflow_manifest.json` for total transcript, position, and variant counts
- `final/run_summary.json` for merged-output validation

