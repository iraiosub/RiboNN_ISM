## SCN2A 5′UTR ISM analysis

End-to-end workflow for running in silico mutagenesis (ISM) on the SCN2A 5′UTR and
testing whether the altAUG positions are ΔTE outliers.

### Prerequisites

```bash
conda activate ribonn
# reference files (CAMP paths – adjust for your system)
FASTA=/camp/lab/ulej/home/shared/oscar_ira_riboloco/ref/human/GRCh38.primary_assembly.genome.fa
GTF=/camp/lab/ulej/home/shared/oscar_ira_riboloco/ref/human/gencode.v44.primary_assembly.annotation.longest_cds_transcripts.gtf.gz
```

### Step 1 – Generate ISM variant sequences

```bash
python prepare_scn2a_ism.py \
    --fasta "${FASTA}" \
    --gtf   "${GTF}" \
    --upstream-bases 9999 \
    --max-deletion   15 \
    --truncate-utr3 \
    --output data/prediction_input.txt
```

`--upstream-bases 9999` is automatically capped to the actual UTR length, producing a
full saturation scan of every base in the 5′UTR (3 SNVs + 1 deletion per position).
Use `--audit` first to inspect the UTR length before committing to the full run.

### Step 2 – Run RiboNN predictions (GPU)

```bash
# interactive / local
python run_ribonn_predict.py \
    --input  data/prediction_input.txt \
    --output results/human/prediction_output.txt

# on the SLURM cluster (steps 1–3 in one job)
sbatch submit_ism_scn2a.sh
```

Output: `results/human/prediction_output.txt` (one row per variant, TE per cell type).

### Step 3 – Plot ISM results

```bash
python plot_te_changes.py \
    --input          results/human/prediction_output.txt \
    --outdir         plots_ism_scn2a \
    --upstream-bases 15
```

Saves SNV heatmap, deletion heatmap, waterfall chart, and neuronal cell-type panel
to `plots_ism_scn2a/`.

### Step 4 – altAUG position-matched null test

Tests whether ΔTE at the altAUG positions (−8/−7/−6 by default) are outliers relative
to the background distribution of all other scanned positions.

```bash
# interactive / local
python analyze_altaug_null.py \
    --input            results/human/prediction_output.txt \
    --altaug-positions -8 -7 -6 \
    --outdir           plots_ism_scn2a

# on the SLURM cluster (run after step 2, or chain with --dependency)
sbatch submit_altaug_null.sh
sbatch --dependency=afterok:<ISM_JOB_ID> submit_altaug_null.sh
```

Outputs:
- `plots_ism_scn2a/altaug_null_summary.tsv` – per-SNV z-scores and empirical p-values
- `plots_ism_scn2a/altaug_null_plot.png`    – background distribution with altAUG signal overlaid

---

# RiboNN: A deep learning model to predict translation efficiency from mRNA sequence

For more information, please see our [RiboNN paper](https://www.nature.com/articles/s41587-025-02712-x).

- System requirements:

  This code has been tested on a system with 4 CPUs, 16 Gb RAM, and 1 NViDIA 10A GPU, with Ubuntu 20.04 as the OS (with CUDA Toolkit 11.3 installed). The required softwares are listed in environment.yml.

- To install project requirements:
  ```bash
  sudo apt install make

  # install mamba (https://github.com/conda-forge/miniforge) into "miniforge3/" in the home directory. 
  curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
  bash Miniforge3-$(uname)-$(uname -m).sh -b 
  ~/miniforge3/bin/mamba shell init 
  source ~/.bashrc

  # clone the repo
  git clone https://github.com/Sanofi-Public/RiboNN.git && cd RiboNN
  
  # install the RiboNN environment
  make install

  # activate the riboNN environment
  mamba activate RiboNN

  ```
  Note: Depending on your network speed, it may take a few minutes to install the required packages.


- To train the RiboNN model from scratch:
   1. Put the training data in a tab-separated text file in the "data" folder, which already contain an example training data file with **fake** TEs. The tab-separated text file should have columns named "tx_id" (unique transcript IDs), "utr5_sequence", "cds_sequence" (including start and stop codons), and "utr3_sequence". Alternatively, the file may have columns named "tx_id", "tx_sequence" (full transcript seuquences containing 5'UTR, CDS, and 3'UTR), "utr5_size" (lengths of the 5'UTRs), and "cds_size" (lengths of the CDSs). The published human and mouse models were trained on data in the Supplementary Tables published in the RiboNN paper.
   2. Edit the path to the training data ("tx_info_path") and other hyperparameters defined in the config/conf.yml file. 
   3. Edit the code below line 18 of src/main.py to control how the model will be trained.
   4. Run `make train` at the terminal to start the training process.
  
- To do transfer learning (using pretrained human multi-task models automatically downloaded from https://zenodo.org/records/17258709):
   1. Put the training data in a tab-separated text file in the "data" folder, which already contain an example training data file. The tab-separated text file should have columns named "tx_id" (unique transcript IDs), "utr5_sequence", "cds_sequence" (including start and stop codons), and "utr3_sequence". Alternatively, the file may have columns named "tx_id", "tx_sequence" (full transcript seuquences containing 5'UTR, CDS, and 3'UTR), "utr5_size" (lengths of the 5'UTRs), and "cds_size" (lengths of the CDSs). 
   2. Edit the path to the training data ("tx_info_path") and other hyperparameters defined in the config/conf.yml file. 
   3. Edit the code below line 118 of src/main.py to control how the model will be trained.
   4. Run `make transfer_learning` at the terminal to start the training process.
  
- To make predictions using pretrained multi-task models automatically downloaded from https://zenodo.org/records/17258709:
  1. Create a tab-separated text file with columns named "tx_id" (unique transcript IDs), "utr5_sequence", "cds_sequence" (including start and stop codons), and "utr3_sequence". Alternatively, the file may have columns named "tx_id", "tx_sequence" (full transcript seuquences containing 5'UTR, CDS, and 3'UTR), "utr5_size" (lengths of the 5'UTRs), and "cds_size" (lengths of the CDSs). 
  2. Save the text file as "prediction_input.txt" in the "data" folder. An example input file can be found in the "data" folder.
  3. (Optional) Edit the code below line 163 of src/main.py to control how the model will be used for prediction.
  4. To use human models for prediction, run `make predict_human` at the terminal. To use mouse models for prediction, run `make predict_mouse`.
  5. The predictions will be automatically written to a tab-separated file named "prediction_output.txt" in the "results/human" or "results/mouse" folder. Pre-existing files with the same name will be overwritten.
  
  **Note:** Input transcripts with 5'UTRs longer than 1,381 nt or combined CDS and 3'UTR sizes larger than 11,937 nt will be excluded in the output.   
