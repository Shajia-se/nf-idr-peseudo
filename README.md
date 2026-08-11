# nf-idr-pseudo

Legacy/optional IDR module for explicit replicate peak pairs and pseudo-replicate IDR.

For the current main ChIP-seq workflow, prefer `nf-idr`. Use this module only when you already have manual CSV files listing peak pairs or BAMs for pseudo-replicate analysis.

## Standard IDR Input

Create `idr_pairs.csv` or pass `--idr_pairs_csv`:

```csv
pair_name,rep1_peaks,rep2_peaks
WT,/path/to/rep1_peaks.narrowPeak,/path/to/rep2_peaks.narrowPeak
```

Both peak files must exist before the module starts.

## Optional Pseudo-IDR Input

Set `--do_pseudo_idr true` and provide `pseudo_idr_bam.csv`:

```csv
rep_name,bam
WT_rep1,/path/to/sample.nomulti.bam
```

The module splits read names into two pseudo-replicate BAMs, calls peaks, then runs IDR.

## Common Run

```bash
nextflow run main.nf -profile hpc \
  --project_folder /path/to/run_root \
  --idr_pairs_csv /path/to/idr_pairs.csv \
  -resume
```

Pseudo-IDR example:

```bash
nextflow run main.nf -profile hpc \
  --project_folder /path/to/run_root \
  --idr_pairs_csv /path/to/idr_pairs.csv \
  --do_pseudo_idr true \
  --pseudo_idr_bam_csv /path/to/pseudo_idr_bam.csv \
  -resume
```

Actual execution should be validated on the target HPC environment before delivery.
