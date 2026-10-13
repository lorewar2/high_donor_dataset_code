# high_donor_dataset_code

Code for generating synthetic high-donor pooled single-cell RNA-seq datasets, with a configurable number of donors, doublet rate and donor imbalance, and for running the clustering pipeline on them.

## Requirements

- [samtools](http://www.htslib.org/) ≥ 1.10
- Python 3 
- The 10x Genomics GRCh38 reference (see below)

## 1. Download the reference genome

The CRAM file was created against the 10x Genomics Cell Ranger GRCh38 reference (`refdata-gex-GRCh38-2024-A`). You need this exact reference to decode it, 

Download it from the [10x Genomics Cell Ranger downloads page](https://www.10xgenomics.com/support/software/cell-ranger/downloads), or directly:

```bash
wget https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2024-A.tar.gz
tar -xzf refdata-gex-GRCh38-2024-A.tar.gz
```

The genome FASTA is then at `refdata-gex-GRCh38-2024-A/fasta/genome.fa`.

## 2. Download the data

Download the CRAM file from Zenodo: link*

```bash
wget -O final.cram <ZENODO_URL>
```

## 3. Convert CRAM to BAM

```bash
samtools view -b \
    -T refdata-gex-GRCh38-2024-A/fasta/genome.fa \
    -o final.bam \
    final.cram

samtools index final.bam
```

## 4. Generate a synthetic dataset

```bash
python datamaker.py -k <num_donors> -d <doublet_rate> -i <imbalance>
```

| Flag | Description | Allowed values |
|------|-------------|----------------|
| `-k` | Number of donors in the pool | 10–250 |
| `-d` | Doublet rate (%) | 0, 4, 8 |
| `-i` | Donor imbalance parameter | 0–1 |

Example, a 250-donor pool with an 8% doublet rate and imbalance 1:

```bash
python datamaker.py -k 250 -d 8 -i 1
```


## 5. Run the clustering pipeline


```bash
python cluster.py --input <datamaker_output>
```

## Citation
