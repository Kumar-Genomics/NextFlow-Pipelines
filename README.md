`# 🧬 RNA-seq Nextflow Pipeline

> **FASTQ → QC → Trimming → Alignment → Quantification → Count Matrix**

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A522.10.1-23aa62.svg)](https://www.nextflow.io/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg)](https://sylabs.io/docs/)
[![run with conda](https://img.shields.io/badge/run%20with-conda-3EB049?logo=anaconda)](https://conda.io/miniconda.html)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![CI](https://github.com/YOUR_USERNAME/rnaseq-nextflow/actions/workflows/ci.yml/badge.svg)](https://github.com/YOUR_USERNAME/rnaseq-nextflow/actions)

A **fully reproducible, portable RNA-seq analysis pipeline** built with Nextflow DSL2. Supports multiple aligners, trimmers, and quantifiers with automatic strandedness detection.

---

## 📋 Table of Contents

- [Pipeline Overview](#-pipeline-overview)
- [Quick Start](#-quick-start)
- [Installation](#-installation)
- [Usage](#-usage)
- [Parameters](#-parameters)
- [Output Structure](#-output-structure)
- [Downstream Analysis](#-downstream-analysis)
- [Running on HPC](#-running-on-hpc)
- [Hands-on Tutorial](#-hands-on-tutorial)
- [Contributing](#-contributing)

---

## 🔬 Pipeline Overview

```
                    ┌─────────────────────────────────────────┐
                    │           INPUT: FASTQ Files            │
                    │    (paired-end or single-end reads)     │
                    └──────────────────┬──────────────────────┘
                                       │
                    ┌──────────────────▼──────────────────────┐
                    │           STEP 1: FastQC                │
                    │     Raw read quality assessment         │
                    └──────────────────┬──────────────────────┘
                                       │
                    ┌──────────────────▼──────────────────────┐
                    │      STEP 2: Read Trimming               │
                    │   Trimmomatic  ──or──  fastp             │
                    │  (adapter removal + quality trimming)    │
                    └──────────────────┬──────────────────────┘
                                       │
                    ┌──────────────────▼──────────────────────┐
                    │       STEP 3: Post-trim FastQC           │
                    └──────────────────┬──────────────────────┘
                                       │
                    ┌──────────────────▼──────────────────────┐
                    │      STEP 4: Genome Indexing             │
                    │  STAR index  /  HISAT2 index  /         │
                    │  Salmon index  (skipped if pre-built)   │
                    └──────────────────┬──────────────────────┘
                                       │
                    ┌──────────────────▼──────────────────────┐
                    │      STEP 5: Alignment / Mapping         │
                    │  ┌─────────┐ ┌────────┐ ┌────────────┐  │
                    │  │  STAR   │ │ HISAT2 │ │   Salmon   │  │
                    │  │ 2-pass  │ │ splice │ │   quasi    │  │
                    │  └────┬────┘ └───┬────┘ └─────┬──────┘  │
                    └───────┼──────────┼─────────────┼─────────┘
                            │  (BAM)   │             │ (quant)
                    ┌───────▼──────────▼─────────────┘
                    │      STEP 6: SAMtools Sort/Index/Flagstat
                    └──────────────────┬──────────────────────
                                       │
                    ┌──────────────────▼──────────────────────┐
                    │   STEP 7: Strandedness Inference         │
                    │        (RSeQC infer_experiment)          │
                    └──────────────────┬──────────────────────┘
                                       │
                    ┌──────────────────▼──────────────────────┐
                    │      STEP 8: Read Quantification         │
                    │  featureCounts  ──or──  HTSeq-count      │
                    └──────────────────┬──────────────────────┘
                                       │
                    ┌──────────────────▼──────────────────────┐
                    │      STEP 9: Merge Count Matrix          │
                    │     genes × samples TSV matrix          │
                    └──────────────────┬──────────────────────┘
                                       │
                    ┌──────────────────▼──────────────────────┐
                    │      STEP 10: MultiQC Report             │
                    │  Aggregates all QC metrics → HTML report │
                    └─────────────────────────────────────────┘
```

### Tools Supported

| Step | Options |
|------|---------|
| **Trimming** | Trimmomatic, fastp |
| **Alignment** | STAR (2-pass), HISAT2, Salmon (quasi-mapping) |
| **Quantification** | featureCounts, HTSeq-count, Salmon |
| **QC** | FastQC, RSeQC, MultiQC |
| **Containers** | Docker, Singularity, Conda |

---

## ⚡ Quick Start

```bash
# 1. Clone the repository
git clone https://github.com/YOUR_USERNAME/rnaseq-nextflow.git
cd rnaseq-nextflow

# 2. Install Nextflow (requires Java 11+)
curl -s https://get.nextflow.io | bash
./nextflow self-update

# 3. Run with Docker (simplest — no manual installs needed)
nextflow run main.nf \
  --reads   'data/*_{1,2}.fastq.gz' \
  --genome  data/genome.fa \
  --gtf     data/genes.gtf \
  --outdir  results \
  -profile  docker

# 4. Run stub test (no data needed — validates syntax)
nextflow run main.nf -stub-run --profile test
```

---

## 📦 Installation

### Prerequisites

- **Java 11+** (Nextflow requirement)
- **Nextflow ≥ 22.10.1**
- One of: **Docker**, **Singularity**, or **Conda**

```bash
# Install Java (if not present)
sudo apt-get install -y default-jdk    # Ubuntu/Debian
brew install openjdk@17               # macOS

# Install Nextflow
curl -s https://get.nextflow.io | bash
sudo mv nextflow /usr/local/bin/

# Verify
nextflow -version
```

### Execution Backends

```bash
# Option A: Docker
docker pull biocontainers/fastqc:0.12.1--hdfd78af_0

# Option B: Singularity (HPC clusters)
singularity pull docker://biocontainers/fastqc:0.12.1--hdfd78af_0

# Option C: Conda (create env automatically)
conda install -c bioconda nextflow
```

---

## 🚀 Usage

### Basic run (STAR + featureCounts)

```bash
nextflow run main.nf \
  --reads   'path/to/data/*_{1,2}.fastq.gz' \
  --genome  /ref/GRCh38.fa \
  --gtf     /ref/gencode.v44.annotation.gtf \
  --outdir  results \
  -profile  docker
```

### Use pre-built STAR index (skips indexing step)

```bash
nextflow run main.nf \
  --reads      'data/*_{1,2}.fastq.gz' \
  --star_index /ref/STAR_index_GRCh38 \
  --gtf        /ref/genes.gtf \
  --outdir     results \
  -profile     singularity
```

### Salmon quasi-mapping (fast, alignment-free)

```bash
nextflow run main.nf \
  --reads   'data/*_{1,2}.fastq.gz' \
  --genome  /ref/GRCh38.fa \
  --gtf     /ref/genes.gtf \
  --aligner salmon \
  --outdir  results_salmon \
  -profile  docker
```

### Single-end reads

```bash
nextflow run main.nf \
  --reads      'data/*.fastq.gz' \
  --single_end true \
  --genome     /ref/genome.fa \
  --gtf        /ref/genes.gtf \
  --outdir     results_se
```

### Resume a failed run

```bash
nextflow run main.nf [same params] -resume
```

---

## ⚙️ Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--reads` | `null` | Glob pattern to input FASTQs (e.g. `'data/*_{1,2}.fastq.gz'`) |
| `--single_end` | `false` | Set `true` for single-end reads |
| `--genome` | `null` | Path to reference genome FASTA |
| `--gtf` | `null` | Path to genome annotation GTF |
| `--star_index` | `null` | Path to pre-built STAR index (skips indexing) |
| `--aligner` | `star` | Aligner: `star`, `hisat2`, `salmon` |
| `--trimmer` | `trimmomatic` | Trimmer: `trimmomatic`, `fastp` |
| `--quantifier` | `featurecounts` | Quantifier: `featurecounts`, `htseq` |
| `--strandedness` | `auto` | Library strandedness: `auto`, `forward`, `reverse`, `unstranded` |
| `--star_twopass` | `true` | Enable STAR 2-pass alignment |
| `--outdir` | `./results` | Output directory |
| `--skip_fastqc` | `false` | Skip FastQC steps |
| `--skip_trimming` | `false` | Skip trimming step |
| `--skip_multiqc` | `false` | Skip MultiQC report generation |
| `--max_cpus` | `16` | Maximum CPUs per process |
| `--max_memory` | `128.GB` | Maximum memory per process |

---

## 📁 Output Structure

```
results/
├── fastqc/
│   ├── raw/                        # Pre-trimming FastQC reports
│   └── trimmed/                    # Post-trimming FastQC reports
│
├── trimmomatic/  (or fastp/)
│   ├── sample1_paired_1.fastq.gz
│   └── sample1_paired_2.fastq.gz
│
├── star/  (or hisat2/ / salmon/)
│   ├── sample1/
│   │   ├── sample1_Aligned.sortedByCoord.out.bam
│   │   ├── sample1_Log.final.out
│   │   └── sample1_ReadsPerGene.out.tab
│   └── ...
│
├── bam/
│   ├── sample1.sorted.bam
│   ├── sample1.sorted.bam.bai
│   └── flagstat/
│       └── sample1.flagstat
│
├── featurecounts/  (or htseq/)
│   ├── sample1.counts.txt
│   └── ...
│
├── counts/
│   └── count_matrix.tsv            ⭐ Main output: gene × sample matrix
│
├── multiqc/
│   ├── multiqc_report.html         ⭐ Aggregated QC report
│   └── multiqc_data/
│
└── pipeline_info/
    ├── execution_report.html
    ├── execution_timeline.html
    ├── execution_trace.txt
    └── pipeline_dag.svg
```

---

## 📊 Downstream Analysis

After the pipeline finishes, run DESeq2 differential expression analysis:

```bash
# Install R packages (once)
Rscript -e "install.packages(c('BiocManager')); BiocManager::install(c('DESeq2','ggplot2','pheatmap','ggrepel'))"

# Edit contrast in bin/downstream_deseq2.R, then run:
Rscript bin/downstream_deseq2.R
```

Outputs:
- `downstream/deseq2/deseq2_results_all.tsv` — All genes with LFC, p-value, padj
- `downstream/deseq2/deseq2_results_significant.tsv` — Filtered DEGs
- `downstream/deseq2/pca_plot.pdf` — PCA of samples
- `downstream/deseq2/volcano_plot.pdf` — Volcano plot
- `downstream/deseq2/heatmap_top50_DEGs.pdf` — Heatmap of top DEGs

---

## 🖥️ Running on HPC (SLURM)

```bash
# Edit conf/slurm.config with your account name, then:
nextflow run main.nf \
  --reads   'data/*_{1,2}.fastq.gz' \
  --genome  /shared/ref/GRCh38.fa \
  --gtf     /shared/ref/genes.gtf \
  --outdir  /scratch/results \
  -profile  singularity,slurm \
  -resume
```

For very large cohorts (100+ samples), submit the Nextflow driver as a SLURM job:

```bash
#!/bin/bash
#SBATCH --job-name=rnaseq-nf
#SBATCH --time=48:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=2

module load nextflow singularity

nextflow run main.nf \
  --reads   'data/*_{1,2}.fastq.gz' \
  --genome  /shared/ref/GRCh38.fa \
  --gtf     /shared/ref/genes.gtf \
  -profile  singularity,slurm \
  -resume
```

---

## 🧪 Hands-on Tutorial

### Tutorial: Analysing a Public RNA-seq Dataset

**Goal:** Take public SRA data from scratch → produce a count matrix.

#### Step 1: Download test data

```bash
# Install SRA toolkit
conda install -c bioconda sra-tools

# Download 6 samples (≈2GB, takes ~15 min)
bash bin/download_test_data.sh
```

#### Step 2: Download reference files (human GRCh38)

```bash
# Genome FASTA (chromosome 22 only for speed)
wget -O data/chr22.fa.gz \
  "https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.22.fa.gz"
gunzip data/chr22.fa.gz

# GTF annotation
wget -O data/genes.gtf.gz \
  "https://ftp.ensembl.org/pub/release-110/gtf/homo_sapiens/Homo_sapiens.GRCh38.110.chr.gtf.gz"
gunzip data/genes.gtf.gz
```

#### Step 3: Run the pipeline

```bash
nextflow run main.nf \
  --reads       'test/data/fastq/*_{1,2}.fastq.gz' \
  --genome      data/chr22.fa \
  --gtf         data/genes.gtf \
  --aligner     star \
  --trimmer     trimmomatic \
  --quantifier  featurecounts \
  --strandedness reverse \
  --outdir      tutorial_results \
  -profile      docker \
  -resume
```

#### Step 4: Inspect outputs

```bash
# Check alignment rates
cat tutorial_results/star/*/Log.final.out | grep "Uniquely mapped"

# Preview count matrix
head tutorial_results/counts/count_matrix.tsv

# Open MultiQC report
open tutorial_results/multiqc/multiqc_report.html
```

#### Step 5: Differential expression

```bash
Rscript bin/downstream_deseq2.R
```

### Interpreting Results

| Metric | Good Range | Warning |
|--------|-----------|---------|
| Uniquely mapped reads | > 75% | < 60% |
| Reads in genes | > 60% | < 40% |
| Duplication rate | < 50% | > 70% |
| Genes detected | > 15,000 (human) | < 10,000 |

---

## 🤝 Contributing

Contributions are welcome! Please:

1. Fork the repository
2. Create a feature branch: `git checkout -b feature/new-aligner`
3. Add tests for new modules
4. Submit a pull request

### Adding a new module

```bash
# 1. Create the module file
vim modules/local/my_tool.nf

# 2. Import it in the workflow
# Add to workflows/rnaseq.nf:
# include { MY_TOOL } from '../modules/local/my_tool'

# 3. Test with stub mode
nextflow run main.nf -stub-run --profile test
```

---

## 📚 Citations

If you use this pipeline, please cite:

- **Nextflow**: Di Tommaso P, et al. *Nat Biotechnol* (2017). doi:10.1038/nbt.3820
- **STAR**: Dobin A, et al. *Bioinformatics* (2013). doi:10.1093/bioinformatics/bts635
- **HISAT2**: Kim D, et al. *Nat Methods* (2019). doi:10.1038/s41592-019-0344-8
- **Salmon**: Patro R, et al. *Nat Methods* (2017). doi:10.1038/nmeth.4197
- **featureCounts**: Liao Y, et al. *Bioinformatics* (2014). doi:10.1093/bioinformatics/btt656
- **DESeq2**: Love MI, et al. *Genome Biol* (2014). doi:10.1186/s13059-014-0550-8
- **FastQC**: Andrews S. https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
- **MultiQC**: Ewels P, et al. *Bioinformatics* (2016). doi:10.1093/bioinformatics/btw354
- **Trimmomatic**: Bolger AM, et al. *Bioinformatics* (2014). doi:10.1093/bioinformatics/btu170

---

## 📄 License

This project is licensed under the MIT License. See [LICENSE](LICENSE) for details.

---

## 👤 Author

**Dr. Suhirthakumar Puvanendran**  
Bioinformatics | Genomics | Nextflow Pipelines  
🔗 [LinkedIn](https://www.linkedin.com/in/suhirthakumar/) · [GitHub](https://github.com/Suhirthakumar)

---

*Built with ❤️ for reproducible genomics*
