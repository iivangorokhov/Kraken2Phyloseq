# Python script for creating OTU and taxonomy tables in Phyloseq format from Kraken2 reports.

This script accepts folder with Kraken2 reports and transforms content of all reports into 2 tables: 
1. **OTU table**: contains column for unique OTU ID and columns with counts of this OTU in each Kraken2 report (sample).
2. **Taxonomy table**: contains column for unique OTU ID and taxonomy lineage of this OTU in next columns.

This allows for further data metagenomic data analysis, samples comparison and statistical testing. 

**Note**: OTUs with underassigned taxonomy (for example, terminated at Genus level) are added to resulting tables as separate OTUs, though their taxonomic lineage is defined only down to Genus level (Species = None). The same is true for missing taxonomic levels (for example, when there is no Family-level taxonomic assignment).

## Requirements

* Python >= 3.12
* external packages: `pandas`, `tqdm` (see `requirements.txt`)

## Installation
    
1. Download git repository:
```bash
git clone https://github.com/iivangorokhov/Kraken2Phyloseq
cd Kraken2Phyloseq
```

2. Create conda environment to install dependencies and to avoid conflicts:
```bash
conda create --name kraken2phyloseq python=3.12 -y
conda activate kraken2phyloseq
pip install -r requirements.txt
```

## Usage

Make the script executable *(if required)*:
```bash
chmod +x Kraken2Phyloseq.py
```

Example command:
```bash
./Kraken2Phyloseq.py \
--inputFolder /path/to/folder/with/Kraken2/reports/ \
--otu /path/to/saved/OTU_table.csv \
--taxonomy /path/to/saved/Taxonomy_table.csv \
--format default \
--extension txt
```

## Arguments help

* `--inputFolder`: Folder with Kraken2 reports (or Bracken, considering they share the same format).
* `--otu`: Path to save the OTU table. *Default:* `./otu_table.csv` (Saved as a comma-separated text file).
* `--taxonomy`: Path to save the taxonomy table. *Default:* `./taxonomy_table.csv` (Saved as a comma-separated text file).
* `--format`: Format of the created Kraken2 report. Options: `default`, `minimizer` (with `--report-minimizer-data` option), or `mpa` (with `--use-mpa-style` option). *Default:* `default`.
* `--extension`: File extension of Kraken2 reports. *Default:* `txt`.

## Citations

1. **Kraken2**:
Wood DE, Lu J, Langmead B. Improved metagenomic analysis with Kraken 2. Genome Biol. 2019 Nov 28;20(1):257. doi: 10.1186/s13059-019-1891-0. PMID: 31779668; PMCID: PMC6883579.

2. **Phyloseq**:
McMurdie PJ, Holmes S. Phyloseq: a bioconductor package for handling and analysis of high-throughput phylogenetic sequence data. Pac Symp Biocomput. 2012:235-46. PMID: 22174279; PMCID: PMC3357092.