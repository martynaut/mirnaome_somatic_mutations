# miRNAome somatic mutations

Project with scripts to analyse somatic mutations found in cancer (LUAD and LUSC) patients 
using somatic mutation data from TCGA

## Description

Python scripts may be reused for other data sources with input data prepared as 
in TCGA (https://cancergenome.nih.gov/).

## Pre-run preparation

All needed Python libraries are gathered in `requirements.txt` (Python 3 is needed).

### Preparing ViennaRNA

Download ViennaRNA distribution from https://www.tbi.univie.ac.at/RNA/#download 
to chosen folder
and install using instructions:
```bash
tar -zxvf ViennaRNA-2.4.11.tar.gz
cd ViennaRNA-2.4.11
./configure
make
sudo make install
```

### Input data description

1) Coordinates

    `\t` separated file without header with values like:

    ```chromosome\tstart\tstop\tsequence_name\n```

    example:

    ```
    chr1	17343	17456	hsa-mir-6859-1
    chr1	30382	30483	hsa-mir-1302-2
    chr1	632306	632428	hsa-mir-6723
    ```
2) Confidence file

    Confidence file may be prepared with a script `prepare_confidence_file.py`
    based on four files downloaded from miRBase `confidence.txt`,
    `confidence_score.txt`,  `aliases.txt` and `mirna_chromosome_build`.
    
    To run this script use:
    ```
    python3 prepare_confidence_file.py ~/Documents/files_for_mirnaome/confidence.txt \
    ~/Documents/files_for_mirnaome/confidence_score.txt \
    ~/Documents/files_for_mirnaome/aliases.txt ~/Documents/files_for_mirnaome/mirna_chromosome_build.txt \
    ~/Documents/files_for_mirnaome/confidence_file.xlsx
    ```
    
3) Localization file

    Localization file may be prepared with a script `prepare_localization_file.py`
    based on two files downloaded from miRBase `hairpin.fa`,
    `hsa.gff` (here saved as `hsa.gff.txt`) and coordinates file. ViennaRNA is needed
    (for installation see above). Add **absolute** path to Vienna package 
    
    To run this script use:
    ```
    python3 prepare_localization_file.py /home/<user>/Documents/ViennaRNA-2.4.11 ~/Documents/files_for_mirnaome/hairpin.fa \
    ~/Documents/files_for_mirnaome/hsa.gff3.txt\
    ~/Documents/files_for_mirnaome/new_coordinates_all_02.bed ~/Documents/output/Data_LUAD
    ```
    
4) mirgenedb file
    
    Genomic coordinates from mirgenedb (http://mirgenedb.org/download) for human.
    `hsa.gff` file
    
5) Cancer exons
    Text file with names of exons that should be included in 
    the first steps of the analysis (first mutations extraction
    from vcf files).
    
    Example:
    ```bash
    EGFR_chr7.20
    EGFR_chr7.21
    EGFR_chr7.22
    ``` 
    
### Output data description

1) 

### How to use it

Example run is prepared in `run_mirnaome.sh` bash script.

To prepare confidence data run

```bash
python3 prepare_confidence_file.py ~/Documents/files_for_mirnaome/confidence.txt \
~/Documents/files_for_mirnaome/confidence_score.txt \
~/Documents/files_for_mirnaome/aliases.txt ~/Documents/files_for_mirnaome/mirna_chromosome_build.txt \
~/Documents/files_for_mirnaome/confidence_file.xlsx
```

Confidence file can be found in defined directory under defined filename.

To prepare localization data run

```bash
python3 prepare_localization_file.py ~/Documents/ViennaRNA-2.4.11 ~/Documents/files_for_mirnaome/hairpin.fa \
~/Documents/files_for_mirnaome/hsa.gff3.txt\
~/Documents/files_for_mirnaome/new_coordinates_all_02.bed ~/Documents/output/Data_LUAD
```

Localization file can be found in output folder.

To run analysis run

```
python3 run_mirnaome_analysis.py ~/dane/HNC/DATA_HNC ~/dane/HNC/RESULTS_HNC ./Reference/new_coordinates_all_02.bed \
./Reference/confidence.xlsx ./Reference/hsa.gff ./Reference/cancer_exons.txt
```

## Authors

Paulina Galka-Marciniak, Martyna O. Urbanek-Trzeciak, Piotr Kozlowski

Institute of Bioorganic Chemistry, Polish Academy of Sciences, Noskowskiego 12/14, 61-704, 
Poznan, Poland

## Citation

TBD


Biorxiv link to preprint: 
https://www.biorxiv.org/

## Contact

For any issues, please create a GitHub Issue.
