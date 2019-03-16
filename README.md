# miRNome somatic mutations

Project containing Python scripts used to analyse somatic mutations in cancer (LUAD and LUSC) patients 
using somatic mutation data from TCGA

Python scripts may be reused for other data sources with input data prepared as somatic mutation data 
in TCGA (https://cancergenome.nih.gov/). 

Results of all four algorithms available in TCGA database (muse, mutect2, somaticsniper, varscan2) were used.

For conditions to reuse of these scripts please refer to `LICENSE` file.

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

    Excel file with columns: `name, score, start, stop, Strand, id, confidence`
    where name is `name` of miRNA, `score` is mirBase confidence score,
    `start` and `stop` are coordinates, `Strand` is strand with `+` and `-` values,
    `id` id mirBase_ID and `confidence` is miRBase confidence label with `High` and `Low` values.
    
    Example row: `hsa-mir-1234, 0, 144400086, 144400165, -, MI0006324, Low`

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
    
3) mirgenedb file
    
    Genomic coordinates from mirgenedb (http://mirgenedb.org/download) for human.
    `hsa.gff` file
    
4) Cancer exons

    Text file with names of exons that should be included in 
    the first steps of the analysis (first mutations extraction
    from vcf files, are included in coordinates file) but are not miRNA genes.
    Single exon name in line.
    
    Example:
    ```bash
    EGFR_chr7.20
    EGFR_chr7.21
    EGFR_chr7.22
    ``` 
     
5) Localization file

    Excel file with columns: `chrom, name, start, stop, orientation, based_on_coordinates, arm, type`
    where name is `chrom` is chromosome id, `name` of miRNA localization, 
    `start` and `stop` are coordinates, `orientation` is strand with `+` and `-` values,
    `based_on_coordinates` states if localization is based on miRBase coordinates or
     structure prediction, `arm` is which arm of miRNA this sequence is on and `type` is type
     of localization within miRNA precursor.

    Example row: `chr1, hsa-mir-6859-1-3p_post-seed, 17369, 17383, -, yes, 3p, post-seed`

    Localization file may be prepared with a script `prepare_localization_file.py`
    based on two files downloaded from miRBase `hairpin.fa`,
    `hsa.gff` (here saved as `hsa.gff.txt` to differentiate from `hsa.gff` from mirgendb) and coordinates file. ViennaRNA is needed
    (for installation see above). Important: add **absolute** path to Vienna package 
    
    To run this script use:
    ```
    python3 prepare_localization_file.py /home/<user>/Documents/ViennaRNA-2.4.11 ~/Documents/files_for_mirnaome/hairpin.fa \
    ~/Documents/files_for_mirnaome/hsa.gff3.txt\
    ~/Documents/files_for_mirnaome/new_coordinates_all_02.bed ~/Documents/output/Data_LUAD
    ```
    
6) (optional) Chromosome file

    Tab-separated file with regions covered by NGS probes with columns 
    `TargetID, Interval, Regions, Size, Databases, Coverage, HighCoverage, LowCoverage`
    
    Example row: 
    `HSA-LET-7A-1, chr9:96938239-96938318, 1, 80, CustomRegion, 100.0, 1, 0`
    
### Output data description

1) **temp folder**

    Contains merged vcf files if there were multiple samples available per single patient.
    
2) **temp_reference folder**
    
    Temporary files created in create localization file script.
    
3) **files_summary_count_per_patient.csv**

    File with information how many files there are per patient per algorithm. Sanity check: if the merging
    of vcf files was successful, we should have only ones.
    
4) **files_summary.csv**

    Files summary including user_id, file localization, sample id name and aliQ,
    and algorithm used.
    
5) **files_count_per_type.csv**

    Count of files per algorithm used.
    
6) **not_unique_patients.csv**

    Patients for which we had multiple files (multiple samples) that were combined in a single vcf.
    
7) **do_not_use.txt**

    Files that were replaced with merged vcf files (stored in temp folder)
    that will not be used in next steps.
    
8) **results_muse.csv, results_mutect2.csv, results_somaticsniper.csv, varscan2.csv**

    Mutations found in vcf files obtained in each of four algorithms.

9) **results_muse_eval.csv, results_mutect2_eval.csv, results_somaticsniper_eval.csv,
varscan2_eval.csv**
    
    Mutations found in vcf files obtained in each of four algorithms after additional evaluation methods
    based on read counts, SSC, BQ and QSS.

10) **all_mutations_filtered.csv**

    Mutations found by each of four algorithms concatenated.

11) **all_mutations_filtered_merge_programs.csv**

    All mutations grouped to deduplicate mutations found by multiple algorithms.

12) **all_mutations_filtered_mut_type_gene.csv**

    All mutations within miRNA genes with gene information, localization and mutation type.

13) **complex.csv**

    Complex mutations are multiple mutations in single miRNA in single patient.
    Column `complex` is `1` if mutation is treated as complex.

14) **miRNA_per_chromosome.csv**

    How many miRNAs were mutated on single chromosome and how many mutations were found
    in total on each chromosome.

15) **occur.csv**

    How many total mutations, unique mutations and patients with mutation found per
    gene.

16) **distinct.csv**

    Unique mutations description with information how many patients had unique mutations.
    
17) **distinct_with_loc.csv**

    Unique mutations description with localization within miRNA precursor.
    
18) (optional) **mirnas_outside_probes.csv**

    Mutations in what miRNAs were detected in vcf files outside probes-defined regions.
    
19) (optional) **high_confidence_mirnas_per_chrom.csv** 

    miRNAs count per chromosome.
    
20) (optional) **mirnas_per_chrom.csv** 

    miRNAs mutated (found in vscf) count per chromosome.

21) (optional) **patients_per_chrom.csv** 

    Patients that had mutations in each chromosome.
    
    
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

See additional features running  

```bash
python3 run_mirnaome_analysis.py --help
```

To skip steps of the analysis (if first steps were already completed) use `-s` argument adding step
from which script should start.

To include chromosome analysis use `-c` argument adding chromosome file path.

## Authors

Martyna O. Urbanek-Trzeciak, Paulina Galka-Marciniak, Piotr Kozlowski

Institute of Bioorganic Chemistry, Polish Academy of Sciences, Noskowskiego 12/14, 61-704, 
Poznan, Poland

## Citation

**Somatic mutations in miRNA genes in lung cancer – potential functional consequences of non-coding sequence variants**


Paulina Galka-Marciniak<sup>1</sup>, Martyna Olga Urbanek-Trzeciak<sup>1</sup>, Paulina Maria Nawrocka<sup>1</sup>, 
Agata Dutkiewicz<sup>1</sup>, Maciej Giefing<sup>2</sup>, Marzena Anna Lewandowska<sup>3,4</sup>, 
and Piotr Kozlowski<sup>1</sup>


<sup>1</sup> Institute of Bioorganic Chemistry, Polish Academy of Sciences,  Poznan, Poland

<sup>2</sup> Institute of Human Genetics, Polish Academy of Sciences, Poznan, Poland

<sup>3</sup> The F. Lukaszczyk Oncology Center, Department of Molecular Oncology and Genetics, Bydgoszcz, Poland

<sup>4</sup> The Ludwik Rydygier Collegium Medicum, Department of Thoracic Surgery and Tumours, Nicolaus Copernicus University, Bydgoszcz, Poland

**Biorxiv link to preprint:** 

http://biorxiv.org/cgi/content/short/579011v1

## Contact

For any issues, please create a GitHub Issue.
