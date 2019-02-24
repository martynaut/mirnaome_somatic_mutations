python3 prepare_confidence_file.py ~/Documents/files_for_mirnaome/confidence.txt \
~/Documents/files_for_mirnaome/confidence_score.txt ~/Documents/files_for_mirnaome/aliases.txt \
~/Documents/files_for_mirnaome/mirna_chromosome_build.txt ~/Documents/files_for_mirnaome/confidence_file.xlsx

python3 prepare_localization_file.py /home/martyna/Documents/ViennaRNA-2.4.11 \
~/Documents/files_for_mirnaome/hairpin3.fa ~/Documents/files_for_mirnaome/hsa2.gff3.txt \
~/Documents/files_for_mirnaome/new_coordinates_all.bed ~/Documents/output/Data_LUAD

python3 prepare_reads_file.py ~/Documents/files_for_mirnaome/hsa2.gff3.txt \
~/Documents/files_for_mirnaome/mirbase_count.csv \
~/Documents/output/Data_LUAD

python3 run_mirnaome_analysis.py ~/dane/HNC/DATA_HNC/ \
 ~/Documents/output/Data_LUAD ~/Documents/files_for_mirnaome/new_coordinates_all_02.bed \
~/Documents/files_for_mirnaome/confidence_file.xlsx ~/Documents/files_for_mirnaome/hsa.gff \
~/Documents/files_for_mirnaome/cancer_exons.txt ~/Documents/output/Data_LUAD/localizations_test.csv \
~/Documents/files_for_mirnaome/mirbase_count.csv
