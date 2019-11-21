#!/bin/bash

#script for pre-processing fungal ITS amplicon sequencing using Qiime 2 v2.0
#originally written by M. Amine Hassani - ahassani@bot.uni-kiel.de - Jan2019
#these scripts were run on the Kiel University Computing Centre rzcluster 

#SBATCH --job-name=q2.1sp.ITS
#SBATCH --mail-user=ahassani@bot.uni-kiel.de
#SBATCH --mail-type=ALL
#SBATCH --output=q2_1sp_ITS.out
#SBATCH --error=q2_1sp_ITS.err
#SBATCH --nodes=2
#SBATCH --tasks-per-node=6
#SBATCH --cpus-per-task=2
#SBATCH --mem=50000
#SBATCH --time=240:00:00
#SBATCH --partition=angus

# loading the configuration file 
config_file=$1
source $config_file

# loading modules
module load qiime2-2018.11
source activate qiime2-2018.11
module load intelmpi16.0.0 
export OMP_NUM_THREADS="$nbr_thr"

# clearing files
rm -fr "$path_to_data/output_*" "$path_to_data/mapping_files" "$path_to_data/export"

# making new directories
mkdir -p "$path_to_data/mapping_files"
mkdir -p "$path_to_data/output_sp001"   	 

# moving mapping file
cp "$path_to_raw/sp001/sample_metadata.tsv" "$path_to_data/mapping_files"
mv "$path_to_raw/sp001/sample_metadata.tsv" "$path_to_data/output_sp001"
	
	
echo -e "\n 1-importing sequences..."
	
	qiime tools import --type EMPPairedEndSequences \
   	 --input-path "$path_to_raw/sp001" \
   	 --output-path "$path_to_data/output_sp001/emp-paired-end-sequences_sp001.qza"
	
echo -e "\n emp-paired-end-sequences_sp001.qza saved in $path_to_data/output_sp001" 
	
mv "$path_to_data/output_sp001/sample_metadata.tsv" "$path_to_raw/sp001"

echo -e "\n 2-demultiplexing sequences..."

	qiime demux emp-paired \
	  --m-barcodes-file "$path_to_raw/sp001/sample_metadata.tsv" \
	  --m-barcodes-column BarcodeSequence \
	  --i-seqs "$path_to_data/output_sp001/emp-paired-end-sequences_sp001.qza" \
	  --o-per-sample-sequences "$path_to_data/output_sp001/demux_sp001" \
	  --p-rev-comp-mapping-barcodes

echo -e "\n 3-Trimming reads with ITSexpress..."	

# did not run
#	qiime itsxpress trim-pair-output-unmerged\
#	  --i-per-sample-sequences "$path_to_data/output_sp001/demux_sp001.qza" \
#	  --p-region ITS1 \
#	  --p-taxa F \
#	  --p-threads 2 \
#	  --o-trimmed "$path_to_data/output_sp001/demux_sp001_trimmed.qza"	

echo -e "\n 4-sequence QC and generating feature table and representative sequences - sp001 ..."

	mpirun -np 12 qiime dada2 denoise-paired \
         --i-demultiplexed-seqs "$path_to_data/output_sp001/demux_sp001.qza" \
         --p-trunc-len-f "$len" \
	 --p-trunc-len-r "$len" \
         --p-n-threads "$nbr_thr" \
	 --p-n-reads-learn 9999  \
	 --o-denoising-stats "$path_to_data/output_sp001/denoise_stat" \
         --o-table "$path_to_data/output_sp001/table_sp001" \
         --o-representative-sequences "$path_to_data/output_sp001/rep-seqs_sp001"

echo -e "\n 5-taxonomic classification of representative sequences ..."

	mpirun -np 4 qiime feature-classifier classify-sklearn  --i-classifier "$class_path" \
 	  --i-reads "$path_to_data/output_sp001/rep-seqs_sp001.qza" \
 	  --o-classification "$path_to_data/output_sp001/taxonomy.qza" --verbose \
	  --p-reads-per-batch 0  --p-n-jobs 4 --p-pre-dispatch 3*n_jobs

mkdir -p "$path_to_data/output_summary"

echo -e "\n 6-summarizing demux file ..."

	qiime demux summarize \
	  --i-data "$path_to_data/output_sp001/demux_sp001.qza" \
	  --o-visualization "$path_to_data/output_summary/demux_sp001.qzv"	

echo -e "\n 7-summarizing feature table ..."

	qiime feature-table summarize \
	  --i-table "$path_to_data/output_sp001/table_sp001.qza" \
	  --o-visualization "$path_to_data/output_summary/table_sp001.qzv" \
	  --m-sample-metadata-file "$path_to_data/mapping_files/sample_metadata.tsv"

echo -e "\n 8-summarizing sequences ..."

	qiime feature-table tabulate-seqs \
	  --i-data "$path_to_data/output_sp001/rep-seqs_sp001.qza" \
	  --o-visualization "$path_to_data/output_summary/rep-seqs.qzv"

echo -e "\n 9-summarizing taxonomy ..."

	qiime metadata tabulate \
	  --m-input-file "$path_to_data/output_sp001/taxonomy.qza" \
	  --o-visualization "$path_to_data/output_summary/taxonomy.qzv"


echo -e "\n 10-exporting .qza files..."

mkdir -p "path_to_data/export"
qiime tools export --input-path "$path_to_data/output_sp001/taxonomy.qza" --output-path "$path_to_data/export"
qiime tools export --input-path "$path_to_data/output_sp001/table_sp001.qza" --output-path "$path_to_data/export"
sed -i '1s/^/\#/' $path_to_data/output_sp001/export/taxonomy.tsv
biom add-metadata -i "$path_to_data/export/feature-table.biom" --output-fp "$path_to_data/export/table-with-taxonomy.biom" --observation-metadata-fp "$path_to_data/export/taxonomy.tsv" --sc-separated taxonomy
biom convert -i "$path_to_data/export/table-with-taxonomy.biom" -o "$path_to_data/export/table-with-taxonomy.txt" --to-tsv


echo -e "\n 11-archiving export & summary ..."

tar czvf  "$path_to_data/archive_export.tar.gz" export/
tar czvf  "$path_to_data/archive_summary.tar.gz" output_summary/

echo -e "\n*** done ***\n"
