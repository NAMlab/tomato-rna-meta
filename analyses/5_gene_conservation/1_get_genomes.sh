#!/bin/bash
# Download and prepare the 14 de-novo assemblies from the Lippman 100 genomes project
# ( https://doi.org/10.1016/j.cell.2020.05.021 ) as well as the ITAG4.1 reference genome.
# The script downloads the assemblies and structural annotations and extracts transcript, cds
# and protein sequences for each gene model, creating one respective fasta file each.
# Usage: bash 1_get_genomes.sh

requiredCommands=(
  wget
  gzip
  gunzip
  gffread
)

for c in ${requiredCommands[@]}; do
  if ! command -v $c &> /dev/null
  then
      echo "Error: $c could not be found."
      exit 1
  fi
done

accessions=(
  BGV006775
  BGV006865
  BGV007931
  BGV007989
  Brandywine
  EA00371
  EA00990
  Fla.8924
  Floradade
  LYC1410
  M82
  PAS014479
  PI169588
  PI303721
)

mkdir -p proteins/ cds/ transcripts/

for a in ${accessions[@]}; do
  echo "Downloading accession $a ..."
  wget "https://solgenomics.net/ftp/genomes/tomato100/March_02_2020_sv_landscape/${a}_MAS2.0.fasta.gz"
  wget "https://solgenomics.net/ftp/genomes/tomato100/March_02_2020_sv_landscape/${a}_MAS2.0_gene_models.gff.gz"
  echo "Generating sequences for accession $a ..."
  gunzip *.gz
  fa_file="${a}_MAS2.0.fasta"
  gff_file="${a}_MAS2.0_gene_models.gff"
  gffread -y proteins/${a}.fa -g $fa_file $gff_file
  gzip proteins/${a}.fa
  gffread -x cds/${a}.fa -g $fa_file $gff_file
  gzip cds/${a}.fa
  gffread -w transcripts/${a}.fa -g $fa_file $gff_file
  gzip transcripts/${a}.fa
  rm $fa_file $gff_file ${fa_file}.fai
done

# Now we do the same with ITAG4.1
echo "Downloading ITAG4.1..."
wget "https://solgenomics.net/ftp/tomato_genome/assembly/build_4.00/S_lycopersicum_chromosomes.4.00.fa.gz"
wget "https://solgenomics.net/ftp/tomato_genome/annotation/ITAG4.1_release/ITAG4.1_gene_models.gff"
echo "Generating sequences for accession $a ..."
gunzip *.gz
fa_file="S_lycopersicum_chromosomes.4.00.fa"
gff_file="ITAG4.1_gene_models.gff"
a="ITAG4.1"
gffread -y proteins/${a}.fa -g $fa_file $gff_file
gzip proteins/${a}.fa
gffread -x cds/${a}.fa -g $fa_file $gff_file
gzip cds/${a}.fa
gffread -w transcripts/${a}.fa -g $fa_file $gff_file
gzip transcripts/${a}.fa

rm $fa_file $gff_file ${fa_file}.fai
