#!/bin/bash
# Get all sequences matching a specific ID from a given folder of FASTA files
# Usage:   bash 2_get_seqs.sh <transcript_id> <folder_with_fastas>
# Example: bash 2_get_seqs.sh "Solyc00g500372.1" proteins/

requiredCommands=(
  seqkit
  sed
)

for c in ${requiredCommands[@]}; do
  if ! command -v $c &> /dev/null
  then
      echo "Error: $c could not be found."
      exit 1
  fi
done

if [ ! -d "$2" ]; then
  echo "Error: Folder '$2' does not exist."
  exit 2
fi

for f in $2/*.fa.gz; do
  bn=$(basename $f)
  accession="${bn%.*.*}"
  seqkit grep -p ".*${1}.*" -r $f | sed "s/>/>${accession}:/"
done
