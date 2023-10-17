module load interproscan/5.55-88.0
interproscan.sh -i ITAG4.1_proteins_without_asterisks.fasta -appl TIGRFAM,SFLD,SUPERFAMILY,PANTHER,ProSiteProfiles,CDD,PRINTS,PIRSR,ProSitePatterns,AntiFam,Pfam -f json -cpu 20 -dp
