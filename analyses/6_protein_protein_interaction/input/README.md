# Network construction
For this analysis, some pre-processing was done outside this repository because the files involved are rather large in size.
The result of this pre-processing is the file `string_edges.lfs.txt.gz` which contains all not-just-predicted protein-protein associations from the [STRING database](https://string-db.org/) in tomato, mapped to their ITAG4.1 protein ids.
The steps for producing this file are described below:

## 1. Map STRING IDs to ITAG4.1 IDs
`string_itag_mapping.txt.gz` contains a mapping which gene models in ITAG4.1 the STRING IDs correspond to.
It was created by blasting each *Solanum lycopersicum* STRING v12.0 protein against the ITAG4.1 proteome (using blast+ version 2.13):

```bash
#!/bin/bash
wget https://stringdb-downloads.org/download/protein.sequences.v12.0/4081.protein.sequences.v12.0.fa.gz
wget https://stringdb-downloads.org/download/protein.sequences.v12.0/4081.protein.links.detailed.v12.0.txt.gz # for later
wget https://solgenomics.net/ftp/tomato_genome/annotation/ITAG4.1_release/ITAG4.1_proteins.fasta
makeblastdb -in ITAG4.1_proteins.fasta -title ITAG4.1 -dbtype prot -out ITAG4.1.blastDB -parse_seqids
gunzip 4081.protein.sequences.v12.0.fa.gz
blastp -query 4081.protein.sequences.v12.0.fa -db ITAG4.1.blastDB -out blast_results.txt -outfmt 6 -evalue 1e-50
```

Blast results were then filtered to return the single best hit for each query so we get a 1:1 mapping (or none if there was no hit with at least 95% sequence identity):

```ruby
#!/bin/ruby
require 'csv'

current_query = nil
current_best_hit = nil
current_score = 0

puts "string itag"
CSV.foreach("blast_results.txt", col_sep: "\t") do |row|
    if current_query != row[0]
        if not current_best_hit.nil? and current_score >= 95.0
            puts "#{current_query} #{current_best_hit}"
        end
        current_query = row[0]
        current_score = 0
    end
    if row[2].to_f > current_score
        current_best_hit = row[1]
        current_score = row[2].to_f
    end
end
if current_score >= 95.0
    puts "#{current_query} #{current_best_hit}"
end
```

## 2. Filter Network and replace IDs
Protein-protein interactions from STRING were filtered to only keep those that had some not-just-predictive evidence, i.e. a non-zero score in at least one of the categories experimental evidence, curated database, or abstract-based text-mining.
The summed evidence score of these three categories is used as the association strength in our analysis.

```R
edges = read.csv("4081.protein.links.detailed.v12.0.txt.gz", sep="") # downloaded above
edges$score = edges$experimental + edges$database + edges$textmining
edges = edges[edges$score > 0, c("protein1", "protein2", "score")]

string_itag_mapping = read.csv("string_itag_mapping.txt.gz", sep="")
m = string_itag_mapping$itag
names(m) = string_itag_mapping$string
# Replace protein IDs with corresponding ITAG4.1 ID if there is one
edges$p1 = ifelse(!is.na(m[edges$protein1]), m[edges$protein1], edges$protein1)
edges$p2 = ifelse(!is.na(m[edges$protein2]), m[edges$protein2], edges$protein2)
edges = edges[c("p1", "p2", "score")]
write.table(edges, "string_edges.txt", row.names=F, quote=F, sep=" ")
```

# Validated and random genes
`validated_heat_proteins.txt` contains a list of proteins which are known to be involved in heat stress tolerance in tomato according to https://doi.org/10.1016/j.scienta.2023.112435 .
The gene IDs in the paper were mapped to the respective canonical protein IDs in ITAG4.1.
One further adjustment was made: The paper lists Solyc01g095320 as a confirmed gene and https://doi.org/10.3390/antiox11081467 as the source but when checking the reference, it was actually confirming Solyc10g084170 (see supplemmentary table). So this was corrected here.

`random_proteins.txt` contains 200 random protein IDs from the ITAG4.1 genome to be used as a comparison in our analyses.