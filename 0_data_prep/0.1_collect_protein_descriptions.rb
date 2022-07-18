# This script combines the ITAG4.1 descriptions for each gene with Arabidopsis homologues from OMA (SL3.0!!) 
# and respective descriptions for each of the Arabidopsis homologues from TAIR
require "pry"
require "csv"
require "zlib"

arath_descriptions = Zlib::GzipReader.open("input/TAIR10_functional_descriptions.tsv.gz").readlines.map do |l|
  p = l.split("\t")
  [p[0].split(".")[0], p[2]]
end.to_h

oma_orthologues = {}
Zlib::GzipReader.open("input/oma_orthologues.tsv.gz").each_line do |l|
  p = l.split("\t")
  g = p[0].split(".")[0]
  a = p[1]
  if oma_orthologues.has_key? g
    oma_orthologues[g] << a
  else
    oma_orthologues[g] = [a]
  end
end

Zlib::GzipWriter.open('output/protein_descriptions.lfs.csv.gz') do |gz|
  csv_content = CSV.generate do |csv|
    csv << ["gene", "ITAG4.1_description", "OMA_orthologues"]
    Zlib::GzipReader.open("input/ITAG4.1_descriptions.txt.gz").each_line do |l|
      p = l.split(" ")
      g = p[0].split(".")[0]
      d = p[1..-1].join(" ")
      if oma_orthologues.has_key? g
        orthos = oma_orthologues[g].map do |o|
          "#{o}(#{arath_descriptions[o]})"
        end.join(", ")
      else
        orthos = ""
      end
      csv << [g, d, orthos]
    end
  end
  gz.write csv_content
end


