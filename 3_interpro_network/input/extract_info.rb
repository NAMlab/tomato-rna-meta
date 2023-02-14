# This script extracts only family and superfamily entries from an InterPro sequence search result and 
# outputs them in JSON format.
require "json"

obj = JSON.parse(File.read(ARGV.first))
families = {}
obj["results"].each do |s|
  entries = s["matches"].map{|m| m["signature"]}.select{|m| !m["entry"].nil? and ["FAMILY", "HOMOLOGOUS_SUPERFAMILY"].include?(m["entry"]["type"])}.map{|m| m["entry"]}
  families[s["xref"].first["id"]] = entries.map{|e| [e["accession"], e["name"], e["description"], e["type"]]}.uniq
end

puts JSON.pretty_generate(families)
