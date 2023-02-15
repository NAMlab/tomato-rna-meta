# This script determines how many proteins are present in tomato for each of the families listed in
# any of the JSON files in input/
require "json"
require "csv"
require "net/http"

families = Dir.glob("input/*.json").map do |f|
  json = JSON.parse(File.read(f))
  json.map {|k, v| v.map{|a| a.first }}.flatten
end
families = families.flatten.sort.uniq

CSV.open("output/family_sizes.csv", "wb") do |csv|
  csv << ["interpro.id", "n.tomato.proteins"]
  families.each do |f|
    puts f
    response = Net::HTTP.get(URI("https://www.ebi.ac.uk/interpro/api/entry/InterPro/#{f}/protein/taxonomy/uniprot/4081/?format=json"))
    if response
      json = JSON.parse(response)
      csv << [f, json["metadata"]["counters"]["proteins"]]
    else
      # "204 No Content"
      csv << [f, 0]
    end
  end
end
