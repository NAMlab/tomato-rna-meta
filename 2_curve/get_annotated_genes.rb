require_relative "ontology_class"
require "set"
require "json"
require "zlib"

# As true positives we use GO:0009408	response to heat
# As a similar but "wrong" term we use GO:0048868	pollen tube development
target_term = ARGV.first.to_sym

ontology = Ontology.from_json_file("GO.json")
Zlib::GzipReader.open("../data/GOMAP_annotation.map.lfs.gz") do |gz|
  gz.readlines.each do |line|
    gene, terms = line.chomp.split("\t")
    terms = Set.new(terms.split(",").map(&:to_sym))
    terms = ontology.set_with_ancestors(terms)
    if terms.include?(target_term)
      puts gene
    end
  end
end
