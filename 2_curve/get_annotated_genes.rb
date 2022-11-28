require_relative "ontology_class"
require "set"
require "json"
require "zlib"

# As true positives we use GO:0009408	response to heat
# As a similar but "wrong" term we use GO:0009203	ribonucleoside triphosphate catabolic pr...
# For drought stress the true/false controls are GO:0009414 response to water deprivation and GO:0042454 ribonucleoside catabolic process
# For salt stress they are GO:0009651 response to salt stress and GO:0098542 defense response to other organism
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
