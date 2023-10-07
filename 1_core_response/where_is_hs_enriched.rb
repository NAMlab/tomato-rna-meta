require_relative "../2_curve/ontology_class"
require "set"
require "json"
require "csv"

# As true positives we use GO:0009408	response to heat
# As a similar but "wrong" term we use GO:0009203	ribonucleoside triphosphate catabolic pr...
# For drought stress the true/false controls are GO:0009414 response to water deprivation and GO:0042454 ribonucleoside catabolic process
# For salt stress they are GO:0009651 response to salt stress and GO:0098542 defense response to other organism
target_term = "GO:0009408".to_sym


ontology = Ontology.from_json_file("../2_curve/GO.json")
["e38.3", "e38.2", "e1.4", "e38.1", "e1.3", "e1.2", "e1.1", "e17.1", "e17.2", "e10.2", "e5.1", "e5.2", "e18.1", "robert.5.3", "robert.1.1", "robert.2.1", "robert.3.1", "robert.4.1", "robert.5.1", "e37.3", "e36.1.1", "huji.1", "robert.1.2", "robert.2.2", "robert.3.2", "robert.4.2", "robert.5.2", "e42.1.1"].each do |contrast|
  terms = CSV.open("output/individual_contrasts_supplement/#{contrast}/GO_upregulated.csv").map{|r| r.first.to_sym }[1..-1] # (remove header)
  terms = ontology.set_with_ancestors(terms)
 if terms.include?(target_term)
   puts "#{contrast},YES"
 else
   puts "#{contrast},NO"
 end
end
# Zlib::GzipReader.open("../data/GOMAP_annotation.map.lfs.gz") do |gz|
#   gz.readlines.each do |line|
#     gene, terms = line.chomp.split("\t")
#     terms = Set.new(terms.split(",").map(&:to_sym))
#     if terms.include?(target_term)
#       puts gene
#     end
#   end
# end
