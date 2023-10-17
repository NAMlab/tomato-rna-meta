require "json"

puts "target_id","family_id"
obj = JSON.parse(File.read(ARGV.first))
obj.each do |gene, fams|
  fams.each do |fam|
    puts(gene + "," + fam.first)
  end
end
