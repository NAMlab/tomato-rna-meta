require "csv"

levels = IO.readlines("clades_hierarchy.txt").map(&:chomp)

puts "target,deepest_level,set"
CSV.foreach("results/hs_core_hogs.csv", headers: true) do |row|
  if row[1]
    level = row[1].split(" ").map{ |h| h.split("->")[1] }.map{ |h| levels.index(h) }.min
  else
    level = levels.length - 1
  end
  puts "#{row[0]},#{level},hs_core"
end

CSV.foreach("results/random_hogs.csv", headers: true) do |row|
  if row[1]
    level = row[1].split(" ").map{ |h| h.split("->")[1] }.map{ |h| levels.index(h) }.min
  else
    level = levels.length - 1
  end
  puts "#{row[0]},#{level},random"
end
