# Pass this script the path to a protein FASTA file and it will print the deepest HOG for each protein on the screen.
# This requires the FASTA file to contain each protein in a single line.
require "csv"
require "net/http"
require "json"

puts "target,hogs_with_deepest_levels"

current_name = nil
current_sequence = ""

File.foreach(ARGV.first) do |line|
  if line.start_with?(">")
    unless current_name.nil?
      print current_name + ","
      begin
        json = JSON.parse(Net::HTTP.get(URI("https://omabrowser.org/api/sequence/?format=json&query=" + current_sequence)))
        if json.has_key?("targets")
          hogs = json["targets"].select { |t| t["species"]["taxon_id"] == 4081 }.map { |t| t["oma_hog_id"].split(".").first }.uniq.compact

          hogs.each do |hog|
            json = JSON.parse(Net::HTTP.get(URI("https://omabrowser.org/api/hog/#{hog}/members/?format=json")))
            print(hog + "->" + json["level"] + " ")
          end
        else
          print("ERROR: #{json['detail']}")
        end
      rescue Net::ReadTimeout => e
        print "TIMEOUT"
      end
      puts
    end
    current_name = line[1..-1].split.first
    current_sequence = ""
  else
    current_sequence += line.chomp
  end
end
