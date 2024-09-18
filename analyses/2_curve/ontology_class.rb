class Ontology
  def initialize(parents, main_id, obsolete)
    @parents  = parents
    @main_id  = main_id
    @obsolete = obsolete
    @ancestors_of = {}
  end

  def self.from_json(json_string)
    json_object = JSON.parse(json_string)
    self.new(
      json_object["parents"].each_with_object({}) {| (k, v), a| a[k.to_sym] = v.map(&:to_sym) },  
      json_object["main_id"].each_with_object({}) {| (k, v), a| a[k.to_sym] = v.to_sym }, 
      json_object["obsolete"].map(&:to_sym)
    )
  end

  def self.from_json_file(filepath)
    if filepath.end_with? (".gz")
      Zlib::GzipReader.open(filepath) do |gz|
        return self.from_json(gz.read)
      end
    else
      return self.from_json(File.read(filepath))
    end
  end

  def is_obsolete?(go_id)
    @obsolete.include? go_id
  end

  def main_id_of(go_id)
    @main_id.include?(go_id) ? @main_id[go_id] : go_id
  end

  def ancestors_of(go_id)
    return @ancestors_of[go_id] if @ancestors_of.keys.include? go_id # We don't need to calculate it if we already have done so before
    ancestors = Set.new
    unless @parents.keys.include? go_id
      puts "Warning: not in Ontology: #{go_id}"
      return ancestors
    end
    @parents[go_id].each do |p|
      if ancestors.add?(p)
        ancestors.merge(ancestors_of(p))
      end
    end
    @ancestors_of[go_id] = ancestors
    ancestors
  end

  def set_with_ancestors(set)
    set.inject(Set.new) do |s, term| 
      if s.add? term # returns nil if term already present, then also all ancestors are already present
        s.merge(ancestors_of(term))
      else
        s
      end
    end
  end

end
