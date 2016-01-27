#coding: utf-8

#reference = "../../7_prod1/npt.gro"

unless ARGV.length == 1
    usage = <<-EOS
    Script for converting GRO file to PDB with centering the protein.
    usage: ruby #{__FILE__} BASE_DIR
    e.g., ruby #{__FILE__} ../prod1/
    EOS
    puts usage
    exit 1
end

base_dir = ARGV.first

# 1 => all protein atoms
# 2 => non-hydrogen atoms
# 3 => CA atoms
# 4 => backbone atoms etc ...
selection_cen = 2
selection_out = 2

Dir.glob(File.join(base_dir, "*.gro")) do |f|
    basename = File.basename(f, ".gro")
    outfile = "./" + basename + ".pdb"
    tpr_file = f.sub(/\.gro/, ".tpr")
    `gmx trjconv -f #{f} -o #{outfile} -s #{tpr_file} -pbc mol -center <<EOS \n #{selection_cen} \n #{selection_out} \n EOS`
end
