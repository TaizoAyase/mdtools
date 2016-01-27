#coding: utf-8

##########################################
# parameters
##########################################

outformat = ".xtc" # .xtc or .trr
base_name = "../../5_prod1/prod1*.trr"
minimum_file_size = 130000 # 134kb
trr_prefix = "prod1"

reference = "../../5_prod1/npt.gro"
ref_tpr = "../../4_eq2/npt.tpr"

# 0 => all system
# 1 => all protein atoms
# 2 => non-hydrogen atoms
# 3 => CA atoms
# 4 => backbone atoms etc...
selection_cen = 2
selection_out = 2

##########################################

### main ###
extension = File.extname(base_name).sub(/^\./, "")

Dir.glob(base_name) do |f|
    f =~ /#{trr_prefix}_(.+)\.#{extension}/
    num = $1.to_i
    outname = "./#{trr_prefix}_%05d#{outformat}"%num

    if File.exists? outname
        # skip convert if file size is more than minimum_file_size
        if File.size?(f).to_i > minimum_file_size
            puts "#{outname} already exists. skip."
            next
        end
    end

    tpr_file = f.sub(/\.#{extension}/, ".tpr")

    #`gmx trjconv -f #{f} -s #{ref} -o #{outname} -pbc nojump <<EOS \n #{selection} \n EOS` # old convert format
    `gmx trjconv -f #{f} -s #{tpr_file} -o #{outname} -pbc mol -center <<EOS \n #{selection_cen} \n #{selection_out} \n EOS`
end

# make ref.gro
`gmx trjconv -f #{reference} -s #{ref_tpr} -o ref.gro -pbc mol -center <<EOS \n #{selection_cen} \n #{selection_out} \n EOS`
