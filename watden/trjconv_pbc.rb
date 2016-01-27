#coding: utf-8

#ref = "../../7_prod1/npt.gro"
outformat = ".xtc" # .xtc or .trr
#minimum_file_size = 130000 # 134kb
minimum_file_size = 0

# 1 => all protein atoms
# 2 => non-hydrogen atoms
# 3 => CA atoms
# 4 => backbone atoms
selection_cen = 2
selection_out = 24

Dir.glob("../../5_prod1/prod1*.trr") do |f|
  f =~ /prod1_(.+)\.trr/
  num = $1.to_i
  outname = "./prod1_%05d#{outformat}"%num
  if File.exists? outname
    # skip convert only if file size is more than minimum_file_size
    if File.size?(f).to_i > minimum_file_size
      puts "#{outname} already exists. skip."
      next
    end
  end
  tpr_file = f.sub(/\.trr/, ".tpr")
  puts tpr_file

  `gmx trjconv -f #{f} -s #{tpr_file} -n index.ndx -o #{outname} -pbc mol -center <<EOS \n #{selection_cen} \n #{selection_out} \n EOS`
end
