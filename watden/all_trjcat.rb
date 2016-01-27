#coding: utf-8

#start_time = 0
#end_time = 250

input = "./prod1_*.xtc"
output = "all.xtc"

file_num = Dir.glob(input).length
setime_com = "c \n"*file_num

command = "gmx_mpi trjcat -f #{input} -o #{output} -overwrite -cat -settime <<EOS\n#{setime_com} EOS"
#puts command
system command
