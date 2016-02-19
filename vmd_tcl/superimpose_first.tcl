#################################################################  
# Generate superimposed frames  
# written by 
# http://oinil.blogspot.jp/2012/07/superimpose-conformations-from-md.html
# 
# modified by @Taizo_Ayase
#################################################################  
  
set nf [molinfo top get numframes]  
   
set domain1  [atomselect top "protein and name CA"]  
set all  [atomselect top "all"]  
set ref [atomselect top "protein and name CA" frame 0]  
  
# superimpose fit loop  
for { set i 0 } { $i <= $nf } { incr i } {  
	$domain1 frame $i  
	$all frame $i  
	$all move [measure fit $domain1 $ref]  
}  
