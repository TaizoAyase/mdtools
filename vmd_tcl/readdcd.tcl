package require pbctools
mol new ../psfgen/ionized.psf type psf

for {set i 0} {$i <= 99} {incr i} {
    set dcd_file [format p1_%02d.dcd $i]
    mol addfile $dcd_file step 20 waitfor all
}

#pbc wrap -centersel "protein" -center com -compound residue -all 
#animate write dcd {test.dcd} beg 0 end -1  

#quit
