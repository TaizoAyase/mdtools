#coding: utf-8

require 'matrix'
require 'bio'

#max_pricomp = 100
natom = 233

### edit Matrix class & functions ###

class Matrix
    def []=(i,j,x)
        @rows[i][j]=x
    end
end

def format_vector_xml(v1, v2)
    %Q|\t\t<line pos1="(%3.5f, %3.5f, %3.5f)" pos2="(%3.5f, %3.5f, %3.5f)"/>\n| % [v1.to_a, v2.to_a].flatten
end

#####################################

### set input files ### 
file = "./eigcomp.xvg"
ref_pdb = "./ref.pdb"
qsc_part = "./show_vec.qsc.part"
eigen_val_file = "./eigenval.xvg"

### parse the eigcomp.xvg ###
ary = File.read(file).split("\n")
ary.reject! { |e| e =~ /^@/ }

# variables
pricomp = 0
n = 0 # section counter: 0=total, 1=x, 2=y, 3=z
m_emp = Matrix.empty(natom, 0) * Matrix.empty(0, 3)
m_tmp = m_emp.clone
vectors = {} # to set the matrices

while line = ary.shift
    # define each section
    if line =~ /&/
        if n == 3
            # if the next is total, 
            # store matrix to hash table
            vectors[pricomp] = m_tmp
            # reset n and increment pc No.
            n = 0; pricomp += 1
            # reset matrix
            m_tmp = m_emp.clone
            next
        end
        n += 1; next
    end

    # skip if I'm in the 'total' section
    next if n == 0

    # other, add to dataframe 
    a = line.split
    atom_idx = a.first.to_i
    coord = a.last.to_f
    m_tmp[atom_idx - 1, n - 1] = coord
end

str_dump = Marshal.dump(vectors)
File.write("evecs_dump.tmp", str_dump)
#ap vectors


### get CA coordinates from reference PDB ###
pdb = Bio::PDB.new(File.read(ref_pdb))
# select CA atoms
ref_ca = pdb.atoms.select { |atom| atom.name == "CA" }


### get eigen values from file ###
eigval_f = File.read(eigen_val_file).split("\n").reject! { |e| e =~ /^@/ or e =~ /^#/ }
eigen_vals = eigval_f.map { |e| e.split }
eigen_vals.map! { |e| [e[0].to_i, e[1].to_f] }
eval_tot = eigen_vals.inject(0) { |sum, e| sum += e[1] }
sum = 0
eigen_vals.each do |e|
    contribute = e[1] / eval_tot * 100
    puts "pc No.#{e[0]}: #{"%3.3f" % contribute}% contribute"
    sum += contribute
    break if sum > 90
end
puts "|--------90% total--------|"


##### OUTPUT RENDER SECTION OF QSC #####
vect_str = File.read(qsc_part)
vect_str << %Q|<object type="MolCoord" name="pdb" sel="*" srctype="pdb" src="ref.pdb">\n|
vect_str << %Q|\t<renderer type="*selection" visible="false"/>\n|
vect_str << %Q|\t<renderer type="*namelabel" style="DefaultAtomLabel"/>\n|

# puts out for each vector point
vectors.each_key do |pc_no|
    vect_str << %Q|\t<renderer type="atomintr" color="#BF0000" mode="fancy" showlabel="false" stipple0="1000.0" name="dom%d" visible="true" width="0.3">\n| % pc_no
    eigen_val = eigen_vals[pc_no][1]
    vect_ary = vectors[pc_no].row_vectors
    ref_ca.zip vect_ary do |coord_ary|
        vec_minus = coord_ary[0].xyz - coord_ary[1] * eigen_val * 5 
        vec_plus    = coord_ary[0].xyz + coord_ary[1] * eigen_val * 5
        vect_str << format_vector_xml(vec_minus, vec_plus) 
    end
    vect_str << %Q|\t</renderer>\n|
end

# tail
vect_str << %Q|\t</object>\n|
vect_str << %Q|<camera center="(62.745417,61.944431,60.189757)" centerMark="crosshair" distance="200" name="__current" perspec="false" rotation="(-0.402869,-0.653961,0.447021,0.458479)" slab="130.747442" stereoDist="1" stereoMode="none" zoom="139.012709"/>\n|
vect_str << %Q|</scene>|

# write out to qsc file
File.write("show_vec.qsc", vect_str)
puts "QSC output done."
