import os
import sys

big_sam_file = open(sys.argv[1], "r")
bc_file = open(sys.argv[2], "r")
protocol = sys.argv[3]
sample = sys.argv[4]
output_path = sys.argv[5]


bc_list = []
dict_demux = {}
for line in bc_file:
    bc = line.split()[0]
    if bc != 'XC':
        bc_list.append(bc)
        dict_demux[bc] = []

header_list = []
for line in big_sam_file:
    if line.startswith("@"):
        header_list.append(line)
    else:
        line_words = line.split()
        BC_tag = [i for i in line_words if i.startswith('BC')][0]
        if len(BC_tag) != 3:
            print(BC_tag)
        cellBC = BC_tag.split(":")[2]
        if cellBC in dict_demux.keys():
            dict_demux[cellBC].append(line)
            
for cBC in bc_list:
    filenames = output_path + "/" + protocol + "." + "mixed." + sample + "." + cBC + ".sam"
    fout = open(filenames, 'w')
    for hline in header_list:
        fout.write(hline)
    for mapped_line in dict_demux[cBC]:
        fout.write(mapped_line)
    fout.close()
        
        