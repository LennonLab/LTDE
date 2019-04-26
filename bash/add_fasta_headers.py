import sys

in_file_name = sys.argv[1]
out_file_name = sys.argv[2]
out_file=open(out_file_name, 'w')

nucs = set(["A","C","G","T"])

for i, line in enumerate(open(in_file_name, "r")):
    line=line.strip()
    line_set=set(line)
    if len(line_set.union(nucs)) > 4:
        continue
    out_file.write(">Barcode_"+str(i) + '\n')
    out_file.write(line + '\n')

out_file.close()
