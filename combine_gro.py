# save as combine_gro.py in ~/Desktop/dengue_sim/
def read_gro(filename):
    with open(filename) as f:
        lines = f.readlines()
    title = lines[0]
    natoms = int(lines[1].strip())
    atoms = lines[2:2+natoms]
    box = lines[2+natoms]
    return title, natoms, atoms, box

prot_title, prot_n, prot_atoms, prot_box = read_gro('protein_processed.gro')
rna_title,  rna_n,  rna_atoms,  rna_box  = read_gro('rna_processed.gro')

total = prot_n + rna_n

with open('complex.gro', 'w') as f:
    f.write('Dengue NS5 + antisense RNA complex\n')
    f.write(f'{total:5d}\n')
    for line in prot_atoms:
        f.write(line)
    for line in rna_atoms:
        f.write(line)
    f.write(prot_box)  # use protein box line

print(f"Combined: {prot_n} protein atoms + {rna_n} RNA atoms = {total} total")