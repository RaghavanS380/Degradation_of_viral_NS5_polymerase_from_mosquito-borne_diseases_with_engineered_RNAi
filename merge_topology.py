# save as merge_topology.py (overwrite the old one)
import re

with open('#topol.top.3#') as f:
    protein_top = f.read()

with open('#topol.top.4#') as f:
    rna_top = f.read()

# Extract RNA moleculetype block
rna_section = re.search(r'(\[ moleculetype \].*?)(?=\[ system \])', rna_top, re.DOTALL)
rna_block = rna_section.group(1) if rna_section else ""

# Remove #include lines — those are already in the protein topology
rna_lines = [l for l in rna_block.split('\n') if not l.strip().startswith('#include')]
rna_block = '\n'.join(rna_lines)

# Take protein topology up to [ system ]
protein_base = re.search(r'^(.*?)(?=\[ system \])', protein_top, re.DOTALL)
prot_block = protein_base.group(1) if protein_base else ""

combined = prot_block
combined += rna_block
combined += "\n[ system ]\nDengue NS5 + antisense RNA in water\n\n"
combined += "[ molecules ]\n"
combined += "; Compound        #mols\n"
combined += "Protein_chain_A     1\n"
combined += "Protein_chain_B     1\n"
combined += "RNA_chain_A         1\n"
combined += "SOL            108487\n"

with open('topol_combined.top', 'w') as f:
    f.write(combined)

print("Done — topol_combined.top written")