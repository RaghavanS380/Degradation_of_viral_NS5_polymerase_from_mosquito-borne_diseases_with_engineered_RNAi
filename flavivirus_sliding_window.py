import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
from Bio import SeqIO

# Load alignment
seqs = {rec.id: str(rec.seq) for rec in
        SeqIO.parse('/Users/smyan/Desktop/flavivirus_aligned.fasta', 'fasta')}

print("Loaded:", list(seqs.keys()))

# Pairs to compare
pairs = [
    ('NC_001477.1', 'NC_012532.1',  'Dengue-1 vs Zika',         '#e63946'),
    ('NC_001477.1', 'NC_002031.1',  'Dengue-1 vs Yellow Fever',  '#f77f00'),
    ('NC_012532.1', 'NC_002031.1',  'Zika vs Yellow Fever',      '#4361ee'),
]

window = 300   # smaller window since genomes are ~10kb not 30kb
step   = 30

def sliding_identity(seq1, seq2, window, step):
    positions, identities = [], []
    for start in range(0, len(seq1) - window, step):
        end = start + window
        s1, s2 = seq1[start:end], seq2[start:end]
        matches = sum(a == b and a != '-' and b != '-' for a, b in zip(s1, s2))
        valid   = sum(a != '-' and b != '-' for a, b in zip(s1, s2))
        if valid > 0:
            positions.append(start + window // 2)
            identities.append(matches / valid * 100)
    return np.array(positions), np.array(identities)

# Flavivirus genome annotation (approximate, based on ~10.8kb genome)
genes = [
    ("5'UTR", 1,    96,    '#dde8f0'),
    ('C',     97,   474,   '#ffd6d6'),
    ('prM',   475,  993,   '#ffe4b5'),
    ('E',     994,  2469,  '#d4edda'),   # Envelope — key for immunity
    ('NS1',   2470, 3525,  '#d6e4ff'),
    ('NS2A',  3526, 4218,  '#f0d6ff'),
    ('NS2B',  4219, 4611,  '#ffd6f0'),
    ('NS3',   4612, 6468,  '#d6fff0'),   # Helicase/protease
    ('NS4A',  6469, 6936,  '#fff0d6'),
    ('NS4B',  6937, 7695,  '#e8d6ff'),
    ('NS5',   7696, 10271, '#d6f5d6'),   # RNA polymerase
    ("3'UTR", 10272,10862, '#dde8f0'),
]

genome_length = 10862

fig, axes = plt.subplots(2, 1, figsize=(14, 8),
    gridspec_kw={'height_ratios': [0.4, 3.6]}, facecolor='white')

# Genome map
ax_map = axes[0]
ax_map.set_xlim(0, genome_length)
ax_map.set_ylim(0, 1)
ax_map.axis('off')
ax_map.axhline(0.5, color='#999', linewidth=1.5)
for name, gs, ge, color in genes:
    rect = mpatches.FancyArrow(gs, 0.5, ge - gs, 0,
        width=0.42, head_width=0.42, head_length=80,
        length_includes_head=True, color=color,
        linewidth=0.5, edgecolor='#888')
    ax_map.add_patch(rect)
    if ge - gs > 200:
        ax_map.text((gs + ge) / 2, 0.5, name,
            ha='center', va='center', fontsize=6.5, fontweight='bold', color='#222')

ax_map.set_title('Flavivirus Sliding Window Pairwise Identity (window=300 nt, step=30 nt)',
    fontsize=12, fontweight='bold', pad=8)

# Identity curves
ax = axes[1]
ax.set_facecolor('#fafafa')

for id1, id2, label, color in pairs:
    if id1 not in seqs or id2 not in seqs:
        print(f"Skipping {label} — ID not found. Available: {list(seqs.keys())}")
        continue
    pos, ident = sliding_identity(seqs[id1], seqs[id2], window, step)
    ax.plot(pos, ident, color=color, linewidth=1.5, label=label, alpha=0.85)

# Shade envelope gene — most variable, most interesting
ax.axvspan(994, 2469, alpha=0.08, color='green', label='Envelope (E) gene')
ax.axvspan(7696, 10271, alpha=0.06, color='purple', label='NS5 (RNA polymerase)')

ax.axhline(50, color='#aaa', linestyle='--', linewidth=0.8, alpha=0.6)
ax.text(100, 50.8, '50%', fontsize=7, color='#888')

ax.set_xlim(0, genome_length)
ax.set_ylim(20, 101)
ax.set_xlabel('Genome position (nt)', fontsize=11)
ax.set_ylabel('Pairwise nucleotide identity (%)', fontsize=11)
ax.legend(loc='lower right', fontsize=9, framealpha=0.95)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

xticks = np.arange(0, genome_length + 1, 2000)
ax.set_xticks(xticks)
ax.set_xticklabels([f'{int(x/1000)}k' if x > 0 else '0' for x in xticks])
axes[0].set_xticks(xticks)

plt.tight_layout(h_pad=0.3)
out = '/Users/smyan/Desktop/flavivirus_sliding_window.png'
plt.savefig(out, dpi=180, bbox_inches='tight', facecolor='white')
print(f"Saved: {out}")