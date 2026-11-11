from collections import Counter
from pathlib import Path
import math
import matplotlib.pyplot as plt
from Bio import SeqIO

BASE = Path(__file__).resolve().parent
FILES = [
    ("DNA 1 (pUC19c)", "dna_sequence.fasta", True),
    ("DNA 2 (Porcine circovirus 2)", "dna_sequence2.fasta", True),
]
COMP = str.maketrans("ACGT", "TGCA")


def complement(seq):
    return seq.translate(COMP)


def reverse_complement(seq):
    return complement(seq)[::-1]


def read_dna(filename):
    path = BASE / filename
    if not path.exists():
        raise FileNotFoundError(f"Missing {filename}. Put it beside exercise1.py")
    record = SeqIO.read(path, "fasta")
    seq = str(record.seq).upper()
    if not seq or set(seq) - set("ACGT"):
        raise ValueError(f"{filename} contains missing or ambiguous DNA bases")
    if not 1000 <= len(seq) <= 3000:
        raise ValueError(f"{filename} must contain 1000-3000 bases; got {len(seq)}")
    return record.id, seq


def positions_of_cuts(dna, motif, cut, circular):
    L, m = len(dna), len(motif)
    rc = reverse_complement(motif)
    extended = dna + dna[:m-1] if circular else dna
    count = L if circular else max(0, L-m+1)
    sites = set()
    for i in range(count):
        found = extended[i:i+m]
        if found == motif:
            sites.add((i+cut) % L if circular else i+cut)
        if rc != motif and found == rc:
            sites.add((i+m-cut) % L if circular else i+m-cut)
    return sorted(sites)


def fragment_lengths(L, cuts, circular):
    if not circular:
        bounds = [0] + [p for p in cuts if 0 < p < L] + [L]
        return [b-a for a,b in zip(bounds, bounds[1:])]
    if not cuts:
        return [L]  
    return [b-a for a,b in zip(cuts, cuts[1:])] + [L-cuts[-1]+cuts[0]]


def cleavage_diagram(motif, cut):
    bottom = complement(motif)
    other_cut = len(motif)-cut
    print("  5' " + motif[:cut] + "|" + motif[cut:] + " 3'")
    print("  3' " + bottom[:other_cut] + "|" + bottom[other_cut:] + " 5'")
    print("  End type:", "blunt" if cut == other_cut else "sticky")


def gel(results):
    ladder = [3000, 2000, 1500, 1000, 750, 500, 250, 100]
    lanes = [("Ladder", ladder, False)] + results
    fig, ax = plt.subplots(figsize=(max(11, len(lanes)*1.4), 7.5))
    fig.patch.set_facecolor('white')
    ax.set_facecolor('#151525')

    max_bp = 3000
    def migration(bp):
        bp = max(50, min(bp, max_bp))
        return 0.12 + 0.78 * math.log10(max_bp/bp) / math.log10(max_bp/50)

    for ix, (label, fragments, uncut_circular) in enumerate(lanes, 1):
        ax.add_patch(plt.Rectangle((ix-.32, .06), .64, .89,
                                   color='#27273e', alpha=.85))
        ax.add_patch(plt.Rectangle((ix-.23, .065), .46, .025,
                                   color='#080815'))
        if uncut_circular:
            ax.text(ix, .50, "uncut\ncircular", ha='center',
                    va='center', fontsize=8, color='white')
            continue
        for bp, copies in Counter(fragments).items():
            y = migration(bp)
            ax.hlines(y, ix-.25, ix+.25, color='#75ffe3',
                      linewidth=min(3 + copies*1.5, 9),
                      alpha=min(.50 + copies*.15, 1))
            if ix == 1:
                ax.text(ix-.36, y, str(bp), ha='right', va='center',
                        fontsize=8, color='white')

    ax.set_xlim(.15, len(lanes)+.6)
    ax.set_ylim(1, 0)
    ax.set_xticks(range(1, len(lanes)+1))
    ax.set_xticklabels([name for name, _, _ in lanes], rotation=25,
                       ha='left', fontsize=9)
    ax.xaxis.tick_top()
    ax.tick_params(axis='x', length=0, pad=7)
    ax.set_yticks([])
    ax.set_title('Simulated gel electrophoresis — separate enzyme digestions',
                 fontsize=13, pad=55)
    plt.tight_layout()
    fig.savefig(BASE / 'electrophoresis_gel.png', dpi=180, bbox_inches='tight')
    plt.show()


def main():
    genomes = []
    print('NCBI DNA INPUTS')
    for label, file, circular in FILES:
        accession, dna = read_dna(file)
        genomes.append((label, dna, circular))
        print(f'{label}: {accession}, {len(dna)} bp, '
              f'{"circular" if circular else "linear"}')

    n = int(input('\nNumber of restriction enzymes (3 for this lab): '))
    if n < 1:
        raise ValueError('Enter at least one restriction enzyme')
    enzymes = []
    for i in range(n):
        print(f'\nRestriction enzyme {i+1}/{n}')
        name = input('Name: ').strip() or f'Enzyme {i+1}'
        motif = input('Recognition sequence (A/C/G/T): ').strip().upper()
        if not motif or set(motif)-set('ACGT'):
            raise ValueError('Recognition sequence must use A/C/G/T only')
        cut = int(input('Cut position (e.g., GAATTC 1): '))
        if not 0 <= cut <= len(motif):
            raise ValueError('Cut position must be between 0 and motif length')
        if motif != reverse_complement(motif):
            print('WARNING: Non-palindromic motif. The complementary-strand '
                  'cleavage offset is assumed to be mirrored. Real enzymes '
                  'may require a separately specified offset.')
        enzymes.append((name, motif, cut))

    lanes = []
    for label, dna, circular in genomes:
        print('\n'+'='*64)
        print(f'{label}: {len(dna)} bp')
        print('='*64)
        for name, motif, cut in enzymes:
            cuts = positions_of_cuts(dna, motif, cut, circular)
            effective_cuts = cuts if circular else [p for p in cuts if 0 < p < len(dna)]
            fragments = fragment_lengths(len(dna), effective_cuts, circular)
            print(f'\n{name} | recognition: {motif} | top-strand cut: {cut}')
            cleavage_diagram(motif, cut)
            print('Number of cleavage sites:', len(effective_cuts))
            print('Cleavage positions (0-based boundaries):', effective_cuts)
            print('Fragment lengths (bp):', fragments)
            print('Number of fragments:', len(fragments))
            print('Length check:', sum(fragments), '/', len(dna), 'bp')
            lanes.append((f'{label.split(" (")[0]}\n{name}', fragments,
                          circular and not effective_cuts))

    gel(lanes)
    print('\nGel saved to electrophoresis_gel.png')


if __name__ == '__main__':
    main()
