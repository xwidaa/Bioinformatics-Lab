import argparse
import csv
from pathlib import Path

import matplotlib
matplotlib.use('Agg')  
import matplotlib.pyplot as plt

WINDOW_SIZE = 30
TEST_SEQUENCE = (
    'CGGACTGATCTATCT' + 'A' * 27 + 'CGTAGCATCTATCGATCTATCTAGCGATCTA'
)


def read_fasta(path):
    records = []
    name, fragments = None, []
    with open(path, encoding='utf-8') as source:
        for raw_line in source:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if name is not None:
                    records.append((name, ''.join(fragments).upper()))
                name = line[1:].strip() or f'promoter_{len(records) + 1}'
                fragments = []
            else:
                if name is None:
                    name = Path(path).stem
                fragments.append(line)
    if name is not None:
        records.append((name, ''.join(fragments).upper()))
    return records


def validate(sequence):
    sequence = ''.join(sequence.upper().split())
    if not sequence or set(sequence) - set('ACGT'):
        raise ValueError('DNA must contain only A, C, G, and T.')
    return sequence


def gc_percent(sequence):
    sequence = validate(sequence)
    return 100 * (sequence.count('C') + sequence.count('G')) / len(sequence)


def kappa_ic(window):
    """Mean percentage of matching bases across offsets 1..n-1."""
    window = validate(window)
    n = len(window)
    if n < 2:
        raise ValueError('Kappa IC needs at least two bases.')
    values = []
    for shift in range(1, n):
        matches = sum(window[i] == window[i + shift] for i in range(n - shift))
        values.append(100.0 * matches / (n - shift))
    return sum(values) / len(values)


def analyze(name, sequence, window_size=WINDOW_SIZE, relative_gc=True):
    sequence = validate(sequence)
    if len(sequence) < window_size:
        raise ValueError(f'{name}: fewer than {window_size} bases')

    whole_gc = gc_percent(sequence)
    windows = [sequence[i:i + window_size]
               for i in range(len(sequence) - window_size + 1)]
    gc_values = [gc_percent(window) for window in windows]
    if relative_gc:
        gc_values = [(whole_gc / 100) * value for value in gc_values]
    ic_values = [kappa_ic(window) for window in windows]

    center = (sum(gc_values) / len(gc_values), sum(ic_values) / len(ic_values))
    return {
        'name': name,
        'length': len(sequence),
        'whole_gc': whole_gc,
        'windows': windows,
        'gc': gc_values,
        'ic': ic_values,
        'center': center,
        'relative_gc': relative_gc,
    }


def plot_pattern(result, output):
    fig, ax = plt.subplots(figsize=(9, 6))
    ax.scatter(
    result['gc'],
    result['ic'],
    s=3,
    marker='s',
    color='blue',
    alpha=0.85
)
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 100)

    ax.axhline(50, color='black', linestyle='--', linewidth=0.6)
    ax.axvline(50, color='black', linestyle='--', linewidth=0.6)
    x, y = result['center']
    ax.scatter([x], [y], marker='X', s=140, color='crimson', zorder=5,
               label=f'Center ({x:.2f}, {y:.2f})')
    xlabel = ('Relative C+G (whole-sequence GC % × window GC fraction)'
              if result['relative_gc'] else 'Window C+G (%)')
    ax.set(xlabel=xlabel, ylabel='Kappa IC (%)',
           title=f"Promoter pattern: {result['name'][:75]}")
    ax.grid(alpha=.25)
    ax.legend()
    fig.tight_layout()
    fig.savefig(output, dpi=120)
    plt.close(fig)


def plot_distribution(results, output):
    fig, ax = plt.subplots(figsize=(10, 7))
    promoter_results = [r for r in results if r['name'] != 'Artificial test sequence']
    if promoter_results:
        ax.scatter([r['center'][0] for r in promoter_results],
                   [r['center'][1] for r in promoter_results],
                   s=28, alpha=.55, color='steelblue', label='FASTA promoters')
    artificial_results = [r for r in results if r['name'] == 'Artificial test sequence']
    for result in artificial_results:
        x, y = result['center']
        ax.scatter(x, y, s=85, marker='X', color='crimson', label='Artificial test')
        ax.annotate('Artificial test', (x, y), xytext=(5, 5),
                    textcoords='offset points', fontsize=8)
    for result in promoter_results[:9]:
        x, y = result['center']
        ax.annotate(result['name'].split(';')[0][:23], (x, y), xytext=(4, 4),
                    textcoords='offset points', fontsize=7)
    xlabel = ('Center of relative C+G' if results[0]['relative_gc']
              else 'Center of window C+G (%)')
    ax.set(xlabel=xlabel, ylabel='Center of Kappa IC (%)',
           title='General distribution of promoter pattern centers')
    ax.grid(alpha=.25)
    ax.legend()
    fig.tight_layout()
    fig.savefig(output, dpi=120)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description='DNA promoter pattern analysis')
    parser.add_argument('--fasta', help='FASTA file containing promoter sequences')
    parser.add_argument('--window', type=int, default=WINDOW_SIZE)
    parser.add_argument('--max-charts', type=int, default=10,
                        help='Maximum individual pattern PNGs (default: 10)')
    parser.add_argument('--absolute-gc', action='store_true',
                        help='Use ordinary GC percentage on the X-axis')
    parser.add_argument('--output', default='dna_pattern_results')
    args = parser.parse_args()

    if args.window < 2:
        parser.error('--window must be >= 2')
    if args.max_charts < 0:
        parser.error('--max-charts must be >= 0')

    output = Path(args.output)
    output.mkdir(parents=True, exist_ok=True)
    sequences = [('Artificial test sequence', TEST_SEQUENCE)]
    if args.fasta:
        sequences.extend(read_fasta(args.fasta)[:100])

    results = []
    for index, (name, sequence) in enumerate(sequences, start=1):
        try:
            result = analyze(name, sequence, args.window, not args.absolute_gc)
        except ValueError as error:
            print(f'Skipping {name}: {error}')
            continue
        results.append(result)
        if len(results) <= args.max_charts:
            path = output / f'pattern_{len(results)}.png'
            plot_pattern(result, path)
            print(f'Saved {path.name} — {name[:55]}')

        if len(results) == 1:
            print(f"Artificial test: {result['length']} bases, "
                  f"{len(result['windows'])} windows")
            print(f"Whole-sequence GC: {result['whole_gc']:.2f}%")
            print(f"First-window ordinary GC: {gc_percent(result['windows'][0]):.2f}%")
            print(f"First-window plotted C+G: {result['gc'][0]:.2f}")
            print(f"First-window Kappa IC: {result['ic'][0]:.2f}%")
        if index % 100 == 0:
            print(f'Analyzed {index} sequences...')

    if not results:
        raise SystemExit('No valid DNA sequences to analyze.')
    plot_distribution(results, output / 'general_distribution.png')

    with (output / 'pattern_centers.csv').open('w', newline='', encoding='utf-8') as file:
        writer = csv.writer(file)
        writer.writerow(['name', 'bases', 'windows', 'whole_gc_percent',
                         'center_cg', 'center_ic_percent'])
        for r in results:
            writer.writerow([r['name'], r['length'], len(r['windows']),
                             f"{r['whole_gc']:.4f}", f"{r['center'][0]:.4f}",
                             f"{r['center'][1]:.4f}"])

    print(f'\nSequences analyzed: {len(results)}')
    print(f'Individual pattern charts saved: {min(len(results), args.max_charts)}')
    print('Distribution chart: general_distribution.png')
    print('All centers: pattern_centers.csv')
    print(f'Output folder: {output.resolve()}')
    print('Note: benchmark CG≈29.27 and IC≈27.53 are not yet reproduced.')


if __name__ == '__main__':
    main()
