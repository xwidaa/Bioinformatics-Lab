
import math
import random
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

FASTA_FILE = "pUC19c_L09137.2.fasta"


def read_fasta(filename):
    with open(filename, encoding="utf-8") as file:
        return "".join(
            line.strip().upper()
            for line in file
            if not line.startswith(">")
        )


def take_samples(dna):
    samples = []

    for _ in range(10):
        length = random.randint(100, min(3000, len(dna)))
        start = random.randint(0, len(dna) - length)

        samples.append(dna[start:start + length])

    return samples


def migration_distance(length):
    return 1 + 8 * (
        (math.log10(3000) - math.log10(length))
        / (math.log10(3000) - math.log10(100))
    )


def main():
    dna = read_fasta(FASTA_FILE)

    samples = take_samples(dna)
    lengths = [len(sample) for sample in samples]
    distances = [migration_distance(length) for length in lengths]

    print(f"DNA length: {len(dna)} bp")

    for i, (length, distance) in enumerate(
        zip(lengths, distances), 1
    ):
        print(
            f"Sample {i}: {length} bp, "
            f"migration distance = {distance:.2f} cm"
        )

    fig, ax = plt.subplots(figsize=(12, 8))

    ax.set_facecolor("black")
    ax.set_xlim(-0.6, 10.6)
    ax.set_ylim(10, -0.5)

    ax.set_xticks(
        range(11),
        ["L"] + [str(i) for i in range(1, 11)]
    )

    ax.set_ylabel("Distance (cm)")
    ax.set_title("Gel electrophoresis")

    ladder = [3000, 2500, 2000, 1500, 1000,
              750, 500, 250, 100]

    for size in ladder:
        y = migration_distance(size)

        ax.hlines(
            y, -0.35, 0.35,
            color="deepskyblue",
            linewidth=5
        )

        ax.text(
            -0.42, y, str(size),
            color="white",
            ha="right",
            va="center",
            fontsize=9
        )

    bands = [
        ax.hlines(
            0, i + 0.65, i + 1.35,
            color="orange",
            linewidth=5
        )
        for i in range(10)
    ]

    labels = [
        ax.text(
            i + 1, 0.25, str(lengths[i]),
            color="white",
            ha="center",
            fontsize=9
        )
        for i in range(10)
    ]

    def update(frame):
        progress = frame / 50

        for i, band in enumerate(bands):
            y = distances[i] * progress

            band.set_segments([
                [(i + 0.65, y), (i + 1.35, y)]
            ])

            labels[i].set_position((i + 1, y + 0.25))

        return bands + labels

    animation = FuncAnimation(
        fig,
        update,
        frames=51,
        interval=50,
        repeat=False
    )

    plt.show()

    update(50)

    fig.savefig(
        "gel_electrophoresis.png",
        dpi=150,
        bbox_inches="tight"
    )


if __name__ == "__main__":
    main()
