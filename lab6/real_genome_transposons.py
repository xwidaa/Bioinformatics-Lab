
import csv
import os
from urllib.request import urlopen


GENOME_ID = "NC_001416.1"
FASTA_FILE = "viral_genome.fasta"

MIN_TIR = 6
MAX_TIR = 8
MIN_INTERNAL = 10
MAX_INTERNAL = 150
MIN_TSD = 2
MAX_TSD = 3


def download_genome():
    if os.path.exists(FASTA_FILE):
        print("Using existing genome file:", FASTA_FILE)
        return

    url = (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
        f"?db=nuccore&id={GENOME_ID}&rettype=fasta&retmode=text"
    )

    print("Downloading genome from NCBI...")

    with urlopen(url, timeout=30) as response:
        data = response.read().decode("utf-8")

    if not data.startswith(">"):
        raise ValueError("Invalid FASTA response")

    with open(FASTA_FILE, "w") as file:
        file.write(data)

    print("Genome downloaded successfully!")


def read_genome(filename):
    with open(filename, "r") as file:
        lines = file.readlines()

    dna = "".join(
        line.strip().upper()
        for line in lines
        if not line.startswith(">")
    )

    if not dna:
        raise ValueError("The FASTA file is empty")

    return dna


def reverse_complement(sequence):
    complement = {"A": "T", "T": "A", "G": "C", "C": "G"}

    return "".join(
        complement[base]
        for base in reversed(sequence)
    )


def detect_transposons(dna):
    results = {}
    n = len(dna)

    for tir_length in range(MAX_TIR, MIN_TIR - 1, -1):
        print(f"Searching TIR length {tir_length}...")

        for i in range(n - 2 * tir_length - MIN_INTERNAL + 1):
            left_tir = dna[i:i + tir_length]

            if any(base not in "ATGC" for base in left_tir):
                continue

            right_tir = reverse_complement(left_tir)

            search_start = i + tir_length + MIN_INTERNAL
            search_end = min(
                i + tir_length + MAX_INTERNAL + tir_length,
                n
            )

            j = dna.find(right_tir, search_start, search_end)

            while j != -1:
                element_end = j + tir_length
                internal_length = j - i - tir_length

                for tsd_length in range(MAX_TSD, MIN_TSD - 1, -1):
                    if i < tsd_length:
                        continue

                    if element_end + tsd_length > n:
                        continue

                    left_tsd = dna[i - tsd_length:i]
                    right_tsd = dna[element_end:element_end + tsd_length]

                    if left_tsd != right_tsd:
                        continue

                    if any(base not in "ATGC" for base in left_tsd):
                        continue

                    start = i + 1
                    end = element_end
                    key = (start, end)

                    candidate = {
                        "start": start,
                        "end": end,
                        "length": end - start + 1,
                        "tir_length": tir_length,
                        "left_tir": left_tir,
                        "right_tir": right_tir,
                        "internal_length": internal_length,
                        "tsd": left_tsd,
                        "tsd_length": tsd_length
                    }

                    if key not in results:
                        results[key] = candidate
                    else:
                        old = results[key]

                        if (
                            tir_length,
                            tsd_length
                        ) > (
                            old["tir_length"],
                            old["tsd_length"]
                        ):
                            results[key] = candidate

                    break

                j = dna.find(right_tir, j + 1, search_end)

    return sorted(
        results.values(),
        key=lambda x: (
            -x["tir_length"],
            -x["tsd_length"],
            x["start"]
        )
    )


def check_overlaps(results):
    overlaps = []

    ordered = sorted(results, key=lambda x: x["start"])

    for i in range(len(ordered)):
        a = ordered[i]

        for j in range(i + 1, len(ordered)):
            b = ordered[j]

            if b["start"] > a["end"]:
                break

            if b["end"] <= a["end"]:
                relation = "Nested"
            else:
                relation = "Partial overlap"

            overlaps.append((
                a["start"],
                a["end"],
                b["start"],
                b["end"],
                relation
            ))

    return overlaps


def save_results(results):
    with open("possible_transposons.csv", "w", newline="") as file:
        writer = csv.DictWriter(
            file,
            fieldnames=[
                "start",
                "end",
                "length",
                "tir_length",
                "left_tir",
                "right_tir",
                "internal_length",
                "tsd",
                "tsd_length"
            ]
        )

        writer.writeheader()
        writer.writerows(results)


def main():
    try:
        download_genome()

    except Exception as error:
        print("Download failed:", error)
        print("Download the FASTA manually and save it as", FASTA_FILE)
        return

    dna = read_genome(FASTA_FILE)

    print("\nGENOME INFORMATION")
    print("Genome:", GENOME_ID)
    print("DNA length:", len(dna))

    print("\nSTARTING TRANSPOSON DETECTION")

    results = detect_transposons(dna)

    print("\nTOTAL POSSIBLE TRANSPOSONS:", len(results))

    print("\nTOP 20 CANDIDATES")

    for index, candidate in enumerate(results[:20], 1):
        print(
            f"{index}. "
            f"Start: {candidate['start']}, "
            f"End: {candidate['end']}, "
            f"Length: {candidate['length']}, "
            f"TIR: {candidate['left_tir']} / "
            f"{candidate['right_tir']}, "
            f"Internal: {candidate['internal_length']}, "
            f"TSD: {candidate['tsd']}"
        )

    save_results(results)

    print("\nOVERLAPPING CANDIDATES")

    overlaps = check_overlaps(results)

    print("Total overlapping pairs:", len(overlaps))

    for overlap in overlaps[:10]:
        print(overlap)

    print("\nResults saved to possible_transposons.csv")


if __name__ == "__main__":
    main()
