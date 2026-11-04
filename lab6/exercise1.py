"""make an artificial DNA seq of 200 to 400 in length in while to simulate 3 to 4 transposable elements.
Implement a software app which can detect the positions of these transposable elements (start,end) within the created DNA sequence.
"""


import random


def random_dna(length):
    return "".join(random.choices("ATGC", k=length))


def reverse_complement(sequence):
    complement = {"A": "T", "T": "A", "G": "C", "C": "G"}
    return "".join(complement[base] for base in reversed(sequence))


def create_dna(length, positions, transposon_types):
    dna = list(random_dna(length))
    occupied = set()

    for name, tir, tsd in transposon_types:
        start, end = positions[name]

        left = tsd + tir
        right = reverse_complement(tir) + tsd

        if end - start + 1 < len(left) + len(right):
            raise ValueError("Transposon is too short")

        for index, marker in [
            (start - 1, left),
            (end - len(right), right)
        ]:
            cells = set(range(index, index + len(marker)))

            if occupied & cells:
                raise ValueError("Transposon markers overlap")

            occupied.update(cells)
            dna[index:index + len(marker)] = marker

    return "".join(dna)


def detect_transposons(dna, transposon_types):
    results = []

    for name, tir, tsd in transposon_types:
        right_tir = reverse_complement(tir)

        for i in range(len(dna) - len(tir) + 1):
            if dna[i:i + len(tir)] != tir:
                continue

            start = i - len(tsd)

            if start < 0 or dna[start:i] != tsd:
                continue

            for j in range(i + len(tir), len(dna) - len(tir) + 1):
                if dna[j:j + len(tir)] != right_tir:
                    continue

                end = j + len(tir) + len(tsd)

                if dna[j + len(tir):end] == tsd:
                    results.append((name, start + 1, end))

    return sorted(set(results), key=lambda x: x[1])


def check_overlaps(results):
    found = False

    for i in range(len(results)):
        name1, start1, end1 = results[i]

        for j in range(i + 1, len(results)):
            name2, start2, end2 = results[j]

            if start1 <= end2 and start2 <= end1:
                found = True

                if start1 <= start2 and end1 >= end2:
                    print(f"{name2} is nested inside {name1}")

                elif start2 <= start1 and end2 >= end1:
                    print(f"{name1} is nested inside {name2}")

                else:
                    print(f"{name1} and {name2} partially overlap")

    if not found:
        print("No overlapping transposons")


random.seed(42)

transposon_types = [
    ("T1", "ATGCCGA", "TAC"),
    ("T2", "CGTATGC", "GGA"),
    ("T3", "GCTACGA", "ACT"),
    ("T4", "TGCAGAT", "CCA")
]


cases = {
    "1. SEPARATED": {
        "T1": (15, 55),
        "T2": (80, 120),
        "T3": (150, 190),
        "T4": (225, 265)
    },

    "2. OVERLAPPING": {
        "T1": (30, 110),
        "T2": (85, 165),
        "T3": (185, 225),
        "T4": (250, 290)
    },

    "3. NESTED": {
        "T1": (30, 185),
        "T2": (65, 155),
        "T3": (95, 130),
        "T4": (225, 265)
    },

    "4. THREE INTERSECTING": {
        "T1": (25, 120),
        "T2": (60, 155),
        "T3": (95, 190),
        "T4": (230, 270)
    }
}


for case_name, positions in cases.items():

    print("\n" + "=" * 60)
    print("CASE:", case_name)
    print("=" * 60)

    dna = create_dna(301, positions, transposon_types)

    print("\nARTIFICIAL DNA SEQUENCE")
    print(dna)
    print("\nDNA length:", len(dna))

    results = detect_transposons(dna, transposon_types)

    print("\nDETECTED TRANSPOSONS")

    for name, start, end in results:
        print(
            f"{name}: Start = {start}, "
            f"End = {end}, Length = {end - start + 1}"
        )

    print("\nOVERLAPPING TRANSPOSONS")
    check_overlaps(results)

    print("\nTRANSPOSON POSITIONS")

    for name, start, end in results:
        print(f"{name}: [{start} - {end}]")
        print(dna[start - 1:end])
        print()

    expected = {
        name: (start, end)
        for name, (start, end) in positions.items()
    }

    detected = {
        name: (start, end)
        for name, start, end in results
    }

    if detected == expected:
        print("RESULT: All transposons detected correctly")
    else:
        print("RESULT: Detection does not match expected positions")
