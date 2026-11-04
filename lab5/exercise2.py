
FASTA_FILE = "pUC19c_L09137.2.fasta"


def read_fasta(filename):
    with open(filename, encoding="utf-8") as file:
        return "".join(
            line.strip().upper()
            for line in file
            if not line.startswith(">")
        )


def is_primitive(unit):
    for size in range(1, len(unit)):
        if len(unit) % size == 0:
            if unit[:size] * (len(unit) // size) == unit:
                return False

    return True


def find_repeats(dna):
    results = []

    for unit_size in range(3, 6):
        i = 0

        while i + 2 * unit_size <= len(dna):
            unit = dna[i:i + unit_size]

            if not is_primitive(unit):
                i += 1
                continue

            copies = 1

            while (
                i + (copies + 1) * unit_size <= len(dna)
                and dna[
                    i + copies * unit_size:
                    i + (copies + 1) * unit_size
                ] == unit
            ):
                copies += 1

            if copies >= 2:
                allowed_copies = min(
                    copies, 10 // unit_size
                )

                if allowed_copies >= 2:
                    repeat = unit * allowed_copies

                    results.append((
                        i + 1,
                        i + len(repeat),
                        unit,
                        allowed_copies,
                        repeat
                    ))

                    i += copies * unit_size
                    continue

            i += 1

    return sorted(results)


def main():
    dna = read_fasta(FASTA_FILE)
    repeats = find_repeats(dna)

    print(f"DNA length: {len(dna)} bp")
    print(f"Total tandem repeats found: {len(repeats)}")

    for start, end, unit, copies, repeat in repeats:
        print(
            f"{start}-{end}: {repeat} "
            f"(unit={unit}, "
            f"copies={copies}, "
            f"length={len(repeat)})"
        )


if __name__ == "__main__":
    main()
