from collections import Counter
import matplotlib.pyplot as plt
from pathlib import Path


filename = Path(__file__).resolve().parent / "GCF_000909615.1_ViralProj212954_genomic.fna"

with open(filename, "r") as file:
    next(file)
    genome = "".join(line.strip() for line in file).upper()


codons = []

for i in range(len(genome) - 2):
    codon = genome[i:i+3]

    if all(base in "ACGT" for base in codon):
        codons.append(codon)


frequencies = Counter(codons)


top_10 = frequencies.most_common(10)

print("Top 10 most frequent codons:")

for codon, count in top_10:
    print(codon, ":", count)


names = [codon for codon, count in top_10]
counts = [count for codon, count in top_10]

plt.bar(names, counts)
plt.xlabel("Codons")
plt.ylabel("Frequency")
plt.title("Top 10 Most Frequent Codons")
plt.show()