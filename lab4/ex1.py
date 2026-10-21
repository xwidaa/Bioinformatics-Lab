#Take a DNA sequence from NCBA entire genoms and take 2000 random samples from the entire sequence of about 100-120 bases and store this samples in an array or list. Rebuild the DNA sequence using this random samples and find out what will be the main problem without the algorithm approach.


import random

filename = "genome.fna"

with open(filename, "r") as file:
    genome = ""

    for line in file:
        if not line.startswith(">"):
            genome += line.strip().upper()

print("Original genome length:", len(genome))


samples = []

for i in range(2000):
    sample_length = random.randint(100, 120)

    start = random.randint(0, len(genome) - sample_length)

    sample = genome[start:start + sample_length]

    samples.append(sample)

print("Total samples:", len(samples))


print("\nFirst 5 random samples:")

for i in range(5):
    print(f"Sample {i + 1}: {samples[i]}")


reconstructed = "".join(samples)

print("\nOriginal genome length:", len(genome))
print("Reconstructed sequence length:", len(reconstructed))


print("\nOriginal genome (first 100 bases):")
print(genome[:100])

print("\nReconstructed sequence (first 100 bases):")
print(reconstructed[:100])

if reconstructed == genome:
    print("\nGenome reconstructed successfully!")
else:
    print("\nGenome reconstruction failed!")
    print("The random samples are not in the correct order.")