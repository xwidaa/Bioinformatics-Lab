S = "ACGTTGCATTACGTTGCAGT"

for sequence in "ACGT":
    frequency = S.count(sequence) / len(S) * 100
    print(sequence, ":", frequency, "%")