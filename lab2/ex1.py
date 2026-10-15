# Bioinformatics - Lab 2 - 7 October
# The melting temperature (Tm) is the temperature at which one half
# of a DNA duplex dissociates and becomes a single strand of DNA.
# The Tm increases with the length of the DNA sequence and with a higher
# G and C content.
# Formula:
# Tm = 4 * (G + C) + 2 * (A + T)
# Task: Make a simple application that calculates the Tm
# based on an input DNA sequence.
# Example sequence: S = ATCGCGTA

sequence = input("Enter DNA sequence: ").upper()

A = sequence.count("A")
T = sequence.count("T")
C = sequence.count("C")
G = sequence.count("G")

Tm = 4 * (G + C) + 2 * (A + T)

print("Melting temperature =", Tm, "°C")