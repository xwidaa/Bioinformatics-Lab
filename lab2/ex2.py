# Bioinformatics - Lab 2 - Exercise 2
# The actual melting temperature (Tm) is influenced by ion concentration.
# Formula:
# Tm = 81.5 + 16.6 * log10([Na+]) + 41 * (%GC) - 600 / length
# Task: Implement an application that calculates the Tm of a DNA sequence.
# Input sequence: ACGCGTGCCA

import math

sequence = input("Enter DNA sequence: ").upper()
na = float(input("Enter Na+ concentration: "))

G = sequence.count("G")
C = sequence.count("C")

length = len(sequence)
gc_content = (G + C) / length

Tm = 81.5 + 16.6 * math.log10(na) + 41 * gc_content - 600 / length

print("Melting temperature =", round(Tm, 2), "°C")