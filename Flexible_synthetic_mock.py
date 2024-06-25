#!/usr/bin/env python

import random
from natsort import natsorted
from Bio.SeqUtils import gc_fraction as GC

def random_seq(length, gc_content, max_homopolymer=4):
    """Generates random DNA sequence with a specified length, GC content, and max homopolymer"""
    # check if length is odd. If so, add 1 for calculations, then truncate the last base at the end
    odd_length = False
    if length % 2 != 0:
        length += 1
        odd_length = True

    gc_count = int(length * gc_content / 2)
    at_count = int(length / 2) - gc_count
    # Create a list of bases with the calculated AT/GC counts, uses a list to allow shuffling function
    bases = ["A"] * at_count + ["T"] * at_count + ["G"] * gc_count + ["C"] * gc_count
    # Shuffle the bases to add randomness
    random.shuffle(bases)
    # Join the bases from list into a string
    seq = "".join(bases)
    # Check for homopolymers, reshuffle if needed
    while any(base * max_homopolymer in seq for base in "ATGC"):
        seq = list(seq)
        random.shuffle(seq)
        seq = "".join(seq)
    
    # Remove extra base if odd length
    if odd_length:
        seq = seq[:-1]

    return seq

def generate_sequence(num_seqs, conserved_regions, total_length, gc_content=0.50):
    """Builds simulated sequences with specified conserved regions and total length."""
    
    # Perform couple checks for the input values
    # Check if total length is long enough for all conserved regions
    conserved_length = sum(len(seq) for _, seq in conserved_regions)
    random_length = total_length - conserved_length
    if random_length < 0:
        raise ValueError("Total length is too short to accommodate conserved regions.")
    
    # Check if conserved regions overlap
    for i in range(len(conserved_regions)-1):
        if conserved_regions[i][0] + len(conserved_regions[i][1]) > conserved_regions[i+1][0]:
            raise ValueError("Conserved regions overlap, please check the positions.")
    
    # Initialize sequences dictionary
    sequences_dict = {}
    # Loop for desired number of sequences
    for i in range(num_seqs):
        # Initialize the new sequence as an empty string
        sequence = ""
        # Track the current position in the sequence
        current_pos = 0

        for position, conserved_seq in conserved_regions:
            # Generate random sequence to fill the gap before the next conserved region
            gap_length = position - current_pos
            if gap_length > 0:
                sequence += random_seq(gap_length, gc_content)
            # Insert the conserved sequence
            sequence += conserved_seq
            # Update current position
            current_pos = position + len(conserved_seq)  

        # Generate random sequence to fill the remaining space after the last conserved region
        if current_pos < total_length:
            sequence += random_seq((total_length - current_pos), gc_content)

        # Calculate total GC content
        total_gc = round(GC(sequence) * 100)

        name = f"SynMock_{i+1};length={len(sequence)};GC={total_gc}"
        sequences_dict[name] = sequence

    return sequences_dict


def main():
    # Define conserved regions with positions (ie, primer sequences to be included)
    # Format: (position, sequence)
    conserved_regions = [
        (0, "TTCTCAACCAACCACAAAGATATTGGATGAACTGTATATCCTCC"),  # COI_1F
        (300, "GGWACWGGWTGAACWGTWTAYCCYCC"),                # COI_375F
        (600, "AAGCCAATCTGATTCTTCGGTCACCCAGAAGTTTA")  # COI_750R
    ]

    total_length = 700  # adjust as needed

    sequences = generate_sequence(10, conserved_regions, total_length, 0.35)

    # Output sequences in FASTA format
    for name, sequence in natsorted(sequences.items()):
        print(f">{name}\n{sequence}")

if __name__ == "__main__":
    main()
