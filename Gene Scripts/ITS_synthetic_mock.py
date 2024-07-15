# Synthetic mock sequence generator for ITS barcoding gene
# Author: Zach Geurin, 2024-07-03
# Modified script from Amptk to generate simulated gene sequences for mock community

import random
from natsort import natsorted
from Bio.SeqUtils import gc_fraction as GC

def random_seq(length, gc_content, max_homopolymer=3):
    """Generates random DNA sequence with a specified length, GC content, and max homopolymer"""
    gc_count = int(length * gc_content / 2)
    at_count = int(length / 2) - gc_count
    # Create a list of bases with the calculated AT/GC counts, uses a list to allow shuffling function
    bases = ["A"] * at_count + ["T"] * at_count + ["G"] * gc_count + ["C"] * gc_count
    # Shuffle the bases to add randomness
    random.shuffle(bases)
    # Join the bases from list into a string
    seq = "".join(bases)
    # Check and break up long homopolymers
    while any(base * max_homopolymer in seq for base in "ATGC"):
        seq = list(seq)
        random.shuffle(seq)
        seq = "".join(seq)

    return seq

def main():
    # Conserved sequences including primers
    SSU = "CTTGGTCATTTAGAGGAAGTAAAAGTCGTAACAAGGTTTCCGTAGGTGAACCTGCGG"
    FES = "CAACAACGGATCTCTTGGTTCTCGCATCGATGAAGAACGCAGCNNNNNNNNNNNNNNNNNNGTGAATCATCGAATCTTTGAACGCACATTGCGCTCCTTGGTAT"
    # LSU_plus includes extended range for LR3 & LR5 primers
    LSU_plus = "TGACCTCAAATCAGGTAGGAGTACCCGCTGAACTTAAGCATATCAATAAGCGGAGGATGACAATTAACCACCGTGTATTCGTCCCGTCTTGAAACACGGTATAACATCAGGCAGTTTAAGTCGGCGAAGTTTCCCTCAGGA"
    
    ''' Extra primers:
    ITS2 = "GCTGCGTTCTTCATCGATGC"
    fITS7 = "GTGARTCATCGAATCTTTG" 
    ITS3 = "GCATCGATGAAGAACGCAGC"
    ITS3_KYO1 = "AHCGATGAAGAACRYAG"
    fITS9 = "GAACGCAGCRAAIIGYGA"
    58A1F	= "GCATCGATGAAGAACGC"
    58A2F	= "ATCGATGAAGAACGCAG"
    58A2R	= "CTGCGTTCTTCATCGAT"

    LR3 = "CCGTGTTTCAAGACGGG"
    LR3_RC = "CCCGTCTTGAAACACGG" # reverse complement to be integrated

    LR5 = "TCCTGAGGGAAACTTCG"
    LR5_RC = "CGAAGTTTCCCTCAGGA" # reverse complement to be integrated
    
    LSU = "TGACCTCAAATCAGGTAGGAGTACCCGCTGAACTTAAGCATATCAATAAGCGGAGGA" # LSU without extended range
    '''

    #Create single qPCR probe 18 bp for use in all sequences
    probe = random_seq(18, 0.65)
       
    # Construct full sequences with conserved regions and probe
    sequences = {}
    random_length = 100
    for i in range(1,16): # change range if more sequences are needed
        if i < 6: # first 5 sequences have 200 bp for ITS1 and ITS2
            random_its1 = random_seq(200, 0.55)
            random_its2 = random_seq(200, 0.55)
        else:
            random_its1 = random_seq(random_length, 0.55)
            random_its2 = random_seq(500-random_length, 0.55)
            random_length += (500//15)

        its1_len = len(random_its1)
        its2_len = len(random_its2)
        its1_gc = round(GC(random_its1) * 100)
        its2_gc = round(GC(random_its2) * 100)
        # build the full sequence
        sequence = (
            SSU # SSU with ITS1F primer 
            + random_its1
            + FES.replace("NNNNNNNNNNNNNNNNNN", probe) # 5.8S with qPCR probe, ITS2_R, and ITS7_F primers 
            + random_its2
            + LSU_plus # LSU extended with LR5 and LR3 primers
        )
        name = f"SynMock_{i};qPCR_probe={probe};Len_ITS1={its1_len};GC_ITS1={its1_gc};Len_ITS2={its2_len};GC_ITS2={its2_gc}"
        sequences[name] = sequence

    # Output sequences in FASTA format
    for name, sequence in natsorted(sequences.items()):
        print(f">{name}\n{sequence}")

if __name__ == "__main__":
    main()