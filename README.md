# SynMockPlus
Scripts to generate simulated barcoding gene sequences (AKA "synthetic mock" sequences).

**The flexible script can be copied and modified for additional genes. This works by editting the file to specify the main features to be included**:
1. Total gene sequence length
2. Percent GC content
3. Conserved regions, typically primer/probe sequences to be inserted, and their respective base positions within the total sequence

The resulting gene-specific scripts are found in the Gene Scripts folder.  

## Instructions for use:
1. Install Python dependencies: [Biopython](https://pypi.org/project/biopython/), [Natsort](https://pypi.org/project/natsort/)    
`pip install biopython natsort`
2. Open a terminal within the script folder, then run the desired gene script. By default, the sequences are output to the terminal in fasta syntax (stdout). Each execution of the script generates new random sequences. The intended use is to redirect stdout to a proper fasta file.  
`python ITS_synthetic_mock.py > SynMock_ITS.fasta`

### In progress:
1. Finalizing conserved primers for CO1, 18S, 16S, rbcL, and trnL.
2. Generating sequence maps that show the structure and conserved primers within the random sequences.