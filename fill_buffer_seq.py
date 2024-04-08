#!/usr/bin/env python3

import os, sys
from itertools import product
import random
import multiprocessing
from Bio import SeqIO
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import subprocess

def rvcomp_re(re_sites):
    re_avoid = set()
    for r in re_sites:
        r_seq = Seq(r)
        r_seq_rc = str(r_seq.reverse_complement())
        re_avoid.add(str(r))
        re_avoid.add(str(r_seq_rc))
    re_avoid = list(re_avoid)
    return re_avoid

def check_full_sequences(sequence, re_avoid):
    # Check if any element in the re_avoid is found in the sequence
    for element in re_avoid:
        if element in sequence:
            return False
    # ensure the full sequence maintains gc content 40-60%
    gc_count = 0
    for i in sequence:
        if i == "G" or i == "C":
            gc_count += 1
    gc_content = gc_count / float(len(sequence))
    if gc_content > 0.6 or gc_content < 0.4:
        return False
    return True

def generate_sequences(dna_sequence):
    n_count = dna_sequence.count('N')
    bases = 'ATGC'
    
    sequences = []
    n=1000
    for _ in range(n):  # Generate n random buffer sequences
        filled_sequence = ''
        for char in dna_sequence:
            if char == 'N':
                filled_sequence += random.choice(bases)
            else:
                filled_sequence += char
        sequences.append(filled_sequence)
    return sequences

def get_n_indices(sequence):
    indices = []
    for i in range(len(sequence)):
        if sequence[i] == "N":
            indices.append(i)
    return indices

def extract_buffers(indices, sequences):
    reduced_sequences = []
    for sequence in sequences:
        reduced_sequence = ''.join([sequence[i] for i in indices])
        reduced_sequences.append(reduced_sequence)
    return reduced_sequences

def perform_blast_search(sequence):
    result_handle = NCBIWWW.qblast("blastn", "nt", sequence, entrez_query="Homo sapiens[organism] OR Mus musculus[organism]")
    blast_records = NCBIXML.parse(result_handle)
    result_handle.close()
    hits_count = 0
    for blast_record in blast_records:
        for alignment in blast_record.alignments:
            for hsp in alignment.hsps:
                if hsp.align_length > 18:
                    hits_count += 1
                    return ""
    return sequence 

def get_seqs_by_name(fasta,names):
    selected_seqs = []
    for record in SeqIO.parse(fasta, "fasta"):
        if record.id in names:
            selected_seqs.append(str(record.seq))
    return selected_seqs

def local_blast_search(sequences,ref,outname):
    db_cmd = f"makeblastdb -dbtype nucl -in {ref}.fasta -out {ref}"
    # use below to build ref database once, then comment out for speed
    # subprocess.run(db_cmd, shell=True)
    cmd = f"blastn -query {sequences} -db {ref} -out {outname}.xml -evalue 0.001 -outfmt 5"
    subprocess.run(cmd, shell=True)
    result_handle = open(f"{outname}.xml")
    blast_records = NCBIXML.parse(result_handle)
    hits_per_query = {}
    for blast_record in blast_records:
        hits_count = 0
        for alignment in blast_record.alignments:
            for hsp in alignment.hsps:
                if hsp.align_length > 10:
                    hits_count += 1
        hits_per_query[blast_record.query] = hits_count
    unmatched = []
    for query,hits in hits_per_query.items():
        if hits == 0:
            unmatched.append(query)
        else:
            print(f"Removed buffer: {query}")
    unmatched_seqs = get_seqs_by_name("results/buffer_sequences.fasta",unmatched)
    result_handle.close()
    return unmatched_seqs



# user set input sequence
dna_sequence = sys.argv[1]
# dna_sequence = "GCAATGNNNNNNNNNNNNNNCTATACGTTGNNNNNNNNNNNNNCGTAGC"

# restriction sites to avoid
re_sites = ["CACCTGC", "GAATTC", "CGTCTC", "CTCGAG", "GCTAGC", "GGTACC", "TCTAGA", "ACCGGT"]
re_avoid = rvcomp_re(re_sites)

# find indices of N characters
n_indices = get_n_indices(dna_sequence)

# generate potential sequence options randomly
sequences = generate_sequences(dna_sequence)

# filter potential full sequences pool based on RE sites and GC content
valid_sequences = []
for sequence in sequences:
    if check_full_sequences(sequence, re_avoid):
        valid_sequences.append(sequence)
if len(valid_sequences) == 0:
    print("Error: no randomly generated sequences meet restriction site and GC content criteria, please re-run and consider increasing the number of initial randomly generated sequences.")
    sys.exit()

# extract filler buffer sequence from full sequences
buffers = extract_buffers(n_indices, valid_sequences)

with open("results/buffer_sequences.fasta", "w") as outfile:
    for i,b in enumerate(buffers):
        outfile.write(">" + "seq" + str(i) + "\n" + b + "\n")
hg = "genomes/GCF_000001405.40_GRCh38.p14_genomic"
outname_hg = "results/GCF_000001405.40_GRCh38.p14_genomic"
mm = "genomes/GCA_000001635.9_GRCm39_genomic"
outname_mm = "results/GCA_000001635.9_GRCm39_genomic"
valid_buffers_hg = local_blast_search("results/buffer_sequences.fasta",hg,outname_hg)
valid_buffers_mm = local_blast_search("results/buffer_sequences.fasta",mm,outname_mm)
valid_buffers = list(set(valid_buffers_hg).intersection(valid_buffers_mm))

if len(valid_buffers) == 0:
    print("Error: no buffer sequences remain after alignment filtering with the human or mouse genome, please try again.")
    sys.exit()

# recombine buffers with the original template sequences in the N positions
final_sequences = []
for buffer in valid_buffers:
    seq = []
    i = 0
    for position in range(len(dna_sequence)):
        if dna_sequence[position] == "N":
            seq.append(buffer[i])
            i+=1
        else:
            seq.append(dna_sequence[position])
    final_sequences.append("".join(seq))

# output first n options to FASTA file
n=100
subset_final_sequences = final_sequences[0:n]
records = [SeqRecord(Seq(seq), id=f"seq{i+1}", description="") for i, seq in enumerate(subset_final_sequences)]
SeqIO.write(records, "results/valid_final_sequences.fasta", "fasta")