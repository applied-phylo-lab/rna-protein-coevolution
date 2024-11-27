import sys
import RNA

from Bio import SeqIO
from numpy import zeros
from ast import literal_eval
import re

WINDOW = 35
FASTA_FILE ="human_transcripts_for_genes_in_ba_etal.fasta"
SEQS = SeqIO.index(FASTA_FILE, "fasta")

five_prime_cap = False

def fold_five_prime_cap(seq, WINDOW):
  energy = RNA.pf_fold(seq[0:WINDOW])[1]
  return energy

def fold_transcript(seq, WINDOW):
    energies = zeros(len(seq) - WINDOW)
    for i in range(len(seq) - WINDOW):
        energies[i] = RNA.pf_fold(seq[i:i+WINDOW])[1]
    return energies

if five_prime_cap:
  with open("human_ss_five_prime_cap_70nt.tsv", "w+") as f:
    f.write("\t".join(("TranscriptID", "ss_efe_profile")) + "\n")
    for transcript in SEQS:
      seq = str(SEQS[transcript].seq)
      if len(seq) < 100: continue
      f.write(transcript + "\t" + str(fold_five_prime_cap(seq, WINDOW)) + "\n")
else :
  pat = re.compile(r'(?<=CDS:)[0-9]+')
  with open("human_ss_around_start_site.tsv", "w+") as f:
    f.write("TranscriptID\t"+"\t".join([str(i) for i in range(-50,51)]) + "\n")
    for transcript in SEQS:
      seq = str(SEQS[transcript].seq)
      desc = SEQS[transcript].description
      cds_start = int(re.findall(pat,desc)[0]) - 1 #Python is 0-indexed, CDS coordinates are not
      seq_start = seq[(cds_start - 50):(cds_start + 86)]
      if seq_start[50:53] != "ATG":
        print(desc)
        print(seq[cds_start:(cds_start+3)])
      if len(seq_start) < 136: continue
      f.write(transcript + "\t" + "\t".join([str(i) for i in fold_transcript(seq_start, WINDOW)]) + "\n")
    
