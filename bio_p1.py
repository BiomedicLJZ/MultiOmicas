from Bio.Seq import Seq
from Bio.SeqIO import parse
from Bio import SeqIO
from Bio.Data import CodonTable
from tkinter import Tk as tk, filedialog as fd
from matplotlib import pyplot as plt
import matplotlib as mpl
from tqdm import tqdm

mpl.use("QT5Agg")
def get_path():
    '''Opens a file dialog and returns the path of the selected file'''
    root = tk()
    root.withdraw()
    path = fd.askopenfilename(title="Select a file", filetypes=(
        ("fasta files", "*.fasta"),("FASTA Nucleotide Files","*.fna"), ("all files", "*.*")))
    return path

def get_sequences(path):
    '''Reads a fasta file and returns a list of sequences'''
    sequences = []
    for record in parse(path, "fasta"):
        sequences.append(record.seq)
    return sequences


my_seq =Seq("AGTACACTGGACATGCATCGCACGCTTGGCTACACGCTCACGACTACGACTACTGACGTACGTAC")
print(my_seq.complement())
print(my_seq.reverse_complement())
print(my_seq.transcribe())
print(my_seq.translate())


path=get_path()
sequences = get_sequences(path)

print(sequences)

#Calculate GC content as a percentage
gc_content = (sequences[0].count("G") + sequences[0].count("C")) / len(sequences[0]) * 100

#Read all posible kmers in the sequence
kmers=[]
for i in tqdm(range(0,len(sequences))):
    for j in range(0,len(sequences[i]),3):
        kmers.append(sequences[i][j:j+3])
    for k in range(1,len(sequences[i]),3):
        kmers.append(sequences[i][k:k+3])
    for l in range(2,len(sequences[i]),3):
        kmers.append(sequences[i][l:l+3])
    #Read from end to start
    for m in range(len(sequences[i]),0,-3):
        kmers.append(sequences[i][m-3:m])
    for n in range(len(sequences[i])-1,0,-3):
        kmers.append(sequences[i][n-3:n])
    for o in range(len(sequences[i])-2,0,-3):
        kmers.append(sequences[i][o-3:o])

#Convert all kmers to uppercase
kmers=[kmer.upper() for kmer in kmers]
#convert all sequences to uppercase
sequences=[sequence.upper() for sequence in sequences]
#remove all duplicates
kmers=list(set(kmers))
#Count the number of times each kmer appears in the sequence
kmer_counts=[]
for kmer in tqdm(kmers):
    kmer_counts.append(sequences[0].count(kmer))
#Create a dictionary with the kmer and its count
kmer_dict=dict(zip(kmers,kmer_counts))
#Eliminate empty kmers
kmer_dict={kmer:count for kmer,count in kmer_dict.items() if kmer != ""}
#Eliminate kmer with count 0
kmer_dict={kmer:count for kmer,count in kmer_dict.items() if count != 0}
#Sort the dictionary by the number of times each kmer appears
kmer_dict={kmer:count for kmer,count in sorted(kmer_dict.items(),key=lambda item: item[1],reverse=True)}
#convert all keys to string
kmer_dict={str(kmer):count for kmer,count in kmer_dict.items()}
#eliminate kmers with less than 3 nucleotides
kmer_dict={kmer:count for kmer,count in kmer_dict.items() if len(kmer) >= 3}
#Plot all kmers and their counts
fig=plt.figure(figsize=(30,10))
ax=fig.add_axes([0,0,1,1])
ax.bar(kmer_dict.keys(),kmer_dict.values(),color="red")
ax.set_title("Kmer distribution")
ax.set_xlabel("Kmers")
ax.set_ylabel("Counts")
ax.set_xticklabels(kmer_dict.keys(),rotation=90)

#Add the count of each bar
ax.bar_label(ax.containers[0],label_type="edge",fontsize=10,rotation=90,padding=3)
plt.show()
#Get all possible ORFs in the sequence from an ATG start codon to a stop codon (TAA, TAG, TGA) 
ORFs1=[]
ORFs2=[]
ORFs3=[]
ORFs4=[]
ORFs5=[]
ORFs6=[]
for i in tqdm(range(0,len(sequences))):
    for j in range(0,len(sequences[i]),3):
        if sequences[i][j:j+3] == "ATG":
            for k in range(j,len(sequences[i]),3):
                if sequences[i][k:k+3] in ["TAA","TAG","TGA"]:
                    ORFs1.append(sequences[i][j:k+3])
                    break
    for l in range(1,len(sequences[i]),3):
        if sequences[i][l:l+3] == "ATG":
            for m in range(l,len(sequences[i]),3):
                if sequences[i][m:m+3] in ["TAA","TAG","TGA"]:
                    ORFs2.append(sequences[i][l:m+3])
                    break
    for n in range(2,len(sequences[i]),3):
        if sequences[i][n:n+3] == "ATG":
            for o in range(n,len(sequences[i]),3):
                if sequences[i][o:o+3] in ["TAA","TAG","TGA"]:
                    ORFs3.append(sequences[i][n:o+3])
                    break
    #Read from end to start
    for p in range(len(sequences[i]),0,-3):
        if sequences[i][p-3:p] == "ATG":
            for q in range(p,0,-3):
                if sequences[i][q-3:q] in ["TAA","TAG","TGA"]:
                    ORFs4.append(sequences[i][q-3:p])
                    break
    for r in range(len(sequences[i])-1,0,-3):
        if sequences[i][r-3:r] == "ATG":
            for s in range(r,0,-3):
                if sequences[i][s-3:s] in ["TAA","TAG","TGA"]:
                    ORFs5.append(sequences[i][s-3:r])
                    break
    for t in range(len(sequences[i])-2,0,-3):
        if sequences[i][t-3:t] == "ATG":
            for u in range(t,0,-3):
                if sequences[i][u-3:u] in ["TAA","TAG","TGA"]:
                    ORFs6.append(sequences[i][u-3:t])
                    break
#Get 5 longest ORFs from each reading frame
ORFs1=sorted(ORFs1,key=len,reverse=True)[:5]
ORFs2=sorted(ORFs2,key=len,reverse=True)[:5]
ORFs3=sorted(ORFs3,key=len,reverse=True)[:5]
ORFs4=sorted(ORFs4,key=len,reverse=True)[:5]
ORFs5=sorted(ORFs5,key=len,reverse=True)[:5]
ORFs6=sorted(ORFs6,key=len,reverse=True)[:5]
#Combine all ORFs
ORFs=ORFs1+ORFs2+ORFs3+ORFs4+ORFs5+ORFs6
#convert all ORFs to uppercase
ORFs=[ORF.upper() for ORF in ORFs]
#Print 10 longest ORFs
print(sorted(ORFs,key=len,reverse=True)[:10])
#Eliminate ORFs with less than 60 nucleotides
ORFs=[ORF for ORF in ORFs if len(ORF) >= 60]
#Convert all ORFs to SeqIO records
ORFs_records=[SeqIO.SeqRecord(Seq(ORF),id=f"Ophiocordyceps unilateralis Seq: {i+1}",description=f"Possible Seq obtained from NCBI FNA file lenght of the ORF: {len(ORF)} nucleotides or {(len(ORF)/3)} codons") for i,ORF in  (ORFs)]
#Sort ORFs by length
ORFs_records=sorted(ORFs_records,key=lambda x: len(x.seq),reverse=True)
#Write all ORFs to a fasta file
SeqIO.write(ORFs_records,"ORFs.fasta","fasta")