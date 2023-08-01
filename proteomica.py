from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import Entrez
import os
import pandas as pd
import numpy as np
from tkinter import Tk,filedialog

# Importar datos
root = Tk()
root.withdraw()
root.attributes("-topmost", True)
file_path = filedialog.askopenfilename(filetypes=[("FASTA files", "*.fasta"),("FASTA Nucleotide files","*.fna"), ("All files", "*.*")])
root.destroy()

# Leer archivo FASTA
data= SeqIO.parse(file_path, "fasta")

seqs=[]
for record in data:
    seqs.append(record.seq)

# Importar datos de proteínas
root = Tk()
root.withdraw()
root.attributes("-topmost", True)
prot_file_path = filedialog.askopenfilename(filetypes=[("FASTA files", "*.fasta"),("FASTA Nucleotide files","*.fna"), ("All files", "*.*")])
root.destroy()

# Leer archivo FASTA
prot_data= SeqIO.parse(prot_file_path, "fasta")
# Crear lista de secuencias
prot_seqs=[]
for record in prot_data:
    prot_seqs.append(record.seq)

# Traducir secuencias de nucleótidos
proto_prot_seqs=[]
for seq in seqs:
    proto_prot_seqs.append(seq.translate())

# Comparar secuencias de proteínas reales con las traducidas de nucleótidos
seqs_comp=[]
for seq in proto_prot_seqs:
    seqs_comp.append(seq in prot_seqs)
