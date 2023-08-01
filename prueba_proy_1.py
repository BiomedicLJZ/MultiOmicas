#Importa las diferntes librerias de las extensiones que se necesitan para la practica
from Bio.Seq import Seq #funcion de Biopython para trabajar con secuencias
from Bio.SeqIO import parse #Funcion de Biopython para trabajar con secuencias de archivos y traducirlas a texto
from tkinter import Tk as tk, filedialog as fd #Funcion de tkinter para crear menu de archivos
from matplotlib import pyplot as plt 
from Bio import SeqIO

#from tqdm import tqdm

#Seccion de trabajo con archivos fasta y su conversion a secuencias manejables

#Por medio de tkinter se crea una interfaz con la que el usuario puede acceder a sus archivos Fasta y seleccionar uno 
root=tk()
root.withdraw()
way=fd.askopenfilename(title="Archivos",filetypes=[("Archivos FASTA","*.fasta"),("Todos los archivos","*.*")])

mi_archivo=parse(way,"fasta") #Transformamos nuestro archivo fasta en una sequencia para poder trabajar con el


strings=[] #Realizamos un ciclo for para poder obtener todas las secuencias de nuestro archivo y con ello guardar cada uno
for strg in mi_archivo:
    strings.append(strg.seq)

#Declaracion de las listas en las cuales se guardaran cada uno de los recorridos segun su marco.
MarcoNormal1=[] #Marco de recorrido a partir del Nucleotido 0 (Primer marco de Lectura)
MarcoNormal2=[] #Marco de recorrido a partir del Nucleotido 1 (Segundo marco de Lectura)
MarcoNormal3=[] #Marco de recorrido a partir del Nucleotido 2 (Tercer marco de Lectura)

#Conjunto de ciclos for, los cuales recorren la secuencia para obtener cada una de las posibles proteinas, dependiendo del Marco a recorrer.
for a in (range(0, len(strings))):#Usamos la variable a como el contador principal para cada una de los marcos recorridos
    #Explicaremos el funcionamiento de un solo ciclo, ya que una explicacion en sencillo comprender el porque de los demas 
    #usamos una ciclo for anidado para cada uno de los marcos teoricos; para este grupos de tres for se realizan las lecturas de los tres
    #marcos de lectura comenzando de la primera posicion a la tercera.
    for s in range(0,len(strings[a]),3): #Se declara que se realice una lectura desde la posicion 0 hasta la longitud de la secuencia, esta de tres en tres
        if strings [a][s:s+3]=="ATG":#Si en la lectura encontramos en este recorrido de tres en tres, encontramos un codon "ATG" comenzaremos a guardar estos codones
            for d in range (s,len(strings[a]),3): #Ciclo for para guardar desde el "ATG" codones de tres en tres
                if strings[a][d:d+3] in ["TAA","TAG","TGA"]:#Hasta que encontremos uno de estos tres condones
                    MarcoNormal1.append(strings[a][s:d+3])#Guardamos esta secuencia en la variable marconormal1
                    break#Detendremos estos ciclos para volverlo a repertir hasta encontrar todas las posibles secuencias
                
    #Lo mismo que el ciclo anterior pero comenzando la lectura desde la posicion 1 y guardar lecturas en la variable marconormal2
    for z in range(1,len(strings[a]),3):
        if strings [a][z:z+3]=="ATG":
            for x in range (z,len(strings[a]),3):
                if strings[a][x:x+3] in ["TAA","TAG","TGA"]:
                    MarcoNormal2.append(strings[a][z:x+3])
                    break
    
    #Lo mismo que el ciclo anterior pero comenzando la lectura desde la posicion 2 y guardar lecturas en la variable marconormal3
    for q in range(2,len(strings[a]),3):
        if strings [a][q:q+3]=="ATG":
            for w in range (q,len(strings[a]),3):
                if strings[a][w:w+3] in ["TAA","TAG","TGA"]:
                    MarcoNormal3.append(strings[a][q:w+3])
                    break

#Declaracion de variables en las cuales por medio de la funcion sorted acomodamos de mayor longitud a menor (key=len), acomodamos al reves la lectura
#(porqur sorted acomoda de menor a mayor) y seleccionamos las primeras cinco secuencias (mas largas) de cada lectura
MarcoNormal1=sorted(MarcoNormal1, key=len, reverse=True)[:5]
MarcoNormal2=sorted(MarcoNormal2, key=len, reverse=True)[:5]
MarcoNormal3=sorted(MarcoNormal3, key=len, reverse=True)[:5]

#Declaracion de las listas en las cuales se guardaran cada uno de los recorridos segun su marco (reverso).
Marcoreversa1=[]#Marco de recorrido a partir del Nucleotido 0, pero en orden inverso(Cuarto marco de Lectura)
Marcoreversa2=[]#Marco de recorrido a partir del Nucleotido 1, pero en orden inverso(Quinto marco de Lectura)
Marcoreversa3=[]#Marco de recorrido a partir del Nucleotido 2, pero en orden inverso(Sexto marco de Lectura)

#Para poder realizar la misma logica de programacion de los for anteriores reacomodamos la secuencia del archivo fasta, pero al orden reverso
Rstrings=list(reversed(strings)) #Guardamos la nueva secuencia inversa en Rstrings en una lista 

#Conjunto de ciclos for, los cuales recorren la secuencia para obtener cada una de las posibles proteinas, dependiendo del Marco a recorrer.
for l in (range(0, len(Rstrings))):

     #Lo mismo que los ciclo anterior pero comenzando la lectura desde la posicion 0 y guardar lecturas en la variable Marcoreversa1
    for k in range(0,len(Rstrings[l]),3):
        if Rstrings [l][k:k+3]=="ATG":
            for j in range (k,len(Rstrings[l]),3):
                if Rstrings[l][j:j+3] in ["TAA","TAG","TGA"]:
                    Marcoreversa1.append(Rstrings[l][k:j+3])
                    break

     #Lo mismo que los ciclo anterior pero comenzando la lectura desde la posicion 1 y guardar lecturas en la variable Marcoreversa2      
    for m in range(1,len(Rstrings[l]),3):
        if Rstrings [l][m:m+3]=="ATG":
            for n in range (m,len(Rstrings[l]),3):
                if Rstrings[l][n:n+3] in ["TAA","TAG","TGA"]:
                    Marcoreversa2.append(Rstrings[l][m:n+3])
                    break

    #Lo mismo que los ciclo anterior pero comenzando la lectura desde la posicion 2 y guardar lecturas en la variable Marcoreversa3  
    for p in range(2,len(Rstrings[l]),3):
        if Rstrings [l][p:p+3]=="ATG":
            for o in range (p,len(Rstrings[l]),3):
                if Rstrings[l][o:o+3] in ["TAA","TAG","TGA"]:
                    Marcoreversa3.append(Rstrings[l][p:o+3])
                    break

#Declaracion de variables en las cuales por medio de la funcion sorted acomodamos de mayor longitud a menor (key=len), acomodamos al reves la lectura
#(porqur sorted acomoda de menor a mayor) y seleccionamos las primeras cinco secuencias (mas largas) de cada lectura
Marcoreversa1=sorted(Marcoreversa1, key=len, reverse=True)[:5]
Marcoreversa2=sorted(Marcoreversa2, key=len, reverse=True)[:5]
Marcoreversa3=sorted(Marcoreversa3, key=len, reverse=True)[:5]

#Unimos todas las listas de marcos de lectura en una nueva variable llamada ORF.
ORF=MarcoNormal1+MarcoNormal2+MarcoNormal3+Marcoreversa1+Marcoreversa2+Marcoreversa3

#Imprimimos los ORF obtenidos
print(ORF)

#Para transformar estas secuencias en un archivo fasta, es necesario trasnformarlas en una variable record
#Todo ello por medio de las funciones que otorga la biblioteca SeqIO, siendo Record para transformalo en esta variable y anotando una descriocion y enumeracion a cada orf
ORFs_records=[SeqIO.SeqRecord(Seq(ORF),id="Human Gamma Herpes Virus 4",description=f"Numero de ORF: '{i+1}'") for i,ORF in enumerate(ORF)]
#Tomamos todos los Records de los ORF para transformarlos en una archivo fasta
#Usamos la funcion write para transformar estos records en archivo fasta.
SeqIO.write(ORFs_records,"Prueba_1.fasta","fasta")