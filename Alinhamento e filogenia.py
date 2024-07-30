from Bio import Entrez
from Bio import pairwise2
from Bio import Align
import math

def PrintMatrix(matrix : list):
    for i in matrix:
        print("")
        for j in i:
            print(f"{j}\t", end='')

Entrez.email = 'dezinho_dh@hotmail.com'


handle = Entrez.esearch('nucleotide', term='16S ribosomal RNA gene[Titl] NOT partial sequence[Titl] NOT uncultured[Titl] NOT clone[Titl]', retmax=10)
results = Entrez.read(handle)

idList = results["IdList"]



fetched = Entrez.efetch(db='nucleotide', id=idList, rettype="fasta")
organismsOriginal = fetched.read().split('>')

for i in range(10):
    organismsOriginal[i] = organismsOriginal[i+1]

listSequences = list()
organismsNames = []

for organism in organismsOriginal:
    organismsNames.append(organism[:organism.find('\n')])
    organism = organism[organism.find('\n'):].replace('\n','')   
    listSequences.append(organism)

count = 0
aligner = Align.PairwiseAligner(match_score=5, mismatch_score=-4, open_gap_score=-10, extend_gap_score=-1, mode='global')
pointsMatrix = []

for i in range(10):
    xthList = []
    for j in range(10):
        newAlignment = 0
        count += 1 
        if i != j:
            newAlignment = aligner.score(listSequences[i], listSequences[j])
            xthList.append(newAlignment)
            continue
        xthList.append(0) 
    pointsMatrix.append(xthList)
    
        
PrintMatrix(pointsMatrix)

distMatrix = list()
print()

for i in pointsMatrix:
    newMatrix = list()
    for j in i:
        newMatrix.append(round(j**-1, 5) if j != 0 else 0)
    distMatrix.append(newMatrix)

PrintMatrix(distMatrix)

        


print("bah")