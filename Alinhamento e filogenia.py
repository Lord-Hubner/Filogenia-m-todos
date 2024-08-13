from Bio import Entrez
from Bio import Align
from Bio import Phylo
from io import StringIO
import re
import string
from matplotlib import pyplot as plt
import random

def PrintMatrix(matrix : list):
    for i in matrix:
        print("")
        for j in i:
            print(f"{j}\t", end='')


def SearchMinimum(distMatrix: list):
    currentMinimum: int = 99999
    currentMinimumIndexes = []
    #i=linhas, j=colunas, k= alvo
    i=0
    for k in range(1,10):    
        for j in range(k):   
            if distMatrix[i][j] != 0 and distMatrix[i][j] < currentMinimum:
                currentMinimum = distMatrix[i][j]
                currentMinimumIndexes = [i, j]
        i=i+1
    #Última iteração
    for j in range(k):    
        if distMatrix[i][j] != 0 and distMatrix[i][j] < currentMinimum:
            currentMinimum = distMatrix[i][j]
            currentMinimumIndexes = [i, j]
    

    return currentMinimum, currentMinimumIndexes

def BalanceParenthesis(newParenthesisIndex: int, tree: str):
    balance = 0
    for i, c in enumerate(tree[newParenthesisIndex:]):
        # if(re.match(r',[A-Z],[A-Z]', tree[i:i+4])):
        #     return i
        if c == '(':
            balance+=1
        elif c == ')':
            balance-=1
        if balance==1 and (c == ','):
            return i


def ModifyTree(tree: str, minimumIndexesDict: dict, repeatedIndex: int, notRepeatedIndex: int, organismsCodes: list, counter: int):
        repeatedPositionIndex = tree.find(minimumIndexesDict[repeatedIndex])
        # Se é folha 
        # ehFolha = tree[repeatedIndex-1] and tree[repeatedIndex-3]!= '('
        # if ehFolha:
        #     newString = f",{alphabet[counter]})"
        #     tree=tree[:repeatedIndex]+'('+tree[repeatedIndex:repeatedIndex+1]+newString+tree[repeatedIndex+1:]
        teste = tree[repeatedPositionIndex::-1]
        newParenthesisIndex = tree.find('(', 0, repeatedPositionIndex)
        for i in range(newParenthesisIndex, 0, -1):
            if tree[i] == '(':
                newParenthesisIndex-=1
            elif tree[i] != '(':
                break

        tree=tree[:newParenthesisIndex]+'('+tree[newParenthesisIndex:]
        repeatedPositionIndex+=1
        indexToGo = BalanceParenthesis(newParenthesisIndex, tree) + newParenthesisIndex

        addVirgula = True
        if addVirgula:
            tree=tree[:indexToGo]+','+organismsCodes[notRepeatedIndex]+')'+tree[indexToGo:]
        else:
            tree=tree[:indexToGo]+organismsCodes[notRepeatedIndex]+tree[indexToGo:]
        # newString = f",{alphabet[counter]})"
        # tree=tree[:indexToGo+1]+newString+tree[indexToGo+1:]

        minimumIndexesDict[notRepeatedIndex] = organismsCodes[notRepeatedIndex]

        return tree, minimumIndexesDict

def ConnectBranches(tree: str, indexToConnect: int, mainBranchIndex: int):

    newParenthesisIndex = tree.find('(', 0, mainBranchIndex)
    for i in range(newParenthesisIndex, 0, -1):
        if tree[i] == '(':
            newParenthesisIndex-=1
        elif tree[i] != '(':
            break
    
    tree=tree[:newParenthesisIndex]+'('+tree[newParenthesisIndex:]

    indexToGo = BalanceParenthesis(0, tree)
    
    target = tree[indexToGo:]
    match = re.search(r"\([A-Z]{2}\d{6}\.1,[A-Z]{2}\d{6}\.1\)", target)
    branchToTransfer = target[match.start():match.end()]
    tree=tree[:tree.rfind(branchToTransfer)]

    tree=tree[:indexToGo]+','+branchToTransfer+')'+tree[indexToGo:]
    return tree


def BuildTree(distMatrix: list, organismsCodes: list) -> tuple[str, dict]:
    tree:str = ""
    alphabet:list = list(string.ascii_uppercase)
    minimumIndexesdict:dict = dict()
    currentMinimum:int = 100000
    currentMinimumIndexes:list = []
    secondaryBranchsSet:list = list()


    #Faz primeira iteração fora do loop pois é previsível, e para não atribuir a um galho secundário
    currentMinimum, currentMinimumIndexes = SearchMinimum(distMatrix)

    distMatrix[currentMinimumIndexes[0]][currentMinimumIndexes[1]] = 99999
    tree = tree+'('+organismsCodes[currentMinimumIndexes[0]]+','+organismsCodes[currentMinimumIndexes[1]]+'),'

    minimumIndexesdict[currentMinimumIndexes[0]] = organismsCodes[currentMinimumIndexes[0]]
    minimumIndexesdict[currentMinimumIndexes[1]] = organismsCodes[currentMinimumIndexes[1]]
    
    counter=2
    while currentMinimum != 99999:
        firstIndexIsRepeated=False
        secondIndexIsRepeated=False
        
        currentMinimum, currentMinimumIndexes = SearchMinimum(distMatrix)
        if(currentMinimum) == 99999:
            break
        distMatrix[currentMinimumIndexes[0]][currentMinimumIndexes[1]] = 99999

        if currentMinimumIndexes[0] in minimumIndexesdict:
            firstIndexIsRepeated=True
         
        if currentMinimumIndexes[1] in minimumIndexesdict:
            secondIndexIsRepeated=True

        if firstIndexIsRepeated and secondIndexIsRepeated:
            for item in secondaryBranchsSet:
                for index in item:
                    if index == currentMinimumIndexes[0]:
                        tree = ConnectBranches(tree, currentMinimumIndexes[0], currentMinimumIndexes[1])
                        secondaryBranchsSet.remove(item)
                    if index == currentMinimumIndexes[1]:
                        tree = ConnectBranches(tree, currentMinimumIndexes[1], currentMinimumIndexes[0])
                        secondaryBranchsSet.remove(item)

            if len(minimumIndexesdict) == len(distMatrix):
                if tree[-1] == ',':
                    tree = tree[:-1]
                tree = tree+';'  
                return tree, minimumIndexesdict
            else:
                continue

        elif not firstIndexIsRepeated and not secondIndexIsRepeated:
            tree = tree+'('+organismsCodes[currentMinimumIndexes[0]]+','+organismsCodes[currentMinimumIndexes[1]]+'),'
            minimumIndexesdict[currentMinimumIndexes[0]] = organismsCodes[currentMinimumIndexes[0]]
            minimumIndexesdict[currentMinimumIndexes[1]] = organismsCodes[currentMinimumIndexes[1]]

            secondaryBranchsSet.append([max(minimumIndexesdict), min(minimumIndexesdict)])

        elif (firstIndexIsRepeated and not secondIndexIsRepeated):
            tree, minimumIndexesdict = ModifyTree(tree, minimumIndexesdict, currentMinimumIndexes[0], currentMinimumIndexes[1], organismsCodes, counter)
       

        elif (not firstIndexIsRepeated and secondIndexIsRepeated):
            tree, minimumIndexesdict = ModifyTree(tree, minimumIndexesdict, currentMinimumIndexes[1], currentMinimumIndexes[0], organismsCodes, counter)
       
        
        if firstIndexIsRepeated or secondIndexIsRepeated:
            counter+=1
        else:
            counter+=2

    return tree




Entrez.email = 'dezinho_dh@hotmail.com'

handle = Entrez.esearch('nucleotide', term='16S ribosomal RNA gene[Titl] NOT partial sequence[Titl] NOT uncultured[Titl] NOT clone[Titl]', retmax=100)
results = Entrez.read(handle)

idList = random.sample(results["IdList"], 10)

fetched = Entrez.efetch(db='nucleotide', id=idList, rettype="fasta")
organismsOriginal = fetched.read().split('>')[1:]

listSequences = list()
organismsNames = []
organismsCodes = []

with open('seqs.fasta', 'w') as file:
    for item in organismsOriginal:
        file.write('>'+item[:item.rfind('\n')])

for organism in organismsOriginal:
    organismsNames.append(organism[:organism.find('\n')])
    organismsCodes.append(organism[:organism.find(' ')])
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
        newMatrix.append(0 if j == 0 else (-j) if j < 0 else round(j**-1, 5))
    distMatrix.append(newMatrix)

    

PrintMatrix(distMatrix)

indicesDict: dict
tree, indicesDict = BuildTree(distMatrix, organismsCodes)

indicesToAdd = 0
i=0
newTree=tree

# while i <= len(newTree):
#     current = newTree[i:i+10]
#     if re.match(r'(PP|ON)[0-9]{6}\.1', current):
#         for item in organismsNames:
#             if item[0:10] == current:        
#                 newTree= newTree[:i]+item[:len(item) if len(item) < 50 else 50]+newTree[i+10:]
#                 i += 30
#                 break
#     i += 1

handle = StringIO(newTree)
treeFile = Phylo.read(handle, "newick")

for terminal in treeFile.get_terminals():
    for i, code in enumerate(organismsCodes):
        if terminal.name == code:
            terminal.name = organismsNames[i]
            break

# handle2 = StringIO(tree)
# treeFile2 = Phylo.read(handle2, "newick")



        # busca = Entrez.esearch(db='nucleotide', term=current, rettype="fasta")
        # bah = Entrez.read(busca)
        # id = bah["IdList"]
        # fetched = Entrez.efetch(db='nucletide', id=id, rettype='fasta')
        # result = fetched.read()
        # tree[i:i+10] = bah


# matches = re.search(r'(PP|ON)[0-9]{6}\.1', tree);
# for match in matches.groups():
#     print(match)


# Phylo.draw_ascii(treeFile);
# Phylo.draw(treeFile)    

Phylo.draw_ascii(treeFile);
Phylo.draw(treeFile)    


        


print("bah")