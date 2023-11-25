"""
15-110 Hw6 - Protein Sequencing Project
Name: Sarah Yu
AndrewID: sayu
"""

import hw6_protein_tests as test

project = "Protein" # don't edit this

### WEEK 1 ###

def readFile(filename):
    f = open(filename, 'r')
    dna = f.read()
    dna = dna.replace('\n', '')
    f.close()
    return dna


def dnaToRna(dna, startIndex):
    dna = dna.replace('T', 'U')

    codons = []
    for base in range(startIndex, len(dna), 3):
        codon = dna[base:(base + 3)]

        if codon == 'UAA' or codon == 'UAG' or codon == 'UGA':
            codons.append(codon)
            return codons

        codons.append(codon)

    return codons #["AUG", "GAU", "GGA", "CUC", "UAA"]


def makeCodonDictionary():
    import json
    f = open('data/codon_table.json', 'r')
    aminoD = json.load(f)

    codonD = {}
    for amino in aminoD:
        key = aminoD[amino]
        for value in key:
            value = value.replace('T', 'U')
            codonD[value] = amino

    return codonD #{'AUG': 'Met', 'GAU': 'Asp', ...}
    f.close()


def generateProtein(codons, codonD):
    protein = []

    for i in range(len(codons)):
        if codonD[codons[i]] == 'Met' and i == 0:
            protein.append('Start')

        else:
            protein.append(codonD[codons[i]])

    return protein


def synthesizeProteins(filename):
    dna = readFile(filename)
    codonD = makeCodonDictionary()

    proteinsMade = []
    unusedCount = 0
    base = 0

    while base < len(dna):
        codon = dna[base:(base + 3)]
        if codon == 'ATG':
            codons = dnaToRna(dna, base)
            protein = generateProtein(codons, codonD)
            proteinsMade.append(protein)
            base += (3*len(codons))

        else:
            base += 1
            unusedCount += 1

    print('Total number of bases:', len(dna))
    print('Unused-base count:', unusedCount)
    print('Total number of proteins synthesized:', len(proteinsMade))

    return proteinsMade


def runWeek1():
    print("Human DNA")
    humanProteins = synthesizeProteins("data/human_p53.txt")
    print("Elephant DNA")
    elephantProteins = synthesizeProteins("data/elephant_p53.txt")


### WEEK 2 ###

def commonProteins(proteinList1, proteinList2):
    commonProteins = []
    for proteins in proteinList1:
        if proteins in proteinList2 and proteins not in commonProteins:
            commonProteins += [proteins]
    return commonProteins


def combineProteins(proteinList):
    aaList = []
    for proteins in proteinList:
        aaList += proteins
    return aaList


def aminoAcidDictionary(aaList):
    aaDic = {}
    for aa in aaList:
        if aa not in aaDic:
            aaDic[aa] = 1
        elif aa in aaDic:
            aaDic[aa] += 1
    return aaDic


def findAminoAcidDifferences(proteinList1, proteinList2, cutoff):
    result = []
    aaList1 = combineProteins(proteinList1)
    aaList2 = combineProteins(proteinList2)
    aaDic1 = aminoAcidDictionary(aaList1)
    aaDic2 = aminoAcidDictionary(aaList2)

    aaSum1 = 0
    for aa in aaDic1:
        aaSum1 += aaDic1[aa] #total num of amino acids made from gene 1

    aaSum2 = 0
    for aa in aaDic2:
        aaSum2 += aaDic2[aa] #total num of amino acids made from gene 2

    aaDic1.pop('Start')
    aaDic1.pop('Stop')
    aaDic2.pop('Start')
    aaDic2.pop('Stop')

    aaSeen = []
    for aa in aaDic1:
        if aa in aaDic2:
            aaFreq1 = aaDic1[aa] / aaSum1 #freq of an amino acid from gene 1
            aaFreq2 = aaDic2[aa] / aaSum2 #freq of an amino acid from gene 2
            aaDiff = abs(aaFreq1 - aaFreq2) #diff btwn freqs

            if aaDiff > cutoff:
                result.append([aa, aaFreq1, aaFreq2])
        else: #if aa not in aaDic2, that aa's freq would be 0
            aaFreq1 = aaDic1[aa] / aaSum1
            aaFreq2 = 0
            if aaFreq1 > cutoff:
                result.append([aa, aaFreq1, aaFreq2])
        aaSeen.append(aa)

    for aa in aaDic2:
        if aa not in aaSeen:
            aaFreq1 = 0
            aaFreq2 = aaDic2[aa] / aaSum2
            if aaFreq2 > cutoff:
                result.append([aa, aaFreq1, aaFreq2])
    return result


def displayTextResults(commonalities, differences):
    print('The following proteins occurred in both DNA sequences:')
    for proteins in commonalities:
        proteins.remove('Start')
        proteins.remove('Stop')
        print('-'.join(proteins))

    diffStr = ''
    for lst in differences:
        aaFreq1 = round(float(lst[1]) * 100, 2)
        aaFreq2 = round(float(lst[2]) * 100, 2)
        diffStr = diffStr + str(lst[0]) + ': ' + str(aaFreq1) + '% in Seq1, ' + str(aaFreq2) + '% in Seq2' + '\n'
    print('The following amino acids occurred at very different rates in the two DNA sequences:' + '\n' + diffStr)


def runWeek2():
    humanProteins = synthesizeProteins("data/human_p53.txt")
    elephantProteins = synthesizeProteins("data/elephant_p53.txt")

    commonalities = commonProteins(humanProteins, elephantProteins)
    differences = findAminoAcidDifferences(humanProteins, elephantProteins, 0.005)
    displayTextResults(commonalities, differences)


### WEEK 3 ###

def makeAminoAcidLabels(proteinList1, proteinList2):
    labels = []
    aaList1 = combineProteins(proteinList1)
    aaList2 = combineProteins(proteinList2)

    for aa in aaList1:
        if aa not in labels:
            labels.append(aa)

    for aa in aaList2:
        if aa not in labels:
            labels.append(aa)

    labels.sort()
    return labels


def setupChartData(labels, proteinList):
    aaList = combineProteins(proteinList)
    aaDic = aminoAcidDictionary(aaList)
    aaDicValues = aaDic.values()
    freq = []
    for aa in labels:
        if aa in aaDic:
            freq += [aaDic[aa] / sum(aaDicValues)]
        else:
            freq += [0]
    return freq


def createChart(xLabels, freqList1, label1, freqList2, label2, edgeList=None):
    #this function is taken from lecture 12-3 (Wednesday)
    import matplotlib.pyplot as plt
    import numpy as np
    x = np.arange(len(xLabels))
    width = 0.35

    fig, ax = plt.subplots()
    rects1 = ax.bar(x - width/2, freqList1, width, label = label1, edgecolor = edgeList)
    rects2 = ax.bar(x + width/2, freqList2, width, label = label2, edgecolor = edgeList)

    ax.set_xticks(x)
    ax.set_xticklabels(xLabels)
    ax.legend()

    plt.show()


def makeEdgeList(labels, biggestDiffs):
    edgeList = []
    for aa in labels:
        if aa in (i[0] for i in biggestDiffs):
            edgeList += ['black']
        else:
            edgeList += ['white']
    return edgeList


def runFullProgram():
    humanProteins = synthesizeProteins("data/human_p53.txt")
    elephantProteins = synthesizeProteins("data/elephant_p53.txt")

    commonalities = commonProteins(humanProteins, elephantProteins)
    differences = findAminoAcidDifferences(humanProteins, elephantProteins, 0.005)
    displayTextResults(commonalities, differences)

    xLabels = makeAminoAcidLabels(humanProteins, elephantProteins)
    label1 = 'Human proteins'
    freqList1 = setupChartData(xLabels, humanProteins)
    label2 = 'Elephant proteins'
    freqList2 = setupChartData(xLabels, elephantProteins)
    edgeList = makeEdgeList(xLabels, differences)
    barChart = createChart(xLabels, freqList1, label1, freqList2, label2, edgeList)

    return barChart


### RUN CODE ###

# This code runs the test cases to check your work
if __name__ == "__main__":
    print("\n" + "#"*15 + " WEEK 1 TESTS " +  "#" * 16 + "\n")
    test.week1Tests()
    print("\n" + "#"*15 + " WEEK 1 OUTPUT " + "#" * 15 + "\n")
    runWeek1()

    ## Uncomment these for Week 2 ##

    print("\n" + "#"*15 + " WEEK 2 TESTS " +  "#" * 16 + "\n")
    test.week2Tests()
    print("\n" + "#"*15 + " WEEK 2 OUTPUT " + "#" * 15 + "\n")
    runWeek2()


    ## Uncomment these for Week 3 ##

    print("\n" + "#"*15 + " WEEK 3 TESTS " +  "#" * 16 + "\n")
    test.week3Tests()
    print("\n" + "#"*15 + " WEEK 3 OUTPUT " + "#" * 15 + "\n")
    runFullProgram()

