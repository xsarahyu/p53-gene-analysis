# p53 Gene Analysis 🧬📊

## Overview
This project involves DNA sequence processing, protein synthesis, data analysis, and visualization using Matplotlib and NumPy. The goal is to compare the p53 genes in humans and elephants to identify similarities and differences. The project is divided into three weeks, each focusing on specific tasks and functionalities.

## Project Structure

### Files
- **`protein.py`**: Contains functions for reading DNA, converting DNA to RNA, creating a codon dictionary, generating proteins, and synthesizing proteins from a given DNA sequence file.
- **`hw6_protein_tests.py`**: Includes test cases for functions in `protein.py`. Provides comprehensive testing for different aspects of the protein sequencing process.

### Functions
#### `readFile(filename)`
- Reads the DNA sequence from the specified file and returns it as a string.
#### `dnaToRna(dna, startIndex)`
- Converts a given DNA sequence to RNA, starting from the specified index.
- Returns a list of codons.
#### `makeCodonDictionary()`
- Creates a dictionary mapping codons to amino acids using data from the codon table JSON file.
- Returns the codon dictionary.
#### `generateProtein(codons, codonD)`
- Translates a list of codons into amino acids using the provided codon dictionary.
- Returns the resulting protein sequence.
#### `synthesizeProteins(filename)`
- Synthesizes proteins from a DNA sequence file.
- Utilizes the previously defined functions to achieve protein synthesis.
- The synthesized proteins are returned as a list.
#### `commonProteins(proteinList1, proteinList2)`
- Takes two lists of proteins (where each protein is a list of amino acids) and returns a list of all the unique proteins that occur in both lists.
- Each protein occurs only once in the result list, even if it shows up multiple times in both genes.
#### `combineProteins(proteinList)`
- Collapses the list of proteins into a list of amino acids that occur across all proteins.
- Returns a 1D list of amino acids in their original order.
#### `aminoAcidDictionary(aaList)`
- Generates a dictionary mapping each amino acid in the list to a count of how often it occurs.
- Returns the amino acid dictionary.
#### `findAminoAcidDifferences(proteinList1, proteinList2, cutoff)`
- Compares amino acids between two genes by analyzing their frequencies.
- Returns a list of three-element lists, where each element represents an amino acid, its frequency in proteinList1, and its frequency in proteinList2.
- Only includes amino acids if the difference in frequencies is greater than the provided cutoff.
#### `displayTextResults(commonalities, differences)`
- Prints out a textual report of common proteins and most-different amino acids between two genes.
- Displays common proteins and amino acid frequencies for each sequence.
#### `makeAminoAcidLabels(proteinList1, proteinList2)`
- Finds all the amino acids that occur across both genes.
- Returns a sorted list of amino acids.
#### `setupChartData(labels, proteinList)`
- Processes data into the format expected by Matplotlib for chart plotting.
- Returns a list of frequencies for each amino acid in a particular gene.
#### `createChart(xLabels, freqList1, label1, freqList2, label2, edgeList=None)`
- Draws a side-by-side bar chart comparing the frequency of amino acids in two genes.
- Uses Matplotlib to create the chart with specified labels and frequencies.
- Optionally sets the edge color of amino acid bars based on the provided edgeList.
#### `makeEdgeList(labels, biggestDiffs)`
- Creates a list indicating whether each amino acid is in the biggestDiffs list.
- Returns a list with "black" for amino acids in biggestDiffs and "white" for others.
#### `runFullProgram()`
- Loads DNA data from two p53 files, processes them into protein lists, generates a text report comparing the two genes, and produces a bar chart comparing the two genes.
- Uses a 0.5% cutoff for differences in amino acid frequencies and outlines sufficiently-different amino acids in black on the chart.

### Testing
The file `hw6_protein_tests.py` contains a set of test cases to ensure accuracy of the implemented functions. The tests cover various aspects, including reading files, DNA to RNA conversion, codon dictionary creation, protein generation, overall protein synthesis, and data analysis using Matplotlib and NumPy.

## Running Tests
To run the tests, execute the following command:
```bash
python hw6_protein_tests.py