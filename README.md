# Refoldibility-Tools
Tools to analyze proteome-wide protein folding experimental data from mass spectrometry. Analyzer is a python script written to analyze raw normalized extracted ion intensity data from Label Free Quantification (LFQ) of detected peptides in refolded and native LiP-MS samples. 

## Running Analyzer

Analyzer code can be run with a stock install of anaconda, which should contain all the required dependencies. Raw ion intensity data need to be first exported from the .pdResult file to a three-level hierarchy (protein > peptide group > consensus feature) to excel and then saved as a .txt file. Analyzer code and associated meta-data files should be moved to the same directory as raw ion intensity data .txt files to run Analyzer code. 

##Input and Output for Analyzer 

For example, once in the correct directory, in ipython, the following command below can be used to run the analyzer code. The exported raw ion intensity data .txt file is inputted as the first argument after the Analyzer code.  

run Analyzer_V18.py protein_lfq.txt

Analyzer returns a file listing all the peptides that can be confidently quantified, and provides their effect-size, P-value, refolded CV, proteinase K site (if half-tryptic), and associated protein metadata as .txt files. 

## Citation

## Funding
