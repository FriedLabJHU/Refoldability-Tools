# Refoldability-Tools
Tools to analyze proteome-wide protein folding experimental data from mass spectrometry data. This python script was written to analyze raw normalized extracted ion intensity data from Label Free Quantification (LFQ) from detected peptides in refolded and native LiP-MS samples. 

## Citation
**Nonrefoldability is Pervasive Across the E. coli Proteome**  
Philip To, Briana Whitehead, Haley E. Tarbox, and Stephen D. Fried  
*Journal of the American Chemical Society* **2021** *143* (30), 11435-11448  
DOI: 10.1021/jacs.1c03270  

## Funding
* HFSP Research Grant (RGY0074/2019)  
* NSF Career Grant (MCB-2045844)  
* NIH Training Grant - Program in Molecular Biophysics (T32GM008403)  
* NIH Training Grant - Chemistry-Biology Interface (T32GM080189)  

## Legacy Files  
All code used for the analysis of materials in our [paper](https://www.pubs.acs.org/doi/abs/10.1021/jacs.1c03270) can be found in the **Legacy_Analyzer** directory  

### (Legacy) Running the Analyzer
Analyzer code can be run with a stock install of anaconda, which should contain all the required dependencies. Raw ion intensity data need to be first exported from the .pdResult file to a three-level hierarchy (protein > peptide group > consensus feature) to excel and then saved as a .txt file. Analyzer code and associated meta-data files should be moved to the same directory as raw ion intensity data .txt files to run Analyzer code. 

### (Legacy) Input and Output for the Analyzer 
Once in the correct directory, in ipython, the following command below can be used to run the analyzer code. The exported raw ion intensity data .txt file is inputted as the first argument after the Analyzer code.  

run Analyzer_V18.py protein_lfq.txt

Analyzer returns a file listing all the peptides that can be confidently quantified, and provides their effect-size, P-value, refolded CV, proteinase K site (if half-tryptic), and associated protein metadata as .txt files. 
