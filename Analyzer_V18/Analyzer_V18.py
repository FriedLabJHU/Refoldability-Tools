import numpy as np
import sys
import re
from scipy.stats import ttest_ind
from scipy.stats import combine_pvalues
from scipy.stats import variation
from scipy.stats import chi2
from scipy.stats import rankdata
import pandas as pd
import ast

def isclose(a, b, rel_tol=1e-05, abs_tol=0.0):
    return abs(a-b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)

#Call this program as python analyzer.py proteinControl_file_name.txt peptideLip_file_name.txt numConditions
#Input files should be tab-delimited files, exported from Excel (BE SURE TO CHANGE TO UTF-8/UNIX CARRIAGE RETURNS)

proteinControl_fn = sys.argv[1]
peptideLip_fn = sys.argv[1]
numConditions = 2

proteinControl_output_fn = peptideLip_fn.split('.')[0] + '_ctrl_out18_10+.txt'
peptideLip_output_fn = peptideLip_fn.split('.')[0] + '_out18_10+.txt'

proteinControl_file = open( proteinControl_fn , 'r' ) 
peptideLip_file = open( peptideLip_fn, 'r' )
proteinControl_output_file = open( proteinControl_output_fn, 'w' )
peptideLip_output_file = open( peptideLip_output_fn, 'w' )

#Reads data from the metadata file and creates a working dictionary
metadata_file = open('20200705_Protein_List_Updated_pi_final.txt' , 'r')
next(metadata_file) #skip the header line
metadata_dict = {}
for line in metadata_file:
    ProteinGeneName = line.split('\t')[0] #extract the name from the metadata file
    if line.split('\t')[2] == 'None': #extract potential synonyms from the metadata file, arrange in a list
        ProteinSynonyms = []
    elif ';' not in line.split('\t')[2]:
        ProteinSynonyms = [ line.split('\t')[2].strip('\"') ]
    else:
        ProteinSynonyms = line.split('\t')[2].strip('\"').split('; ')

    if line.split('\t')[3] == 'Unknown': #extract the locations from the metadata file
        ProteinLocations = 'Unknown'
    elif ',' not in line.split('\t')[3]:
        ProteinLocations = [ line.split('\t')[3].strip('\"') ]
    else:
        ProteinLocations =  line.split('\t')[3].strip('\"').split(', ') 

    if line.split('\t')[4].lower() == 'monomer': #extract complexes info from metadata file, arrange in a list
        ProteinComplex = 'Monomer'
    else:
        ProteinComplex = line.split('\t')[4].strip('\'\"][').split(', ')
        ProteinComplex = [x.strip('\'') for x in ProteinComplex]

    if line.split('\t')[5] == 'None': #extract the cofactors from the metadata file
        ProteinCofactors = 'None'
    elif ',' not in line.split('\t')[5]:
        ProteinCofactors = [ line.split('\t')[5].strip('\'\"') ]
        ProteinCofactors = [x.strip() for x in ProteinCofactors]
    else:
        ProteinCofactors = line.split('\t')[5].strip('\'\"').split(', ')
        ProteinCofactors = [x.strip() for x in ProteinCofactors]

    if line.split('\t')[6].startswith('Growth'): #extract essentiality from the metadata file
        Essential = 'N'
    elif line.split('\t')[6].startswith('No Growth'):
        Essential = 'Y'
    else:
        Essential = '?'

    if line.split('\t')[7].startswith('No Info Available'): #extracts copy number information
        CopyNumber = '?'
    else:
        CopyNumber = int(line.split('\t')[7])

    if line.split('\t')[8].startswith('No Info Available'): #extracts molecular weight
        MW = '?'
    else:
        MW = float( line.split('\t')[8] ) 

    if line.split('\t')[9].startswith('No Info Available'): #extracts abundance (femtograms)
        fg = '?'
    else:
        fg = float(line.split('\t')[9])

    numsubunits = int(line.split('\t')[10])
    uniquesubunits = int(line.split('\t')[11])

    if line.split('\t')[12].startswith('No Info Available'): #extracts pI
        pI = '?'
    else:
        pI = float(line.split('\t')[12])

    metadata_dict.update({ProteinGeneName:[ProteinSynonyms, ProteinLocations, ProteinComplex,ProteinCofactors,Essential,CopyNumber,MW,fg,numsubunits,uniquesubunits,pI]})
metadata_file.close()

#Reads domain information from a curated supfam file and creates a working dictionary
domain_dict = {}
domain_info_file = open('SuperFamilyListComplete_pt_1.txt', 'r')
next(domain_info_file)
foldList = []
aa_rangeList = []
counter = 0
for line in domain_info_file:
    info_in_line = line.split('\t')
    if counter == 1:
        if re.match(r'[a-z][a-z][a-z]', info_in_line[0]): #Starting a new protein
            totalFolds = len(foldList)
            foldCount = len(set(foldList))
            numFolds = "%s" %str(foldCount)
            totFolds = "%s" %str(totalFolds)
            domain_dict[ genename ].append( totFolds )
            domain_dict[ genename ].append( numFolds )
            foldList = []
            genename = info_in_line[0]
            fold = info_in_line[4][:-1]
            aa_range = re.sub( "," , "+" , info_in_line[1] )
            foldList.append(fold)
            aa_rangeList.append(aa_range)
            domain_dict[ genename ] = [[],[]]
            domain_dict[ genename ][0].append( fold )
            domain_dict[ genename ][1].append( aa_range )
        else:
            fold = info_in_line[4][:-1]
            aa_range = re.sub( "," , "+" , info_in_line[1] )
            foldList.append(fold)
            aa_rangeList.append(aa_range)
            domain_dict[ genename ][0].append( fold )
            domain_dict[ genename ][1].append( aa_range )
    else:
        genename = info_in_line[0]
        fold = info_in_line[4][:-1]
        aa_range = re.sub( "," , "+" , info_in_line[1] )
        foldList.append(fold)
        aa_rangeList.append(aa_range)
        domain_dict[ genename ] = [[],[]]
        domain_dict[ genename ][0].append( fold )
        domain_dict[ genename ][1].append( aa_range )
        counter += 1
totalFolds = len(foldList)
foldCount = len(set(foldList))
numFolds = foldCount
totFolds = totalFolds
domain_dict[ genename ].append( totFolds )
domain_dict[ genename ].append( numFolds )
domain_info_file.close()

#Put the domain data into the metadata dictionary
for key in domain_dict.keys():
    if key in metadata_dict.keys():
        metadata_dict[key].append( domain_dict[key][:-2])
        metadata_dict[key].append( domain_dict[key][-2])
        metadata_dict[key].append( domain_dict[key][-1])
    else:
        for key2 in metadata_dict.keys():
            if key in metadata_dict[key2][0]:
                metadata_dict[key2].append( domain_dict[key][:-2])
                metadata_dict[key2].append( domain_dict[key][-2])
                metadata_dict[key2].append( domain_dict[key][-1])
                break
for key in metadata_dict.keys():
    if len(metadata_dict[key]) != 14:
        metadata_dict[key].append( 'Unknown' )
        metadata_dict[key].append( 0 )
        metadata_dict[key].append( 0 )

#Loop through Protein Control txt file line by line 
proteinControl_output_file.write( 'GeneName\tAccession\tLocations\tComposition\tCofactors\tEssential?\tCopy Number\tMolecular Weight\tAmount(fg)\tNum Subunits\tUnique Subunits\tpI\tFolds\tTotal Domains\tUnique Domains\tProteinRatio\tProteinPValue\n' )
for line in proteinControl_file:
    info_in_line = line.split('\t')
    if info_in_line[0] == 'FALSE':
        if info_in_line[1] == 'Medium' or info_in_line[1] == 'High':
            Accession = info_in_line[3]
            ProteinGeneName = info_in_line[4].split('=')[3].split(' ')[0] #Takes Protein Gene Name from Description line 
            if info_in_line[26] != '' and info_in_line[27] != '': #If quans are available
                ProteinRatio = float(0.0) 
                ProteinPValue = float(0.0)
            else: #if protein_ratio and protein_p_value was blank 
                ProteinRatio = float(0.0) 
                ProteinPValue = float(0.0)

            if ProteinGeneName in metadata_dict.keys():
                metadata_dict[ProteinGeneName].append( Accession )
                metadata_dict[ProteinGeneName].append( ProteinRatio )
                metadata_dict[ProteinGeneName].append( ProteinPValue )
                proteinControl_output_file.write( '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' %(ProteinGeneName, metadata_dict[ProteinGeneName][14], metadata_dict[ProteinGeneName][1], metadata_dict[ProteinGeneName][2], metadata_dict[ProteinGeneName][3], metadata_dict[ProteinGeneName][4], metadata_dict[ProteinGeneName][5], metadata_dict[ProteinGeneName][6], metadata_dict[ProteinGeneName][7], metadata_dict[ProteinGeneName][8],metadata_dict[ProteinGeneName][9], metadata_dict[ProteinGeneName][10], metadata_dict[ProteinGeneName][11][0], metadata_dict[ProteinGeneName][12], metadata_dict[ProteinGeneName][13], metadata_dict[ProteinGeneName][15], metadata_dict[ProteinGeneName][16] ) )
            else:
                for key2 in metadata_dict.keys():
                    if ProteinGeneName in metadata_dict[key2][0]:
                        metadata_dict[key2].append( Accession )
                        metadata_dict[key2].append( ProteinRatio )
                        metadata_dict[key2].append( ProteinPValue )
                        proteinControl_output_file.write( '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' %(ProteinGeneName, metadata_dict[key2][14], metadata_dict[key2][1], metadata_dict[key2][2], metadata_dict[key2][3], metadata_dict[key2][4], metadata_dict[key2][5], metadata_dict[key2][6], metadata_dict[key2][7], metadata_dict[key2][8], metadata_dict[key2][9], metadata_dict[key2][10], metadata_dict[key2][11][0], metadata_dict[key2][12], metadata_dict[key2][13], metadata_dict[key2][15], metadata_dict[key2][16] ) )
                        break
proteinControl_file.close()
proteinControl_output_file.close()

#Create data structures important for looping through the Peptide LiP file (with consensus features)
Protein_info = []
Peptides_table = [] #Create Empty Peptide Data Table from Peptides LiP File
quanArrayR = [[] for i in range(numConditions-1)]
quanArrayN = [[] for i in range(numConditions-1)]
quanArrayRfilt = [[] for i in range(numConditions-1)]
quanArrayNfilt = [[] for i in range(numConditions-1)]
Pep_Seq = ''
charge_states_seen = []

#Proteins_with_halftryptics = []
peptideLip_output_file.write( 'Accession\tGeneName\tPeptide Sequence\tproteinaseKsite\tLocation\tComposition\tCofactors\tEssential?\tCopy Number\tMolecular Weight\tAmount (fg)\tNum Subunits\tUnique Subunits\tpI\tFolds\tDomain Ranges\tNum Domains\tUnique Domains\tDomain Location\tPosition In Domain\t' + '\t'.join('PeptideRatio'+str(s) for s in range(1,numConditions)) + '\t' + '\t'.join('NormPeptideRatio'+str(s) for s in range(1,numConditions)) + '\t' + '\t'.join('PeptidePValue'+str(s) for s in range(1,numConditions)) + '\t' + '\t'.join('RefoldedCoeffVar'+str(s) for s in range(1,numConditions))+ '\n' )

#Loop through Peptide LiP File and append AbundanceRatio Data to it 
for line in peptideLip_file:
    info_in_line = line.split('\t')
    if info_in_line[0] == 'FALSE': #Checks if this line is a protein
        Pep_Seq_thisLine = info_in_line[3]
        #Process the quanArray for the previous peptide group and write it out.  skip on the first time (when Pep_seq is not yet defined)
        if Pep_Seq and quanArrayR and Pep_Seq_thisLine != Pep_Seq: #Pep_Seq_thisLine should equal to Pep_Seq unless it is a new peptide sequence that is starting 
            #what to keep and what to reject (on the basis of missing data)
            for c in range(numConditions-1):
                for j in range(len(quanArrayN[c])):
                    if quanArrayN[c][j][1].count( 0 ) == 0 and quanArrayR[c][j][1].count( 0 ) == 0:
                        quanArrayNfilt[c].append( quanArrayN[c][j] )
                        quanArrayRfilt[c].append( quanArrayR[c][j] )
                    elif quanArrayN[c][j][1].count( 0 ) == 3 and quanArrayR[c][j][1].count( 0 ) == 0: #Missing data being used to impute low value for Ns
                        quanArrayNfilt[c].append( [quanArrayN[c][j][0], [1000,1000,1000]] )
                        quanArrayRfilt[c].append( quanArrayR[c][j] )
                    elif quanArrayN[c][j][1].count( 0 ) == 0 and quanArrayR[c][j][1].count( 0 ) == 3: #Missing data being used to impute low value for Rs
                        quanArrayNfilt[c].append( quanArrayN[c][j] )
                        quanArrayRfilt[c].append( [quanArrayR[c][j][0], [1000,1000,1000]] )
                    elif (quanArrayN[c][j][1].count( 0 ) + quanArrayR[c][j][1].count( 0 )) == 1:
                        if quanArrayN[c][j][1].count( 0 ) == 1:
                            quanArrayN[c][j][1].remove( 0 )
                            quanArrayNfilt[c].append( quanArrayN[c][j] )
                            quanArrayRfilt[c].append( quanArrayR[c][j] )
                        if quanArrayR[c][j][1].count( 0 ) == 1:
                            quanArrayR[c][j][1].remove( 0 )
                            quanArrayNfilt[c].append( quanArrayN[c][j] )
                            quanArrayRfilt[c].append( quanArrayR[c][j] )                        
                    else:
                        pass
            PeptideRatioByCondition = []
            VariationByCondition = []
            normalizedPeptideRatioByCondition = []
            PeptidePValueByCondition = []
            for c in range(numConditions-1):
                if quanArrayRfilt[c]:
                    ratioArray = [ np.log2(np.mean( quanArrayRfilt[c][i][1] )/np.mean( quanArrayNfilt[c][i][1] )) for i in range(len(quanArrayRfilt[c])) ]
                    variationArray = [ variation(quanArrayRfilt[c][i][1] ) for i in range(len(quanArrayRfilt[c])) ]
                    pvalueArray = [ ttest_ind( quanArrayRfilt[c][i][1], quanArrayNfilt[c][i][1] , equal_var=False )[1] for i in range(len(quanArrayRfilt[c])) ]
                    teststatArray = [ ttest_ind( quanArrayRfilt[c][i][1], quanArrayNfilt[c][i][1] , equal_var=False )[0] for i in range(len(quanArrayRfilt[c])) ]
                    PeptideRatioByCondition.append( np.median( ratioArray ) )#use the median from all the concensus features for this peptide group
                    VariationByCondition.append( np.min( variationArray ) )
                    if ProteinPValue > 2 and np.abs(ProteinRatio) > 1: #Is the ProteinRatio Significant???
                        normalizedPeptideRatioByCondition.append( PeptideRatioByCondition[c] - ProteinRatio )
                    else:
                        normalizedPeptideRatioByCondition.append( PeptideRatioByCondition[c] )
                
                    if len(pvalueArray) == 1:
                        PeptidePValueByCondition.append( np.abs(np.log10( pvalueArray[0] )) )
                    else: 
                        if all(z > 0 for z in teststatArray) or all(z < 0 for z in teststatArray): #Fisher's Method
                            ChiSquareTeststat = 2*np.abs( np.sum( np.log(pvalueArray)* np.sign(teststatArray) ) )
                            PeptidePValueByCondition.append( np.abs( np.log10 ( 1 - chi2.cdf(ChiSquareTeststat,len(ratioArray)*2) )) )
                        else: #if the different consensus features are inconsistent with each other
                            PeptidePValueByCondition.append( 0 )
                else:
                    PeptideRatioByCondition.append( 0 )
                    VariationByCondition.append( 0 )
                    normalizedPeptideRatioByCondition.append( 0 )
                    PeptidePValueByCondition.append( 0 )
            Peptides_table.append( [Accession, ProteinGeneName, Pep_Seq, proteinaseKsite, ProteinLocation, ProteinComposition, ProteinCofactors, Essential, CopyNumber, ProteinMW, fg, numsubunits, uniquesubunits, pI, ProteinFolds, ProteinDomainRanges, ProteinNumDomains, ProteinNumUniqueDomains, DomainLocation, PositionInDomain, PeptideRatioByCondition, normalizedPeptideRatioByCondition, PeptidePValueByCondition, VariationByCondition] )
            quanArrayR = [[] for i in range(numConditions-1)]
            quanArrayN = [[] for i in range(numConditions-1)]
            quanArrayRfilt = [[] for i in range(numConditions-1)]
            quanArrayNfilt = [[] for i in range(numConditions-1)]
            charge_states_seen = []
            Pep_Seq = ''    

        ProteinGeneName = info_in_line[4].split('=')[3].split(' ')[0]
        if ProteinGeneName in metadata_dict.keys():
            ProteinLocation = metadata_dict[ProteinGeneName][1]
            ProteinComposition = metadata_dict[ProteinGeneName][2]
            ProteinCofactors = metadata_dict[ProteinGeneName][3]
            Essential = metadata_dict[ProteinGeneName][4]
            CopyNumber = metadata_dict[ProteinGeneName][5]
            ProteinMW = metadata_dict[ProteinGeneName][6]
            fg = metadata_dict[ProteinGeneName][7]
            numsubunits = metadata_dict[ProteinGeneName][8]
            uniquesubunits = metadata_dict[ProteinGeneName][9]
            pI = metadata_dict[ProteinGeneName][10]
            if metadata_dict[ProteinGeneName][11] != 'Unknown':
                ProteinFolds = metadata_dict[ProteinGeneName][11][0]
                ProteinDomainRanges = metadata_dict[ProteinGeneName][11][1]
                ProteinDomainLists = []
                for entry in metadata_dict[ProteinGeneName][11][1]: 
                    if '+' not in entry: 
                        start=int(entry.split('-')[0]) 
                        end=int(entry.split('-')[1]) 
                        dr = [x for x in range(start,end+1)] 
                        ProteinDomainLists.append( dr ) 
                    else:
                        dr = []
                        componentranges = entry.split('+')
                        for componentrange in componentranges:
                            start=int(componentrange.split('-')[0]) 
                            end=int(componentrange.split('-')[1]) 
                            dr.append( [x for x in range(start,end+1)] )
                        ProteinDomainLists.append( [item for sublist in dr for item in sublist] )
            else:
                ProteinFolds = 'Unknown'
                ProteinDomainRanges = 'N/A'
            ProteinNumDomains = metadata_dict[ProteinGeneName][12]
            ProteinNumUniqueDomains = metadata_dict[ProteinGeneName][13]
        else:
            for key2 in metadata_dict.keys():
                if ProteinGeneName in metadata_dict[key2][0]:
                    ProteinLocation = metadata_dict[key2][1]
                    ProteinComposition = metadata_dict[key2][2]
                    ProteinCofactors = metadata_dict[key2][3]
                    Essential = metadata_dict[key2][4]
                    CopyNumber = metadata_dict[key2][5]
                    ProteinMW = metadata_dict[key2][6]
                    fg = metadata_dict[key2][7]
                    numsubunits = metadata_dict[key2][8]
                    uniquesubunits = metadata_dict[key2][9]
                    pI = metadata_dict[key2][10]
                    if metadata_dict[key2][11] != 'Unknown':
                        ProteinFolds = metadata_dict[key2][11][0]
                        ProteinDomainRanges = metadata_dict[key2][11][1]
                        ProteinDomainLists = []
                        for entry in metadata_dict[key2][11][1]: 
                            if '+' not in entry: 
                                start=int(entry.split('-')[0]) 
                                end=int(entry.split('-')[1]) 
                                dr = [x for x in range(start,end+1)] 
                                ProteinDomainLists.append( dr ) 
                            else:
                                dr = []
                                componentranges = entry.split('+')
                                for componentrange in componentranges:
                                    start=int(componentrange.split('-')[0]) 
                                    end=int(componentrange.split('-')[1]) 
                                    dr.append( [x for x in range(start,end+1)] )
                                ProteinDomainLists.append( [item for sublist in dr for item in sublist] )
                    else:
                        ProteinFolds = 'Unknown'
                        ProteinDomainRanges = 'N/A'
                    ProteinNumDomains = metadata_dict[key2][12]
                    ProteinNumUniqueDomains = metadata_dict[key2][13]
                    break

    if info_in_line[1] == 'FALSE': #Checks if this line is a peptide
        Pep_Seq_thisLine = info_in_line[3]
        if '1xOxidation' in info_in_line[4]:
            MetOxidation = 1
        elif '2xOxidation' in info_in_line[4]:
            MetOxidation = 2
        else:
            MetOxidation = 0

        if '1xLabel' in info_in_line[4]:
            HeavyResidue = 1 
        elif '2xLabel' in info_in_line[4]:
            HeavyResidue = 2 
        elif '3xLabel' in info_in_line[4]: 
            HeavyResidue = 3 
        else:   
            HeavyResidue = 0 
    
        Missing_cleavages = int(info_in_line[11])

        #Process the quanArray for the previous peptide group and write it out.  skip on the first time (when Pep_seq is not yet defined)
        if Pep_Seq and quanArrayR and Pep_Seq_thisLine != Pep_Seq:
            #what to keep and what to reject (on the basis of missing data)
            for c in range(numConditions-1):
                for j in range(len(quanArrayN[c])):
                    if quanArrayN[c][j][1].count( 0 ) == 0 and quanArrayR[c][j][1].count( 0 ) == 0:
                        quanArrayNfilt[c].append( quanArrayN[c][j] )
                        quanArrayRfilt[c].append( quanArrayR[c][j] )
                    elif quanArrayN[c][j][1].count( 0 ) == 3 and quanArrayR[c][j][1].count( 0 ) == 0: #Missing data being used to impute low value for Ns
                        quanArrayNfilt[c].append( [quanArrayN[c][j][0], [1000,1000,1000]] )
                        quanArrayRfilt[c].append( quanArrayR[c][j] )
                    elif quanArrayN[c][j][1].count( 0 ) == 0 and quanArrayR[c][j][1].count( 0 ) == 3: #Missing data being used to impute low value for Rs
                        quanArrayNfilt[c].append( quanArrayN[c][j] )
                        quanArrayRfilt[c].append( [quanArrayR[c][j][0], [1000,1000,1000]] )
                    elif (quanArrayN[c][j][1].count( 0 ) + quanArrayR[c][j][1].count( 0 )) == 1:
                        if quanArrayN[c][j][1].count( 0 ) == 1:
                            quanArrayN[c][j][1].remove( 0 )
                            quanArrayNfilt[c].append( quanArrayN[c][j] )
                            quanArrayRfilt[c].append( quanArrayR[c][j] )
                        if quanArrayR[c][j][1].count( 0 ) == 1:
                            quanArrayR[c][j][1].remove( 0 )
                            quanArrayNfilt[c].append( quanArrayN[c][j] )
                            quanArrayRfilt[c].append( quanArrayR[c][j] )                        
                    else:
                        pass
            PeptideRatioByCondition = []
            VariationByCondition = []
            normalizedPeptideRatioByCondition = []
            PeptidePValueByCondition = []
            for c in range(numConditions-1):
                if quanArrayRfilt[c]:
                    ratioArray = [ np.log2(np.mean( quanArrayRfilt[c][i][1] )/np.mean( quanArrayNfilt[c][i][1] )) for i in range(len(quanArrayRfilt[c])) ]
                    variationArray = [ variation(quanArrayRfilt[c][i][1] ) for i in range(len(quanArrayRfilt[c])) ]
                    pvalueArray = [ ttest_ind( quanArrayRfilt[c][i][1], quanArrayNfilt[c][i][1] , equal_var=False )[1] for i in range(len(quanArrayRfilt[c])) ]
                    teststatArray = [ ttest_ind( quanArrayRfilt[c][i][1], quanArrayNfilt[c][i][1] , equal_var=False )[0] for i in range(len(quanArrayRfilt[c])) ]
                    PeptideRatioByCondition.append( np.median( ratioArray ) )#use the median from all the concensus features for this peptide group
                    VariationByCondition.append( np.min( variationArray ) )
                    if ProteinPValue > 2 and np.abs(ProteinRatio) > 1: #Is the ProteinRatio Significant???
                        normalizedPeptideRatioByCondition.append( PeptideRatioByCondition[c] - ProteinRatio )
                    else:
                        normalizedPeptideRatioByCondition.append( PeptideRatioByCondition[c] )
                
                    if len(pvalueArray) == 1:
                        PeptidePValueByCondition.append( np.abs(np.log10( pvalueArray[0] )) )
                    else: 
                        if all(z > 0 for z in teststatArray) or all(z < 0 for z in teststatArray): #Fisher's Method
                            ChiSquareTeststat = 2*np.abs( np.sum( np.log(pvalueArray)* np.sign(teststatArray) ) )
                            PeptidePValueByCondition.append( np.abs( np.log10 ( 1 - chi2.cdf(ChiSquareTeststat,len(ratioArray)*2) )) )
                        else: #if the different consensus features are inconsistent with each other
                            PeptidePValueByCondition.append( 0 )
                else:
                    PeptideRatioByCondition.append( 0 )
                    VariationByCondition.append( 0 )
                    normalizedPeptideRatioByCondition.append( 0 )
                    PeptidePValueByCondition.append( 0 )
            Peptides_table.append( [Accession, ProteinGeneName, Pep_Seq, proteinaseKsite, ProteinLocation, ProteinComposition, ProteinCofactors, Essential, CopyNumber, ProteinMW, fg, numsubunits, uniquesubunits, pI, ProteinFolds, ProteinDomainRanges, ProteinNumDomains, ProteinNumUniqueDomains, DomainLocation, PositionInDomain, PeptideRatioByCondition, normalizedPeptideRatioByCondition, PeptidePValueByCondition, VariationByCondition] )
            quanArrayR = [[] for i in range(numConditions-1)]
            quanArrayN = [[] for i in range(numConditions-1)]
            quanArrayRfilt = [[] for i in range(numConditions-1)]
            quanArrayNfilt = [[] for i in range(numConditions-1)]
            charge_states_seen = []
            Pep_Seq = ''    
        #Start up the next peptide group
        if info_in_line[2] == 'Medium' or info_in_line[2] == 'High':
            if not re.search(';',info_in_line[8]) : #Checks that the peptide maps to a single protein
                if (info_in_line[13 + (numConditions-1)*3 + numConditions*2]) == '': #Checks that there is quan information (i.e., not 'No Quan Value' or 'Not Reliable')
                    Accession = info_in_line[8]
                    FoundProteinRecord = False
                    for key in metadata_dict.keys(): #Get Protein Info if the protein exists in the dictionary (i.e., it was found in the Control sample)
                        if len(metadata_dict[key]) == 17:
                            if Accession == metadata_dict[key][14]:
                                ProteinRatio = metadata_dict[key][15]
                                ProteinPValue = metadata_dict[key][16]
                                FoundProteinRecord = True
                                break
                    if not FoundProteinRecord: #if protein was not identifid in the control sample, then generate filler values
                        ProteinRatio = float(1.0)
                        ProteinPValue = float(0.0)
                    # Check if it's half-tryptic, and if so, locate the site of the proteinase K cut
                    Pep_Seq = info_in_line[3]
                    points = 0
                    if (Pep_Seq.split('.')[1][-1] == 'R' or Pep_Seq.split('.')[1][-1] == 'K') and (Pep_Seq.split('.')[2][1] != '-'): 
                        points += 1
                        proteinaseKsite = 'Front'
                    if Pep_Seq.split('.')[0][1] == 'R' or Pep_Seq.split('.')[0][1] == 'K' or Pep_Seq.split('.')[0][1] == '-':
                        points += 1 
                        proteinaseKsite = 'Back'
                    if Pep_Seq.split('.')[2][1] == '-':
                        points += 1
                        proteinaseKsite = 'Front'
                    if points == 1: #This is the marker of whether it is a 'half-tryptic' peptide
                        PositionInProtein = re.sub ( '[^0-9\[\]-]' ,'' , info_in_line[9].split(' ')[1] )
                        PositionInProtein = PositionInProtein.split('-')
                        PositionInProtein[0] = PositionInProtein[0][1:]
                        PositionInProtein[1] = PositionInProtein[1][:-1]
                        if proteinaseKsite == 'Front':
                            proteinaseKsite = Pep_Seq.split('.')[1][0]+PositionInProtein[0]
                            #Find which domain this site is in.
                            if ProteinFolds != 'Unknown':
                                numDomainsFoundIn = 0
                                for (ProteinDomainList,ProteinDomainRange) in zip(ProteinDomainLists,ProteinDomainRanges):
                                    if int(PositionInProtein[0]) in ProteinDomainList:
                                        DomainLocation = ProteinGeneName + ProteinDomainRange
                                        numDomainsFoundIn += 1
                                        PositionInDomain = ProteinDomainList.index(int(PositionInProtein[0])) + 1
                                if numDomainsFoundIn != 1:
                                    DomainLocation = ProteinGeneName + 'LNK'
                                    PositionInDomain = 'N/A'
                            else:
                                DomainLocation = ProteinGeneName
                                PositionInDomain = 'N/A'
                        if proteinaseKsite == 'Back':
                            proteinaseKsite = Pep_Seq.split('.')[1][-1]+PositionInProtein[1]
                            #Find which domain this site is in.
                            if ProteinFolds != 'Unknown':
                                numDomainsFoundIn = 0
                                for (ProteinDomainList,ProteinDomainRange) in zip(ProteinDomainLists,ProteinDomainRanges):
                                    if int(PositionInProtein[1]) in ProteinDomainList:
                                        DomainLocation = ProteinGeneName + ProteinDomainRange
                                        numDomainsFoundIn += 1
                                        PositionInDomain = ProteinDomainList.index(int(PositionInProtein[1])) + 1
                                if numDomainsFoundIn != 1:
                                    DomainLocation = ProteinGeneName + 'LNK'
                                    PositionInDomain = 'N/A'
                            else:
                               DomainLocation = ProteinGeneName 
                               PositionInDomain = 'N/A'
                    else: #This is NOT a half-tryptic peptide
                        proteinaseKsite = info_in_line[9].split(' ')[1]
                        trypticstart = int(  re.sub( '[^0-9]' , '' ,proteinaseKsite.split('-')[0] ) )
                        trypticend = int(  re.sub( '[^0-9]' , '' ,proteinaseKsite.split('-')[1] ) )
                        trypticList = [x for x in range(trypticstart,trypticend+1)] 
                        if ProteinFolds != 'Unknown':
                            numDomainsFoundIn = 0
                            for (ProteinDomainList,ProteinDomainRange) in zip(ProteinDomainLists,ProteinDomainRanges):
                                if all(site in ProteinDomainList for site in trypticList):
                                    DomainLocation = ProteinGeneName + ProteinDomainRange
                                    numDomainsFoundIn += 1
                                    PositionInDomain = 'N/A'
                            if numDomainsFoundIn != 1:
                                DomainLocation = ProteinGeneName + 'LNK'
                                PositionInDomain = 'N/A'
                        else:
                            DomainLocation = ProteinGeneName
                            PositionInDomain = 'N/A' 
    if info_in_line[2] == 'FALSE' and Pep_Seq: #Checks if this line is a consensus feature
        massPerCharge = float(info_in_line[3])
        if MetOxidation == 1:
            charge_state = info_in_line[4] + 'O'
        elif MetOxidation == 2:
            charge_state = info_in_line[4] + 'O2'
        else:
            charge_state = info_in_line[4]
        
        if HeavyResidue == 1:
            charge_state = charge_state + 'H'
        elif HeavyResidue == 2: 
            charge_state = charge_state + 'H2'
        elif HeavyResidue == 3: 
            charge_state = charge_state + 'H3'
        else:
            pass


        averageRT = float(info_in_line[6])
        PSMs = int(info_in_line[7])
        
        currentConcensusFeature = [z for z in info_in_line[12:12+numConditions*3]]
        currentConcensusFeatureFilled = []
        for z in currentConcensusFeature:
            try:
                currentConcensusFeatureFilled.append( float (z) )
            except ValueError:
                currentConcensusFeatureFilled.append( float(0) )
        if HeavyResidue == Missing_cleavages + 1 or HeavyResidue == 0: 
        #if info_in_line[8] == "1": 
            if charge_state not in charge_states_seen: 
                for c in range(numConditions-1):
                    quanArrayN[c].append( [ [massPerCharge, charge_state, averageRT, PSMs] , currentConcensusFeatureFilled[0:3] ] )
                    quanArrayR[c].append( [ [massPerCharge, charge_state, averageRT, PSMs] , currentConcensusFeatureFilled[3*c+3:3*c+6] ] )
                charge_states_seen.append(charge_state)
            else:
                for c in range(numConditions-1):
                    for j in range(len(quanArrayN[c])):
                        if quanArrayN[c][j][0][1] == charge_state:
                            quanArrayN[c][j][1][0] = quanArrayN[c][j][1][0] + currentConcensusFeatureFilled[0]
                            quanArrayN[c][j][1][1] = quanArrayN[c][j][1][1] + currentConcensusFeatureFilled[1]
                            quanArrayN[c][j][1][2] = quanArrayN[c][j][1][2] + currentConcensusFeatureFilled[2]
                            quanArrayN[c][j][0][3] += PSMs
                
                    for j in range(len(quanArrayR[c])):
                        if quanArrayR[c][j][0][1] == charge_state:
                            quanArrayR[c][j][1][0] = quanArrayR[c][j][1][0] + currentConcensusFeatureFilled[3+3*c]
                            quanArrayR[c][j][1][1] = quanArrayR[c][j][1][1] + currentConcensusFeatureFilled[4+3*c]
                            quanArrayR[c][j][1][2] = quanArrayR[c][j][1][2] + currentConcensusFeatureFilled[5+3*c]
                            quanArrayR[c][j][0][3] += PSMs


if Pep_Seq and quanArrayR and Pep_Seq_thisLine != Pep_Seq:
    #what to keep and what to reject (on the basis of missing data)
    for c in range(numConditions-1):
        for j in range(len(quanArrayN[c])):
            if quanArrayN[c][j][1].count( 0 ) == 0 and quanArrayR[c][j][1].count( 0 ) == 0:
                quanArrayNfilt[c].append( quanArrayN[c][j] )
                quanArrayRfilt[c].append( quanArrayR[c][j] )
            elif quanArrayN[c][j][1].count( 0 ) == 3 and quanArrayR[c][j][1].count( 0 ) == 0: #Missing data being used to impute low value for Ns
                quanArrayNfilt[c].append( [quanArrayN[c][j][0], [1000,1000,1000]] )
                quanArrayRfilt[c].append( quanArrayR[c][j] )
            elif quanArrayN[c][j][1].count( 0 ) == 0 and quanArrayR[c][j][1].count( 0 ) == 3: #Missing data being used to impute low value for Rs
                quanArrayNfilt[c].append( quanArrayN[c][j] )
                quanArrayRfilt[c].append( [quanArrayR[c][j][0], [1000,1000,1000]] )
            elif (quanArrayN[c][j][1].count( 0 ) + quanArrayR[c][j][1].count( 0 )) == 1:
                if quanArrayN[c][j][1].count( 0 ) == 1:
                    quanArrayN[c][j][1].remove( 0 )
                    quanArrayNfilt[c].append( quanArrayN[c][j] )
                    quanArrayRfilt[c].append( quanArrayR[c][j] )
                if quanArrayR[c][j][1].count( 0 ) == 1:
                    quanArrayR[c][j][1].remove( 0 )
                    quanArrayNfilt[c].append( quanArrayN[c][j] )
                    quanArrayRfilt[c].append( quanArrayR[c][j] )                        
            else:
                pass
    PeptideRatioByCondition = []
    VariationByCondition = []
    normalizedPeptideRatioByCondition = []
    PeptidePValueByCondition = []
    for c in range(numConditions-1):
        if quanArrayRfilt[c]:
            ratioArray = [ np.log2(np.mean( quanArrayRfilt[c][i][1] )/np.mean( quanArrayNfilt[c][i][1] )) for i in range(len(quanArrayRfilt[c])) ]
            variationArray = [ variation(quanArrayRfilt[c][i][1] ) for i in range(len(quanArrayRfilt[c])) ]
            pvalueArray = [ ttest_ind( quanArrayRfilt[c][i][1], quanArrayNfilt[c][i][1] , equal_var=False )[1] for i in range(len(quanArrayRfilt[c])) ]
            teststatArray = [ ttest_ind( quanArrayRfilt[c][i][1], quanArrayNfilt[c][i][1] , equal_var=False )[0] for i in range(len(quanArrayRfilt[c])) ]
            PeptideRatioByCondition.append( np.median( ratioArray ) )#use the median from all the concensus features for this peptide group
            VariationByCondition.append( np.min( variationArray ) )
            if ProteinPValue > 2 and np.abs(ProteinRatio) > 1: #Is the ProteinRatio Significant???
                normalizedPeptideRatioByCondition.append( PeptideRatioByCondition[c] - ProteinRatio )
            else:
                normalizedPeptideRatioByCondition.append( PeptideRatioByCondition[c] )
        
            if len(pvalueArray) == 1:
                PeptidePValueByCondition.append( np.abs(np.log10( pvalueArray[0] )) )
            else: 
                if all(z > 0 for z in teststatArray) or all(z < 0 for z in teststatArray): #Fisher's Method
                    ChiSquareTeststat = 2*np.abs( np.sum( np.log(pvalueArray)* np.sign(teststatArray) ) )
                    PeptidePValueByCondition.append( np.abs( np.log10 ( 1 - chi2.cdf(ChiSquareTeststat,len(ratioArray)*2) )) )
                else: #if the different consensus features are inconsistent with each other
                    PeptidePValueByCondition.append( 0 )
        else:
            PeptideRatioByCondition.append( 0 )
            VariationByCondition.append( 0 )
            normalizedPeptideRatioByCondition.append( 0 )
            PeptidePValueByCondition.append( 0 )
    Peptides_table.append( [Accession, ProteinGeneName, Pep_Seq, proteinaseKsite, ProteinLocation, ProteinComposition, ProteinCofactors, Essential, CopyNumber, ProteinMW, fg, numsubunits, uniquesubunits, pI, ProteinFolds, ProteinDomainRanges, ProteinNumDomains, ProteinNumUniqueDomains, DomainLocation, PositionInDomain, PeptideRatioByCondition, normalizedPeptideRatioByCondition, PeptidePValueByCondition, VariationByCondition] )
    quanArrayR = [[] for i in range(numConditions-1)]
    quanArrayN = [[] for i in range(numConditions-1)]
    quanArrayRfilt = [[] for i in range(numConditions-1)]
    quanArrayNfilt = [[] for i in range(numConditions-1)]
    charge_states_seen = []
    Pep_Seq = ''    

#peptideLip_output_file.write('Number of proteins with half-tryptic petides:\t%s' %str(len(Proteins_with_halftryptics)))
for i in range(len(Peptides_table)):
    if np.count_nonzero( Peptides_table[i][22] ) != 0:   
        line_to_write = '\t'.join(str(s) for s in Peptides_table[i][0:20]) + '\t' + '\t'.join(str(s) for s in Peptides_table[i][20][:]) + '\t' + '\t'.join(str(s) for s in Peptides_table[i][21][:]) + '\t' + '\t'.join(str(s) for s in Peptides_table[i][22][:]) + '\t' + '\t'.join(str(s) for s in Peptides_table[i][23][:])+ '\n'
        peptideLip_output_file.write( line_to_write )

peptideLip_file.close()
peptideLip_output_file.close()

peptideLip_output_file = open( peptideLip_output_fn, 'r' )
protein_summary_fn = peptideLip_fn.split('.')[0] + '_summary18Protein.txt'
protein_summary_file = open( protein_summary_fn, 'w')

GeneName = ''
Pvalues = [[] for c in range(numConditions-1)]
Ratios = [[] for c in range(numConditions-1)]
next(peptideLip_output_file)

protein_summary_file.write( 'Accession\tGeneName\tLocation\tComposition\tCofactors\tEssential?\tCopy Number\tMolecular Weight\tAmount (fg)\tNum Subunits\tUnique Subunits\tpI\tFolds\tDomain Ranges\tNum Domains\tUnique Domains\t' + '\t'.join('SigPeptides'+str(s)+'\t'+'TotPeptides'+str(s)+'\t'+'AonPeptides'+str(s)+'\t'+'AllSigPeptides'+str(s) for s in range(1,numConditions)) + '\n' )
for line in peptideLip_output_file:
    info_in_line = line.split('\t')
    if info_in_line[1] == GeneName: #we're on the same gene
        for c in range(numConditions-1):
            Ratios[c].append( float(info_in_line[19 + numConditions + c] ) )
            Pvalues[c].append( float(info_in_line[18 + 2*numConditions + c] ) )
    else: #we've started looking at a new gene
        if GeneName == '': #We've just started; this was the first gene
            Accession = info_in_line[0]
            GeneName = info_in_line[1]
            Location = info_in_line[4]
            Composition = info_in_line[5]
            Cofactors = info_in_line[6]
            Essential = info_in_line[7]
            CopyNumber = info_in_line[8]
            ProteinMW = info_in_line[9]
            fg = info_in_line[10]
            numsubunits = info_in_line[11]
            uniquesubunits = info_in_line[12]
            pI = info_in_line[13]
            ProteinFolds = info_in_line[14]
            ProteinDomainRanges = info_in_line[15]
            ProteinNumDomains = info_in_line[16]
            ProteinNumUniqueDomains = info_in_line[17]
            DomainLocation = info_in_line[18]
            PositionInDomain = info_in_line[19]
            for c in range(numConditions-1):
                Ratios[c].append( float(info_in_line[19 + numConditions + c] ) )
                Pvalues[c].append( float(info_in_line[18 + 2*numConditions + c] ) )
        else: #We've just started looking at a new gene
            sigPeptidesByCondition = []
            totPeptidesByCondition = []
            aonPeptidesByCondition = [] #All or Nothing Peptides, those which are only present (or completely not present) in the refolded sample
            allsigPeptidesByCondition = [] 
            for c in range(numConditions-1):
                sigPeptides = 0
                totPeptides = 0
                aonPeptides = 0
                allsigPeptides = 0
                for (Ratio,Pval) in zip(Ratios[c],Pvalues[c]):
                    if Ratio != 0 and Pval != 0:
                        totPeptides += 1
                        if (6 > np.abs(Ratio) > 1 and Pval > 2):
                            sigPeptides += 1
                            allsigPeptides += 1
                        if np.abs(Ratio) > 6 and Pval > 1.8:
                            aonPeptides += 1
                            allsigPeptides += 1
                sigPeptidesByCondition.append( sigPeptides )
                totPeptidesByCondition.append( totPeptides )
                aonPeptidesByCondition.append( aonPeptides )
                allsigPeptidesByCondition.append( allsigPeptides )
            Pvalues = [[] for c in range(numConditions-1)]
            Ratios = [[] for c in range(numConditions-1)]
            
            protein_summary_file.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t' %(Accession,GeneName,Location,Composition,Cofactors,Essential,CopyNumber,ProteinMW,fg,numsubunits,uniquesubunits,pI,ProteinFolds,ProteinDomainRanges,ProteinNumDomains,ProteinNumUniqueDomains) +  '\t'.join(str(s)+'\t'+str(t)+'\t'+str(a)+'\t'+str(aS) for (s,t,a,aS) in zip(sigPeptidesByCondition,totPeptidesByCondition,aonPeptidesByCondition,allsigPeptidesByCondition)) + '\n' )
            
            Accession = info_in_line[0]
            GeneName = info_in_line[1]
            Location = info_in_line[4]
            Composition = info_in_line[5]
            Cofactors = info_in_line[6]
            Essential = info_in_line[7]
            CopyNumber = info_in_line[8]
            ProteinMW = info_in_line[9]
            fg = info_in_line[10]
            numsubunits = info_in_line[11]
            uniquesubunits = info_in_line[12]
            pI = info_in_line[13]
            ProteinFolds = info_in_line[14]
            ProteinDomainRanges = info_in_line[15]
            ProteinNumDomains = info_in_line[16]
            ProteinNumUniqueDomains = info_in_line[17]
            DomainLocation = info_in_line[18]
            PositionInDomain = info_in_line[19]
            for c in range(numConditions-1):
                Ratios[c].append( float(info_in_line[19 + numConditions + c] ) )
                Pvalues[c].append( float(info_in_line[18 + 2*numConditions + c] ) )

sigPeptidesByCondition = []
totPeptidesByCondition = []
aonPeptidesByCondition = []
allsigPeptidesByCondition = [] 
for c in range(numConditions-1):
    sigPeptides = 0
    totPeptides = 0
    aonPeptides = 0
    allsigPeptides = 0
    for (Ratio,Pval) in zip(Ratios[c],Pvalues[c]):
        if Ratio != 0 and Pval != 0:
            totPeptides += 1
            if (6 > np.abs(Ratio) > 1 and Pval > 2):
                sigPeptides += 1
                allsigPeptides += 1
            if np.abs(Ratio) > 6 and Pval > 1.8:
                aonPeptides += 1 
                allsigPeptides += 1
    sigPeptidesByCondition.append( sigPeptides )
    totPeptidesByCondition.append( totPeptides )
    aonPeptidesByCondition.append( aonPeptides )
    allsigPeptidesByCondition.append( allsigPeptides )
protein_summary_file.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t' %(Accession,GeneName,Location,Composition,Cofactors,Essential,CopyNumber,ProteinMW,fg,numsubunits,uniquesubunits,pI,ProteinFolds,ProteinDomainRanges,ProteinNumDomains,ProteinNumUniqueDomains) +  '\t'.join(str(s)+'\t'+str(t)+'\t'+str(a)+'\t'+str(aS) for (s,t,a,aS) in zip(sigPeptidesByCondition,totPeptidesByCondition,aonPeptidesByCondition,allsigPeptidesByCondition)) + '\n' )


protein_summary_file.close()
peptideLip_output_file.close()

#Sort the peptide file on the Domain Location so that way we can use the same ordering feature to make a domain-level output
peptideLip_output_df = pd.read_csv( peptideLip_output_fn , sep='\t', header='infer')   
peptideLip_output_df = peptideLip_output_df.sort_values(by=['GeneName','Domain Location','NormPeptideRatio1'])
peptideLip_output_df.to_csv( peptideLip_output_fn , sep='\t', index=False)


peptideLip_output_file = open( peptideLip_output_fn, 'r' )
domain_summary_fn = peptideLip_fn.split('.')[0] + '_summary18Domain.txt'
domain_summary_file = open( domain_summary_fn, 'w')

DomainLocation = ''
Pvalues = [[] for c in range(numConditions-1)]
Ratios = [[] for c in range(numConditions-1)]
next(peptideLip_output_file)

domain_summary_file.write( 'Accession\tDomain Location\tFold\tDomain Order\t' + '\t'.join('SigPeptides'+str(s)+'\t'+'TotPeptides'+str(s)+'\t'+'AonPeptides'+str(s)+'\t'+'AllSigPeptides'+str(s) for s in range(1,numConditions)) + '\n' )
for line in peptideLip_output_file:
    info_in_line = line.split('\t')
    if info_in_line[18] == DomainLocation: #we're on the same domain location
        for c in range(numConditions-1):
            Ratios[c].append( float(info_in_line[19 + numConditions + c] ) )
            Pvalues[c].append( float(info_in_line[18 + 2*numConditions + c] ) )
    else: #we've started looking at a new gene
        if DomainLocation == '': #We've just started; this was the first gene
            Accession = info_in_line[0]
            GeneName = info_in_line[1]
            DomainLocation = info_in_line[18]
            ProteinFolds = info_in_line[14]
            if ProteinFolds != 'Unknown':
                ProteinFolds = ast.literal_eval( info_in_line[14] )
                ProteinDomainRanges = ast.literal_eval( info_in_line[15] )
                if type(ProteinFolds) != list:
                    try:
                        ProteinFolds = ast.literal_eval( ProteinFolds )
                    except SyntaxError:
                        ProteinFolds = [ 'domain'+str(x) for x in range(len(ProteinDomainRanges))] #for one edge case, ushA, where the quotation remarks refuse to get formatted into a list no matter what I try!
                #Find the particular fold and order of domain (in sequence)
                domainRange = re.sub( '[^0-9-,]', '', DomainLocation) 
                if domainRange != '':
                    try:
                        thisFold = ProteinFolds[ ProteinDomainRanges.index(domainRange) ]
                        ProteinDomainStartPoints = [ int(z.split('-')[0])  for z in ProteinDomainRanges ]     
                        orderings = rankdata( ProteinDomainStartPoints ) 
                        thisOrder = orderings[ ProteinDomainRanges.index(domainRange) ]  
                    except ValueError:
                        thisFold = ProteinFolds[0]
                        thisOrder = 1
                else:
                    thisFold = 'Linkers'
                    thisOrder = 'Linkers'
            else:
                thisFold = 'Unknown'
                thisOrder = 'None'
            for c in range(numConditions-1):
                Ratios[c].append( float(info_in_line[19 + numConditions + c] ) )
                Pvalues[c].append( float(info_in_line[18 + 2*numConditions + c] ) )
        else: #We've just started looking at a new gene
            sigPeptidesByCondition = []
            totPeptidesByCondition = []
            aonPeptidesByCondition = []
            allsigPeptidesByCondition = []
            for c in range(numConditions-1):
                sigPeptides = 0
                totPeptides = 0
                aonPeptides = 0
                allsigPeptides = 0
                for (Ratio,Pval) in zip(Ratios[c],Pvalues[c]):
                    if Ratio != 0 and Pval != 0:
                        totPeptides += 1
                        if (6 > np.abs(Ratio) > 1) and Pval > 2:
                            sigPeptides += 1
                            allsigPeptides += 1
                        if np.abs(Ratio) > 6 and Pval > 1.8:
                            aonPeptides += 1
                            allsigPeptides += 1
                sigPeptidesByCondition.append( sigPeptides )
                totPeptidesByCondition.append( totPeptides )
                aonPeptidesByCondition.append( aonPeptides )
                allsigPeptidesByCondition.append( allsigPeptides )
            Pvalues = [[] for c in range(numConditions-1)]
            Ratios = [[] for c in range(numConditions-1)]
            
            domain_summary_file.write('%s\t%s\t%s\t%s\t' %(Accession,DomainLocation,thisFold,thisOrder) +  '\t'.join(str(s)+'\t'+str(t)+'\t'+str(a)+'\t'+str(aS) for (s,t,a,aS) in zip(sigPeptidesByCondition,totPeptidesByCondition,aonPeptidesByCondition,allsigPeptidesByCondition)) + '\n' )
            
            Accession = info_in_line[0]
            GeneName = info_in_line[1]
            DomainLocation = info_in_line[18]
            ProteinFolds = info_in_line[14]
            if ProteinFolds != 'Unknown':
                ProteinFolds = ast.literal_eval( info_in_line[14] )
                ProteinDomainRanges = ast.literal_eval( info_in_line[15] )
                if type(ProteinFolds) != list:
                    try:
                        ProteinFolds = ast.literal_eval( ProteinFolds )
                    except SyntaxError:
                        ProteinFolds = [ 'domain'+str(x) for x in range(len(ProteinDomainRanges))] #for one edge case, ushA, where the quotation remarks refuse to get formatted into a list no matter what I try!
                #Find the particular fold and order of domain (in sequence)
                domainRange = re.sub( '[^0-9-,]', '', DomainLocation) 
                if domainRange != '':
                    try:
                        thisFold = ProteinFolds[ ProteinDomainRanges.index(domainRange) ]
                        ProteinDomainStartPoints = [ int(z.split('-')[0])  for z in ProteinDomainRanges ]     
                        orderings = rankdata( ProteinDomainStartPoints ) 
                        thisOrder = orderings[ ProteinDomainRanges.index(domainRange) ]  
                    except ValueError:
                        thisFold = ProteinFolds[0]
                        thisOrder = 1
                else:
                    thisFold = 'Linkers'
                    thisOrder = 'Linkers'
            else:
                thisFold = 'Unknown'
                thisOrder = 'None'
            for c in range(numConditions-1):
                Ratios[c].append( float(info_in_line[19 + numConditions + c] ) )
                Pvalues[c].append( float(info_in_line[18 + 2*numConditions + c] ) )

sigPeptidesByCondition = []
totPeptidesByCondition = []
aonPeptidesByCondition = []
allsigPeptidesByCondition = []
for c in range(numConditions-1):
    sigPeptides = 0
    totPeptides = 0
    aonPeptides = 0
    allsigPeptides = 0
    for (Ratio,Pval) in zip(Ratios[c],Pvalues[c]):
        if Ratio != 0 and Pval != 0:
            totPeptides += 1
            if (6 > np.abs(Ratio) > 1) and Pval > 2:
                sigPeptides += 1
                allsigPeptides +=1
            if np.abs(Ratio) > 6 and Pval > 1.8:
                aonPeptides += 1
                allsigPeptides += 1
    sigPeptidesByCondition.append( sigPeptides )
    totPeptidesByCondition.append( totPeptides )
    aonPeptidesByCondition.append( aonPeptides )
    allsigPeptidesByCondition.append( allsigPeptides )
domain_summary_file.write('%s\t%s\t%s\t%s\t' %(Accession,DomainLocation,thisFold,thisOrder) +  '\t'.join(str(s)+'\t'+str(t)+'\t'+str(a)+'\t'+str(aS) for (s,t,a,aS) in zip(sigPeptidesByCondition,totPeptidesByCondition,aonPeptidesByCondition,allsigPeptidesByCondition)) + '\n' )


domain_summary_file.close()
peptideLip_output_file.close()