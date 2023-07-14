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

proteinControl_fn = '20210722_Y_Ctrl_LFQ_MsFragger.txt'
peptideLip_fn = sys.argv[1]
#consensus_fn = sys.argv[2]
numConditions = 2

proteinControl_output_fn = peptideLip_fn.split('.')[0] + '_ctrl_out18_10+.txt'
peptideLip_output_fn = peptideLip_fn.split('.')[0] + '_out18_10+.txt'

proteinControl_file = open( proteinControl_fn , 'r' )
peptideLip_file = open( peptideLip_fn, 'r' )
#consensus_file = open( consensus_fn, 'r' )
proteinControl_output_file = open( proteinControl_output_fn, 'w' )
peptideLip_output_file = open( peptideLip_output_fn, 'w' )

'''
#Read through the location info file and create a working dictionary
location_info_file = open('/Users/stephenfried/Dropbox/Systems Folding/location_info_pt.txt' , 'r')
location_dict = {}
for line in location_info_file:
    location_dict[line.split('\t')[0].lower()] = [z.lower() for z in line.split('\t')[1].split(',')] + [ line.split('\t')[2][:] ]
location_info_file.close()
'''

#Reads data from the metadata file and creates a working dictionary
metadata_file = open('yeast_metadata_complete_20221214.txt' , 'r')
next(metadata_file) #skip the header line
complex_dict = {}
for line in metadata_file:
    if line.split('\t')[32] != '':
        ProteinComplex = line.split('\t')[32].strip("\"")
        SubunitNumber = line.split('\t')[33].strip()
        complex_dict[ ProteinComplex ] = SubunitNumber
metadata_file.close()

metadata_file = open('yeast_metadata_complete_20221214.txt' , 'r')
next(metadata_file) #skip the header line
metadata_dict = {}
for line in metadata_file:
    ProteinGeneName = line.split('\t')[1] #extract the name from the metadata file
    Accession = line.split('\t')[3]
    YeastLocusTag = line.split('\t')[0]

    if line.split('\t')[8] == "[]": #extract the locations from the metadata file
        ProteinLocations = 'Unknown'
    else:
        ProteinLocations = line.split('\t')[8].strip("\"")
        ProteinLocations = ast.literal_eval(ProteinLocations)

    if line.split('\t')[9] == "['N/A']":
        ProteinComplexes = 'None'
        SubunitCounts = 'None'
    else:
        ProteinComplexes = line.split('\t')[9].strip("\"")
        ProteinComplexes = ast.literal_eval(ProteinComplexes)
        SubunitCounts = [complex_dict[ProteinComplex] for ProteinComplex in ProteinComplexes]

    if line.split('\t')[7] == '': #extracts copy number information
        CopyNumber = 'N/A'
    else:
        CopyNumber = int(line.split('\t')[7])

    if line.split('\t')[5] == '': #extracts molecular weight
        MW = 'N/A'
    else:
        MW = float( line.split('\t')[5] )

    if line.split('\t')[6] == '': #extracts pI
        pI = 'N/A'
    else:
        pI = float(line.split('\t')[6])

    if line.split('\t')[11] != "N/A":
        pSup = line.split('\t')[11]
    else:
        pSup = "N/A"

    if line.split('\t')[12] != '':
        dpSup_30_42 = line.split('\t')[12]
    else:
        dpSup_30_42 = 'N/A'

    if line.split('\t')[13] != '':
        dpSup_42_46 = line.split('\t')[13]
    else:
        dpSup_42_46 = 'N/A'

    if line.split('\t')[16] != '[]':
        # aa_ranges = re.sub( "," , "+" , line.split('\t')[12] )
        aa_ranges = ast.literal_eval(line.split('\t')[15])

        if type(aa_ranges) == str:
            aa_ranges = ast.literal_eval(aa_ranges)

        architectures = line.split('\t')[16].strip("\"")
        architectures = ast.literal_eval(architectures)
        xgroups = line.split('\t')[17] # .strip("\"")
        xgroups = ast.literal_eval(xgroups)
        tgroups = line.split('\t')[18] # .strip("\"")
        tgroups = ast.literal_eval(tgroups)
    else:
        aa_ranges = 'NA'
        architectures = 'NA'
        xgroups = 'NA'
        tgroups = 'NA'

    numDomains = int( line.split('\t')[20] )

    if line.split('\t')[21] == '':
        percentWithinDomain = 'NA'
    else:
        percentWithinDomain = float( line.split('\t')[21] )

    if line.split('\t')[22] == '':
        sequence = 'NA'
    else:
        sequence = line.split('\t')[22]

    if line.split('\t')[23] == ' ':
        percentDisorder = 'NA'
    else:
        percentDisorder = float( line.split('\t')[23] )

    if line.split('\t')[24] == '[]' and line.split('\t')[25].strip('\n') == '[]' :
        disorderedRegions = 'NA'
        orderedRegions = 'NA'
    else:
        if line.split('\t')[24] != '[]':
            disorderedRegions = ast.literal_eval(ast.literal_eval(line.split('\t')[24]))
        else:
            disorderedRegions = []

        if line.split('\t')[25].strip('\n') != '[]':
            orderedRegions = ast.literal_eval(ast.literal_eval(line.split('\t')[25]))
        else:
            orderedRegions = []

    if line.split('\t')[26] == '':
        CCT_sites = 'NA'
    else:
        CCT_sites = int( line.split('\t')[26] )

    if line.split('\t')[27] == '':
        Hsp70_sites = 'NA'
    else:
        Hsp70_sites = int( line.split('\t')[27] )

    if 'X' in line.split('\t')[28]:
        SGMayorLab = 'Y'
    else:
        SGMayorLab = 'NA'

    if 'X' in line.split('\t')[29]:
        SGParkerLab = 'Y'
    else:
        SGParkerLab = 'NA'

    if 'X' in line.split('\t')[30]:
        SGBukauLab = 'Y'
    else:
        SGBukauLab = 'NA'

    if 'X' in line.split('\t')[31]:
        SGDrummondLab = 'Y'
    else:
        SGDrummondLab = 'NA'

    metadata_dict.update({Accession:[ProteinGeneName, YeastLocusTag, ProteinLocations, ProteinComplexes, SubunitCounts, CopyNumber, MW, pI, pSup, dpSup_30_42, dpSup_42_46, aa_ranges, architectures, xgroups, tgroups, numDomains, percentWithinDomain, sequence, percentDisorder, disorderedRegions, orderedRegions, CCT_sites, Hsp70_sites, SGMayorLab, SGParkerLab, SGBukauLab, SGDrummondLab]})
metadata_file.close()

#Loop through Protein Control txt file line by line
proteinControl_output_file.write( 'GeneName\tLocus Tag\tAccession\tLocations\tComplexes\tSubunitCounts\tCopy Number\tMolecular Weight\tpI\tpSup\tdpSup_30_42\tdpSup_42_46\tDomainRanges\tArchitectures\tXgroups\tTgroups\tNumDomains\tpercentWithinDomain\tpercentDisordered\tdisorderedRegions\torderedRegions\tCCT_sites\tHsp70_sites\tSGMayorLab\tSGParkerLab\tSGBukauLab\tSGDrummondLab\tProteinRatio\tProteinPValue\n' )
for line in proteinControl_file:
    info_in_line = line.split('\t')
    if info_in_line[0] == 'FALSE':
        if info_in_line[1] == 'Medium' or info_in_line[1] == 'High':
            Accession = info_in_line[3]
            #GeneName = info_in_line[4] #We will adjust this to look for 'GN=' once we are using the new FASTA
            ProteinGeneName = info_in_line[4].split('=')[3].split(' ')[0] #Takes Protein Gene Name from Description line
            if info_in_line[28] != '' and info_in_line[29] != '': #If quans are available
                #ProteinRatio = np.log2(float(info_in_line[26]))
                #ProteinPValue = np.abs(np.log10(float(info_in_line[27])))
                ProteinRatio = np.log2(float(info_in_line[28]))
                ProteinPValue = -np.log10(float(info_in_line[29]))
            else: #if protein_ratio and protein_p_value was blank
                ProteinRatio = float(0.0)
                ProteinPValue = float(0.0)
            metadata_dict[Accession].append( ProteinRatio )
            metadata_dict[Accession].append( ProteinPValue )
            proteinControl_output_file.write( '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' %(metadata_dict[Accession][0], metadata_dict[Accession][1], Accession, metadata_dict[Accession][2], metadata_dict[Accession][3], metadata_dict[Accession][4],metadata_dict[Accession][5],metadata_dict[Accession][6],metadata_dict[Accession][7],metadata_dict[Accession][8],metadata_dict[Accession][9],metadata_dict[Accession][10],metadata_dict[Accession][11],metadata_dict[Accession][12],metadata_dict[Accession][13],metadata_dict[Accession][14],metadata_dict[Accession][15],metadata_dict[Accession][16],metadata_dict[Accession][18],metadata_dict[Accession][19],metadata_dict[Accession][20],metadata_dict[Accession][21],metadata_dict[Accession][22],metadata_dict[Accession][23],metadata_dict[Accession][24],metadata_dict[Accession][25],metadata_dict[Accession][26],metadata_dict[Accession][27],metadata_dict[Accession][28]) )

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
peptideLip_output_file.write( 'GeneName\tLocus Tag\tAccession\tPeptide Sequence\tproteinaseKsite\tLocations\tComplexes\tSubunitCounts\tCopy Number\tMolecular Weight\tpI\tpSup\tdpSup_30_42\tdpSup_42_46\tArchitectures\tXgroups\tTgroups\tDomainRanges\tNumDomains\tpercentWithinDomain\tpercentDisordered\tdisorderedRegions\torderedRegions\tCCT_sites\tHsp70_sites\tSGMayorLab\tSGParkerLab\tSGBukauLab\tSGDrummondLab\tDomain Location\tPosition in Domain\t' + '\t'.join('PeptideRatio'+str(s) for s in range(1,numConditions)) + '\t' + '\t'.join('NormPeptideRatio'+str(s) for s in range(1,numConditions)) + '\t' + '\t'.join('PeptidePValue'+str(s) for s in range(1,numConditions)) + '\t' + '\t'.join('RefoldedCoeffVar'+str(s) for s in range(1,numConditions))+ '\n' )

#Loop through Peptide LiP File and append AbundanceRatio Data to it
#for line in peptideLip_file:
for x in range(177):
    line = next(peptideLip_file)
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
                        #if quanArrayN[c][j][0][3] > 1:
                        #    for feature in consensus_features: #first check if there were any non-assigned features that could possibly match to this peptide's consensus feature
                        #        if isclose(feature[0],quanArrayN[c][j][0][0]): #m/z is within 10 ppm
                        #            if True: #charge state match
                        #                if abs(feature[2] - quanArrayN[c][j][0][2]) < 5: #retention times within five minutes
                        #                    if feature[3] == 0: #this concensus feature was not already mapped to other peptides; e.g., no PSMs
                        #                        quanArrayN[c][j][1][0] = quanArrayN[c][j][1][0] + feature[4][0]
                        #                        quanArrayN[c][j][1][1] = quanArrayN[c][j][1][1] + feature[4][1]
                        #                        quanArrayN[c][j][1][2] = quanArrayN[c][j][1][2] + feature[4][2]
                        #                        quanArrayR[c][j][1][0] = quanArrayR[c][j][1][0] + feature[4][3*c+3]
                        #                        quanArrayR[c][j][1][1] = quanArrayR[c][j][1][1] + feature[4][3*c+4]
                        #                        quanArrayR[c][j][1][2] = quanArrayR[c][j][1][2] + feature[4][3*c+5]
                        #    if quanArrayN[c][j][1].count( 0 ) == 3 and quanArrayR[c][j][1].count( 0 ) == 0:
                        quanArrayNfilt[c].append( [quanArrayN[c][j][0], [1000,1000,1000]] )
                        quanArrayRfilt[c].append( quanArrayR[c][j] )
                        #    else:
                        #        quanArrayNfilt[c].append( quanArrayN[c][j] )
                        #        quanArrayRfilt[c].append( quanArrayR[c][j] )
                    elif quanArrayN[c][j][1].count( 0 ) == 0 and quanArrayR[c][j][1].count( 0 ) == 3: #Missing data being used to impute low value for Rs
                        #if quanArrayN[c][j][0][3] > 1:
                        #    for feature in consensus_features: #first check if there were any non-assigned features that could possibly match to this peptide's consensus feature
                        #        if isclose(feature[0],quanArrayN[c][j][0][0]): #m/z is within 10 ppm
                        #            if True: #charge state match
                        #                if abs(feature[2] - quanArrayN[c][j][0][2]) < 5: #retention times within five minutes
                        #                    if feature[3] == 0: #this concensus feature was not already mapped to other peptides; e.g., no PSMs
                        #                        quanArrayN[c][j][1][0] = quanArrayN[c][j][1][0] + feature[4][0]
                        #                        quanArrayN[c][j][1][1] = quanArrayN[c][j][1][1] + feature[4][1]
                        #                        quanArrayN[c][j][1][2] = quanArrayN[c][j][1][2] + feature[4][2]
                        #                        quanArrayR[c][j][1][0] = quanArrayR[c][j][1][0] + feature[4][3*c+3]
                        #                        quanArrayR[c][j][1][1] = quanArrayR[c][j][1][1] + feature[4][3*c+4]
                        #                        quanArrayR[c][j][1][2] = quanArrayR[c][j][1][2] + feature[4][3*c+5]
                        #    if quanArrayN[c][j][1].count( 0 ) == 0 and quanArrayR[c][j][1].count( 0 ) == 3:
                        quanArrayNfilt[c].append( quanArrayN[c][j] )
                        quanArrayRfilt[c].append( [quanArrayR[c][j][0], [1000,1000,1000]] )
                        #    else:
                        #        quanArrayNfilt[c].append( quanArrayN[c][j] )
                        #        quanArrayRfilt[c].append( quanArrayR[c][j] )
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

            Peptides_table.append( [ProteinGeneName, ProteinLocusTag, Accession, Pep_Seq, proteinaseKsite, ProteinLocation, ProteinComplexes, SubunitCounts, CopyNumber, ProteinMW, pI, pSup, ProteinDeltapSup_30_42, ProteinDeltapSup_42_46, ProteinDomainsArch, ProteinDomainsX, ProteinDomainsT, ProteinDomainRanges, ProteinNumDomains, ProteinPercentWithinDomain, ProteinPerecentDisordered, ProteinDisorderedRanges, ProteinOrderedRanges, ProteinCCTsites, ProteinHsp70sites, ProteinSGMayor, ProteinSGParker, ProteinSGBukau, ProteinSGDrummond, DomainLocation, PositionInDomain, PeptideRatioByCondition, normalizedPeptideRatioByCondition, PeptidePValueByCondition, VariationByCondition] )
            quanArrayR = [[] for i in range(numConditions-1)]
            quanArrayN = [[] for i in range(numConditions-1)]
            quanArrayRfilt = [[] for i in range(numConditions-1)]
            quanArrayNfilt = [[] for i in range(numConditions-1)]
            charge_states_seen = []
            Pep_Seq = ''
        #ProteinDescription = info_in_line[4]

        Accession = info_in_line[3]

        #   metadata_dict.update({Accession:[ProteinGeneName, YeastLocusTag, ProteinLocations, ProteinComplexes, SubunitCounts, CopyNumber, MW, pI, pSup, aa_ranges, architectures, xgroups, tgroups, numDomains, percentWithinDomain]})
        if Accession in metadata_dict.keys():
            ProteinGeneName = metadata_dict[Accession][0]
            ProteinLocusTag = metadata_dict[Accession][1]
            ProteinLocation = metadata_dict[Accession][2]
            ProteinComplexes = metadata_dict[Accession][3]
            SubunitCounts = metadata_dict[Accession][4]
            CopyNumber = metadata_dict[Accession][5]
            ProteinMW = metadata_dict[Accession][6]
            pI = metadata_dict[Accession][7]
            pSup = metadata_dict[Accession][8]
            ProteinDeltapSup_30_42 = metadata_dict[Accession][9]
            ProteinDeltapSup_42_46 = metadata_dict[Accession][10]
            ProteinDomainRanges = metadata_dict[Accession][11]
            ProteinDomainsArch = metadata_dict[Accession][12]
            ProteinDomainsX = metadata_dict[Accession][13]
            ProteinDomainsT = metadata_dict[Accession][14]
            ProteinNumDomains = metadata_dict[Accession][15]
            ProteinPercentWithinDomain = metadata_dict[Accession][16]
            ProteinPerecentDisordered = metadata_dict[Accession][18]
            ProteinDisorderedRanges = metadata_dict[Accession][19]
            ProteinOrderedRanges = metadata_dict[Accession][20]
            ProteinCCTsites = metadata_dict[Accession][21]
            ProteinHsp70sites = metadata_dict[Accession][22]
            ProteinSGMayor = metadata_dict[Accession][23]
            ProteinSGParker = metadata_dict[Accession][24]
            ProteinSGBukau = metadata_dict[Accession][25]
            ProteinSGDrummond = metadata_dict[Accession][26]
            ProteinDomainLists = []

            if ProteinDomainRanges == 'NA':
                ProteinDomainLists.append( 'NA' )
            else:
                for entry in ProteinDomainRanges:
                    if '+' not in entry:
                        start=int(entry.split('-')[0])
                        end=int(entry.split('-')[1])
                        dr = [x for x in range(start,end+1)]
                        ProteinDomainLists.append( dr )
                    elif '+' in entry:
                        dr = []
                        componentranges = entry.split('+')
                        for componentrange in componentranges:
                            start=int(componentrange.split('-')[0])
                            end=int(componentrange.split('-')[1])
                            dr.append( [x for x in range(start,end+1)] )
                        ProteinDomainLists.append( [item for sublist in dr for item in sublist] )
        else:
            ProteinGeneName = info_in_line[4].split('=')[3].split(' ')[0]
            ProteinLocusTag = 'NA'
            ProteinLocation = 'NA'
            ProteinComplexes = 'NA'
            SubunitCounts = 'NA'
            CopyNumber = 'NA'
            ProteinMW = 'NA'
            pI = 'NA'
            pSup = 'NA'
            ProteinDeltapSup_30_42 = 'NA'
            ProteinDeltapSup_42_46 = 'NA'
            ProteinDomainRanges = 'NA'
            ProteinDomainsArch = 'NA'
            ProteinDomainsX = 'NA'
            ProteinDomainsT = 'NA'
            ProteinNumDomains = 'NA'
            ProteinPercentWithinDomain = 'NA'
            ProteinPerecentDisordered = 'NA'
            ProteinDisorderedRanges = 'NA'
            ProteinOrderedRanges = 'NA'
            ProteinCCTsites = 'NA'
            ProteinHsp70sites = 'NA'
            ProteinSGMayor = 'NA'
            ProteinSGParker = 'NA'
            ProteinSGBukau = 'NA'
            ProteinSGDrummond = 'NA'


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
                        #if quanArrayN[c][j][0][3] > 1:
                        #    for feature in consensus_features: #first check if there were any non-assigned features that could possibly match to this peptide's consensus feature
                        #        if isclose(feature[0],quanArrayN[c][j][0][0]): #m/z is within 10 ppm
                        #            if True: #charge state match
                        #                if abs(feature[2] - quanArrayN[c][j][0][2]) < 5: #retention times within five minutes
                        #                    if feature[3] == 0: #this concensus feature was not already mapped to other peptides; e.g., no PSMs
                        #                        quanArrayN[c][j][1][0] = quanArrayN[c][j][1][0] + feature[4][0]
                        #                        quanArrayN[c][j][1][1] = quanArrayN[c][j][1][1] + feature[4][1]
                        #                        quanArrayN[c][j][1][2] = quanArrayN[c][j][1][2] + feature[4][2]
                        #                        quanArrayR[c][j][1][0] = quanArrayR[c][j][1][0] + feature[4][3*c+3]
                        #                        quanArrayR[c][j][1][1] = quanArrayR[c][j][1][1] + feature[4][3*c+4]
                        #                        quanArrayR[c][j][1][2] = quanArrayR[c][j][1][2] + feature[4][3*c+5]
                        #    if quanArrayN[c][j][1].count( 0 ) == 3 and quanArrayR[c][j][1].count( 0 ) == 0:
                        quanArrayNfilt[c].append( [quanArrayN[c][j][0], [1000,1000,1000]] )
                        quanArrayRfilt[c].append( quanArrayR[c][j] )
                        #    else:
                        #        quanArrayNfilt[c].append( quanArrayN[c][j] )
                        #        quanArrayRfilt[c].append( quanArrayR[c][j] )
                    elif quanArrayN[c][j][1].count( 0 ) == 0 and quanArrayR[c][j][1].count( 0 ) == 3: #Missing data being used to impute low value for Rs
                        #if quanArrayN[c][j][0][3] > 1:
                        #    for feature in consensus_features: #first check if there were any non-assigned features that could possibly match to this peptide's consensus feature
                        #        if isclose(feature[0],quanArrayN[c][j][0][0]): #m/z is within 10 ppm
                        #            if True: #charge state match
                        #                if abs(feature[2] - quanArrayN[c][j][0][2]) < 5: #retention times within five minutes
                        #                    if feature[3] == 0: #this concensus feature was not already mapped to other peptides; e.g., no PSMs
                        #                        quanArrayN[c][j][1][0] = quanArrayN[c][j][1][0] + feature[4][0]
                        #                        quanArrayN[c][j][1][1] = quanArrayN[c][j][1][1] + feature[4][1]
                        #                        quanArrayN[c][j][1][2] = quanArrayN[c][j][1][2] + feature[4][2]
                        #                        quanArrayR[c][j][1][0] = quanArrayR[c][j][1][0] + feature[4][3*c+3]
                        #                        quanArrayR[c][j][1][1] = quanArrayR[c][j][1][1] + feature[4][3*c+4]
                        #                        quanArrayR[c][j][1][2] = quanArrayR[c][j][1][2] + feature[4][3*c+5]
                        #    if quanArrayN[c][j][1].count( 0 ) == 0 and quanArrayR[c][j][1].count( 0 ) == 3:
                        quanArrayNfilt[c].append( quanArrayN[c][j] )
                        quanArrayRfilt[c].append( [quanArrayR[c][j][0], [1000,1000,1000]] )
                        #    else:
                        #        quanArrayNfilt[c].append( quanArrayN[c][j] )
                        #        quanArrayRfilt[c].append( quanArrayR[c][j] )
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
            Peptides_table.append( [ProteinGeneName, ProteinLocusTag, Accession, Pep_Seq, proteinaseKsite, ProteinLocation, ProteinComplexes, SubunitCounts, CopyNumber, ProteinMW, pI, pSup, ProteinDeltapSup_30_42, ProteinDeltapSup_42_46, ProteinDomainsArch, ProteinDomainsX, ProteinDomainsT, ProteinDomainRanges, ProteinNumDomains, ProteinPercentWithinDomain, ProteinPerecentDisordered, ProteinDisorderedRanges, ProteinOrderedRanges, ProteinCCTsites, ProteinHsp70sites, ProteinSGMayor, ProteinSGParker, ProteinSGBukau, ProteinSGDrummond, DomainLocation, PositionInDomain, PeptideRatioByCondition, normalizedPeptideRatioByCondition, PeptidePValueByCondition, VariationByCondition] )
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
                    #PeptideRatio = np.log2(float(info_in_line[16]))
                    #PeptidePValue = np.abs(np.log10(float(info_in_line[17])))
                    if len(metadata_dict[Accession]) == 29:
                        ProteinRatio = metadata_dict[Accession][27]
                        ProteinPValue = metadata_dict[Accession][28]
                    else: #if protein was not identifid in the control sample, then generate filler values
                        ProteinRatio = float(0.0)
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
                            if ProteinDomainsT != 'NA':
                                numDomainsFoundIn = 0
                                for (ProteinDomainList,ProteinDomainRange) in zip(ProteinDomainLists,ProteinDomainRanges):
                                    if int(PositionInProtein[0]) in ProteinDomainList:
                                        DomainLocation = ProteinGeneName + ',' + ProteinDomainRange
                                        numDomainsFoundIn += 1
                                        PositionInDomain = ProteinDomainList.index(int(PositionInProtein[0])) + 1
                                if numDomainsFoundIn != 1:
                                    DomainLocation = ProteinGeneName + ',LNK'
                                    PositionInDomain = 'NA'
                            else:
                                DomainLocation = ProteinGeneName
                                PositionInDomain = 'NA'
                        if proteinaseKsite == 'Back':
                            proteinaseKsite = Pep_Seq.split('.')[1][-1]+PositionInProtein[1]
                            #Find which domain this site is in.
                            if ProteinDomainsT != 'NA':
                                numDomainsFoundIn = 0
                                for (ProteinDomainList,ProteinDomainRange) in zip(ProteinDomainLists,ProteinDomainRanges):
                                    if int(PositionInProtein[1]) in ProteinDomainList:
                                        DomainLocation = ProteinGeneName + ',' + ProteinDomainRange
                                        numDomainsFoundIn += 1
                                        PositionInDomain = ProteinDomainList.index(int(PositionInProtein[1])) + 1
                                if numDomainsFoundIn != 1:
                                    DomainLocation = ProteinGeneName + ',LNK'
                                    PositionInDomain = 'NA'
                            else:
                               DomainLocation = ProteinGeneName
                               PositionInDomain = 'NA'
                    else: #This is NOT a half-tryptic peptide
                        proteinaseKsite = info_in_line[9].split(' ')[1]
                        trypticstart = int(  re.sub( '[^0-9]' , '' ,proteinaseKsite.split('-')[0] ) )
                        trypticend = int(  re.sub( '[^0-9]' , '' ,proteinaseKsite.split('-')[1] ) )
                        trypticList = [x for x in range(trypticstart,trypticend+1)]
                        if ProteinDomainsT != 'NA':
                            numDomainsFoundIn = 0
                            for (ProteinDomainList,ProteinDomainRange) in zip(ProteinDomainLists,ProteinDomainRanges):
                                if all(site in ProteinDomainList for site in trypticList):
                                    DomainLocation = ProteinGeneName + ',' + ProteinDomainRange
                                    numDomainsFoundIn += 1
                                    PositionInDomain = 'NA'
                            if numDomainsFoundIn != 1:
                                DomainLocation = ProteinGeneName + ',LNK'
                                PositionInDomain = 'NA'
                        else:
                            DomainLocation = ProteinGeneName
                            PositionInDomain = 'NA'
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

        currentConcensusFeature = [z for z in info_in_line[12:12+numConditions*3]] # based on normalized conditions
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
                #if quanArrayN[c][j][0][3] > 1:
                #    for feature in consensus_features: #first check if there were any non-assigned features that could possibly match to this peptide's consensus feature
                #        if isclose(feature[0],quanArrayN[c][j][0][0]): #m/z is within 10 ppm
                #            if True: #charge state match
                #                if abs(feature[2] - quanArrayN[c][j][0][2]) < 5: #retention times within five minutes
                #                    if feature[3] == 0: #this concensus feature was not already mapped to other peptides; e.g., no PSMs
                #                        quanArrayN[c][j][1][0] = quanArrayN[c][j][1][0] + feature[4][0]
                #                        quanArrayN[c][j][1][1] = quanArrayN[c][j][1][1] + feature[4][1]
                #                        quanArrayN[c][j][1][2] = quanArrayN[c][j][1][2] + feature[4][2]
                #                        quanArrayR[c][j][1][0] = quanArrayR[c][j][1][0] + feature[4][3*c+3]
                #                        quanArrayR[c][j][1][1] = quanArrayR[c][j][1][1] + feature[4][3*c+4]
                #                        quanArrayR[c][j][1][2] = quanArrayR[c][j][1][2] + feature[4][3*c+5]
                #    if quanArrayN[c][j][1].count( 0 ) == 3 and quanArrayR[c][j][1].count( 0 ) == 0:
                quanArrayNfilt[c].append( [quanArrayN[c][j][0], [1000,1000,1000]] )
                quanArrayRfilt[c].append( quanArrayR[c][j] )
                #    else:
                #        quanArrayNfilt[c].append( quanArrayN[c][j] )
                #        quanArrayRfilt[c].append( quanArrayR[c][j] )
            elif quanArrayN[c][j][1].count( 0 ) == 0 and quanArrayR[c][j][1].count( 0 ) == 3: #Missing data being used to impute low value for Rs
                #if quanArrayN[c][j][0][3] > 1:
                #    for feature in consensus_features: #first check if there were any non-assigned features that could possibly match to this peptide's consensus feature
                #        if isclose(feature[0],quanArrayN[c][j][0][0]): #m/z is within 10 ppm
                #            if True: #charge state match
                #                if abs(feature[2] - quanArrayN[c][j][0][2]) < 5: #retention times within five minutes
                #                    if feature[3] == 0: #this concensus feature was not already mapped to other peptides; e.g., no PSMs
                #                        quanArrayN[c][j][1][0] = quanArrayN[c][j][1][0] + feature[4][0]
                #                        quanArrayN[c][j][1][1] = quanArrayN[c][j][1][1] + feature[4][1]
                #                        quanArrayN[c][j][1][2] = quanArrayN[c][j][1][2] + feature[4][2]
                #                        quanArrayR[c][j][1][0] = quanArrayR[c][j][1][0] + feature[4][3*c+3]
                #                        quanArrayR[c][j][1][1] = quanArrayR[c][j][1][1] + feature[4][3*c+4]
                #                        quanArrayR[c][j][1][2] = quanArrayR[c][j][1][2] + feature[4][3*c+5]
                #    if quanArrayN[c][j][1].count( 0 ) == 0 and quanArrayR[c][j][1].count( 0 ) == 3:
                quanArrayNfilt[c].append( quanArrayN[c][j] )
                quanArrayRfilt[c].append( [quanArrayR[c][j][0], [1000,1000,1000]] )
                #    else:
                #        quanArrayNfilt[c].append( quanArrayN[c][j] )
                #        quanArrayRfilt[c].append( quanArrayR[c][j] )
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
    Peptides_table.append( [ProteinGeneName, ProteinLocusTag, Accession, Pep_Seq, proteinaseKsite, ProteinLocation, ProteinComplexes, SubunitCounts, CopyNumber, ProteinMW, pI, pSup, ProteinDeltapSup_30_42, ProteinDeltapSup_42_46, ProteinDomainsArch, ProteinDomainsX, ProteinDomainsT, ProteinDomainRanges, ProteinNumDomains, ProteinPercentWithinDomain, ProteinPerecentDisordered, ProteinDisorderedRanges, ProteinOrderedRanges, ProteinCCTsites, ProteinHsp70sites, ProteinSGMayor, ProteinSGParker, ProteinSGBukau, ProteinSGDrummond, DomainLocation, PositionInDomain, PeptideRatioByCondition, normalizedPeptideRatioByCondition, PeptidePValueByCondition, VariationByCondition] )
    quanArrayR = [[] for i in range(numConditions-1)]
    quanArrayN = [[] for i in range(numConditions-1)]
    quanArrayRfilt = [[] for i in range(numConditions-1)]
    quanArrayNfilt = [[] for i in range(numConditions-1)]
    charge_states_seen = []
    Pep_Seq = ''

#peptideLip_output_file.write('Number of proteins with half-tryptic petides:\t%s' %str(len(Proteins_with_halftryptics)))
for i in range(len(Peptides_table)):
    if np.count_nonzero( Peptides_table[i][33] ) != 0:   #This would apply when there are multiple conditions.
        line_to_write = '\t'.join(str(s) for s in Peptides_table[i][0:31]) + '\t'+ '\t'.join(str(s) for s in Peptides_table[i][31][:]) + '\t' + '\t'.join(str(s) for s in Peptides_table[i][32][:]) + '\t' + '\t'.join(str(s) for s in Peptides_table[i][33][:]) + '\t' + '\t'.join(str(s) for s in Peptides_table[i][34][:])+ '\n'
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

protein_summary_file.write( 'GeneName\tLocus Tag\tAccession\tLocations\tComplexes\tSubunitCounts\tCopy Number\tMolecular Weight\tpI\tpSup\tdpSup_30_42\tdpSup_42_46\tArchitectures\tXgroups\tTgroups\tDomainRanges\tNumDomains\tPercent within Domain\tPerecent Disordered\tDisordered Regions\tOrdered Regions\tCCT sites\tHsp70 sites\tSG MayorLab\tSG ParkerLab\tSG BukauLab\tSG DrummondLab\t' + '\t'.join('SigPeptides'+str(s)+'\t'+'TotPeptides'+str(s)+'\t'+'AonPeptides'+str(s)+'\t'+'AllSigPeptides'+str(s)+'\t'+'NonrefoldabilityScore'+str(s) for s in range(1,numConditions)) + '\n' )
for line in peptideLip_output_file:
    info_in_line = line.split('\t')
    if info_in_line[0] == GeneName: #we're on the same gene
        for c in range(numConditions-1):
            # edited values to test, revert to ensure consideration of multiple conditions
            # Ratios[c].append( float(info_in_line[19 + numConditions + c] ) ) # normalized
            # Pvalues[c].append( float(info_in_line[20 + 2*numConditions + c] ) )

            Ratios[c].append( float(info_in_line[32] ) ) # normalized
            Pvalues[c].append( float(info_in_line[33] ) )
    else: #we've started looking at a new gene
        if GeneName == '': #We've just started; this was the first gene
            GeneName = info_in_line[0]
            YeastLocusTag = info_in_line[1]
            Accession = info_in_line[2]
            Location = info_in_line[5]
            Complexes = info_in_line[6]
            SubunitCounts = info_in_line[7]
            CopyNumber = info_in_line[8]
            ProteinMW = info_in_line[9]
            pI = info_in_line[10]
            pSup = info_in_line[11]
            ProteinDeltapSup_30_42 = info_in_line[12]
            ProteinDeltapSup_42_46 = info_in_line[13]
            ProteinDomainsArch = info_in_line[14]
            ProteinDomainsX = info_in_line[15]
            ProteinDomainsT = info_in_line[16]
            ProteinDomainRanges = info_in_line[17]
            ProteinNumDomains = info_in_line[18]
            ProteinPercentWithinDomain = info_in_line[19]
            ProteinPerecentDisordered = info_in_line[20]
            ProteinDisorderedRanges = info_in_line[21]
            ProteinOrderedRanges = info_in_line[22]
            ProteinCCTsites = info_in_line[23]
            ProteinHsp70sites = info_in_line[24]
            ProteinSGMayor = info_in_line[25]
            ProteinSGParker = info_in_line[26]
            ProteinSGBukau = info_in_line[27]
            ProteinSGDrummond = info_in_line[28]

            for c in range(numConditions-1):
                # edited values to test, revert to ensure consideration of multiple conditions
                # Ratios[c].append( float(info_in_line[19 + numConditions + c] ) ) # normalized
                # Pvalues[c].append( float(info_in_line[20 + 2*numConditions + c] ) )

                Ratios[c].append( float(info_in_line[32] ) ) # normalized
                Pvalues[c].append( float(info_in_line[33] ) )
        else: #We've just started looking at a new gene
            sigPeptidesByCondition = []
            totPeptidesByCondition = []
            aonPeptidesByCondition = [] #All or Nothing Peptides, those which are only present (or completely not present) in the refolded sample
            allsigPeptidesByCondition = []
            nonrefoldabilityScoreByCondition = []
            for c in range(numConditions-1):
                sigPeptides = 0
                totPeptides = 0
                aonPeptides = 0
                allsigPeptides = 0
                nonrefoldabilityScore = float(0)
                for (Ratio,Pval) in zip(Ratios[c],Pvalues[c]):
                    if Ratio != 0 and Pval != 0:
                        totPeptides += 1
                        if (6 > np.abs(Ratio) > 1 and Pval > 2):
                            sigPeptides += 1
                            allsigPeptides += 1
                            nonrefoldabilityScore += np.abs(Ratio)
                        if np.abs(Ratio) > 6 and Pval > 1.8:
                            aonPeptides += 1
                            allsigPeptides += 1
                            nonrefoldabilityScore += np.abs(Ratio)
                sigPeptidesByCondition.append( sigPeptides )
                totPeptidesByCondition.append( totPeptides )
                aonPeptidesByCondition.append( aonPeptides )
                allsigPeptidesByCondition.append( allsigPeptides )
                nonrefoldabilityScoreByCondition.append( nonrefoldabilityScore / totPeptides )
            Pvalues = [[] for c in range(numConditions-1)]
            Ratios = [[] for c in range(numConditions-1)]

            protein_summary_file.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t' %(GeneName, YeastLocusTag, Accession, Location, Complexes, SubunitCounts, CopyNumber, ProteinMW, pI, pSup, ProteinDeltapSup_30_42, ProteinDeltapSup_42_46, ProteinDomainsArch, ProteinDomainsX, ProteinDomainsT, ProteinDomainRanges, ProteinNumDomains, ProteinPercentWithinDomain, ProteinPerecentDisordered, ProteinDisorderedRanges, ProteinOrderedRanges, ProteinCCTsites, ProteinHsp70sites, ProteinSGMayor, ProteinSGParker, ProteinSGBukau, ProteinSGDrummond) +  '\t'.join(str(s)+'\t'+str(t)+'\t'+str(a)+'\t'+str(aS)+'\t'+str(nrs) for (s,t,a,aS,nrs) in zip(sigPeptidesByCondition,totPeptidesByCondition,aonPeptidesByCondition,allsigPeptidesByCondition,nonrefoldabilityScoreByCondition)) + '\n' )

            GeneName = info_in_line[0]
            YeastLocusTag = info_in_line[1]
            Accession = info_in_line[2]
            Location = info_in_line[5]
            Complexes = info_in_line[6]
            SubunitCounts = info_in_line[7]
            CopyNumber = info_in_line[8]
            ProteinMW = info_in_line[9]
            pI = info_in_line[10]
            pSup = info_in_line[11]
            ProteinDeltapSup_30_42 = info_in_line[12]
            ProteinDeltapSup_42_46 = info_in_line[13]
            ProteinDomainsArch = info_in_line[14]
            ProteinDomainsX = info_in_line[15]
            ProteinDomainsT = info_in_line[16]
            ProteinDomainRanges = info_in_line[17]
            ProteinNumDomains = info_in_line[18]
            ProteinPercentWithinDomain = info_in_line[19]
            ProteinPerecentDisordered = info_in_line[20]
            ProteinDisorderedRanges = info_in_line[21]
            ProteinOrderedRanges = info_in_line[22]
            ProteinCCTsites = info_in_line[23]
            ProteinHsp70sites = info_in_line[24]
            ProteinSGMayor = info_in_line[25]
            ProteinSGParker = info_in_line[26]
            ProteinSGBukau = info_in_line[27]
            ProteinSGDrummond = info_in_line[28]
            #DomainLocation = info_in_line[18]
            for c in range(numConditions-1):
                # edited values to test, revert to ensure consideration of multiple conditions
                # Ratios[c].append( float(info_in_line[19 + numConditions + c] ) )
                # Pvalues[c].append( float(info_in_line[20 + 2*numConditions + c] ) )

                Ratios[c].append( float(info_in_line[32] ) )
                Pvalues[c].append( float(info_in_line[33] ) )

sigPeptidesByCondition = []
totPeptidesByCondition = []
aonPeptidesByCondition = []
allsigPeptidesByCondition = []
nonrefoldabilityScoreByCondition = []
for c in range(numConditions-1):
    sigPeptides = 0
    totPeptides = 0
    aonPeptides = 0
    allsigPeptides = 0
    nonrefoldabilityScore = float(0)
    for (Ratio,Pval) in zip(Ratios[c],Pvalues[c]):
        if Ratio != 0 and Pval != 0:
            totPeptides += 1
            if (6 > np.abs(Ratio) > 1 and Pval > 2):
                sigPeptides += 1
                allsigPeptides += 1
                nonrefoldabilityScore += np.abs(Ratio)
            if np.abs(Ratio) > 6 and Pval > 1.8:
                aonPeptides += 1
                allsigPeptides += 1
                nonrefoldabilityScore += np.abs(Ratio)
    sigPeptidesByCondition.append( sigPeptides )
    totPeptidesByCondition.append( totPeptides )
    aonPeptidesByCondition.append( aonPeptides )
    allsigPeptidesByCondition.append( allsigPeptides )
    nonrefoldabilityScoreByCondition.append( nonrefoldabilityScore / totPeptides )

protein_summary_file.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t' %(GeneName, YeastLocusTag, Accession, Location, Complexes, SubunitCounts, CopyNumber, ProteinMW, pI, pSup, ProteinDeltapSup_30_42, ProteinDeltapSup_42_46, ProteinDomainsArch, ProteinDomainsX, ProteinDomainsT, ProteinDomainRanges, ProteinNumDomains, ProteinPercentWithinDomain, ProteinPerecentDisordered, ProteinDisorderedRanges, ProteinOrderedRanges, ProteinCCTsites, ProteinHsp70sites, ProteinSGMayor, ProteinSGParker, ProteinSGBukau, ProteinSGDrummond) +  '\t'.join(str(s)+'\t'+str(t)+'\t'+str(a)+'\t'+str(aS)+'\t'+str(nrs) for (s,t,a,aS,nrs) in zip(sigPeptidesByCondition,totPeptidesByCondition,aonPeptidesByCondition,allsigPeptidesByCondition,nonrefoldabilityScoreByCondition)) + '\n' )
protein_summary_file.close()
peptideLip_output_file.close()

#########################################################################
peptideLip_output_file = open( peptideLip_output_fn, 'r' )
protein_summary_fn = peptideLip_fn.split('.')[0] + '_summary18ProteinNoNorm.txt'
protein_summary_file = open( protein_summary_fn, 'w')

GeneName = ''
Pvalues = [[] for c in range(numConditions-1)]
Ratios = [[] for c in range(numConditions-1)]
next(peptideLip_output_file)

protein_summary_file.write( 'GeneName\tLocus Tag\tAccession\tLocations\tComplexes\tSubunitCounts\tCopy Number\tMolecular Weight\tpI\tpSup\tdpSup_30_42\tdpSup_42_46\tArchitectures\tXgroups\tTgroups\tDomainRanges\tNumDomains\tPercent within Domain\tPerecent Disordered\tDisordered Regions\tOrdered Regions\tCCT sites\tHsp70 sites\tSG MayorLab\tSG ParkerLab\tSG BukauLab\tSG DrummondLab\t' + '\t'.join('SigPeptides'+str(s)+'\t'+'TotPeptides'+str(s)+'\t'+'AonPeptides'+str(s)+'\t'+'AllSigPeptides'+str(s)+'\t'+'NonrefoldabilityScore'+str(s) for s in range(1,numConditions)) + '\n' )
for line in peptideLip_output_file:
    info_in_line = line.split('\t')
    if info_in_line[0] == GeneName: #we're on the same gene
        for c in range(numConditions-1):
            # edited values to test, revert to ensure consideration of multiple conditions
            # Ratios[c].append( float(info_in_line[19 + numConditions + c] ) ) # normalized
            # Pvalues[c].append( float(info_in_line[20 + 2*numConditions + c] ) )

            Ratios[c].append( float(info_in_line[31] ) ) # not normalized
            Pvalues[c].append( float(info_in_line[33] ) )
    else: #we've started looking at a new gene
        if GeneName == '': #We've just started; this was the first gene
            GeneName = info_in_line[0]
            YeastLocusTag = info_in_line[1]
            Accession = info_in_line[2]
            Location = info_in_line[5]
            Complexes = info_in_line[6]
            SubunitCounts = info_in_line[7]
            CopyNumber = info_in_line[8]
            ProteinMW = info_in_line[9]
            pI = info_in_line[10]
            pSup = info_in_line[11]
            ProteinDeltapSup_30_42 = info_in_line[12]
            ProteinDeltapSup_42_46 = info_in_line[13]
            ProteinDomainsArch = info_in_line[14]
            ProteinDomainsX = info_in_line[15]
            ProteinDomainsT = info_in_line[16]
            ProteinDomainRanges = info_in_line[17]
            ProteinNumDomains = info_in_line[18]
            ProteinPercentWithinDomain = info_in_line[19]
            ProteinPerecentDisordered = info_in_line[20]
            ProteinDisorderedRanges = info_in_line[21]
            ProteinOrderedRanges = info_in_line[22]
            ProteinCCTsites = info_in_line[23]
            ProteinHsp70sites = info_in_line[24]
            ProteinSGMayor = info_in_line[25]
            ProteinSGParker = info_in_line[26]
            ProteinSGBukau = info_in_line[27]
            ProteinSGDrummond = info_in_line[28]

            for c in range(numConditions-1):
                # edited values to test, revert to ensure consideration of multiple conditions
                # Ratios[c].append( float(info_in_line[19 + numConditions + c] ) ) # normalized
                # Pvalues[c].append( float(info_in_line[20 + 2*numConditions + c] ) )

                Ratios[c].append( float(info_in_line[31] ) ) # not normalized
                Pvalues[c].append( float(info_in_line[33] ) )
        else: #We've just started looking at a new gene
            sigPeptidesByCondition = []
            totPeptidesByCondition = []
            aonPeptidesByCondition = [] #All or Nothing Peptides, those which are only present (or completely not present) in the refolded sample
            allsigPeptidesByCondition = []
            nonrefoldabilityScoreByCondition = []
            for c in range(numConditions-1):
                sigPeptides = 0
                totPeptides = 0
                aonPeptides = 0
                allsigPeptides = 0
                nonrefoldabilityScore = float(0)
                for (Ratio,Pval) in zip(Ratios[c],Pvalues[c]):
                    if Ratio != 0 and Pval != 0:
                        totPeptides += 1
                        if (6 > np.abs(Ratio) > 1 and Pval > 2):
                            sigPeptides += 1
                            allsigPeptides += 1
                            nonrefoldabilityScore += np.abs(Ratio)
                        if np.abs(Ratio) > 6 and Pval > 1.8:
                            aonPeptides += 1
                            allsigPeptides += 1
                            nonrefoldabilityScore += np.abs(Ratio)
                sigPeptidesByCondition.append( sigPeptides )
                totPeptidesByCondition.append( totPeptides )
                aonPeptidesByCondition.append( aonPeptides )
                allsigPeptidesByCondition.append( allsigPeptides )
                nonrefoldabilityScoreByCondition.append( nonrefoldabilityScore / totPeptides )
            Pvalues = [[] for c in range(numConditions-1)]
            Ratios = [[] for c in range(numConditions-1)]

            protein_summary_file.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t' %(GeneName, YeastLocusTag, Accession, Location, Complexes, SubunitCounts, CopyNumber, ProteinMW, pI, pSup, ProteinDeltapSup_30_42, ProteinDeltapSup_42_46, ProteinDomainsArch, ProteinDomainsX, ProteinDomainsT, ProteinDomainRanges, ProteinNumDomains, ProteinPercentWithinDomain, ProteinPerecentDisordered, ProteinDisorderedRanges, ProteinOrderedRanges, ProteinCCTsites, ProteinHsp70sites, ProteinSGMayor, ProteinSGParker, ProteinSGBukau, ProteinSGDrummond) +  '\t'.join(str(s)+'\t'+str(t)+'\t'+str(a)+'\t'+str(aS)+'\t'+str(nrs) for (s,t,a,aS,nrs) in zip(sigPeptidesByCondition,totPeptidesByCondition,aonPeptidesByCondition,allsigPeptidesByCondition,nonrefoldabilityScoreByCondition)) + '\n' )

            GeneName = info_in_line[0]
            YeastLocusTag = info_in_line[1]
            Accession = info_in_line[2]
            Location = info_in_line[5]
            Complexes = info_in_line[6]
            SubunitCounts = info_in_line[7]
            CopyNumber = info_in_line[8]
            ProteinMW = info_in_line[9]
            pI = info_in_line[10]
            pSup = info_in_line[11]
            ProteinDeltapSup_30_42 = info_in_line[12]
            ProteinDeltapSup_42_46 = info_in_line[13]
            ProteinDomainsArch = info_in_line[14]
            ProteinDomainsX = info_in_line[15]
            ProteinDomainsT = info_in_line[16]
            ProteinDomainRanges = info_in_line[17]
            ProteinNumDomains = info_in_line[18]
            ProteinPercentWithinDomain = info_in_line[19]
            ProteinPerecentDisordered = info_in_line[20]
            ProteinDisorderedRanges = info_in_line[21]
            ProteinOrderedRanges = info_in_line[22]
            ProteinCCTsites = info_in_line[23]
            ProteinHsp70sites = info_in_line[24]
            ProteinSGMayor = info_in_line[25]
            ProteinSGParker = info_in_line[26]
            ProteinSGBukau = info_in_line[27]
            ProteinSGDrummond = info_in_line[28]
            #DomainLocation = info_in_line[18]
            for c in range(numConditions-1):
                # edited values to test, revert to ensure consideration of multiple conditions
                # Ratios[c].append( float(info_in_line[19 + numConditions + c] ) )
                # Pvalues[c].append( float(info_in_line[20 + 2*numConditions + c] ) )

                Ratios[c].append( float(info_in_line[31] ) )
                Pvalues[c].append( float(info_in_line[33] ) )

sigPeptidesByCondition = []
totPeptidesByCondition = []
aonPeptidesByCondition = []
allsigPeptidesByCondition = []
nonrefoldabilityScoreByCondition = []
for c in range(numConditions-1):
    sigPeptides = 0
    totPeptides = 0
    aonPeptides = 0
    allsigPeptides = 0
    nonrefoldabilityScore = float(0)
    for (Ratio,Pval) in zip(Ratios[c],Pvalues[c]):
        if Ratio != 0 and Pval != 0:
            totPeptides += 1
            if (6 > np.abs(Ratio) > 1 and Pval > 2):
                sigPeptides += 1
                allsigPeptides += 1
                nonrefoldabilityScore += np.abs(Ratio)
            if np.abs(Ratio) > 6 and Pval > 1.8:
                aonPeptides += 1
                allsigPeptides += 1
                nonrefoldabilityScore += np.abs(Ratio)
    sigPeptidesByCondition.append( sigPeptides )
    totPeptidesByCondition.append( totPeptides )
    aonPeptidesByCondition.append( aonPeptides )
    allsigPeptidesByCondition.append( allsigPeptides )
    nonrefoldabilityScoreByCondition.append( nonrefoldabilityScore / totPeptides )

protein_summary_file.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t' %(GeneName, YeastLocusTag, Accession, Location, Complexes, SubunitCounts, CopyNumber, ProteinMW, pI, pSup, ProteinDeltapSup_30_42, ProteinDeltapSup_42_46, ProteinDomainsArch, ProteinDomainsX, ProteinDomainsT, ProteinDomainRanges, ProteinNumDomains, ProteinPercentWithinDomain, ProteinPerecentDisordered, ProteinDisorderedRanges, ProteinOrderedRanges, ProteinCCTsites, ProteinHsp70sites, ProteinSGMayor, ProteinSGParker, ProteinSGBukau, ProteinSGDrummond) +  '\t'.join(str(s)+'\t'+str(t)+'\t'+str(a)+'\t'+str(aS)+'\t'+str(nrs) for (s,t,a,aS,nrs) in zip(sigPeptidesByCondition,totPeptidesByCondition,aonPeptidesByCondition,allsigPeptidesByCondition,nonrefoldabilityScoreByCondition)) + '\n' )
protein_summary_file.close()
peptideLip_output_file.close()
######################################################


#Sort the peptide file on the Domain Location so that way we can use the same ordering feature to make a domain-level output
peptideLip_output_df = pd.read_csv( peptideLip_output_fn , sep='\t', header='infer')
peptideLip_output_df = peptideLip_output_df.sort_values(by=['GeneName','Domain Location','NormPeptideRatio1'])
peptideLip_output_df.to_csv( peptideLip_output_fn , sep='\t', index=False ) # Change from False to True to keep GeneName, otherwise was lost and indices were shifted

peptideLip_output_file = open( peptideLip_output_fn, 'r' )
domain_summary_fn = peptideLip_fn.split('.')[0] + '_summary18Domain.txt'
domain_summary_file = open( domain_summary_fn, 'w')

DomainLocation = ''
Pvalues = [[] for c in range(numConditions-1)]
Ratios = [[] for c in range(numConditions-1)]
next(peptideLip_output_file)

# domain_summary_file.write( 'Accession\tDomain Location\tDomain Tgroup\tDomain Xgroup\tDomain Arch\tDomain Order\t' + '\t'.join('SigPeptides'+str(s)+'\t'+'TotPeptides'+str(s)+'\t'+'AonPeptides'+str(s)+'\t'+'AllSigPeptides'+str(s) for s in range(1,numConditions)) + '\n' )
domain_summary_file.write( 'Accession\tDomain Location\tDomain Arch\tDomain Xgroup\tDomain Tgroup\tDomain Order\t' + '\t'.join('SigPeptides'+str(s)+'\t'+'TotPeptides'+str(s)+'\t'+'AonPeptides'+str(s)+'\t'+'AllSigPeptides'+str(s) for s in range(1,numConditions)) + '\n' )
for line in peptideLip_output_file:
    info_in_line = line.split('\t')
    if info_in_line[29] == DomainLocation: #we're on the same domain location
        for c in range(numConditions-1):
            # Ratios[c].append( float(info_in_line[21 + numConditions + c] ) )
            # Pvalues[c].append( float(info_in_line[20 + 2*numConditions + c] ) )

            Ratios[c].append( float(info_in_line[32] ) )
            Pvalues[c].append( float(info_in_line[33] ) )
    else: #we've started looking at a new gene
        if DomainLocation == '': #We've just started; this was the first gene
            Accession = info_in_line[2]
            GeneName = info_in_line[0]
            DomainLocation = info_in_line[29]
            ProteinDomainsArch =  info_in_line[14]
            ProteinDomainsX = info_in_line[15]
            ProteinDomainsT = info_in_line[16]
            if ProteinDomainsT != '':
                ProteinDomainsArch = ast.literal_eval( info_in_line[14] )
                ProteinDomainsX = ast.literal_eval( info_in_line[15] )
                ProteinDomainsT = ast.literal_eval( info_in_line[16] )
                ProteinDomainRanges = ast.literal_eval( info_in_line[17] )
                if type(ProteinDomainsT) != list:
                    try:
                        ProteinDomainsT = ast.literal_eval( ProteinDomainsT )
                    except SyntaxError:
                        ProteinDomainsT = [ 'domain'+str(x) for x in range(len(ProteinDomainRanges))] #for one edge case, ushA, where the quotation remarks refuse to get formatted into a list no matter what I try!
                #Find the particular fold and order of domain (in sequence)
                if ',' in GeneName:
                    domainRange = DomainLocation.split(',')[2]
                else:
                    domainRange = DomainLocation.split(",")[1]
                if domainRange != '' and domainRange != 'LNK':
                    try:
                        thisDomainArch = ProteinDomainsArch[ ProteinDomainRanges.index(domainRange) ]
                        thisDomainX = ProteinDomainsX[ ProteinDomainRanges.index(domainRange) ]
                        thisDomainT = ProteinDomainsT[ ProteinDomainRanges.index(domainRange) ]
                        ProteinDomainStartPoints = [ int(z.split('-')[0])  for z in ProteinDomainRanges ]
                        orderings = rankdata( ProteinDomainStartPoints )
                        thisOrder = orderings[ ProteinDomainRanges.index(domainRange) ]
                    except ValueError:
                        thisDomainArch = ProteinDomainsArch[0]
                        thisDomainX = ProteinDomainsX[0]
                        thisDomainT = ProteinDomainsT[0]
                        thisOrder = 1
                elif domainRange == 'LNK':
                    thisDomainArch = 'Linkers'
                    thisDomainX = 'Linkers'
                    thisDomainT = 'Linkers'
                    thisOrder = 'Linkers'
            else:
                thisDomainArch = 'None'
                thisDomainX = 'None'
                thisDomainT = 'None'
                thisOrder = 'None'
            for c in range(numConditions-1):
                # Ratios[c].append( float(info_in_line[21 + numConditions + c] ) )
                # Pvalues[c].append( float(info_in_line[20 + 2*numConditions + c] ) )

                Ratios[c].append( float(info_in_line[32] ) )
                Pvalues[c].append( float(info_in_line[33] ) )
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

            domain_summary_file.write('%s\t%s\t%s\t%s\t%s\t%s\t' %(Accession,DomainLocation,thisDomainArch,thisDomainX,thisDomainT,thisOrder) +  '\t'.join(str(s)+'\t'+str(t)+'\t'+str(a)+'\t'+str(aS) for (s,t,a,aS) in zip(sigPeptidesByCondition,totPeptidesByCondition,aonPeptidesByCondition,allsigPeptidesByCondition)) + '\n' )

            Accession = info_in_line[2]
            GeneName = info_in_line[0]
            DomainLocation = info_in_line[29]
            ProteinDomainsArch =  info_in_line[14]
            ProteinDomainsX = info_in_line[15]
            ProteinDomainsT = info_in_line[16]
            if ProteinDomainsT != '':
                ProteinDomainsArch = ast.literal_eval( info_in_line[14] )
                ProteinDomainsX = ast.literal_eval( info_in_line[15] )
                ProteinDomainsT = ast.literal_eval( info_in_line[16] )
                ProteinDomainRanges = ast.literal_eval( info_in_line[17] )
                if type(ProteinDomainsT) != list:
                    try:
                        ProteinDomainsT = ast.literal_eval( ProteinDomainsT )
                    except SyntaxError:
                        ProteinDomainsT = [ 'domain'+str(x) for x in range(len(ProteinDomainRanges))] #for one edge case, ushA, where the quotation remarks refuse to get formatted into a list no matter what I try!
                #Find the particular fold and order of domain (in sequence)
                if ',' in GeneName:
                    domainRange = DomainLocation.split(',')[2]
                else:
                    domainRange = DomainLocation.split(",")[1]
                if domainRange != '' and domainRange != 'LNK':
                    try:
                        thisDomainArch = ProteinDomainsArch[ ProteinDomainRanges.index(domainRange) ]
                        thisDomainX = ProteinDomainsX[ ProteinDomainRanges.index(domainRange) ]
                        thisDomainT = ProteinDomainsT[ ProteinDomainRanges.index(domainRange) ]
                        ProteinDomainStartPoints = [ int(z.split('-')[0])  for z in ProteinDomainRanges ]
                        orderings = rankdata( ProteinDomainStartPoints )
                        thisOrder = orderings[ ProteinDomainRanges.index(domainRange) ]
                    except ValueError:
                        thisDomainArch = ProteinDomainsArch[0]
                        thisDomainX = ProteinDomainsX[0]
                        thisDomainT = ProteinDomainsT[0]
                        thisOrder = 1
                elif domainRange == 'LNK':
                    thisDomainArch = 'Linkers'
                    thisDomainX = 'Linkers'
                    thisDomainT = 'Linkers'
                    thisOrder = 'Linkers'
            else:
                thisDomainArch = 'None'
                thisDomainX = 'None'
                thisDomainT = 'None'
                thisOrder = 'None'
            for c in range(numConditions-1):
                # Ratios[c].append( float(info_in_line[21 + numConditions + c] ) )
                # Pvalues[c].append( float(info_in_line[20 + 2*numConditions + c] ) )

                Ratios[c].append( float(info_in_line[32] ) )
                Pvalues[c].append( float(info_in_line[33] ) )

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

domain_summary_file.write('%s\t%s\t%s\t%s\t%s\t%s\t' %(Accession,DomainLocation,thisDomainArch,thisDomainX,thisDomainT,thisOrder) +  '\t'.join(str(s)+'\t'+str(t)+'\t'+str(a)+'\t'+str(aS) for (s,t,a,aS) in zip(sigPeptidesByCondition,totPeptidesByCondition,aonPeptidesByCondition,allsigPeptidesByCondition)) + '\n' )

domain_summary_file.close()
peptideLip_output_file.close()


##############################################################
peptideLip_output_file = open( peptideLip_output_fn, 'r' )
domain_summary_fn = peptideLip_fn.split('.')[0] + '_summary18DomainNoNorm.txt'
domain_summary_file = open( domain_summary_fn, 'w')

DomainLocation = ''
Pvalues = [[] for c in range(numConditions-1)]
Ratios = [[] for c in range(numConditions-1)]
next(peptideLip_output_file)

# domain_summary_file.write( 'Accession\tDomain Location\tDomain Tgroup\tDomain Xgroup\tDomain Arch\tDomain Order\t' + '\t'.join('SigPeptides'+str(s)+'\t'+'TotPeptides'+str(s)+'\t'+'AonPeptides'+str(s)+'\t'+'AllSigPeptides'+str(s) for s in range(1,numConditions)) + '\n' )
domain_summary_file.write( 'Accession\tDomain Location\tDomain Arch\tDomain Xgroup\tDomain Tgroup\tDomain Order\t' + '\t'.join('SigPeptides'+str(s)+'\t'+'TotPeptides'+str(s)+'\t'+'AonPeptides'+str(s)+'\t'+'AllSigPeptides'+str(s) for s in range(1,numConditions)) + '\n' )
for line in peptideLip_output_file:
    info_in_line = line.split('\t')
    if info_in_line[29] == DomainLocation: #we're on the same domain location
        for c in range(numConditions-1):
            # Ratios[c].append( float(info_in_line[21 + numConditions + c] ) )
            # Pvalues[c].append( float(info_in_line[20 + 2*numConditions + c] ) )

            Ratios[c].append( float(info_in_line[31] ) )
            Pvalues[c].append( float(info_in_line[33] ) )
    else: #we've started looking at a new gene
        if DomainLocation == '': #We've just started; this was the first gene
            Accession = info_in_line[2]
            GeneName = info_in_line[0]
            DomainLocation = info_in_line[29]
            ProteinDomainsArch =  info_in_line[14]
            ProteinDomainsX = info_in_line[15]
            ProteinDomainsT = info_in_line[16]
            if ProteinDomainsT != '':
                ProteinDomainsArch = ast.literal_eval( info_in_line[14] )
                ProteinDomainsX = ast.literal_eval( info_in_line[15] )
                ProteinDomainsT = ast.literal_eval( info_in_line[16] )
                ProteinDomainRanges = ast.literal_eval( info_in_line[17] )
                if type(ProteinDomainsT) != list:
                    try:
                        ProteinDomainsT = ast.literal_eval( ProteinDomainsT )
                    except SyntaxError:
                        ProteinDomainsT = [ 'domain'+str(x) for x in range(len(ProteinDomainRanges))] #for one edge case, ushA, where the quotation remarks refuse to get formatted into a list no matter what I try!
                #Find the particular fold and order of domain (in sequence)
                if ',' in GeneName:
                    domainRange = DomainLocation.split(',')[2]
                else:
                    domainRange = DomainLocation.split(",")[1]
                if domainRange != '' and domainRange != 'LNK':
                    try:
                        thisDomainArch = ProteinDomainsArch[ ProteinDomainRanges.index(domainRange) ]
                        thisDomainX = ProteinDomainsX[ ProteinDomainRanges.index(domainRange) ]
                        thisDomainT = ProteinDomainsT[ ProteinDomainRanges.index(domainRange) ]
                        ProteinDomainStartPoints = [ int(z.split('-')[0])  for z in ProteinDomainRanges ]
                        orderings = rankdata( ProteinDomainStartPoints )
                        thisOrder = orderings[ ProteinDomainRanges.index(domainRange) ]
                    except ValueError:
                        thisDomainArch = ProteinDomainsArch[0]
                        thisDomainX = ProteinDomainsX[0]
                        thisDomainT = ProteinDomainsT[0]
                        thisOrder = 1
                elif domainRange == 'LNK':
                    thisDomainArch = 'Linkers'
                    thisDomainX = 'Linkers'
                    thisDomainT = 'Linkers'
                    thisOrder = 'Linkers'
            else:
                thisDomainArch = 'None'
                thisDomainX = 'None'
                thisDomainT = 'None'
                thisOrder = 'None'
            for c in range(numConditions-1):
                # Ratios[c].append( float(info_in_line[21 + numConditions + c] ) )
                # Pvalues[c].append( float(info_in_line[20 + 2*numConditions + c] ) )

                Ratios[c].append( float(info_in_line[31] ) )
                Pvalues[c].append( float(info_in_line[33] ) )
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

            domain_summary_file.write('%s\t%s\t%s\t%s\t%s\t%s\t' %(Accession,DomainLocation,thisDomainArch,thisDomainX,thisDomainT,thisOrder) +  '\t'.join(str(s)+'\t'+str(t)+'\t'+str(a)+'\t'+str(aS) for (s,t,a,aS) in zip(sigPeptidesByCondition,totPeptidesByCondition,aonPeptidesByCondition,allsigPeptidesByCondition)) + '\n' )

            Accession = info_in_line[2]
            GeneName = info_in_line[0]
            DomainLocation = info_in_line[29]
            ProteinDomainsArch =  info_in_line[14]
            ProteinDomainsX = info_in_line[15]
            ProteinDomainsT = info_in_line[16]
            if ProteinDomainsT != '':
                ProteinDomainsArch = ast.literal_eval( info_in_line[14] )
                ProteinDomainsX = ast.literal_eval( info_in_line[15] )
                ProteinDomainsT = ast.literal_eval( info_in_line[16] )
                ProteinDomainRanges = ast.literal_eval( info_in_line[17] )
                if type(ProteinDomainsT) != list:
                    try:
                        ProteinDomainsT = ast.literal_eval( ProteinDomainsT )
                    except SyntaxError:
                        ProteinDomainsT = [ 'domain'+str(x) for x in range(len(ProteinDomainRanges))] #for one edge case, ushA, where the quotation remarks refuse to get formatted into a list no matter what I try!
                #Find the particular fold and order of domain (in sequence)
                if ',' in GeneName:
                    domainRange = DomainLocation.split(',')[2]
                else:
                    domainRange = DomainLocation.split(",")[1]
                if domainRange != '' and domainRange != 'LNK':
                    try:
                        thisDomainArch = ProteinDomainsArch[ ProteinDomainRanges.index(domainRange) ]
                        thisDomainX = ProteinDomainsX[ ProteinDomainRanges.index(domainRange) ]
                        thisDomainT = ProteinDomainsT[ ProteinDomainRanges.index(domainRange) ]
                        ProteinDomainStartPoints = [ int(z.split('-')[0])  for z in ProteinDomainRanges ]
                        orderings = rankdata( ProteinDomainStartPoints )
                        thisOrder = orderings[ ProteinDomainRanges.index(domainRange) ]
                    except ValueError:
                        thisDomainArch = ProteinDomainsArch[0]
                        thisDomainX = ProteinDomainsX[0]
                        thisDomainT = ProteinDomainsT[0]
                        thisOrder = 1
                elif domainRange == 'LNK':
                    thisDomainArch = 'Linkers'
                    thisDomainX = 'Linkers'
                    thisDomainT = 'Linkers'
                    thisOrder = 'Linkers'
            else:
                thisDomainArch = 'None'
                thisDomainX = 'None'
                thisDomainT = 'None'
                thisOrder = 'None'
            for c in range(numConditions-1):
                # Ratios[c].append( float(info_in_line[21 + numConditions + c] ) )
                # Pvalues[c].append( float(info_in_line[20 + 2*numConditions + c] ) )

                Ratios[c].append( float(info_in_line[31] ) )
                Pvalues[c].append( float(info_in_line[33] ) )

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

domain_summary_file.write('%s\t%s\t%s\t%s\t%s\t%s\t' %(Accession,DomainLocation,thisDomainArch,thisDomainX,thisDomainT,thisOrder) +  '\t'.join(str(s)+'\t'+str(t)+'\t'+str(a)+'\t'+str(aS) for (s,t,a,aS) in zip(sigPeptidesByCondition,totPeptidesByCondition,aonPeptidesByCondition,allsigPeptidesByCondition)) + '\n' )

domain_summary_file.close()
peptideLip_output_file.close()
