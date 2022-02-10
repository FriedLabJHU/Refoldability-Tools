import Bio
from Bio import SeqIO   
from Bio import Seq
import pandas as pd
from numpy import nan
import math
import re 
from Bio.SeqUtils import molecular_weight
import requests
from bs4 import BeautifulSoup

### Create dictionary keyed on accession number from original fasta sequence with pI. Specifically want accession number, gene name, sequence, and final avg pI value. 
input_filename_step = 'uniprot-proteome-UP000000592_iso.fasta'
input_file_step = open( input_filename_step , 'r' )

thermusECOD = pd.read_csv("20220120_Thermus_Domains.txt", index_col=False, delimiter= "\t", skiprows=31, names = ["Identifier", "E-Value", "Residue Range", "Property", "Architecture", "X-Group", "T-Group", "F-Group", "F-ID"])
### Okay why is it not printing out the right stuff????

thermusHomologues = pd.read_csv("20211213_thermus_homologues_updated.txt", index_col=False, delimiter = "\t")
metadata_dict={}
with open("uniprot-proteome-UP000000592_iso.fasta") as handle: 
    for record in SeqIO.parse(handle, "fasta"): 
        Accession = record.id.split('|')[1]
        FASTAgene = record.description.split('=')[3].split()[0] 
        FASTAseq = record.seq   
        pIaverage = record.description.split('=')[5].split(',')[-1]
        MW=molecular_weight(FASTAseq, "protein")
        Length=len(FASTAseq)
        metadata_dict[Accession]={}
        metadata_dict[Accession]["Gene Name"]=FASTAgene
        metadata_dict[Accession]["Sequence"]=FASTAseq
        metadata_dict[Accession]["pI"]=pIaverage
        metadata_dict[Accession]["Molecular Weight"]=MW
        metadata_dict[Accession]["Length"]=Length
        annotatedRegions = 0
        proteinECODinfo = thermusECOD[thermusECOD["Identifier"].str.contains(str(Accession))]
        metadata_dict[Accession]["Domain E-Value"]=(list(proteinECODinfo["E-Value"]))
        metadata_dict[Accession]["Residue Range"]=(list(proteinECODinfo["Residue Range"]))
        metadata_dict[Accession]["Residue Range"] = [re.sub("," , "+", residueRange) for residueRange in metadata_dict[Accession]["Residue Range"]]
        metadata_dict[Accession]["Architecture"]=(list(proteinECODinfo["Architecture"]))
        metadata_dict[Accession]["X-Group"]=(list(proteinECODinfo["X-Group"]))
        metadata_dict[Accession]["T-Group"]=(list(proteinECODinfo["T-Group"]))
        metadata_dict[Accession]["F-Group"]=(list(proteinECODinfo["F-Group"]))
        metadata_dict[Accession]["Number of Domains"]=len(list(proteinECODinfo["X-Group"]))
        for range in list(proteinECODinfo["Residue Range"]):
            range2 = range.replace("," , "+")
            annotatedRegions += eval(range2)
            percentECOD = abs(annotatedRegions)/metadata_dict[Accession]["Length"]
            metadata_dict[Accession]["Percent within Domain"]=percentECOD
        
        proteinHit = thermusHomologues.loc[thermusHomologues["Unnamed: 0"] == Accession]
        
        if Accession in list(proteinHit["Unnamed: 0"]):
            metadata_dict[Accession]["E. coli Homolog Accession"] = list(proteinHit["E. coli Homolog Accession"])[0]
            metadata_dict[Accession]["E. coli Homolog Gene Name"] = list(proteinHit["E. coli Homolog Gene Name"])[0]
            metadata_dict[Accession]["Hmmer E-Value"] = list(proteinHit["E-Value"])[0]
        else:
            metadata_dict[Accession]["E. coli Homolog Accession"] = "N/A"
            metadata_dict[Accession]["E. coli Homolog Gene Name"] = "N/A"
            metadata_dict[Accession]["Hmmer E-Value"] = "N/A"
#Now let's make a dictionary from our E. coli scrape
# #print(metadata_dict)
eColi_scrape_file = open('20220124_GeneSubLocSyn_Total.txt' , encoding="utf8")

next(eColi_scrape_file) #skip the header line
eColi_scrape_dict = {}
for i,line in enumerate(eColi_scrape_file):

    line = line.replace(u'\xa0', u'')
    line = line.split("\t")
    line[-1] = line[-1].strip()

    ecGeneName, ecAccession, Frame_id, ecSynonyms, ecLocation, ecCofactor, ecComplex, ecOverallSubunits, ecUniqueSubunits = line

    eColi_scrape_dict[ecAccession]={}
    eColi_scrape_dict[ecAccession]["E Coli Gene Name"]=ecGeneName
    eColi_scrape_dict[ecAccession]["E Coli Synonyms"]=ecSynonyms
    eColi_scrape_dict[ecAccession]["E Coli Location"]=ecLocation
    eColi_scrape_dict[ecAccession]["E Coli Cofactors"]=ecCofactor
    eColi_scrape_dict[ecAccession]["E Coli Complex"]=ecComplex
    eColi_scrape_dict[ecAccession]["E Coli Overall Subunits"]=ecOverallSubunits
    eColi_scrape_dict[ecAccession]["E Coli Unique Subunits"]=ecUniqueSubunits
    
eColi_scrape_file.close()

#Matchup Time
for Accession in metadata_dict.keys():
    # matching with same gene name
    try:
        gene_name = metadata_dict[Accession]["Gene Name"]
        for ecAccession in eColi_scrape_dict.keys():
            ec_gene_name = eColi_scrape_dict[ecAccession]["E Coli Gene Name"]
            if ec_gene_name == gene_name:
                metadata_dict[Accession]["E. Coli Accession Number"] = ecAccession
                metadata_dict[Accession]["E. Coli Gene Name"] = eColi_scrape_dict[ecAccession]["E Coli Gene Name"]
                metadata_dict[Accession]["E. Coli Synonyms"] = eColi_scrape_dict[ecAccession]["E Coli Synonyms"]
                metadata_dict[Accession]["E. Coli Location"] = eColi_scrape_dict[ecAccession]["E Coli Location"]
                metadata_dict[Accession]["E. Coli Cofactors"]= eColi_scrape_dict[ecAccession]["E Coli Cofactors"]
                metadata_dict[Accession]["E. Coli Complex"] = eColi_scrape_dict[ecAccession]["E Coli Complex"]
                metadata_dict[Accession]["E. Coli Overall Subunits"]= eColi_scrape_dict[ecAccession]["E Coli Overall Subunits"]
                metadata_dict[Accession]["E. Coli Unique Subunits"]= eColi_scrape_dict[ecAccession]["E Coli Unique Subunits"]
                break
    except:
        pass

    # matching with same synonyms
    try:
        gene_name = metadata_dict[Accession]["Gene Name"]
        for ecAccession in eColi_scrape_dict.keys():
            ec_synonym_names = eColi_scrape_dict[ecAccession]["E Coli Synonyms"]
            if gene_name in ec_synonym_names:
                metadata_dict[Accession]["E. Coli Accession Number"] = ecAccession
                metadata_dict[Accession]["E. Coli Gene Name"] = eColi_scrape_dict[ecAccession]["E Coli Gene Name"]
                metadata_dict[Accession]["E. Coli Synonyms"] = eColi_scrape_dict[ecAccession]["E Coli Synonyms"]
                metadata_dict[Accession]["E. Coli Location"] = eColi_scrape_dict[ecAccession]["E Coli Location"]
                metadata_dict[Accession]["E. Coli Cofactors"]= eColi_scrape_dict[ecAccession]["E Coli Cofactors"]
                metadata_dict[Accession]["E. Coli Complex"] = eColi_scrape_dict[ecAccession]["E Coli Complex"]
                metadata_dict[Accession]["E. Coli Overall Subunits"]= eColi_scrape_dict[ecAccession]["E Coli Overall Subunits"]
                metadata_dict[Accession]["E. Coli Unique Subunits"]= eColi_scrape_dict[ecAccession]["E Coli Unique Subunits"]
                break
    except:
        pass

    # if all else fails... 
    # matching with same accession number // redundant to matching with same gene name
    if float(metadata_dict[Accession]["Hmmer E-Value"]) < 1e-5:
        pot_ecAccession = metadata_dict[Accession]["E. coli Homolog Accession"]
        try:
            if metadata_dict[Accession]["E. Coli Gene Name"] != '': 
                pass
        except KeyError: 
            if pot_ecAccession != "N/A" and pot_ecAccession in eColi_scrape_dict.keys():
                metadata_dict[Accession]["E. Coli Accession Number"] = pot_ecAccession
                metadata_dict[Accession]["E. Coli Gene Name"] = eColi_scrape_dict[pot_ecAccession]["E Coli Gene Name"]
                metadata_dict[Accession]["E. Coli Synonyms"] = eColi_scrape_dict[pot_ecAccession]["E Coli Synonyms"]
                metadata_dict[Accession]["E. Coli Location"] = eColi_scrape_dict[pot_ecAccession]["E Coli Location"]
                metadata_dict[Accession]["E. Coli Cofactors"]= eColi_scrape_dict[pot_ecAccession]["E Coli Cofactors"]
                metadata_dict[Accession]["E. Coli Complex"] = eColi_scrape_dict[pot_ecAccession]["E Coli Complex"]
                metadata_dict[Accession]["E. Coli Overall Subunits"]= eColi_scrape_dict[pot_ecAccession]["E Coli Overall Subunits"]
                metadata_dict[Accession]["E. Coli Unique Subunits"]= eColi_scrape_dict[pot_ecAccession]["E Coli Unique Subunits"]

### Otherwise we don't have enough metadata
with open("20220203_thermus_final_metadata.txt", 'w', encoding="utf8") as file: 
    file.write(('{}\t'*25 +"\n").format("Accession Number", "Gene Name", "Sequence", "pI", "Molecular Weight", "Length", "Domain E-Value", "Residue Range", "Architecture", "X-Group", "T-Group", "F-Group", "Number of Domains", "Percent within Domain", "E. Coli Homolog Accession Number", "E. Coli Homolog Gene Name", "Hmmer E-Value", "E. Coli Scrape Accession Number", "E. Coli Scrape Gene Name", "E. Coli Scrape Synonyms", "E. Coli Scrape Locations", "E. Coli Scrape Cofactors", "E. Coli Scrape Complexes","E. Coli Scrape Overall Subunits","E. Coli Scrape Unique Subunits"))
    for item in metadata_dict.keys():
        gene = metadata_dict[item]["Gene Name"] 
        FASTAseq = metadata_dict[item]["Sequence"]  
        pI = metadata_dict[item]["pI"] 
        MW=molecular=metadata_dict[item]["Molecular Weight"]
        Length=metadata_dict[item]["Length"]
        domainE=metadata_dict[item]["Domain E-Value"] 
        resrange=metadata_dict[item]["Residue Range"]    
        arch=metadata_dict[item]["Architecture"]
        xgro=metadata_dict[item]["X-Group"]
        tgro=metadata_dict[item]["T-Group"]
        fgro=metadata_dict[item]["F-Group"]
        numdom=metadata_dict[item]["Number of Domains"]
        ecolihomacc=metadata_dict[item]["E. coli Homolog Accession"]
        ecolihomgene=metadata_dict[item]["E. coli Homolog Gene Name"]
        hmmerE=metadata_dict[item]["Hmmer E-Value"]

        try:
            percdom=metadata_dict[item]["Percent within Domain"]
        except: 
            percdom="0"
        try: 
            ecoliscrapeacc=metadata_dict[item]["E. Coli Accession Number"]
        except:
            ecoliscrapeacc="N/A"
        try:
            ecoliscrapegene=metadata_dict[item]["E. Coli Gene Name"]
        except: 
            ecoliscrapegene="N/A"
        try:
            ecoliscrapesyn=metadata_dict[item]["E. Coli Synonyms"]
        except:
            ecoliscrapesyn="N/A"
        try:
            ecoliscrapeloc=metadata_dict[item]["E. Coli Location"]
        except:
            ecoliscrapeloc="N/A"
        try:
            ecoliscrapecofactor= metadata_dict[item]["E. Coli Cofactors"]
        except:
            ecoliscrapecofactor="N/A"
        try:
            ecoliscrapecomplex=metadata_dict[item]["E. Coli Complex"]
        except:
            ecoliscrapecomplex="N/A"
        try:
            ecoliscrapeoverallsubunit=metadata_dict[item]["E. Coli Overall Subunits"]
        except:
            ecoliscrapeoverallsubunit="N/A"
        try:
            ecoliscrapeuniquesubunit=metadata_dict[item]["E. Coli Unique Subunits"]
        except:
            ecoliscrapeuniquesubunit="N/A"   

        file.write(('{}\t'*25).format(item, gene, FASTAseq, pI, MW, Length, domainE, resrange, arch, xgro, tgro, fgro, numdom, percdom, ecolihomacc, ecolihomgene, hmmerE, ecoliscrapeacc, ecoliscrapegene, ecoliscrapesyn, ecoliscrapeloc, ecoliscrapecofactor, ecoliscrapecomplex, ecoliscrapeoverallsubunit, ecoliscrapeuniquesubunit)+"\n")
