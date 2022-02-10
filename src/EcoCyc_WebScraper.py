#These are important for parsing the websites
import requests
from bs4 import BeautifulSoup
import re
import time
from requests.packages.urllib3.exceptions import InsecureRequestWarning
requests.packages.urllib3.disable_warnings(InsecureRequestWarning)

def fileWrite(writeStuff, genes, Accession, frameid, Synonyms, locations, cofactors, complexes, overall_subunits, unique_subunits): #Simple function to write into a file
    lineToWrite = ('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' %(genes, Accession, frameid, Synonyms, locations, cofactors, complexes, overall_subunits, unique_subunits))
    return writeStuff.write(lineToWrite)
def locationSyn(synLoc): #Returns the location or synonyms should there be any
    locSynText = synLoc[1].split('\n')
    return [x for x in locSynText if x != ''][0]
def cof(cof): #Helps identify the cofactors 
    letters = "\]+[A-Za-z]"
    cof = cof[1]
    cofText = re.split(letters, cof)
    return [x for x in cofText if x != ''][0]


#in_fn = '20192030_ProteinIDListMinusTheGook.txt'
in_fn = '20220124_ProteinIDList_Complete.txt'
#in_fn = '20192030_Last1500ishProteinIDList.txt'
in_file = open(in_fn, 'r')
#subunitFile = open('proteinCompositionPractice.txt', 'w')
subunitFile = open("20220124_GeneSubLocSyn_Total.txt", "w", encoding="utf-8")
#subunitFile = open("2019114_GeneSubLocSyn_CompleteList_v3.txt", "w") #Creates a file containing the gene name and subunits for each protein
#subunitFile = open("2019114_GeneSubLocSyn_ptIII.txt", "w")
subunitFile.write('Gene'+'\t'+'Accession'+'\t'+'Frame ID'+'\t'+'Synonyms'+'\t'+'Location'+'\t'+'Cofactors'+'\t'+'Complex'+'\t'+'Overall Subunits'+'\t'+'UniqueSubunits'+'\n')
for line_num, line in enumerate(in_file):
    #time.sleep(1.1)
    print(line_num)
    gene = line.split('\t')[0]
    frameid = line.split('\t')[1][:-1]
    url = 'https://ecocyc.org//gene-tab?id=' + frameid + '&orgid=ECOLI&tab=SUMMARY' #URL to find gene name and complex state
    url2 = 'https://ecocyc.org//gene?id=' + frameid #URL to determine alternate names and location in cell
    url3 = 'https://ecocyc.org/gene-tab?id='  + frameid + '&origid=ECOLI&tab=RXNS'   #URL to determine cofactors in cell
    site = requests.get(url, verify = False)
    time.sleep(1.1)
    site2 = requests.get(url2, verify = False)
    time.sleep(1.1)
    site3 = requests.get(url3, verify = False)
    time.sleep(1.1)
    content = BeautifulSoup(site.content, "html.parser")
    synLocContent = BeautifulSoup(site2.content, "html.parser")
    cofContent = BeautifulSoup(site3.content, "html.parser")
    comp = content.find(class_='componentOf') #Finds the CSS class 'componentOf' and returns either the HTML code of the class or None (if the class does not exist)
    accessionText = synLocContent.get_text().split('(UniProt)')
    #<+TEST+>
    #if line_num >= 142:
    #    print(line_num, accessionText)
    #    print(url,url2,url3,end="\n")
    #<+TEST+>
    synonymsText = synLocContent.get_text().split('Synonyms') #Splits the code by the word Synonyms
    synonymsTextEdge = synLocContent.get_text().split('Synonym') #Edge case where 'Synonym' is used in place of 'Synonyms' 
    locationText = synLocContent.get_text().split('Locations') #Splits the code by the the keyword Location
    locationTextEdge = synLocContent.get_text().split('Location') #Just in case they do the same thing with location they did with synonyms
    cofactors = cofContent.get_text().split('Cofactors or Prosthetic Groups') #Identify Cofactors

    #Returns the location and gene synonyms
    if len(synonymsText) == 1:
        if len(synonymsTextEdge) == 1:
            synonyms = 'None'
        else:
            synonyms = locationSyn(synonymsTextEdge)
    elif len(synonymsText) > 1:
        synonyms = locationSyn(synonymsText)
    if len(locationText) == 1:
        if len(locationTextEdge) == 1:
            location = 'Unknown'
        else:
            location = locationSyn(locationTextEdge)
    elif len(locationText) > 1:
        location = locationSyn(locationText)
    #Records Accession Number
    if len(accessionText) == 1:
        accession= 'NA'
    else:
        accession=accessionText[0][-7:]

    #Use this code to identify cofactors
    if len(cofactors) == 1:
        cofactors = 'None'
    else:
        cofactors = cof(cofactors)
        cofactors = cofactors + ']'
        cofactors = re.sub("\["+"[A-Za-z]"+"[^\]]*"+"[a-z0-9]+\]", '', cofactors)
        cofactors = re.split("\+" + "[A-Za-z]", cofactors)
        cofactors = cofactors[0]

    #This chunk of code finds the subunit compositions and complexes of the proteins
    counter = 0
    subunit = []
    for line in content.stripped_strings:
        if counter == 0:
            if line != "Subunit Composition":
                pass
            elif line == "Subunit Composition":
                counter = 1
        elif counter == 1:
            if line ==   "Component of": # The following elif lines make sure list only has the subunit composition
                counter = 0
            elif line == "Gene-Reaction Schematic":
                counter = 0
            elif line == "In Paralogous Gene Groups":
                counter = 0
            elif line == "In Paralogous Gene Group":
                counter = 0
            elif line == "Alternative form":
                counter = 0
            elif line == "Alternative forms":
                counter = 0
            else:
                subunit.append(line)

    #If there is no subunit composition, looks through component of to determine the protein complex. If component of does not exist, it is a monomer. 
    if not subunit:
        if comp == None:

            cistronic = 'No'
            overall_subunits = '1'
            unique_subunits  = '1'

            fileWrite(subunitFile, gene, accession, frameid, synonyms, location, cofactors, "Monomer", overall_subunits, unique_subunits)
            pass

        else:
            component = comp.find('dd').get_text()

            print (component)
            info_in_string = re.split(r'[([]', component)
            info_in_string = list(filter(None, info_in_string))

            print (info_in_string)

            info_in_string_reduced = []

            overall_subunits = 0

            for proteins in info_in_string:
                if re.match( r'[A-Z][a-z][a-z][A-Z]', proteins):
                    info_in_string_reduced.append(proteins)
                    split_p = re.split(r'[)\]]', proteins)
                    protein = split_p[0]
                    num = 1

                    print (split_p)
                    for int_ in split_p[1:]:
                        print ('4 letter match Int Value is: ' + str(int_))
                        if int_.isdigit():
                            num *= int(int_)
                        else:
                            pass

                elif re.match( r'[A-Z][a-z][a-z]', proteins):
                    info_in_string_reduced.append(proteins)
                    split_p = re.split(r'[)\]]', proteins)
                    protein = split_p[0]
                    num = 1

                    print (split_p)
                    for int_ in split_p[1:]:
                        print ('3 letter match Int Value is: ' + str(int_))
                        if int_.isdigit():
                            num *= int(int_)
                        else:
                            pass
                print(protein, num)
                overall_subunits += num

            print ('Overall Subunit is: ' + str(overall_subunits))

            if len(info_in_string_reduced) == 1:
                overall_subunits = '1'
                unique_subunits = '1'
                fileWrite(subunitFile, gene, accession, frameid, synonyms, location, cofactors, "Monomer", overall_subunits, unique_subunits)

            else:
                #This part of the code is used to identify of proteins are multi-cistronics by checking the first three letters of all items in a list and seeing if they are the same 
                counter = 0
                counter_max = len(info_in_string_reduced)
                print (info_in_string)
                print (info_in_string_reduced[0])

                if len(info_in_string_reduced[0]) > 3:
                    first_entry = info_in_string_reduced[0][:3]

                else:
                    first_entry = info_in_string_reduced[0]
                print (first_entry)

                for entry in info_in_string_reduced:
                    if entry[:3] == first_entry:
                        counter += 1
                    else:
                        pass

                if counter == counter_max:
                    cistronic = 'Yes'

                else:
                    cistronic = 'Double Check'

                print ("Used Component to determine Subunit Info" + str(info_in_string_reduced))
                print ('This protein is multi cistronic?' + cistronic)

                #This part of the code uses info_in_string_reduced to detemrine the number of unique subunits

                unique_subunits = str(len(info_in_string_reduced))
                print ("Number of Unique Subunits is: " + unique_subunits)
                print ("The overall subunit count is: " + str(overall_subunits))

                fileWrite(subunitFile, gene, accession, frameid, synonyms, location, cofactors, info_in_string_reduced, overall_subunits, unique_subunits)

    #Uses subunit to determine Complex if available. 
    else:
        if comp == None:
            cistronic = 'No'
            subunit = ''.join(subunit)
            print ('Subunit is: ' + str(subunit))

            info_in_string = re.split(r'[([]', subunit)
            info_in_string = list(filter(None, info_in_string))

            print (info_in_string)

            info_in_string_reduced = []

            overall_subunits = 0

            for proteins in info_in_string:
                if re.match( r'[A-Z][a-z][a-z][A-Z]', proteins):
                    info_in_string_reduced.append(proteins)
                    split_p = re.split(r'[)\]]', proteins)
                    protein = split_p[0]
                    num = 1

                    print (split_p)
                    for int_ in split_p[1:]:
                        print ('4 letter Int Value is: ' + str(int_))
                        if int_.isdigit():
                            num *= int(int_)
                        else:
                            pass

                elif re.match( r'[A-Z][a-z][a-z]', proteins):
                    info_in_string_reduced.append(proteins)
                    split_p = re.split(r'[)\]]', proteins)
                    protein = split_p[0]
                    num = 1

                    print (split_p)
                    for int_ in split_p[1:]:
                        print ('3 letter Int Value is: ' + str(int_))
                        if int_.isdigit():
                            num *= int(int_)
                        else:
                            pass

                overall_subunits += num
                print(protein, num)

            print ('Overall Subumit is: ' + str(overall_subunits))
            unique_subunits = '1'

            fileWrite(subunitFile, gene, accession, frameid, synonyms, location, cofactors, info_in_string_reduced, overall_subunits, unique_subunits)

        else:
                component = comp.find('dd').get_text()

                print (component)

                info_in_string = re.split(r'[([]', component)
                info_in_string = list(filter(None, info_in_string))

                print (info_in_string)

                info_in_string_reduced = []

                overall_subunits = 0

                for proteins in info_in_string:
                    if re.match( r'[A-Z][a-z][a-z][A-Z]', proteins):
                        info_in_string_reduced.append(proteins)
                        split_p = re.split(r'[)\]]', proteins)
                        protein = split_p[0]

                        num = 1

                        print (split_p)
                        for int_ in split_p[1:]:
                            print ('4 letter Int Value is: ' + str(int_))
                            if int_.isdigit():
                                num *= int(int_)
                            else:
                                pass
                    elif re.match( r'[A-Z][a-z][a-z]', proteins):
                        info_in_string_reduced.append(proteins)
                        split_p = re.split(r'[)\]]', proteins)
                        protein = split_p[0]
                        num = 1

                        print (split_p)
                        for int_ in split_p[1:]:
                            print ('3 letter Int Value is: ' + str(int_))
                            if int_.isdigit():
                                num *= int(int_)
                            else:
                                pass
                    overall_subunits += num
                print(protein, num)
                print (info_in_string_reduced)


                if len(info_in_string_reduced) == 1:
                    cistronic = 'No'
                    overall_subunits = '1'
                    unique_subunits = '1'
                    fileWrite(subunitFile, gene, accession, frameid, synonyms, location, cofactors, "Monomer", overall_subunits, unique_subunits)

                else:
                    #This part of the code is used to identify of proteins are multi-cistronics by checking the first three letters of all items in a list and seeing if they are the same 
                    counter = 0
                    counter_max = len(info_in_string_reduced)
                    print (info_in_string)
                    print (info_in_string_reduced[0])

                    if len(info_in_string_reduced[0]) > 3:
                        first_entry = info_in_string_reduced[0][:3]

                    else:
                        first_entry = info_in_string_reduced[0]
                    print (first_entry)

                    for entry in info_in_string_reduced:
                        if entry[:3] == first_entry:
                            counter += 1
                        else:
                            pass

                    if counter == counter_max:
                        cistronic = 'Yes'

                    else:
                        cistronic = 'Double Check'

                    print ("Used Component to determine Subunit Info" + str(info_in_string_reduced))
                    print ('This protein is multi cistronic?' + cistronic)

                    #This part of the code uses info_in_string_reduced to detemrine the number of unique subunits

                    unique_subunits = str(len(info_in_string_reduced))
                    print ("Number of Unique Subunits is: " + unique_subunits)
                    print ("The overall subunit count is: " + str(overall_subunits))
                    fileWrite(subunitFile, gene, accession, frameid, synonyms, location, cofactors, info_in_string_reduced, overall_subunits, unique_subunits)

in_file.close()
subunitFile.close()
