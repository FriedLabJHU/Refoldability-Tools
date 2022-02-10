import re
import sys

in_fn1 = sys.argv[1]	#provide a filename
in_fn2 = sys.argv[2]	#provide a filename
out_fnI = sys.argv[3]	#provide an output filename for the intersection
#out_fnU = sys.argv[4]	#provide an output filename for the union

inputFile1 = open(in_fn1,'r')
inputFile2 = open(in_fn2,'r')
outputFileI = open(out_fnI,'w')
#outputFileU = open(out_fnU,'w')

data1 = []
data2 = []
listOfPeptides1 = []
listOfPeptides2 = []
dataOutI = []
#dataOutU = []

for line in inputFile1:
    infoInLine = line.split('\t')
    data1.append( infoInLine[:-1] )
    data1[-1].append( infoInLine[-1][:-1])
    listOfPeptides1.append( infoInLine[1] )
    
for line in inputFile2:
    infoInLine = line.split('\t')
    data2.append( infoInLine[:-1] )
    data2[-1].append( infoInLine[-1][:-1] )
    listOfPeptides2.append( infoInLine[1] )

data1rows = len( data1 )
data2rows = len( data2 )

inputFile1.close()
inputFile2.close()

#start by listing the peptides in common    
for indexin1 in range( data1rows ):
    peptide1 = listOfPeptides1[indexin1]
    if peptide1 in listOfPeptides2:
        indexin2 = listOfPeptides2.index( peptide1 )
        dataOutI.append( data1[indexin1] + data2[indexin2][15:19] )
        #dataOutU.append( data1[indexin1] )

#peptides only in set1
#for indexin1 in range( data1rows ):
#    peptide1 = listOfPeptides1[indexin1]
#    if peptide1 not in listOfPeptides2:
#        dataOutU.append( data1[indexin1] )

#peptides only in set2        
#for indexin2 in range( data2rows ):
#    peptide2 = listOfPeptides2[indexin2]
#    if peptide2 not in listOfPeptides1:
#        dataOutU.append( data2[indexin2] )
        
#write output file for the union
#outputLine = 'peptide\tgene name\tposition\tfrom start\tfrom end\tabundance\tback-context\tcodon\tfwd-context\tredundancy\tanticodon1-1\tanticodon1-2\tanticodon1-3\tanticodon2-1\tanticodon2-2\tanticodon2-3\tanticodon3-1\tanticodon3-2\tanticodon3-3\tanticodon1-1N\tanticodon1-2N\tanticodon1-3N\tanticodon2-1N\tanticodon2-2N\tanticodon2-3N\tanticodon3-1N\tanticodon3-2N\tanticodon3-3N \tanticodon1-1NL\tanticodon1-2NL\tanticodon1-3NL\tanticodon2-1NL\tanticodon2-2NL\tanticodon2-3NL\tanticodon3-1NL\tanticodon3-2NL\tanticodon3-3NL\tanticodon1-ave\tanticodon2-ave\tanticodon3-ave\tanticodon1-std\tanticodon2-std\tanticodon3-std\tanticodon1-CV\tanticodon2-CV\tanticodon3-CV\tanticodon2-anticodon1\tanticodon3-anticodon1\tlogpval12\tlogpval13\n' 
#outputFileU.write( outputLine )
#for i in range(len(dataOutU)):
#    outputLine = '\t'.join(str(x) for x in dataOutU[i]) + '\n'
#    outputFileU.write( outputLine )    
#outputFileU.close()         
        
#write output file for the intersection
#outputLine = 'peptide\tgene name\tposition\tfrom start\tfrom end\tabundance\tback-context\tcodon\tfwd-context\tredundancy\tanticodon1-1\tanticodon1-2\tanticodon1-3\tanticodon2-1\tanticodon2-2\tanticodon2-3\tanticodon3-1\tanticodon3-2\tanticodon3-3\tanticodon1-1N\tanticodon1-2N\tanticodon1-3N\tanticodon2-1N\tanticodon2-2N\tanticodon2-3N\tanticodon3-1N\tanticodon3-2N\tanticodon3-3N \tanticodon1-1NL\tanticodon1-2NL\tanticodon1-3NL\tanticodon2-1NL\tanticodon2-2NL\tanticodon2-3NL\tanticodon3-1NL\tanticodon3-2NL\tanticodon3-3NL\tanticodon1-ave\tanticodon2-ave\tanticodon3-ave\tanticodon1-std\tanticodon2-std\tanticodon3-std\tanticodon1-CV\tanticodon2-CV\tanticodon3-CV\tanticodon2-anticodon1\tanticodon3-anticodon1\tlogpval12\tlogpval13\n' 
#outputFileI.write( outputLine )
for i in range(len(dataOutI)):
    outputLine = '\t'.join(str(x) for x in dataOutI[i]) + '\n'
    outputFileI.write( outputLine )    
outputFileI.close()         
