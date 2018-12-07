#Pseudo:
#Subset motifs
#write to file in meme format
#use tomtom compare two motifs
#if e-val is small enough return true
#delete tomtom/motif files -> file operations high overhead... ignore for now
#we can optimize this later
import os
from copy import deepcopy

#TODO: do this better
def merge(motifs):
    newMotif = deepcopy(motifs[0])
    for i,motif in enumerate(motifs):
        if i > 0:
            newMotif['Motif_Locations'].extend(motif['Motif_Locations'])
    return newMotif


def ParseTomTom(path):

    pval = 0.0
    with open(path + '/tomtom.tsv','r') as tomtomfile:
        header = tomtomfile.readline()
        results = tomtomfile.readline()
        elems = results.split('\t')
        pval = float(elems[3])

    return pval

#TODO: use subprocess and handle errors
def RunTomTom(path1,path2,outpath):
    tomtom_command ='/kb/deployment/bin/meme/bin/tomtom -oc ' + outpath + ' ' + path1 + ' ' + path2

    os.system(tomtom_command)
    return

#TODO: put this in motif utils
#TODO: use background frequency..

def WriteMotifAsMEME(motif,path):
    MEMEText = 'MEME version 4\n\n'
    sortedAlph = sorted(MSO['Alphabet'])
    alphStr = ''.join(sortedAlph)
    MEMEText += alphStr + '\n'
    MEMEText += '\n'
    MEMEText += 'strands: + -\n\n'
    MEMEText += 'Background letter frequencies\n'
    for letter in sortedAlph:
        MEMEText += letter + ' ' + str(.25) + ' '
    MEMEText += '\n\n'
    MEMEText += 'MOTIF ' + motif['Iupac_sequence'] + '\n'
    MEMEText += 'letter-probability matrix: alength= 4 w= '
    MEMEText += str(len(motif['Iupac_sequence'])) + ' nsites= '
    if motif['Motif_Locations'] is not None:
        MEMEText += str(len(motif['Motif_Locations']))
    else:
        MEMEText += '0'
    MEMEText += ' E= 0.0\n'
    #TODO: PVAL! -> ADD TO MSO -> PARSE IT HERE
    #MEMEText +=
    #print(len(motif['Iupac_sequence']))
    #print(len(motif['PWM']['T']))
    #exit()
    for i in range(0,len(motif['Iupac_sequence'])):
        for letter in sortedAlph:
            #print(i)
            try:
                MEMEText += str(motif['PWM'][letter][i]) + ' '
            except IndexError:
                print(len(motif['Iupac_sequence']))
                print(len(motif['PWM'][letter]))
                print(letter)
        MEMEText += '\n'
    MEMEText += '\n'
    with open(path,'w') as motifFile:
        motifFile.write(MEMEText)

    return

def CompareMotifs(motif1,motif2):

    pvalThresh = 1.0e-10
    path1 = '/kb/module/tmp/motif1.txt'
    path2 = '/kb/module/tmp/motif2.txt'
    WriteMotifAsMEME(motif1,path1)
    WriteMotifAsMEME(motif2,path2)

    tomtompath = '/kb/module/tmp/tomtom_out'
    RunTomTom(path1,path2,tomtompath)

    pval = ParseTomTom(tomtompath)

    os.remove(path1)
    os.remove(path2)
    if pval <= pvalThresh:
        return True

    return False