# -*- coding: utf-8 -*-

import os
import sqlite3

from functions import pull_sentences, pull_abstracts

def check_keyword(sentence, key_elements, mutation_names, gene_names, ptmod):

    Flag = "green"
    keFlag = False; ptFlag = False; gnFlag = False; mnFlag = False
    
    for x in key_elements:
        if x in sentence:
            keFlag = True
    
    for x in ptmod:
        if x in sentence:
            ptFlag = True
    
    for x in gene_names:
        if x.strip("\"") in sentence:
            gnFlag = True
    
    for x in mutation_names:
        if x in sentence:
            mnFlag = True

    if len([x for x in [keFlag, ptFlag, gnFlag, mnFlag] if x == True]) > 2:
        Flag = "Red"
    elif (keFlag == True and ptFlag == True) or (mnFlag == True and ptFlag == True):
        Flag = "Red"
    elif len([x for x in [keFlag, ptFlag, gnFlag, mnFlag] if x == True]) == 1:
        Flag = "Orange"

    return Flag



def parse_output():
    
    fdata = open("out.txt", "r")
    logProbab = dict()
    fetchFlag = True

    
    for lines in fdata:
        

    
        if "file" in lines[:4]:
            break
        
        if fetchFlag == True and lines.strip() != "":
            key = lines.strip()
            fetchFlag = False
    
        if fetchFlag == False and "logprob=" in lines:
            logProbab[key] = (float(lines.split(" ")[3]), float(lines.split(" ")[5]))
            fetchFlag = True

    return logProbab


def execute_probab(pmid, mutation, gene_name, key_elements, ptmod):

    sentences = pull_abstracts(str(pmid))
    
    
    # Its UV component is the major epidemiologic risk factor for squamous cell carcinoma of the skin
    # Funny output - We found remarkably few mutations.
    
    stopwords = ["analyzed", "tested", "performed", "examine", "determine", "encode", "was obtained", "previous report", "is the", "assess", "reviewed", "were selected", "used to", "studied", "should", "considered", "used for", "was used", "looked into", "were mutated", "comment", "using", "searched", "performed"]


    with open("sentence.txt", "w") as fp:
        for lines in sentences:
            Flag = False
            for words in stopwords:
                if words in lines:
                    Flag = True
                    break
        
            if Flag == True:
                continue
        
            fp.write("%s\n" %lines.encode('utf-8'))
    
    os.system("./srilm-1.7.1/bin/macosx/ngram -lm bigram.lm -ppl sentence.txt -debug 2 > out.txt")

    logprobability = parse_output()

    for key, val in logprobability.items():
        logprobability[key] = (val[0], check_keyword(key, key_elements, mutation, gene_name, ptmod), val[1])

    return logprobability



def parse_sentences(pmid, mutation, gene, gene_name_dict, gene_name_list, ptmod, key_elements):
    
    retSentence = list(); negativeTags = list()
    
    
    try:
        gene_name = gene_name_dict[gene]
    except KeyError:
        for val in gene_name_list:
            if gene in val:
                gene_name = val
                break

    crucial_elements = ["associated", "elevated", "reduced", "increase", "rise", "loss", "decrease", "gain", "change", "damage", "increased", "decreased", "lost", "exhibit", "expression", "revealed", "amplification"]


    logprobability = execute_probab(pmid, mutation, gene_name, key_elements, ptmod)
    for key, val in logprobability.items():
        check_flag = False
        for elements in crucial_elements:
            if (elements in key) or (";" not in key and "(" not in key and ")" not in key and ":" not in key): # second part for eliminating references (if found in the abstracts)
                retSentence.append(key)
                check_flag = True
                break

        if check_flag == True:
            continue

        if ((val[0] > -20 or val[1] == "Red") or (-20 > val[0] > -30 and val[1] == "Orange")) and val[2] > 100:
            retSentence.append(key)
        else:
            negativeTags.append(key)

    return retSentence, negativeTags



def pull_cosmic(geneNames):
    con = sqlite3.connect(r"../data/Cosmic.db")
    
    retDict = dict()
    #geneNames = ["TP53"]
    with con:
        cur = con.cursor()
        for gene in geneNames:
            retDict[gene] = dict()
            try:
                cur.execute("select Pubmed_PMID, Mutation_CDS, Mutation_AA from '%s';" %gene)
                retData = cur.fetchall()
            except sqlite3.OperationalError:
                continue
            for i in range(len(retData)):
                pmid, cds_mut, aa_mut = [str(x) for x in retData[i]]
                if "DataNotPresent" in pmid:
                    continue
                else:
                    mutation = [cds_mut, aa_mut, cds_mut.lstrip(".c"), aa_mut.lstrip(".p"), cds_mut.lstrip(".c")[-3:-2]+cds_mut.lstrip(".c")[:-3]+cds_mut.lstrip(".c")[-1]]
                retDict[gene][pmid] = (mutation)
    return retDict


def two_to_one(lObj):
    retObj = list()
    for val in lObj:
        if (type(val) is list) == True:
            for inval in val:
                retObj.append(inval)
        else:
            retObj.append(val)

    return retObj


def pull_geneNames():
    fdata = open("../data/aspm.txt", "r")
    retList = list()
    retDict = dict()
    for lines in fdata:
        if "Approved Symbol" in lines:
            continue
        
        lineObj = lines.split("\t")[:-3]
        for i, objects in enumerate(lineObj):
            if "," in objects:
                if '"' in objects:
                    lineObj[i] = [x.strip('"') for x in objects.split(", ")]
                else:
                    lineObj[i] = objects.split(", ")
    
        retDict[lineObj[0]] = two_to_one([x for x in lineObj if x != ""])
        retList.append(two_to_one([x for x in lineObj if x != ""]))

    return retDict, retList


key_elements = ["associated", "elevated", "reduced", "increase", "rise", "loss", "decrease", "gain", "change", "damage", "increased", "decreased", "lost", "exhibit", "cause", "expression", "high", "low", "indicate", "showed", "revealed", "resulted"]
ptmod = ["acetylation", "glycosylation", "amidation", "disulphide bond", "methylation", "phosphorylation", "sulfation", "kinase activity"]


gene_name_dict, gene_name_list = pull_geneNames()
print gene_name_dict, gene_name_list
cosmic_data = pull_cosmic()


writer = open("positive_tagged_sequences.txt", "w")
writer_neg = open("negative_tagged_sequences.txt", "w")
for gene, val in cosmic_data.items():
    for pmid, mutation in val.items():
        sentences, negativeTags = parse_sentences(pmid, mutation, gene, gene_name_dict, gene_name_list, ptmod, key_elements)
        print pmid, mutation, "\n\n"
        print sentences, "\n\n"
        if sentences != []:
            writer.write("\n\n%s\t%s\t%s\n\n" %(gene, pmid, mutation))
            for lines in sentences:
                writer.write("%s\n" %lines)

        if negativeTags != []:
            writer_neg.write("\n\n%s\t%s\t%s\n\n" %(gene, pmid, mutation))
            for lines in negativeTags:
                writer_neg.write("%s\n" %lines)

writer.close()












# Isues
#
# In this study, we describe the somatic landscape of pediatric Ewing sarcoma.
# We found remarkably few mutations.
# Pediatric Ewing sarcoma is characterized by the expression of chimeric fusions of EWS and ETS family transcription factors, representing a paradigm for studying cancers driven by transcription factor rearrangements.
# The aberrations were absent in well-differentiated carcinomas.
# word therapeutic
# We identified high frequencies of somatic mutations in CHD4 (17%), EP300 (8%), ARID1A (6%), TSPYL2 (6%), FBXW7 (29%), SPOP (8%), MAP3K4 (6%) and ABCC9 (6%).
# The aim of the present work was to study the prognostic significance of these DNA amplifications.
# It has been shown that it is possible to delineate subgroups of breast tumors according to specific sets of DNA amplifications.
# strand conformational polymorphism, using an automated sequencer.
# Interestingly, stratified analysis according to nodal status confirmed results obtained in the univariate tests: significance of MDM2 amplification and p53 mutations in node-negative and that of CCND1, EMS1, and FGFR1 in node-positive patients.
# Multivariate analysis on an unselected set of patients retained significance for the amplification of EMS1, FGFR1, and MDM2 with DFS, of CCND1 with OVS, and of RMC20C001 with both DFS and OVS.
# Three cases were positive for both HPV DNA and p53 antibody.
# To perform simultaneous analyses of oncogene or tumor suppressor gene mutations and related protein expression in single histologic sections, we have developed a novel method using an antigen-retrieval solution for a polymerase chain reaction template before immunohistochemical staining.
# DNA analysis of antigen-retrieval solutions was possible in all 20 cases and revealed completely consistent results (100%) with fresh cancer tissue and microdissected cancer tissue of paraffin-embedded histologic sections.














