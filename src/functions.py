# -*- coding: utf-8 -*-

from __future__ import print_function

import os
import re
import string
import nltk
import nltk.data
import sqlite3

from Bio import Entrez
from Bio.Entrez import efetch, read

from nltk.collocations import *
from nltk.model import NgramModel
from nltk.probability import LidstoneProbDist
from xgoogle.search import GoogleSearch, SearchError
import xml.etree.ElementTree as ET
from bs4 import BeautifulSoup


################################################################################################

class google_search:
    
    """
        google search class for pulling gene associated sentences
        @search_object - the search keyword
        
        """
    
    def __init__(self):
        self.object = search_object
    
    def perform_search(self):
        url_list = list()
        
        try:
            gs = GoogleSearch(self.object)
            gs.results_per_page = 50
            results = gs.get_results()
            
            for res in results:
                url_list.append(res.url.encode("utf8"))
        
            return url_list
        
        except SearchError, e:
            print("Search failed: %s" %e)

    def pull_xml(xmlFile):
        tree = ET.parse(xmlFile)
        root = tree.getroot()
        for child in root:
            print(child.tag, child.attrib)

    def pull_content(url):
        xmlData = urllib2.urlopen(url)
        return xmlData

    def clean_xml(xmlFile, keyObject):
        soup = BeautifulSoup(open(xmlFile, "r"), 'html.parser')
        pull_data = soup.get_text()
        filtered_object = list()
        for lines in pull_data.split("\n"):
            try:
                if lines !=  u'' and keyObject in lines and "{" not in lines and "http" not in lines:
                    filtered_object.append(str(lines))
            except UnicodeEncodeError:
                continue
        return filtered_object
        

################################################################################################

"""
        
def clean_xml(xmlFile, keyObject):
    soup = BeautifulSoup(open(xmlFile, "r"), 'html.parser')
    pull_data = soup.get_text()
    filtered_object = list()
    for lines in pull_data.split("\n"):
        try:
            if lines !=  u'' and keyObject in lines and "{" not in lines and "http" not in lines:
                filtered_object.append(str(lines))
        except UnicodeEncodeError:
            continue
    return filtered_object
    
    
"""

################################################################################################


def find_mutation_line(line, mutation_text): return True if mutation_text in line else False

def lexical_diversity(text): return float(len(set(text))) / len(text)

def percentage(count, total): return 100 * (float(count) / total)

def rem_punctuation(content): return content.translate(string.maketrans(string.punctuation, ' '*len(string.punctuation)))

def frequency_cal(listObj): return nltk.FreqDist(w for w in listObj)

def perplexity_cal(lm, test_data): return lm.perplexity(test_data)

def hasNumbers(inputString): return any(char.isdigit() for char in inputString)




################################################################################################



def clean_abstract(abstract_data):
    
    abstract_data = abstract_data.split("\n")
    Flag = False
    parserFlag = False
    abstract_lines = list()
    
    for lines in abstract_data:
        if "Author information:" in lines:
            Flag = True
        
        if Flag == True and lines == "":
            parserFlag = True
            Flag = False
            continue
        
        if Flag == False and lines == "":
            parserFlag = False
        
        if parserFlag == True:
            abstract_lines.append(lines.strip().decode('utf-8'))
    
    return " ".join(abstract_lines)



def fetch_abstract1(pmid):
    handle = efetch(db='pubmed', id=pmid, retmode='text', rettype='abstract')
    data = handle.read()
    return data


def pull_sentences(filename):
    """
        Breaks abstract into sentences
        """
    tokenizer = nltk.data.load('tokenizers/punkt/english.pickle')
    fp = open(filename)
    data = fp.read()
    return tokenizer.tokenize(data.decode('utf-8'))


def pull_abstracts(pmid):

    fq = open("abstracts.txt", "w")
    abstract_para = clean_abstract(fetch_abstract1(pmid))

    fq.write(abstract_para.encode('utf-8'))
    fq.close()
    retSentence = pull_sentences("abstracts.txt")
    return retSentence


####################################################################################################




def two_one(lObj):
    retList = list()
    for inL in lObj:
        for x in inL:
            retList.append(x)
    
    return retList



def train_model(fdist, listObj, n):
    
    """ 
        @n - size of ngram
        @fdist - frequency distribution
        @listObj - ngram data list
        
        """
    
    
    estimator = lambda fdist, bins: LidstoneProbDist(fdist, 0.2)
    lm = NgramModel(n, fdist, estimator=estimator)
    return lm


def remove_rare(fdist, listObj):
    vocabulary = set(map(lambda x: x[0], filter(lambda x: x[1] >= 5, fdist.iteritems())))
    return map(lambda x: x if x in vocabulary else "*unknown*", listObj)


def fetch_abstract(pmid):
    
    Entrez.email = "sendambuj@gmail.com"
    
    handle = efetch(db='pubmed', id=pmid, retmode='xml')
    xml_data = read(handle)[0]
    
    try:
        article = xml_data['MedlineCitation']['Article']
        abstract = article['Abstract']['AbstractText'][0]
        return abstract
    except IndexError:
        return None



def getMultigrams(sentence, multigramObj):
    
    sentence = rem_punctuation(sentence) # Remove punctuations
    wordList = nltk.wordpunct_tokenize(sentence) #Tokenize
    
    if multigramObj == "bigram":
        bigram_measures = nltk.collocations.BigramAssocMeasures()
        finder = BigramCollocationFinder.from_words(wordList)
        finder.apply_word_filter(lambda w: w in ('I', 'me', 'is', 'a', 'the'))
        scored = finder.score_ngrams(bigram_measures.raw_freq)
        return sorted(bigram for bigram, score in scored if hasNumbers(bigram[0]) != True and hasNumbers(bigram[1]) != True)
    elif multigramObj == "trigram":
        trigram_measures = nltk.collocations.TrigramAssocMeasures()
        finder = TrigramCollocationFinder.from_words(wordList)
        finder.apply_word_filter(lambda w: w in ('I', 'me', 'is', 'a', 'the'))
        scored = finder.score_ngrams(trigram_measures.raw_freq)
        return sorted(trigram for trigram, score in scored if hasNumbers(trigram[0]) != True and hasNumbers(trigram[1]) != True)



def classify(model, input_file):
    if model == "bigram":
        os.system("./srilm-1.7.1/bin/macosx/ngram -order 2 -lm bigram.lm -ppl %s -debug 2 > output_bi.txt" %input_file)
    elif model == "trigram":
        os.system("./srilm-1.7.1/bin/macosx/ngram -order 3 -lm trigram.lm -ppl %s -debug 2 > output_tri.txt" %input_file)




def create_model(senetence_file):
    sentence_obj = open(senetence_file, "rU")
    store_sent = list()
    sentList = [x.strip().rstrip(".").lower() for x in sentence_obj]
    
    bigrams = list()
    trigrams = list()

    for lines in sentList:
        bigrams.append(getMultigrams(lines, "bigram"))
        trigrams.append(getMultigrams(lines, "trigram"))

    bigrams = two_one(bigrams)
    trigrams = two_one(trigrams)

    bigram_fdist = frequency_cal(bigrams)
    #print(bigram_fdist.most_common())
    trigram_fdist = frequency_cal(trigrams)

    #bigram_clean = remove_rare(bigram_fdist, bigrams)
    #trigram_clean = remove_rare(trigram_fdist, trigrams)

    #bigram_model = train_model(bigram_fdist, bigrams, 2)
    #trigram_model = train_model(trigram_fdist, trigrams, 3)

    #print(bigram_model, trigram_model)
    with open("bigram_freq.txt", "w") as fp:
        for key, val in bigram_fdist.items():
            fp.write("%s\t%s\t%s\n" %(key[0], key[1], val))

    os.system("./srilm-1.7.1/bin/macosx/ngram-count -read bigram_freq.txt -order 2 -lm bigram.lm")

    with open("trigram_freq.txt", "w") as fp:
        for key, val in trigram_fdist.items():
            fp.write("%s\t%s\t%s\n" %(key[0], key[1], val))

    os.system("./srilm-1.7.1/bin/macosx/ngram-count -read trigram_freq.txt -order 3 -lm trigram.lm")




####################################################################################################



create_model("../data/sentences.txt")




















