import sqlite3


def fragmentation(string_obj): return [x.strip(".").strip(";") for x in string_obj[7:].split(" ") if x != "{ECO" and x != ['(in', 'dbSNP']]


def two_one(lObj):
    retList = list()
    for inL in lObj:
        for x in inL:
            retList.append(x)

    return retList


def test_sentence(sentence, omit_sentence):
    retFlag = True
    for object in omit_sentence:
        if object in sentence:
            retFlag = False
            return retFlag

    return retFlag


con = sqlite3.connect(r"../data/Uniprot_Human.db") # Locally available
cur = con.cursor()
fragments = list()


idData = open("../data/IDs.txt", "r")

idList = list()

for lines in idData:
    objects = [x for x in lines.strip().split(" ") if x != ""]
    idList.append(objects)


omit_sentence = ["allele DRB3", "sporadic cancers somatic mutation", "MMAM mu", "allele DQB1", "-EBS dbSNP", "allele DRB1", "sporadic cancer somatic mutation", "sample somatic mutation", "unknown pathological significance", "uncertain pathological significance", "->", "no significant activity", "a primary colorectal cancer", "confirmed at protein level", "a colorectal cancer patient"]


idList = two_one(idList)

for ids in idList:
    #print ids
    try:
        cur.execute("select Mutation from '%s';" %ids)
        score_collect = cur.fetchall()
        functional_strings = [str(x[0].strip()) for x in score_collect if str(x[0].strip()) != "NA"]
        #print functional_strings
        if functional_strings == []:
            continue
        for string_obj in functional_strings:
            
            fr_out = fragmentation(string_obj)
            if fr_out != ['(in', 'dbSNP']:
                sentence = " ".join(fr_out[1:]).rstrip(")").lstrip("(")
                if sentence.split(" ")[0].isupper() == True:
                    sentence = " ".join(sentence.split()[1:])
                
                if " " in sentence:
                    if test_sentence(sentence, omit_sentence) == True and sentence.count(" ") > 2 and sentence.split(" ")[0] not in ["allele", "in", "NSHA", "allotype"]:
                        if sentence.count(" ") < 5 and "cell line" in sentence:
                            continue
                        else:
                            fragments.append(sentence)
                        #print sentence
    except sqlite3.OperationalError:
        print "Panga", ids
        continue




fragments = set(fragments)

print len(fragments)

with open("sentences.txt", "w") as fp:
    for lines in fragments:
        fp.write("%s\n" %lines)










