from mechanize import Browser
from BeautifulSoup import BeautifulSoup
import urllib2, urllib, os
import requests



def _extract(soup, year):
    table = soup.find("table", border=1)
    for row in table.findAll('tr')[1:]:
        col = row.findAll('td')
        rank = col[0].string
        artist = col[1].string
        album = col[2].string
        cover_link = col[3].img['src']
        record = (str(year), rank, artist, album, cover_link)
        print >> outfile, "|".join(record)
        save_as = os.path.join("./", album + ".jpg")
        urllib.urlretrieve("http://www.palewire.com" + cover_link, save_as)
        print "Downloaded %s album cover" % album


def _pull_from_web:
    outfile = open("albums.txt", "w")
    mech = Browser()
    url = "http://www.j-biomed-inform.com/article/S1532-0464(07)00145-1/fulltext"
    page1 = mech.open(url)
    html1 = page1.read()
    soup1 = BeautifulSoup(html1)
    extract(soup1, 2007)
    page2 = mech.follow_link(text_regex="Next")
    html2 = page2.read()
    soup2 = BeautifulSoup(html2)
    extract(soup2, 2006)
    outfile.close()


def _pull_from_pdf(fileName):
    files = {'f': (fileName, open(fileName, 'rb'))}
    response = requests.post("https://pdftables.com/api?key=mzx18gwkft5r&format=xml", files=files)
    response.raise_for_status() # ensure we notice bad responses
    print(response.text)


def _pull_xml(url):
    xmlData = urllib2.urlopen(url)
    soup = BeautifulSoup(xmlData)
    return soup


def _pull_pmc_data(ID):
    url = "http://www.ncbi.nlm.nih.gov/pmc/oai/oai.cgi?verb=GetRecord&identifier=oai:pubmedcentral.nih.gov:" + str(ID) + "&metadataPrefix=pmc"
    refined_xml = _pull_xml(url)
    return refined_xml








