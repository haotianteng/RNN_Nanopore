# -*- coding: utf-8 -*-
"""
Spyder Editor

Downloading fast5 file
"""
import urllib,urllib2,re,os,sys

DataURL = "http://data.genomicsresearch.org/internal/Nanopore/GN_003_R9_2D_20160928/reads/downloads/pass/"
FileSaveLoc = "/home/haotianteng/UQ/BINF7000/Nanopore/GN_003_R9_2D_20160928/fast5/"

####################################For tidy output
def convert_bytes(num):
    """
    this function will convert bytes to MB.... GB... etc
    """
    for x in ['bytes', 'KB', 'MB', 'GB', 'TB']:
        if num < 1024.0:
            return "%3.1f %s" % (num, x)
        num /= 1024.0



####################################Extract fast5 file list
def FileExtractor()
if not os.path.exists(FileSaveLoc):
    os.makedirs(FileSaveLoc)
print "Connect to server..."
Page = urllib2.urlopen(DataURL).read()
Links = re.findall('<a href=(.*?)>.*?</a>',Page)
fast5Links = []
for Link in Links:
    if "fast5" in Link:
        fast5Links.append(Link)
print "Totally %d fast5 files" %len(fast5Links)


####################################Download the files
fileNo = 100 #len(fast5Links) #Max number of files need to be download
TryMax = 5#Max download attemps
Try = 0#Current attempt times
fileCount = 1
for fileName in fast5Links[0:fileNo]:
    fileName = fileName[1:-1]#Get rid of the double quotation
    Try = 0
    while Try<TryMax :
        try:
            urllib.urlretrieve(DataURL+fileName,FileSaveLoc+fileName)
            fileSize = os.path.getsize(FileSaveLoc+fileName)
            print "Download %s (%d/%d) %s" %(fileName,fileCount,fileNo,convert_bytes(fileSize))
            break
        except IOError:
            print "Download data files,try again (%d/%d)" %(Try,TryMax)
            Try+=1
        except:
            print "Unexpected Error when downloading(%d/%d):"%(fileCount,fileNo),sys.exc_info()[0]
            Try+=1
    fileCount+=1