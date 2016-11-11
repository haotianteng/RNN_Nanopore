#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 13 06:27:43 2016

@author: haotianteng
"""
import h5py as h5, numpy as np, os, sys, urllib, lxml.html

path = "http://data.genomicsresearch.org/internal/Nanopore/Ecoli/porecamp/"
connection = urllib.urlopen(path)
dom = lxml.html.fromstring(connection.read())
f5_files = []
for link in dom.xpath('//a/@href')[5:]: # first 5 links are not files on this page
    f5_files.append(link)
for filename in f5_files:
    print filename
f5_urls = [path + f for f in f5_files]
for i, f in enumerate(f5_urls):
    opener = urllib.URLopener()
    opener.retrieve(f, f5_files[i])
filename = f5_files[1]
# a list to sotre the names of the "groups" - similar to keys in a dictionary
lst_names = []
# open the fast5 file - closes automatically when loop is finished.
with h5.File(filename, 'r') as f:
    # gets the name of all of the groups and stores them in our list
    f.visit(lst_names.append)
    
    # get the names of the groups we are interested in
    groups = []
    for name in lst_names:
        if name.endswith("Events") or name.endswith("template/Fastq"):
            if not name.endswith("complement/Events") and "2D" not in name:
                groups.append(name)
    
    # catch if there are more than 3 elements in the groups that have been extracted
    if len(groups) != 3:
        raise IndexError("{0} groups detected. Need 3 groups. Speak to Michael if you see this.".format(len(groups)))

    
    # extract the information from the basecalled events
    events = f[groups[0]].value
    called_events = np.array([x for xs in events for x in xs]).reshape(len(events), 
                                                                       len(events[0]))
    
    # extract the FASTQ file
    fastq = f[groups[1]].value
    
    # extract the information from the detected events - before basecalling
    dataset = f[groups[2]].value
    event_detection = np.array([x for xs in dataset for x in xs]).reshape(len(dataset), 
                                                                          len(dataset[0]))
     