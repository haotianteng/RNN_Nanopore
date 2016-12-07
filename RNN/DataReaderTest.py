#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 25 02:24:03 2016

@author: haotianteng
"""

from DataReader import DataReader

#############################super parameter
TRANNING_READS = 30000  #Total reads used
TRAIN_WEIGHT = 0.4  #Proportion of reads used to train
TEST_WEIGHT = 0.4   #Proportion of reads used to test
VALID_WEIGHT = 0.2  #Proportion of reads used to validate

#Structure
EVENT_LENGTH = 20 #Length of each sentence
HIDDEN_UNIT_NUM = 24 #Length of the hidden state of each hidden layer
CellLayer = 3 #Number of the hidden layers

#Training
STEP_RATE = 0.5
BATCH_SIZE = 20
EPOCH = 5000
#############################

###############################Read the data
data = DataReader(TRAIN_WEIGHT,TEST_WEIGHT,VALID_WEIGHT,EVENT_LENGTH,TRANNING_READS,event_total = 10000,file_list = '/home/haotianteng/UQ/BINF7000/Nanopore/GN_003/event_pass.dat')