#!/usr/bin/python
#bionetworks.py
#last update : 14 Aug 2014

__author__ = '''Hyunju Kim'''

import os
import sys
import numpy as np
import networkx as nx
from math import log
from collections import defaultdict
#import random as ran
#from math import log
#import itertools

import input_net as inet
import updating_rule as ur
import time_evol as ev



################# BEGIN: build_historyList (aList, historyLength) ############################
def build_historyList(aList, historyLength, Nbr_States=2):
    '''
    Description: 
                generate decimal state's index for k-previous states with a given list, aList
    Arguments:
                1. list of time series data
                2. history length
    Return:
                1. list of history
    '''
    historyList = []
    historyUnit = aList[:historyLength]
    aList[:historyLength] = []
    
    historyState = 0
    for s in range(historyLength):
        historyState += historyUnit[s] * power(Nbr_States, s)
    historyList.append(historyState)

    for x in aList:
        historyState = historyState / int(Nbr_States) + x * power(Nbr_States, historyLength - 1)
        historyList.append(historyState)

    return historyList
################# END: build_historyList (aList, historyLength) ########################


################## BEGIN: com_AI_over_ensemble(timeSeriesNode, historyLength, nodes_list) ########################
def com_AI_over_ensemble(timeSeriesNode, historyLength, Nbr_Nodes, Nbr_States=2):

    count_currState_hiState = defaultdict(int)
    count_hiState = defaultdict(int)
    count_currState = defaultdict(int)
    Nbr_All_Initial_States = np.power(Nbr_States, Nbr_Nodes)
    
    for initState in range(0, Nbr_All_Initial_States):
        aList = list(timeSeriesNode[initState])
        historyList = build_historyList(aList, historyLength) # aList becomes aList[historyLength:] after historyList function

        for s in range(len(aList)):  # obtain the distribution for each pattern to compute Active Information
            #print s
            count_currState_hiState[(aList[s], historyList[s])] += 1
            count_hiState[historyList[s]] += 1
            count_currState[aList[s]] += 1

    AI = 0
    for initState in range(0, Nbr_All_Initial_States):
        aList = list(timeSeriesNode[initState])
        historyList = build_historyList(aList, historyLength) # aList becomes aList[historyLength:] after historyList function
        sampleLength = len(aList) * pow(2,len(nodes_list))
        for s in range(len(aList)):
            prob_currState_hiState = float(count_currState_hiState[(aList[s], historyList[s])]) / float(sampleLength)
            prob_hiState = float(count_hiState[historyList[s]]) / float(sampleLength)
            prob_currState = float(count_currState[aList[s]]) / float(sampleLength)
            AI = AI + log( prob_currState_hiState / ( prob_currState * prob_hiState)) / log(2.0) # since the summation is over not all possible pattern of currState_hiState
    AI = AI / float(sampleLength)
    return AI
################## END: com_AI_over_ensemble(timeSeriesNode, historyLength, nodes_list) ########################


################## begin : comTE ########################
def com_TE_over_ensemble(timeSeriesNodeA, timeSeriesNodeB, historyLength, nodes_list):
    
    # declare dic for distribution to compute Transfer Entropy
    count_tarCurrState_tarHiState_sourPrevState = defaultdict(int)
    count_tarCurrState_tarHiState = defaultdict(int)
    count_tarHiState_sourPrevState = defaultdict(int)
    count_tarHiState = defaultdict(int)
    
    for si in range(0, pow(2,len(nodes_list))):
        sourList = list(timeSeriesNodeA[si])
        tarList = list(timeSeriesNodeB[si])
        historyList = build_historyList(tarList, historyLength) # tarList becomes tarList[historyLength:] after historyList function
        sourList[:historyLength - 1] = [] # sourList becomes sourList[historyLength-1:]
        #*** obtain the distribution for each pattern to compute Transfer Entropy ***#
        for s in range(len(tarList)):
            count_tarCurrState_tarHiState_sourPrevState[(tarList[s], historyList[s], sourList[s])] += 1
            count_tarCurrState_tarHiState[(tarList[s], historyList[s])] += 1
            count_tarHiState_sourPrevState[(historyList[s], sourList[s])] += 1
            count_tarHiState[historyList[s]] += 1

#    print count_tarCurrState_tarHiState_sourPrevState
#    print count_tarCurrState_tarHiState
#    print count_tarHiState_sourPrevState
#    print count_tarHiState
    #*** obtain the distribution for each pattern to compute Active Information ***#
    TE = 0
    for si in range(0, pow(2,len(nodes_list))):
        sourList = list(timeSeriesNodeA[si])
        tarList = list(timeSeriesNodeB[si])
        historyList = build_historyList(tarList, historyLength) # tarList becomes tarList[historyLength:] after historyList function
        sourList[:historyLength - 1] = [] # sourList becomes sourList[historyLength-1:]
        sampleLength = len(tarList) * pow(2,len(nodes_list))
        for s in range(len(tarList)):
            prob_tarCurrState_tarHiState_sourPrevState = float(count_tarCurrState_tarHiState_sourPrevState[(tarList[s], historyList[s], sourList[s])]) / float(sampleLength)
            prob_tarCurrState_tarHiState = float(count_tarCurrState_tarHiState[(tarList[s], historyList[s])]) / float(sampleLength)
            prob_tarHiState_sourPrevState = float(count_tarHiState_sourPrevState[(historyList[s], sourList[s])]) / float(sampleLength)
            prob_tarHiState = float(count_tarHiState[historyList[s]]) / float(sampleLength)
            TE = TE + log( (prob_tarCurrState_tarHiState_sourPrevState * prob_tarHiState) / ( prob_tarHiState_sourPrevState * prob_tarCurrState_tarHiState)) / log(2.0) # since the summation is over not all possible pattern of tarCurrState_tarHiState_sourPrevState but tarList, there is no prob_tarCurrState_tarHiState_sourPrevState multiplied by the log term.
    TE = TE / float(sampleLength)
    return TE
################## end : comAI ########################


################## begin : comAI ########################
def com_AI_over_trajectory(timeSeries, historyLength):
    
    aList = list(timeSeries)
    historyList = build_historyList(aList, historyLength) # aList becomes aList[historyLength:] after historyList function
    
    #*** declare dic for distribution to compute Active Information ***#
    count_currState_hiState = defaultdict(int)
    count_hiState = defaultdict(int)
    count_currState = defaultdict(int)
    #print "test", len(aList)
    #*** obtain the distribution for each pattern to compute Active Information ***#
    for s in range(len(aList)):
        #print s
        count_currState_hiState[(aList[s], historyList[s])] += 1
        count_hiState[historyList[s]] += 1
        count_currState[aList[s]] += 1
    
    #*** obtain the distribution for each pattern to compute Active Information ***#
    AI = 0
    for s in range(len(aList)):
        prob_currState_hiState = float(count_currState_hiState[(aList[s], historyList[s])]) / float(len(aList))
        #print "joint", prob_currState_hiState
        prob_hiState = float(count_hiState[historyList[s]]) / float(len(aList))
        #print "his", prob_hiState, historyList[s]
        prob_currState = float(count_currState[aList[s]]) / float(len(aList))
        #print "curr", prob_currState
        AI = AI + log( prob_currState_hiState / ( prob_currState * prob_hiState)) / log(2.0) # since the summation is over not all possible pattern of currState_hiState but aList, there is no prob_currState_hiState multiplied by the log term.
    
    AI = AI / float(len(aList))
    return AI
################## end : comAI ########################


################## begin : comTE ########################
def com_TE_over_trajectory(timeSeriesA, timeSeriesB, historyLength):
    
    sourList = list(timeSeriesA)
    tarList = list(timeSeriesB)
    historyList = build_historyList(tarList, historyLength) # tarList becomes tarList[historyLength:] after historyList function
    sourList[:historyLength - 1] = [] # sourList becomes sourList[historyLength-1:]
    #print timeSeriesA
    #print timeSeriesB
    #print tarList
    #print sourList
    #*** declare dic for distribution to compute Transfer Entropy ***#
    count_tarCurrState_tarHiState_sourPrevState = defaultdict(int)
    count_tarCurrState_tarHiState = defaultdict(int)
    count_tarHiState_sourPrevState = defaultdict(int)
    count_tarHiState = defaultdict(int)
    
    #*** obtain the distribution for each pattern to compute Transfer Entropy ***#
    for s in range(len(tarList)):
        #        print "len s", len(sourList)
        #        print "len t", len(tarList)
        #        print s
        #        print "s", sourList[s]
        #        print "h", historyList[s]
        #        print "t", tarList[s]
        count_tarCurrState_tarHiState_sourPrevState[(tarList[s], historyList[s], sourList[s])] += 1
        count_tarCurrState_tarHiState[(tarList[s], historyList[s])] += 1
        count_tarHiState_sourPrevState[(historyList[s], sourList[s])] += 1
        count_tarHiState[historyList[s]] += 1
    
    #*** obtain the distribution for each pattern to compute Active Information ***#
    TE = 0
    for s in range(len(tarList)):
        prob_tarCurrState_tarHiState_sourPrevState = float(count_tarCurrState_tarHiState_sourPrevState[(tarList[s], historyList[s], sourList[s])]) / float(len(tarList))
        prob_tarCurrState_tarHiState = float(count_tarCurrState_tarHiState[(tarList[s], historyList[s])]) / float(len(tarList))
        prob_tarHiState_sourPrevState = float(count_tarHiState_sourPrevState[(historyList[s], sourList[s])]) / float(len(tarList))
        prob_tarHiState = float(count_tarHiState[historyList[s]]) / float(len(tarList))
        TE = TE + log( (prob_tarCurrState_tarHiState_sourPrevState * prob_tarHiState) / ( prob_tarHiState_sourPrevState * prob_tarCurrState_tarHiState)) / log(2.0) # since the summation is over not all possible pattern of tarCurrState_tarHiState_sourPrevState but tarList, there is no prob_tarCurrState_tarHiState_sourPrevState multiplied by the log term.
    TE = TE / float(len(tarList))


return TE
################## end : comAI ########################

def main():
    print "info_dyn module is the main code."
    EDGE_FILE = '../data/example/example-net-edges.dat'
    NODE_FILE = '../data/example/example-net-nodes.dat'

    net = read_network_from_file(EDGE_FILE, NODE_FILE)
    nodes_list = build_nodes_list(NODE_FILE)



if __name__=='__main__':
    main()
