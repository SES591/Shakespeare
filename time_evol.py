#!/usr/bin/python
#bioinfo.py

__author__ = '''Hyunju Kim'''

import os
import sys
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt

import input_net as inet
import updating_rule as ur


################# BEGIN: decimal_to_binary(nodes_list, decState, Nbr_States=2) ########################
def decimal_to_binary(nodes_list, decState, Nbr_States=2): # more left in the nodes list means higher order of 2 in binary
    biStates = {}
    x = len(nodes_list) -1
    for u in nodes_list:
        biStates[u] = decState / np.power(Nbr_States, x)
        decState = decState % np.power(Nbr_States, x)
        x = x - 1
    return biStates
################# END: decimal_to_binary(nodes_list, decState, Nbr_States=2) ########################


################# BEGIN: binary_to_decimal(nodes_list, biStates, Nbr_States=2) ########################
def binary_to_decimal(nodes_list, biStates, Nbr_States=2):  # more left in the nodes list means higher order of 2 in binary
    decState = 0
    x = len(nodes_list) -1
    for u in nodes_list:
        decState = decState + biStates[u]  * np.power(Nbr_States, x)
        x = x - 1
    return decState
################# END: binary_to_decimal(nodes_list, biStates, Nbr_States=2) ########################


################# BEGIN: ensemble_time_series(net, nodes_list, Nbr_States=2, MAX_TimeStep=20, Transition_Step=0) ########################
def ensemble_time_series(net, nodes_list, Nbr_States=2, MAX_TimeStep=20):
    
    '''
        Arguments:
        1. net
        2. Nbr_States
        3. MAX_TimeStep
        
        Return:
        1. timeSeriesData
        '''
    
    Nbr_Nodes = len(net.nodes())
    Nbr_All_Initial_States = np.power(Nbr_States, Nbr_Nodes)
    
    timeSeriesData = {}
    for n in net.nodes():
        timeSeriesData[n] = {}
        for initState in range(0, Nbr_All_Initial_States):
            timeSeriesData[n][initState] = []
    
    for initDecState in range(0, Nbr_All_Initial_States):
        currBiState = decimal_to_binary(nodes_list, initDecState, Nbr_States)
        for step in range(0, MAX_TimeStep):
            prevBiState = currBiState.copy()
            for n in nodes_list:
                timeSeriesData[n][initDecState].append(prevBiState[n])
            currBiState = ur.sigmoid_updating(net, prevBiState)

    return timeSeriesData
################# END: ensemble_time_series(net, nodes_list, Nbr_States=2, MAX_TimeStep=20) ########################


################# BEGIN: net_state_transition_map(net, nodes_list, Nbr_States=2) ########################
def net_state_transition(net, nodes_list, Nbr_States=2):

    '''
    Arguments:
               1. net
               2. Nbr_States
    Return:
               1. decStateTransMap
    '''
    
    Nbr_Nodes = len(net.nodes())
    Nbr_All_Initial_States = np.power(Nbr_States, Nbr_Nodes)
    
    decStateTransMap = nx.DiGraph()
    for prevDecState in range(0, Nbr_All_Initial_States):
        prevBiState = decimal_to_binary(nodes_list, prevDecState, Nbr_States)
        currBiState = ur.sigmoid_updating(net, prevBiState)
        currDecState = binary_to_decimal(nodes_list, currBiState, Nbr_States)
        decStateTransMap.add_edge(prevDecState, currDecState)
    return decStateTransMap
################# END: net_state_transition_map(net, nodes_list, Nbr_States=2) ########################


################# BEGIN: attractor_analysis(decStateTransMap) ########################
def find_attractor(decStateTransMap):
    
    '''
        Arguments:
        1. decStateTransMap
        Return:
        1. attractor
    '''
    attractor_list = nx.simple_cycles(decStateTransMap) #in case of deterministic system, any cycle without considering edge direction will be directed cycle.
    attractors = {}
    attractors['fixed'] = []
    attractors['cycle'] = []

    for u in attractor_list:
        if len(u) == 1:
            attractors['fixed'].append(u)
        else:
            attractors['cycle'].append(u)

    return attractors
################# END: attractor_analysis(decStateTransMap) ########################

def main():
    print "time_evol module is the main code."
    EDGE_FILE = '../data/example/example-net-edges.dat'
    NODE_FILE = '../data/example/example-net-nodes.dat'
    
    net = inet.read_network_from_file(EDGE_FILE, NODE_FILE)
    nodes_list = inet.build_nodes_list(NODE_FILE)

    timeSeriesData = ensemble_time_series(net, nodes_list, 2, 10)#, Nbr_States=2, MAX_TimeStep=20)
    

    initState = 1
    biStates = decimal_to_binary(nodes_list, initState)
    print 'initial state', biStates
    print 'a', timeSeriesData['a'][1]
    print 'b', timeSeriesData['b'][1]
    print 'c', timeSeriesData['c'][1]


    decStateTransMap = net_state_transition(net, nodes_list)
    nx.draw(decStateTransMap)
    plt.show()
    

    attractors = find_attractor(decStateTransMap)
    print attractors





if __name__=='__main__':
    main()
