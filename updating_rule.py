#!/usr/bin/python
#attractor-analysis.py
#last update: 17 DEC 2015

__author__ = '''Hyunju Kim'''


import os
import sys
import numpy as np
import networkx as nx
from collections import OrderedDict

import input_net as inet



################# begin: sigmoid_updating ######################
def sigmoid_updating(net, prevState):

    '''
        Arguments:
                1. net
                2. prevState
        Return:
               1. currState
    '''
    
    currState = {}
    #### compute the current states of nodes in the net ####
    for v in net.nodes():
        #### compute weighted sum for node v over its neighbors u ####
        eSum = 0
        for u in net.predecessors_iter(v):
            w_uv = 1.0*net[u][v]['weight']
            eSum += w_uv * prevState[u]
        #### determine the current state for v as a function of eSum and threshold of v ####
        if eSum < net.node[v]['threshold']:
            currState[v] = 0
        if eSum > net.node[v]['threshold']:
            currState[v] = 1
        if eSum == net.node[v]['threshold']:
            currState[v] = prevState[v]

    return currState
    
################# end: sigmoid_updating ########################

def main():
    print "updating_rule module is the main code."
    EDGE_FILE = '../data/example/example-net-edges.dat'
    NODE_FILE = '../data/example/example-net-nodes.dat'
    
    net = inet.read_network_from_file(EDGE_FILE, NODE_FILE)

    #prevState = {'a':0.0, 'b':0.0, 'c':1.0}
    prevState = {}
    prevState['a'] = float(sys.argv[1])
    prevState['b'] = float(sys.argv[2])
    prevState['c'] = float(sys.argv[3])
    print "network state @ previous step", OrderedDict(sorted(prevState.items(), key=lambda t: t[0]))
    
    currState = sigmoid_updating(net, prevState)
    print "network state @ current step", OrderedDict(sorted(currState.items(), key=lambda t: t[0]))

if __name__=='__main__':
    main()
