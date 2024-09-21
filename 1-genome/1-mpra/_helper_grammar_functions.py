import pickle
import js
from itertools import product
import networkx as nx
import pandas as pd
import matplotlib.pyplot as plt

def can_be_converted_to_int(val):
    try: 
        int(val)
        return True
    except ValueError:
        return False

def get_phrase_string(phraseKey):
    return '_'.join([str(i) for i in phraseKey])

def get_phrase_key(phrase_string):
    phraseKey=[]
    for p in phrase_string.split('_'):
        try:
            p=int(p)
            phraseKey.append(p)
        except ValueError:
            phraseKey.append(p)
    phraseKey=tuple(phraseKey)
    return phraseKey

ori2flip={
    'g':'G',
    'G':'g',
    'e':'E',
    'E':'e',
    'S':'s',
    's':'S',
    'T':'t',
    't':'T',
}

def motif_revcomp(phrase,phraseType):

    if   type(phrase)==tuple: 
        pass
        return_tuple=True
    elif type(phrase)==str:   
        phrase=get_phrase_key(phrase)
        return_tuple=False
    else:                     raise ValueError("Input type not supported")
        
    # determine if you want to flip sites
    orientation='ori' in phraseType
    
    # reverse 
    phrase=phrase[::-1]
    
    # complement
    rcPhrase=[]
    for pi in phrase:
        
        # skip if spacing
        if can_be_converted_to_int(pi): 
            rcPhrase.append(pi)
            continue
        # flip if you are considering ori
        if orientation: 
            pi=ori2flip[pi[0]]+pi[1:]
            rcPhrase.append(pi)
            
        elif not orientation:
            rcPhrase.append(pi)
           
    if return_tuple:
        return tuple(rcPhrase)
    else:
        return '_'.join([str(i) for i in rcPhrase])        
    
def motif_revcomp_with_phrasetype(phraseWithPhrasetype):
       
    phraseType,phrase=phraseWithPhrasetype.split('=')

    if   type(phrase)==tuple:
        raise ValueError('Tuple not supported')

    elif type(phrase)==str:
        phrase=get_phrase_key(phrase)
        return_tuple=False

    else:                     raise ValueError("Input type not supported")
    
    # determine if you want to flip sites
    orientation='ori' in phraseType
    
    # reverse
    phrase=phrase[::-1]
    
    # complement
    rcPhrase=[] 
    for pi in phrase:
    
        # skip if spacing
        if can_be_converted_to_int(pi):
            rcPhrase.append(pi)
            continue
        # flip if you are considering ori
        if orientation:
            pi=ori2flip[pi[0]]+pi[1:]
            rcPhrase.append(pi)
            
        elif not orientation:
            rcPhrase.append(pi)
            
    else:
        return phraseType+'='+'_'.join([str(i) for i in rcPhrase])

def apply_permutations(motif,before2after,affinity):
    
    if   type(motif)==str:   pass
    elif type(motif)==tuple: motif=get_phrase_string(motif)
    elif type(motif)==list:  motif=get_phrase_string(motif)
    else:                    raise ValueError('Input type not supported')
    
    # Generate all possible elements at each position
    allPossibleElementsList=[]
    for mi in motif.split('_'):
        
        # If spacing, just add without any options
        if can_be_converted_to_int(mi):
            allPossibleElementsList.append([mi])
            
        else:
            allPossibleElementsList.append(before2after[mi])

    # Create all posibilities 
    allPermutations={}
    for i in range(len(allPossibleElementsList)):
        if i==0: 
            allPermutations[i]=allPossibleElementsList.pop(0)

        else:    
            allPermutations[i]=[]
            startList=allPermutations[i-1]
            nextPossibilityList=allPossibleElementsList.pop(0)
            for progress in startList:
                for possibility in nextPossibilityList:
                    allPermutations[i].append(progress+'_'+possibility)

    completeListIndex=max(allPermutations.keys())
    allPerms=allPermutations[completeListIndex]
    
    # If an affinity permutation, only allow for one of each site (ie, you cannot allow g1 and G1 to be in same seq)
    if affinity:
        selectPerms=[]
        for perm in allPerms:
            phrase2count={'G1':0,'G2':0,'G3':0,'E1':0,'E2':0}
            tfbsUnique=True
            perm=get_phrase_key(perm)
            for phrase in perm:
                phrase=phrase[0].upper()+phrase[1] # collapse on orientation
                phrase2count[phrase]+=1
                if phrase2count[phrase]>1: 
                    tfbsUnique=False # skip if this phrase has more than 2 occurences of a tfbs
                    break
            if tfbsUnique: selectPerms.append(get_phrase_string(perm))
                
        return selectPerms
    else:
        return allPerms
        
def scale_fc(fc):
    return 2*fc-1

def scale_fc(fc):
    return 2*fc-1

F    = 'limegreen'
NF   = 'red'
null = 'darkgrey'
edge2color= {
'Keep F'   : F,
'Keep NF'  : NF,
'Switch NF': NF,
'Switch F' : F,
'Gain NF'  : NF,
'Gain F'   : F,
'Loss'     : null,
'Loss'     : null,
'Null'     : null
}

node2color= {
    'F':F,
    'NF':NF,
    '':null
    
}

def print_graph_all_descendents(G,MotifGraph,node,descendType,size=(7,6),arrow_weight=20,edge_size=5):
    
    # Determine nodes to plot
    nodeList=[node]
    nodesToAdd=list(MotifGraph[node][descendType])
    for ni in nodesToAdd:
        niRC=motif_revcomp_with_phrasetype(ni)
        if niRC not in nodeList:
            nodeList.append(ni)
    
    s=G.subgraph(nodeList)
    nodelabels={node:node.split('=')[-1].replace("_",'\n') for node in s.nodes()}
    edgelabels={edge:s.edges[edge]['motifDelta'] for edge in s.edges()}
    edgecolors = [edge2color[motifDelta] for motifDelta in nx.get_edge_attributes(s,'motifDelta').values()]
    nodecolors = [node2color[motifType]  for motifType  in nx.get_node_attributes(s,'motifType').values()]
    
    nodesizes   = []
    nodesizeadj = 1500
    Node2FuncFc    = nx.get_node_attributes(s,'funcEnrich')
    Node2NonfuncFc = nx.get_node_attributes(s,'nfunEnrich')
    for ni,motifType in nx.get_node_attributes(s,'motifType').items():
        if   motifType=='NF': fc=scale_fc(Node2NonfuncFc[ni])
        elif motifType=='F' : fc=scale_fc(Node2FuncFc[ni])
        else:                 fc=1
        nodesizes.append(fc*nodesizeadj)
    
    # Determine node2pos
    node2pos={node:(0,0),}
    y_last=0
    for ni in nodeList[1:]:
        node2pos[ni]=(1,y_last)
        y_last+=1
        
    pos = nx.spring_layout(s, pos = node2pos,iterations=0)


    fig,ax=plt.subplots(1,figsize=size,dpi=150)
    nx.draw(s,pos,ax=ax,
            labels=nodelabels,
            edge_color=edgecolors,
            node_color=nodecolors,
            node_size=nodesizes,arrowsize=arrow_weight,width=edge_size)
    nx.draw_networkx_edge_labels(s,pos,edge_labels=edgelabels)
    plt.show()
