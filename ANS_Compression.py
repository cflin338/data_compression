# -*- coding: utf-8 -*-
"""
Created on Fri Feb  7 15:06:24 2020

@author: Chris Lin

Primary functions to call are:
    stream_ANS_encoding(to_encode): requires only list of data to be compressed
        outputs objects: value, additional outputs, decoding table, all needed for decoding
    stream_ANS_decoding(val,  decoding table, additional outputs)
        outputs decoded sequence
"""
import numpy as np

def get_min_dict(d):
    minVal = d[-1]
    minVal_Ind = len(d)-1
    
    for i in range(len(d)-1,-1,-1):
        if minVal>d[i]:
            minVal = d[i]
            minVal_Ind = i
    return minVal_Ind


def ANS_table(weights,b,mode,states):
    #n = number of symbols
    n = len(weights)
    l = sum(weights)
    p = [i/l for i in weights]
    symb_counter=np.zeros(n)
    symb_prob_tracker = np.zeros(n)
    for i in range(0,n):        #symbols are integers 0 -> n-1
        symb_prob_tracker[i] = .5/p[i]
        symb_counter[i] = weights[i]
    if mode=='C':
        C = {}
    elif mode =='D':
        D={}
    else:
        raise("Invalid mode")
    l = states
    for x in range(l,b*l):
        s = get_min_dict(symb_prob_tracker)
        symb_prob_tracker[s]+=1/p[s]
        if mode=='D':
            D[x] = (s,symb_counter[s])
        else:
            C[(s,symb_counter[s])] = x
        symb_counter[s]+=1
    if mode=='C':
        list_comp=[]
        for i in weights:
            list_comp.append({})
        for i in C:
            list_comp[i[0]][i[1]] = C[i]
        return list_comp
    else:
        return D
    
    
def stream_ANS_encoding(to_encode):
    b = 2
    interval = range(2048,2048*b) #assuming number of states is 2**11, binary
    subinterval_listing, decoding_table = create_ans_tables(to_encode)
     
    val = interval[-1]  #initial value
    all_additional_output = []
    for i in to_encode:
        max_poss = max(subinterval_listing[i])
        list_encoding_to = subinterval_listing[i]
        if val>max_poss:
            while val > max_poss:
                all_additional_output.append(val%2)
                val = np.int(np.floor(val/2))        
        val = list_encoding_to[val]
    return val,all_additional_output,decoding_table


def stream_ANS_decoding(final_val,decoding_table,extra_outputs,max_bits=12):
    final_sequence = []
    val = final_val
    while (len(extra_outputs)>0): 
        final_sequence.append(decoding_table[val][0])
        val = decoding_table[val][1]
        for i in range(0,max_bits-np.int(np.ceil(np.log2(val))+(np.log2(val)%1==0))):
            val = val*2 + extra_outputs[-1]
            del extra_outputs[-1]
    final_sequence.reverse()        
    return final_sequence

def table_size(tables,weights,b):
    #2048 * (12 + 8)
    sz = 0
    a = np.ceil(np.log(max(tables[0])))
    c = np.ceil(np.log(weights[0]*b-1))
    for i in tables:
        sz = sz+len(i)*(a+c)
    return sz

def create_ans_tables(file_data):
    weights = np.zeros(256)
    number_of_states = 2048 #dont hardcode, adjust this depending on number of actual values in dictionary, 2**11
    unique_vals = np.unique(file_data)
    for i in unique_vals:
        weights[i] = file_data.count(i)
    
    rev_sort = sorted([(e,i) for i,e in enumerate(weights)],reverse = True)     #should be 12 bits at all time
    diction={}
    counter = 0
    for (i,j) in rev_sort:
        if i>0:
            diction[j] = [counter,i]
            counter = counter+1
    for i in range(0,len(file_data)):
        file_data[i] = diction[file_data[i]][0]       
    
    weights = sorted(weights,reverse = True)
    try:
        weights = weights[0:weights.index(0)]
    except ValueError:
        print('no 0 value')
    
    weights = np.array([max(1,np.int(np.round(i*number_of_states/sum(weights)))) for i in weights])   #rescaling weights
    while sum(weights)>number_of_states:
        weights = weights - 1
        weights = np.where(weights>0,weights,weights+1)
# try different ways of rescaling/ rebalancing the weights table            
    if sum(weights)<number_of_states:
        weights[0] = weights[0] + np.abs(sum(weights)-number_of_states)
    b = 2
    encoding_table = ANS_table(weights,b,'C',number_of_states)
    decoding_table = ANS_table(weights,b,'D',number_of_states)
    
    return encoding_table, decoding_table