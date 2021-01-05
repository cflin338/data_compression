#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  7 15:37:42 2020

@author: linchri
"""
from ANS_Compression import create_ans_tables   #, table_size
from arith_single_encode import arith_encr, arith_decr, entrp_enc, entrp_dec, final_encoding
import numpy as np
import math
from bitarray import bitarray


def arith_encode_sequence(data,initial_weights=None,adaptive = True):
#bits_per indicates number of bits to represent each item
#initial_weights = weight of each item; if given, table must be of size 2^N
#if adaptive is false, weights will be reset to original after every item is encoded
    N = 20
    if initial_weights==None:
        initial_weights = np.ones(2**8)
        
    compression_object = arith_encr([], [],0,2**N-1,0)
    compression_object.counts = initial_weights
    for i in range(len(data)):
        compression_object = entrp_enc(data[i], compression_object, N)
        if adaptive==False:
            compression_object.counts = initial_weights
        
    compression_object = final_encoding(compression_object, N)
    encoded_value = compression_object.ncode
    return encoded_value

def arith_decode_sequence(encoded_value, pt_count,initial_weights=None,adaptive = True):
#N indicates number of bits to represent each item
#initial_weights = weight of each item; if given, table must be of size 2^N
#if adaptive is false, weights will be reset to original after every item is encoded
    N = 20
    tag = encoded_value[0:N]
    tagval = np.sum([tag[N-1-i]*2**i for i in range(N-1,-1,-1)])
    if initial_weights == None:
        initial_weights = np.ones(2**8)
    
    decompression_object = arith_decr([], [],0,2**N-1,0,tagval,N-1)
    decompression_object.counts = initial_weights
    
    final_list = np.zeros(pt_count,dtype = 'uint8')
    counter = 0
    for i in range(pt_count):
        decompression_object = entrp_dec(encoded_value, decompression_object, N)
        final_list[i] = decompression_object.dseq[counter]
        counter+=1
        if adaptive==False:
            decompression_object.counts = initial_weights
        
    return final_list

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

#def ans_Compression_size(val, decodingT, extra):
#    return 2048*(12+8)+len(extra)+
    
def RLEncoding(val_list):
    #need to add in a max possible run, will be based on number of bits used to record a run length
    counter = 1
    prev_val = val_list[0]
    val_list = val_list[1:]
    encoded = []
    for i in val_list:
        if i==prev_val:
            counter+=1
        else:
            encoded.append(prev_val)
            encoded.append(counter)
            prev_val = i
            counter = 1
    encoded.append(prev_val)
    encoded.append(counter)
    return encoded
    
def RLDecoding(val_list):
    l = int(len(val_list)/2)
    decoded = []
    for i in range(l):
        symb = val_list[2*i] 
        count = val_list[2*i+1]
        decoded += [symb]*count
    return decoded
        
def zero_encoding(val_list):
    mask = np.array(val_list)>0
    nonzero_vals = np.array(val_list)[mask]
    return np.append(mask,nonzero_vals)

def zero_decoding(mask, non_zero):
    locs = np.nonzero(mask)
    for i in range(len(locs)):
        mask[locs[i]] = non_zero[i]
    return mask

def LNVEnc(val_list,dist, metric):
    if metric=='xor':
        diff = val_list[0:dist]
        for i in range(dist,len(val_list)):
            prev_val = val_list[i-dist]
            diff.append(val_list^prev_val)
        return diff
    elif metric=='diff':
        diff = val_list[0:dist]
        for i in range(dist,len(val_list)):
            prev_val = val_list[i-dist]
            diff.append(val_list-prev_val)
        return diff
"""
def LNVDec(enc_list, dist, metric):
    if metric == 'xor':
        
    elif metric=='diff':
"""             
        
        

class LZ77Compressor:
    MAX_WINDOW_SIZE = 400
    def __init__(self, window_size=20):
        self.window_size = min(window_size, self.MAX_WINDOW_SIZE) 
        self.lookahead_buffer_size = 15 # length of match is at most 4 bits
        
    def compress(self, data,  verbose=False):
        output_buffer=bitarray(endian='big')
#        data = []
        i = 0
        while i < len(data):
			#print i
            match = self.findLongestMatch(data, i)
            if match: 
				# Add 1 bit flag, followed by 12 bit for distance, and 4 bit for the length
				# of the match 
                (bestMatchDistance, bestMatchLength) = match
                output_buffer.append(True)
                output_buffer.frombytes(bytes([bestMatchDistance >> 4]))
                output_buffer.frombytes(bytes([((bestMatchDistance & 0xf) << 4) | bestMatchLength]))
                if verbose:
                    print("<1, %i, %i>" % (bestMatchDistance, bestMatchLength))
                i += bestMatchLength
            else:
				# No useful match was found. Add 0 bit flag, followed by 8 bit for the character
                output_buffer.append(False)
                output_buffer.frombytes(bytes(data[i]))
                if verbose:
                    print("<0, %s>" % data[i])
                i += 1
        output_buffer.fill()

        return output_buffer


    def decompress(self, data):
        output_buffer = []
        while len(data) >= 9:
            flag = data.pop(0)
            if not flag:
                byte = data[0:8].tobytes()
                output_buffer.append(byte)
                del data[0:8]
            else:
                byte1 = ord(data[0:8].tobytes())
                byte2 = ord(data[8:16].tobytes())
                del data[0:16]
                distance = (byte1 << 4) | (byte2 >> 4)
                length = (byte2 & 0xf)
                for i in range(length):
                    output_buffer.append(output_buffer[-distance])
#        out_data =  ''.join(output_buffer)
        return output_buffer

    def findLongestMatch(self, data, current_position):
        end_of_buffer = min(current_position + self.lookahead_buffer_size, len(data) + 1, )
        best_match_distance = -1
        best_match_length = -1
        print('looking for longest match')
#        for j in range(current_position + 2, end_of_buffer):
        for j in range(end_of_buffer-1,current_position + 1,-1):
            start_index = max(0, current_position - self.window_size)
            substring = data[current_position:j]
            print(substring)
            matchFound = False
#            for i in range(start_index, current_position):
            for i in range(current_position-1, start_index-1,-1):
                matched_string = data[i:i+len(substring)]
                """
                repetitions = int(np.floor(len(substring) / (current_position - i)))
                last = len(substring) % (current_position - i)
                matched_string = data[i:current_position] * repetitions + data[i:i+last]
                """
#                print(matched_string)
                if matched_string == substring:# and len(substring) > best_match_length:
                    best_match_distance = current_position - i 
                    best_match_length = len(substring)
                    matchFound = True
                    break
            if matchFound:
                break
        print('end of window search')
        if best_match_distance > 0 and best_match_length > 0:
            print(best_match_distance, best_match_length)
            return (best_match_distance, best_match_length)
        return None


class lz78:

    def compress(data):
        dictionary, word = {0: ''}, 0
        dynamic_dictionary = lambda dictionary, key: dictionary.get(key) or dictionary.__setitem__(key,len(dictionary)) or 0
        return [token for char in data for token in [(word, char)] for word2 in [dynamic_dictionary(dictionary, token)] if
                not word] + [(word, '')]
    
    
    def decompress(data):
        dictionary = {0: ''}#, j = , ''.join
        dynamic_dictionary = lambda dictionary, value: dictionary.__setitem__(len(dictionary), value) or value
        return [dynamic_dictionary(dictionary, dictionary[codeword] + str(char)) for (codeword, char) in data]

def comp(data,lookahead, backwindow):
    comp_list = []
    comp_list.append(data[0])
    index = 1
    while index<len(data):
        print(data)
        print((index,len(data)))
        dist = backwindow*2
        max_len = -1
        back_start = np.max((index-backwindow,0))
        lookahead_end = np.min((index+lookahead+1,len(data)))
        if data[index] in data[back_start:index]:
            indices = np.where(np.array(data[back_start:index])==data[index])[0]
            for i in range(index+1,lookahead_end):
                val_to_check = data[index:i]
                print(val_to_check)
                
                for j in indices:
                    if data[back_start+j:back_start+j+len(val_to_check)]==val_to_check:
                        max_len = len(val_to_check)
                        dist = index - back_start + j

        if max_len>-1:
            index+=max_len
        else:
            index+=1
        comp_list.append((dist,max_len))
    return comp_list

if __name__ == '__main__':
    data = list(np.uint8([1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3]))# 
    compress_data = lz78.compress(data)
#    decompress_data = lz78.decompress(compress_data)
    
    a = LZ77Compressor()
    c = a.compress(data) 
    print(c)
    d = a.decompress(c)
    print(d)
    
    z = comp(data,10,10)
    print(z)
#    print('DATA:', data)
#    print('COMPRESSING:', compress_data)
#    print('DECOMPRESSING:', decompress_data)