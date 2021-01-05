#symbol is an integer
#counts is a list of weights
#dec_low is an integer
#dec_up is an integer
#E3_count is an integer
#code_index is an integer
#N is an integer
import numpy as np

class arith_encr:
    def __init__(self,ncode,counts,dec_low,dec_up,e3_count):
        self.ncode = ncode
        self.counts = counts
        self.dec_low = dec_low
        self.dec_up = dec_up
        self.E3_count = e3_count

class arith_decr:
    def __init__(self,dseq,counts,dec_low,dec_up,dseq_index,dec_tag,k):
        self.dseq = dseq
        self.counts = counts
        self.dec_low = dec_low
        self.dec_up = dec_up
        self.dseq_index = dseq_index
        self.dec_tag = dec_tag
        self.k = k

def entrp_enc(symbol,arith_obj,N):
    counts = arith_obj.counts
    dec_low = arith_obj.dec_low
    dec_up = arith_obj.dec_up
    E3_count = arith_obj.E3_count
    
    cum_counts = np.cumsum(np.append([0],counts))
    total_count = cum_counts[-1]
    
    dec_low_new = dec_low + np.floor( (dec_up-dec_low+1)*cum_counts[symbol]/total_count )
    dec_up = dec_low + np.floor( (dec_up-dec_low+1)*cum_counts[symbol+1]/total_count )-1
    dec_low = dec_low_new
    E12 = (dec_low>=2**(N-1)) | (dec_up<2**(N-1))
    E3a = (np.mod(np.floor(dec_low/(2**(N-2))),2)==1)
    E3b = (np.mod(np.floor(dec_up/(2**(N-2))),2)==0)
    ncode=[]
    
    
    while E12 | (E3a & E3b ):
        
        if E12:
            b = dec_low>=2**(N-1) #get most significant bit  of N bits
            ncode=np.append(ncode,b)
            # Left shifts
            dec_low = dec_low*2#bitshift(dec_low, 1) + 0;
            dec_up = dec_up*2+1#bitshift(dec_up, 1) + 1;
            
            # Check if E3_count is non-zero and transmit appropriate bits
            if (E3_count > 0):
                
                ncode = np.append(ncode,(np.ones(E3_count)*(b!=1))) #encoding opp of b, so yes?
                E3_count = 0

            # Reduce to N for next loop
            dec_low = np.mod(dec_low,2**N)
            dec_up  = np.mod(dec_up, 2**N)
            
        # Else if it is an E3 condition    
        elif E3a & E3b:
            
            # Left shifts
            dec_low = dec_low*2#bitshift(dec_low, 1) + 0;
            dec_up  = dec_up*2+1#bitshift(dec_up, 1) + 1;

            # Reduce to N for next loop
            dec_low = np.mod(dec_low, 2**(N))
            dec_up  = np.mod(dec_up, 2**(N))
            
            # Complement the new MSB of dec_low and dec_up
            if dec_low<2**(N-1):
                dec_low = dec_low+2**(N-1)
            else:
                dec_low = dec_low - 2**(N-1)
                
            if dec_up<2**(N-1):
                dec_up = dec_up+2**(N-1)
            else:
                dec_up = dec_up+2**(N-1)
#            dec_low = np.int(dec_low) ^ 2**(N-1) #************this may not be correct
#            dec_up  = np.int(dec_up) ^ 2**(N-1) 
            
            # Increment E3_count to keep track of number of times E3 condition is hit.
            E3_count = E3_count+1
        
        E12 = (dec_low>=2**(N-1)) | (dec_up<2**(N-1))
        E3a = (np.mod(np.floor(dec_low/(2**(N-2))),2)==1)
        E3b = (np.mod(np.floor(dec_up/(2**(N-2))),2)==0)  
#    update table
    counts[symbol] = counts[symbol]+1
    
    if np.sum(counts) == N:
        counts = np.ceil(counts/2)
        
    new_obj = arith_encr(np.append(arith_obj.ncode,ncode),counts,dec_low,dec_up,E3_count)
    
    return new_obj


#function [dseq, counts, dec_low, dec_up, dseq_index,dec_tag,k] = 
#    decomp_sing(code,counts,dec_low, dec_up, dseq_index,dec_tag,k,N)
def entrp_dec(code,arith_obj,N):
    counts = arith_obj.counts
    dec_low = arith_obj.dec_low
    dec_up = arith_obj.dec_up
    dseq_index = arith_obj.dseq_index
    dec_tag = arith_obj.dec_tag
    k = arith_obj.k
    cum_counts = np.cumsum(np.append([0],counts))
    total_count = cum_counts[-1]
    
    dec_tag_new =np.floor( ((dec_tag-dec_low+1)*total_count-1)/(dec_up-dec_low+1) )
    
    if dec_tag_new >= cum_counts[-1]:
        ptr = len(cum_counts)-2
    else:
        temp = np.array(cum_counts)<=dec_tag_new
        index = 0
        while temp[index]:
            index=index+1
        ptr = index-1

    dseq = ptr
    dseq_index = dseq_index + 1
    
    dec_low_new = dec_low + np.floor( (dec_up-dec_low+1)*cum_counts[ptr]/total_count )
    dec_up = dec_low + np.floor( (dec_up-dec_low+1)*cum_counts[ptr+1]/total_count )-1
    dec_low = dec_low_new

    E12 = (dec_low>=2**(N-1)) | (dec_up<2**(N-1))
    E3a = (np.mod(np.floor(dec_low/(2**(N-2))),2)==1)
    E3b = (np.mod(np.floor(dec_up/(2**(N-2))),2)==0)    
    while (E12 | ( E3a & E3b)):
        if ( k>=len(code)-1 ):
            break
        k = k + 1
        if E12:
            dec_low = dec_low*2 + 0
            dec_up  = dec_up*2 + 1
            dec_tag = dec_tag*2 + code[k]
            dec_low = np.mod(dec_low, 2**(N))
            dec_up  = np.mod(dec_up, 2**(N))
            dec_tag = np.mod(dec_tag, 2**(N))
        elif E3a & E3b:
            #Left shifts and update
            dec_low = dec_low*2 + 0;
            dec_up  = dec_up*2 + 1;

            dec_tag = dec_tag*2 + code[k];
            
            #Reduce to N for next loop
            dec_low = np.mod(dec_low, 2**(N))
            dec_up  = np.mod(dec_up, 2**(N))
            dec_tag = np.mod(dec_tag, 2**(N))
            
            #Complement the new MSB of dec_low, dec_up and dec_tag
            if dec_low<2**(N-1):
                dec_low = dec_low+2**(N-1)
            else:
                dec_low = dec_low - 2**(N-1)
                
            if dec_up<2**(N-1):
                dec_up = dec_up+2**(N-1)
            else:
                dec_up = dec_up-2**(N-1)
                
            if dec_tag<2**(N-1):
                dec_tag = dec_tag+2**(N-1)
            else:
                dec_tag = dec_tag-2**(N-1)
                
        E12 = (dec_low>=2**(N-1)) | (dec_up<2**(N-1))
        E3a = (np.mod(np.floor(dec_low/(2**(N-2))),2)==1)
        E3b = (np.mod(np.floor(dec_up/(2**(N-2))),2)==0)  
    counts[ptr] = counts[ptr]+1;
    cum_counts = np.cumsum(np.append([0],counts))
    total_count = cum_counts[-1]
    
    if total_count == N:
        counts = np.ceil(counts/2)

    new_obj = arith_decr(np.append(arith_obj.dseq,dseq),counts,dec_low,dec_up,dseq_index,dec_tag,k)
    return new_obj

def final_encoding(comp,N):    
    dlow_bin = np.array([])
    if len(bin(np.int(comp.dec_low)))-2<N:
        for i in range(0,N-len(bin(np.int(comp.dec_low)))+2):
            dlow_bin = np.append(dlow_bin,0)
    
    for i in range(2,len(bin(np.int(comp.dec_low)))):
        dlow_bin = np.append(dlow_bin,np.int(bin(np.int(comp.dec_low))[i]))
        
    if comp.E3_count ==0:
        comp.ncode = np.append(comp.ncode,dlow_bin)
    else:
        b = dlow_bin[0]
        comp.ncode = np.append(comp.ncode,b)
        comp.ncode = np.append(comp.ncode,(np.ones(comp.E3_count)*(b!=1))) #transmit complement of b, e3_count times
        comp.ncode = np.append(comp.ncode,dlow_bin[1:-1])
    return comp
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
