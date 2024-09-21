def read_tsv(fn,pc,header,breakBool=False,sep='\t',pc_list=False):
    '''Read a tsv file'''
    with open(fn,'r') as f:
        
        # If printing columns, skip header
        if header:
            if (pc==True): pass
            else:          next(f)
                
        for i,line in enumerate(f):
            a=line.strip().split(sep)
            if pc:
                if pc_list==False:
                    if i==0:
                        for i,c in enumerate(a):
                            print(i,c)
                        print()
                        if breakBool: break
                        continue
                else:
                    if i==0:
                        print(', '.join([i.replace('-','_').replace(' ','_') for i in a]))
                        if breakBool: break
                        continue
            yield a

def revcomp(dna): 
	'''Takes DNA sequence as input and returns reverse complement'''
	inv={'A':'T','T':'A','G':'C','C':'G', 'N':'N','W':'W'}
	revcomp_dna=[]
	for nt in dna:
		revcomp_dna.append(inv[nt])
	return ''.join(revcomp_dna[::-1])

def zipdf(df,cols):
    return zip(*[df[c] for c in cols])

def dprint(d,n=0):
    '''Print a dictionary'''
    for i,(k,v) in enumerate(d.items()):
        if i<=n:
            print(k,v)

cb={}
cb['lightblue']= [i/255 for  i in [86,180,233]]
cb['green']    = [i/255 for  i in [0,158,115]]
cb['red']      = [i/255 for  i in [213,94,0]]
cb['yellow']   = [i/255 for  i in [240,228,66]]
cb['orange']   = [i/255 for  i in [230,159,0]]
cb['blue']     = [i/255 for  i in [0,114,178]]
cb['pink']     = [i/255 for  i in [204,121,167]]
cb['black']    = [i/255 for  i in [0,0,0]]
cb['orangenature'] = [i/255 for  i in [210,58,40]]

def loadAff(ref):
        '''Load an arbitrary affinity dataset. First column should be 8mer DNA 
sequence and second column should be the normalized affinity (between 0-1).'''
        Seq2EtsAff  = {line.split('\t')[0]:float(line.split('\t')[1]) for line in open(ref,'r').readlines()}
        return Seq2EtsAff