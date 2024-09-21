
print('loading libraries...')
import pandas as pd

import numpy as np
from random import choice
import subprocess
import pickle 
# from Bio import SeqIO

import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib_venn import venn2,venn3

print('done!')

####################################################################################
# Helper Objects 
####################################################################################


Iupac2AllNt= {
        'A':['A'],
        'C':['C'],
        'G':['G'],
        'T':['T'],
        'R':['A','G'],
        'Y':['C','T'],
        'S':['G','C'],
        'W':['A','T'],
        'K':['G','T'],
        'M':['A','C'],
        'B':['C','G','T'],
        'D':['A','G','T'],
        'H':['A','C','T'],
        'V':['A','C','G'],
        'N':['A','C','G','T'],
}

ATGC=set(['A','T','G','C'])

def GenerateSingleRandomSequence(template):
        '''Takes a string with letters A,T,C,G,IUPAC and returns a DNA string with rnadom AGTC nucleotides at position(s) with N.'''
        dna=''
        for i,nt in enumerate(template):
                if nt not in ATGC:
                        dna+=choice(Iupac2AllNt[nt])
                else:   
                        dna+=nt
        return dna

def revcomp(dna):
        '''Takes DNA sequence as input and returns reverse complement'''
        inv={'A':'T','T':'A','G':'C','C':'G', 'N':'N','W':'W'}
        revcomp_dna=[]
        for nt in dna:
                revcomp_dna.append(inv[nt])
        return ''.join(revcomp_dna[::-1])
    
####################################################################################
# Input/Output Helpers
####################################################################################

def to_pickle_iter_obj(iterObj,fnBasename,entriesPerFile):
    
    entriesPerFile=int(entriesPerFile)
    
    if entriesPerFile==None:
        with open(fnBasename+'.pickle','wb') as f:
            pickle.dump(iterObj,f,protocol=pickle.HIGHEST_PROTOCOL)
            
    else:
        
        num_chunks = (len(iterObj) + entriesPerFile - 1) // entriesPerFile
        
        for i in range(num_chunks):
            start_idx = i * entriesPerFile
            end_idx = (i + 1) * entriesPerFile
            
            if    type(iterObj)==dict: chunk_obj = {k: iterObj[k] for k in list(iterObj)[start_idx:end_idx]}
            elif  type(iterObj)==list: chunk_obj = iterObj[start_idx:end_idx]
            else:                      raise ValueError('This function only writes lists or dictionaries')

            with open(f'{fnBasename}-part{i}.pickle', 'wb') as f:
                pickle.dump(chunk_obj, f)
                
def get_files_in_directory(pattern):

    matching_files = glob.glob(pattern)
    return sorted(matching_files)    


def load_pickle_iter_obj(fnPattern,objType):
    
    if objType==dict:
        outObj = {}

        for filename in get_files_in_directory(fnPattern):
            with open(os.path.join(filename), 'rb') as f:
                part_dict = pickle.load(f)
                outObj.update(part_dict)
    
    elif objType==list:
        outObj=[]

        for filename in get_files_in_directory(fnPattern):
            with open(os.path.join(filename), 'rb') as f:
                part_list = pickle.load(f)
                outObj+=part_list
        
    else:
        raise ValueError('This function only loads lists or dictionaries')

    return outObj

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
    
def get_fastq_reads(fn):
    '''Iterates over sequences in a fastq file.'''
    with open(fn,'r') as f:
        readFlag=False
        for line in f:
            if line[0]=='@':
                readFlag=True
                continue
            if readFlag:
                yield line.strip()
                readFlag=False
                
def get_sam_reads(fn):
    '''Iterates over sequences in a sam file.'''
    for row in read_tsv(hifiReads,pc=False,header=False):
        # skip headers
        if row[0][0]=='@': continue
        # get read
        read=row[readCol0idx]
        yield read
        
def write_row(rowList,delim='\t'):
    '''Write a single row of a tsv file.'''
    return delim.join([str(i) for i in rowList])+'\n'


        
def percent(number,rounding_digit=1):
    '''Get percent of fraction'''
    if rounding_digit==0:
        return str(int(100*number))+'%'
    else:
        return str(round(100*number,rounding_digit))+'%'
    
def read2count_from_tsv(tsv):
    '''Converts read2count from tsv to python dictionary.'''
    read2count={}
    for line in read_tsv(tsv,pc=False,header=False):
        if len(line)==2:
            read,count=line
            read2count[read]=float(count)
        else:
            print(f'len(line)!=2... skipping {line}')
    return read2count

def  readsWithEnhancerAndBarcode_from_tsv(tsv):
    '''Converts readsWithEnhancerAndBarcode from tsv to list of tuples.'''
    readsWithEnhancerAndBarcode=[]
    for en,bc,count in read_tsv(tsv,pc=False,header=False):
        readsWithEnhancerAndBarcode.append((en,bc,float(count)))
    return readsWithEnhancerAndBarcode

def bc2en2count_from_tsv(bc2en2count_tsv):
    '''Converts bc2en2count from tsv file and returns python dictionary.'''
    bc2en2count={}
    for bc,en,count in read_tsv(bc2en2count_tsv,pc=False,header=False):
        if bc not in bc2en2count:
            bc2en2count[bc]={}
        bc2en2count[bc][en]=float(count)
    return bc2en2count

def ubc2en_from_tsv(ubc2en_tsv):
    '''Converts ubc2en (unique barcode to enhancer) from tsv to python dictionary.'''
    ubc2en={}
    for ubc,en in read_tsv(ubc2en_tsv,pc=False,header=False):
        if ubc not in ubc2en:
            ubc2en[ubc]={}
        ubc2en[ubc]=en
    return ubc2en

def enhancerOrderSheet_from_tsv(enOrderSheet):
    
    EnOrderName2Seq={}
    EnOrderSeq2Name={}
    for row in read_tsv(enOrderSheet,pc=False,header=True):
        name,seq=row
        # seq=seq[25:-25] # delete [25:-25] if no sequencing adapters a re in the order sheet
        EnOrderName2Seq[name]=seq
        EnOrderSeq2Name[seq ]=name
        
    return EnOrderName2Seq,EnOrderSeq2Name
    
def get_str_million(num):
    return str(round(num/1e6,1)).replace('.','-').replace('-0','')+'M'


####################################################################################
# Dictionary analysis
####################################################################################

# Collapse reads
def dict_step1_raw_to_collapsed_reads(rawReadFile,
                                      fileType,
                                      outfn,
                                      readLengthMax=None,
                                      readCol0idx=None,
                                      printProgress=True,
                                      fastqTotalLines=None,
                                      progressReportInterval=1e6,
                                      numberReadsPerFile=False):
    
    '''Converts raw duplicate read file to a tsv with each unique read sequence and associated read count.
    
    Only use readCol0idx if inputting sam file.
    
    Saves a new tsv. Returns read2count.'''
    
    if numberReadsPerFile:
        if '{part}' not in outfn: raise ValueError('If you specify numberReadsPerFile, then you must include "{part}" in the outfn to indicate which output part it is. ')

    # Determine how to read file
    if fileType=='sam':
        readFileFunc=get_sam_reads
    elif fileType=='fastq':
        readFileFunc=get_fastq_reads
    else:
        raise ValueError(f'fileType must be either "fastq" or "sam"')
    
    # Determine size of file
    if printProgress:
        if fileType=='fastq':
            if not fastqTotalLines:
                print(f'Determining file size of {rawReadFile}...')
                fastqTotalReads=int(int(subprocess.getoutput(f"wc -l < {rawReadFile} | bc"))/4)
                print(f'\t{fastqTotalReads:,} reads detected')
            else: 
                fastqTotalLines=int(fastqTotalLines)
                fastqTotalReads=fastqTotalLines
                print(f'\n\nFile size of {rawReadFile} specified by user to be {fastqTotalLines:,}...')
    else:
        print('User does not want to print progress. To report progress, set "printProgress" to True')
        
    # only use readCol0idx if sam file
    if readCol0idx!=None and fileType=='fastq': 
        raise ValueError('readCol0idx not used if fileType=="fastq"')

    print('\nReading fastq...\n\t0M (0%), ',end='')
    read2count={}
    numReadsTooLong=0
    numReadsWithN=0
    readsProcessed=0
    for read in readFileFunc(rawReadFile):
        
        if printProgress:
            readsProcessed+=1
            if readsProcessed%progressReportInterval==0: print(f'{int(readsProcessed/1e6):,}M ({percent(readsProcessed/fastqTotalReads)}), ',end='')

        # if read longer than max, skip
        if readLengthMax:
            if len(read)>readLengthMax: 
                numReadsTooLong+=1
                continue
        
        if 'N' in read: 
            numReadsWithN+=1
            continue
        
        # Try adding fwd
        if read in read2count:
            read2count[read]+=1

        # Try adding rev
        elif revcomp(read) in read2count:
            read2count[revcomp(read)]+=1

        # If fwd/rev absent, add fwd initialized at 1
        else: 
            read2count[read]=1
    
    print('\nWriting to tsv...')
    if not numberReadsPerFile:
        totokuniqreads=0
        totokreadcount=0
        line_out=''
        for read,count in read2count.items():
            totokreadcount+=count
            totokuniqreads+=1
            line_out+=write_row([read,count])
        with open(outfn,'w') as f: f.write(line_out)
    else:
        suffix=outfn.split('/')[-1].split('.')[0]
        part=0
        totokuniqreads=0
        totokreadcount=0
        readcount=0
        line_out=''
        for read,count in read2count.items():
            if read=='': continue
            readcount+=1
            totokreadcount+=count
            totokuniqreads+=1
            
            if readcount>numberReadsPerFile:
                readcount=0
                part+=1
                with open(outfn.replace('{part}',str(int(part))),'w') as f: f.write(line_out)
                line_out=''
                
            line_out+=write_row([read,count])
            
        part+=1
        with open(outfn.replace('{part}',str(int(part))),'w') as f: f.write(line_out)
    
    # print('\nWriting to pickle...')
    # to_pickle(read2count,outfn+'.read2count.pydict.pickle')
    
    print(f'\nNumber reads excluded due to containing N:   {numReadsWithN:,}')
    print(f'Number reads excluded due to being too long: {numReadsTooLong:,}')
    print()
    print(f'Number Unique Reads OK: {totokuniqreads:,}')
    print(f'Number Reads OK:        {totokreadcount:,}')
    print()
    print(f'=> read2count tsv file outputted here \n\t{outfn}')
    
    return read2count
    
# # jsd, em
# def dict_step2SR_parse_en_bc_from_read(enOrderSheet,
#                                        read2count,
#                                        bcLength,
#                                        relPositionsToLookForLinker,
#                                        linkerSeq,
#                                        outfn,
#                                        returnRevCompEn=False,
#                                        returnRevCompBc=False,
#                                        locationOnlyDoNotLookForLinker=False):
    
#     # read enhancer order sheet
#     print(f'Loading order sheet...\n\t{enOrderSheet}\n')
#     EnOrderName2Seq,EnOrderSeq2Name=enhancerOrderSheet_from_tsv(enOrderSheet)
#     if returnRevCompEn:
#         EnOrderSeq2Name={revcomp(seq):name for seq,name in EnOrderSeq2Name.items()}
#         EnOrderName2Seq={name:revcomp(seq) for name,seq in EnOrderName2Seq.items()}
#     orderedEnhancers= set(EnOrderSeq2Name.keys())
    
#     linkerLen=len(linkerSeq)
    
#     # look where you expect to find the linker sequence
#     relPositionsToLookForLinker=[0]+sorted(relPositionsToLookForLinker,key=lambda i: abs(i)) # reorder to look at most common locations first 
#     relPositionsToLookForLinker=[bcLength+i for i in relPositionsToLookForLinker]

#     totalReads = sum(read2count.values())
#     nreads_withOrderedEn = 0
#     nreads_withOrderedEn_withLinkerSeq = 0
#     readsWithEnhancerAndBarcode=[]
#     if locationOnlyDoNotLookForLinker==False:
#         print(f'Looking for enhancers in the following location (bcLen + relPos + linkerLen)...\n\t{relPositionsToLookForLinker}\n')
#     else:
#         print(f'locationOnlyDoNotLookForLinker={locationOnlyDoNotLookForLinker}, so we assume all reads with enhancers at expected position have the barcode at positions 0-{bcLength}.')

#     # case 1 - look for linker sequence, parse barcode behind it
#     if relPositionsToLookForLinker==False:
#         for read,count in read2count.items():
#             for linkerLoc in relPositionsToLookForLinker:
#                 enLoc=linkerLoc+linkerLen
#                 en=read[enLoc:]

#                 # if enhancer is in the order sheet...
#                 if en in orderedEnhancers: 
#                     nreads_withOrderedEn+=count

#                     # if linker seq is right, extract the enhancer,barcode and read
#                     obsLinkerSeq=read[linkerLoc:enLoc]
#                     if obsLinkerSeq==linkerSeq: 
#                         nreads_withOrderedEn_withLinkerSeq+=count

#                         bc=read[:linkerLoc]

#                         if returnRevCompBc: bc=revcomp(bc)
#                         if returnRevCompEn: en=revcomp(en)

#                         readsWithEnhancerAndBarcode.append((en,bc,count))
#     # case 2 - just use location based indexing to parse linker and barcode
#     else:
#     # ------------------------------|------|
#     # bc       1         2         3.link  . 4         5                                                         
#     # 12345678901234567890123456789012345678901234567890123456789012345678901234567890123456789                                                                                  
#     # 012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789
#         for read,count in read2count.items():
#             for linkerLoc in relPositionsToLookForLinker[0]: # relPositionsToLookForLinker[0] enforces that you only loook at the first position
#                 enLoc=linkerLoc+linkerLen
#                 en=read[enLoc:]

#                 # if enhancer is in the order sheet...
#                 if en in orderedEnhancers: 
#                     nreads_withOrderedEn+=count

#                     # case 2 - assume linker is right, get barcode where expected
#                     bc=read[:linkerLoc]

#                     if returnRevCompBc: bc=revcomp(bc)
#                     if returnRevCompEn: en=revcomp(en)

#                     readsWithEnhancerAndBarcode.append((en,bc,count))
                    
#     if locationOnlyDoNotLookForLinker==False:
#         print(f'{nreads_withOrderedEn:,} ({percent(nreads_withOrderedEn/totalReads)}) reads with ordered enhancer...')
#         print(f'{nreads_withOrderedEn_withLinkerSeq:,} ({percent(nreads_withOrderedEn_withLinkerSeq/totalReads)}) reads with ordered enhancer and linker...')
#     else:
#         print(f'{nreads_withOrderedEn:,} ({percent(nreads_withOrderedEn/totalReads)}) reads with ordered enhancer...')
    
#     line_out='\n'.join(['\t'.join([str(li) for li in line]) for line in readsWithEnhancerAndBarcode])
#     with open(outfn,'w') as f: f.write(line_out)
    
#     print()
#     print(f'=> Parsed enhancer-barcode reads outputted here \n\t{outfn}')
    
#     return readsWithEnhancerAndBarcode
                
                
    

def dict_step2SR_parse_en_bc_from_read(enOrderSheet,
                                       read2count,
                                       bcLength,
                                       relPositionsToLookForLinker,
                                       linkerSeq,
                                       outfn,
                                       returnRevCompEn=False,
                                       returnRevCompBc=False,
                                       locationOnlyDoNotLookForLinker=False,
                                       writeOutLinesByChunkSize=500e3):
    
    # read enhancer order sheet
    print(f'Loading order sheet...\n\t{enOrderSheet}\n')
    EnOrderName2Seq,EnOrderSeq2Name=enhancerOrderSheet_from_tsv(enOrderSheet)
    if returnRevCompEn:
        EnOrderSeq2Name={revcomp(seq):name for seq,name in EnOrderSeq2Name.items()}
        EnOrderName2Seq={name:revcomp(seq) for name,seq in EnOrderName2Seq.items()}
    orderedEnhancers= set(EnOrderSeq2Name.keys())
    
    linkerLen=len(linkerSeq)
    
    # look where you expect to find the linker sequence
    relPositionsToLookForLinker=[0]+sorted(relPositionsToLookForLinker,key=lambda i: abs(i)) # reorder to look at most common locations first 
    relPositionsToLookForLinker=[bcLength+i for i in relPositionsToLookForLinker]

    totalReads = sum(read2count.values())
    nreads_withOrderedEn = 0
    nreads_withOrderedEn_withLinkerSeq = 0
    readsWithEnhancerAndBarcode=[]
    if locationOnlyDoNotLookForLinker==False:
        print(f'\nLooking for enhancers in the following location (bcLen + relPos + linkerLen)...\n\t{relPositionsToLookForLinker}\n')
    else:
        print(f'\nlocationOnlyDoNotLookForLinker={locationOnlyDoNotLookForLinker}, so we assume all reads with enhancers at expected position have the barcode at positions 0-{bcLength}.')

    # case 1 - look for linker sequence, parse barcode behind it
    if locationOnlyDoNotLookForLinker==False:
        for read,count in read2count.items():
            for linkerLoc in relPositionsToLookForLinker:
                enLoc=linkerLoc+linkerLen
                en=read[enLoc:]

                # if enhancer is in the order sheet...
                if en in orderedEnhancers: 
                    nreads_withOrderedEn+=count

                    # if linker seq is right, extract the enhancer,barcode and read
                    obsLinkerSeq=read[linkerLoc:enLoc]
                    if obsLinkerSeq==linkerSeq: 
                        nreads_withOrderedEn_withLinkerSeq+=count

                        bc=read[:linkerLoc]

                        if returnRevCompBc: bc=revcomp(bc)
                        if returnRevCompEn: en=revcomp(en)

                        readsWithEnhancerAndBarcode.append((en,bc,count))
    # case 2 - just use location based indexing to parse linker and barcode
    else: 
        for read,count in read2count.items():
            for linkerLoc in [relPositionsToLookForLinker[0]]: # relPositionsToLookForLinker[0] enforces that you only loook at the first position
                enLoc=linkerLoc+linkerLen
                en=read[enLoc:]

                # if enhancer is in the order sheet...
                if en in orderedEnhancers: 
                    nreads_withOrderedEn+=count

                    # case 2 - assume linker is right, get barcode where expected
                    bc=read[:linkerLoc]

                    if returnRevCompBc: bc=revcomp(bc)
                    if returnRevCompEn: en=revcomp(en)

                    readsWithEnhancerAndBarcode.append((en,bc,count))
                    
    if locationOnlyDoNotLookForLinker==False:
        print(f'\n{nreads_withOrderedEn:,} ({percent(nreads_withOrderedEn/totalReads)}) reads with ordered enhancer...')
        print(f'{nreads_withOrderedEn_withLinkerSeq:,} ({percent(nreads_withOrderedEn_withLinkerSeq/totalReads)}) reads with ordered enhancer and linker...')
    else:
        print(f'\n{nreads_withOrderedEn:,} ({percent(nreads_withOrderedEn/totalReads)}) reads with ordered enhancer...')
        
    # write out as pickled object
    # outpicklelist=outfn+'.readsWithEnhancerAndBarcode.pylist.pickle'
    # if '.tsv' in outpicklelist: outpicklelist=outpicklelist.replace('.tsv','')
    # print(f'\nwriting pickled list to {outpicklelist}...')        
    # with open(outpicklelist,'wb') as f: pickle.dump(readsWithEnhancerAndBarcode,f,protocol=pickle.HIGHEST_PROTOCOL)
    
    # write out as tsv
    print(f'\nwriting out reads with enhancer and barcode to {outfn}...')
    line_out=''
    #with open(outfn,'w') as f: f.write(line_out)
    lc=0
    for line in readsWithEnhancerAndBarcode:
        lc+=1

        line='\t'.join([str(li) for li in line])
        line_out+=line+'\n'
        
        if lc==writeOutLinesByChunkSize:
            with open(outfn,'a') as f: f.write(line_out)
            line_out=''
            lc=0

    # write last chunk        
    with open(outfn,'a') as f: f.write(line_out)
    
    # print('\nWriting to pickle...')
    # to_pickle(readsWithEnhancerAndBarcode,outfn+'.readsWithEnhancerAndBarcode.pylist.pickle')
    
    print('\nDone!')    
    return readsWithEnhancerAndBarcode
        
def plot_barcode_length_dist(readsWithEnhancerAndBarcode):
    
    # Check bc length distribution
    data=[len(bc) for en,bc,readcount in readsWithEnhancerAndBarcode]
    data=pd.Series(data).value_counts()
    fig,ax=plt.subplots(1,dpi=150)
    xmin=bcLength-10
    xmax=bcLength+10
    plt.hist(data,bins=range(xmin,xmax))
    plt.xlabel('Barcode Length')
    ax.set_xticks([int(i) for i in range(xmin,xmax,5)])
    


def dict_step3_assemble_dict(readsWithEnhancerAndBarcode,enOrderSheet,outDictBasename):
    
    '''Takes parsed enhancers and barcodes and assembles into pythonic enhancer-barcode dictionary.
    
    Saves several dictionary files formatted as tsvs and returns bc2en2count python dictionary.'''
    
    ####################################################
    # Get Enhancers ordered
    ####################################################
    
    EnOrderName2Seq,EnOrderSeq2Name=enhancerOrderSheet_from_tsv(enOrderSheet)
    
    ####################################################
    # Assemble en-bc dictionaries
    ####################################################
    
    print('Reading readsWithEnhancerAndBarcode...')
    en2bcset={}
    bc2enset={}
    bc2en2count={}
    for en,bc,count in readsWithEnhancerAndBarcode:

        if en not in en2bcset:
            en2bcset[en]=set()

        if bc not in bc2enset:
            bc2enset[bc]=set()
            bc2en2count[bc]={}

        if en not in bc2en2count[bc]:
            bc2en2count[bc][en]=0

        bc2en2count[bc][en]+=count

        en2bcset[en].add(bc)
        bc2enset[bc].add(en)
        
    # Report how many enhancers are missing
    orderedEns=set(EnOrderSeq2Name.keys())
    observedEns=set(en2bcset.keys())
    nMissingEnhancers=len(orderedEns-observedEns)
    
    # Report which enhancers are missing
    if nMissingEnhancers>0:
        # print('\nMissing Enhancers:')
        missingEnhancers=orderedEns-observedEns
        percentmissing=percent(nMissingEnhancers/len(orderedEns))
        print(f'\nMissing {nMissingEnhancers} ({percentmissing}) enhancers...')
        # for missingEn in missingEnhancers:
            # print(f'\t{EnOrderSeq2Name[missingEn]}')
            
    
    ##############################################################################
    # Write out all barcodes
    ##############################################################################


    # outfn=f'{outDictBasename}_en2bcset.tsv'
    # print(f'\nWriting {outfn}...')
    # line_out=''
    # for en,bcset in en2bcset.items():
    #     bcset=','.join(list(bcset))
    #     if len(bcset)>1:
    #         line_out+=write_row([en,bcset])
    #     else:
    #         line_out+=write_row([en,bcset])
    # with open(outfn,'w') as f: f.write(line_out)
    # print('\t',outfn)
        
    # line_out=''
    # outfn=f'{outDictBasename}_bc2enset.tsv'
    # print(f'\nWriting {outfn}...')
    # for bc,enset in bc2enset.items():
    #     enset=','.join(list(enset))
    #     line_out+=write_row([bc,enset])
    # with open(outfn,'w') as f: f.write(line_out)
    # print('\t',outfn)
    
    line_out=''
    outfn=f'{outDictBasename}_bc2en2count.tsv'
    print(f'\nWriting {outfn}...')
    for bc,en2count in bc2en2count.items():
        for en,count in en2count.items():
            line_out+=write_row([bc,en,count])
    with open(outfn,'w') as f: f.write(line_out)
    print('\t',outfn)
    
    # print('\nWriting to pickle...')
    # to_pickle(bc2en2count,outfn+'.pydict.pickle')
        
    return bc2en2count,missingEnhancers

def visualize_barcode_bc2en2count(bc2en2count,yscale='linear'):
    
    ##############################################################################
    # Report enhancer barcode readcount histogram
    # and... Report readcounts of top enhancer vs second second enhancer for mm barcodes
    ##############################################################################
    
    x=[]
    y=[]
    
    readcounts=[]
    
    
    for bc,en2count in bc2en2count.items():

        countEnTuples=[]
        for en,count in en2count.items():
            countEnTuples.append((count,en))
            readcounts.append(count)

        if len(countEnTuples)<2: continue

        countEnTuples=sorted(countEnTuples,reverse=True)

        topCount,secCount= [count for count, en in countEnTuples[:2]]

        x.append(topCount)
        y.append(secCount)
        
    ######################
    # readcount hist
    ######################
    print('\nPlotting enhancer barcode readcount distribution...')
    fig,ax=plt.subplots(1,figsize=(5,5),dpi=150)

    ax.hist(readcounts,bins=50)

    ax.set_xlabel('Enhancer-Barcode Readcount')
    ax.set_ylabel('Frequency')
    ax.set_title('Enhancer-Barcode Readcount Distribution')
    ax.set_yscale(yscale)

    plt.show()
    
    ######################
    # multi match readcounts
    ######################
    
    # if no multiple match barcodes
    if x==[] and y==[]: 
        print('\nNo multiple match barcodes detected.')


    # if multiple match, show discrepency in readcount from top enhancer to second
    else:
        print('\nPlotting multi match readcount discrepansy...')
        fig,ax=plt.subplots(1,figsize=(5,5),dpi=150)

        ax.scatter(x,y,alpha=.7,s=5)

        ax.set_xlabel('Top enhancer read counts')
        ax.set_ylabel('Second enhancer read counts')
        ax.set_title('Multiple Match Barcode\nEnhancer Readcount Discrepansy')

        plt.show()



# def dict_step4_filter_dictionary_for_unique_and_min_reads(bc2en2count,minReads,outDictBasename,bcTrimLength=None):
#     '''Filters bc2en2count python dictionary and returns only unique barcodes. Also allows for a minimum readcount for barcode-enhancer pairs.
    
#     Returns ubc2en (unique barcode to enhancer) python dictionary.'''

#     if bcTrimLength: print(f'Trimming barcodes to {bcTrimLength}...')
    
#     filtbc2enset={}
#     for bc,en2count in bc2en2count.items():
        
#         if bcTrimLength: bc=bc[:bcTrimLength]
        
#         for en,count in en2count.items():
            
#             if count>minReads:
                
#                 if bc not in filtbc2enset:
#                     filtbc2enset[bc]=set()
                
#                 filtbc2enset[bc].add(en)
               
#     ubc2en={bc:list(enset)[0] for bc,enset in filtbc2enset.items() if len(enset)==1}
#     en2ubcset={}
#     for ubc,en in ubc2en.items():
#         if en not in en2ubcset:
#             en2ubcset[en]=set()
#         en2ubcset[en].add(ubc)
        
#     if not bcTrimLength:  print(f'{len(ubc2en):,} Unique barcodes saved with >={minReads} reads')
#     else:                 print(f'{len(ubc2en):,} Unique barcodes saved with >={minReads} reads and trimmed to {bcTrimLength} nt')
#     # plot number of bc per en
#     data=[len(ubcSet) for ubcSet in en2ubcset.values()]
#     mean_bcperen=np.mean(data)
#     fig,ax=plt.subplots(1,dpi=150)
#     ax.hist(data,bins=20)
#     ax.set_xlabel('# bc / en')
#     ax.set_title('Number barcodes per enhancer')
#     ax.set_ylabel('Frequency')
#     ax.axvline(mean_bcperen,ls='--',color='red',alpha=.7)
#     plt.show()
    
#     # write out ubc2en
#     line_out=''
#     for ubc,en in ubc2en.items():
#         line_out+=write_row([ubc,en])
        
#     if not bcTrimLength:
#         outfn=f'{outDictBasename}_ubc2en_filt_minReads={minReads}.tsv'
#     else: 
#         outfn=f'{outDictBasename}_ubc2en_filt_minReads={minReads}_bcTrimLength={bcTrimLength}.tsv'
        
#     with open(outfn,'w') as f: f.write(line_out)
#     print(f'\n=> Outputted ubc2en to \n\t{outfn}')
#     return ubc2en

def dict_step4_filter_dictionary(bc2en2count,
                                 outDictBasename,
                                 minUReadsRequired=0,
                                 maxMmReadsAcceptable=np.inf,
                                 bcTrimLength=None):
    
    '''Filters bc2en2count python dictionary and returns only unique barcodes. Also allows for a minimum readcount for barcode-enhancer pairs.
    
    Returns ubc2en (unique barcode to enhancer) python dictionary.'''

    if bcTrimLength: print(f'Trimming barcodes to {bcTrimLength}...')
    
    plotmm,plotu,plotcolors=[],[],[]
    filtbc2en2count={}
    ubc2en={}
    en2ubcset={}
    for bc,en2count in bc2en2count.items():
        
        if bcTrimLength: bc=bc[:bcTrimLength]
        
        numEnMapped=len(en2count)
        if numEnMapped==1:  bctype='u'
        else:               bctype='mm'
        
        if bc not in filtbc2en2count:
            filtbc2en2count[bc]={}
        
        countEnTuples=sorted([(c,e) for e,c in en2count.items()],reverse=True) # reverse=True makes most abundant enhancers first
        primaryCount,primaryEn=countEnTuples[0]
        
        # either way, the unique count must be satisfied
        if primaryCount>=minUReadsRequired:
            
            # if u
            if bctype=='u':
                filtbc2en2count[bc][primaryEn]=primaryCount
                ubc2en[bc]=primaryEn
                if primaryEn not in en2ubcset:
                    en2ubcset[primaryEn]=set()
                en2ubcset[primaryEn].add(bc)
                plotmm.append(0)
                plotu.append(primaryCount)
                plotcolors.append('green')
                
            # if mm
            if bctype=='mm':
                segundoCount,segundoEn=countEnTuples[1]
                if segundoCount<maxMmReadsAcceptable:
                    filtbc2en2count[bc][primaryEn]=primaryCount
                    ubc2en[bc]=primaryEn
                    if primaryEn not in en2ubcset:
                        en2ubcset[primaryEn]=set()
                    en2ubcset[primaryEn].add(bc)
                
                    plotmm.append(segundoCount)
                    plotu.append(primaryCount)
                    plotcolors.append('blue')
                    
                else:
                    plotmm.append(segundoCount)
                    plotu.append(primaryCount)
                    plotcolors.append('red')
        else:
            # if u
            if bctype=='u':
                plotmm.append(0)
                plotu.append(primaryCount)
                plotcolors.append('grey')
                
            # if mm
            if bctype=='mm':
                segundoCount,segundoEn=countEnTuples[1]
                if segundoCount<=maxMmReadsAcceptable:
                    
                    plotmm.append(segundoCount)
                    plotu.append(primaryCount)
                    plotcolors.append('grey')
                    
                else:
                    plotmm.append(segundoCount)
                    plotu.append(primaryCount)
                    plotcolors.append('grey')
    
    data=[len(ubcSet) for ubcSet in en2ubcset.values()]
    mean_bcperen=np.mean(data)
    
    # plot primary vs secondary en counts for mm

    # write out ubc2en
    line_out=''
    for ubc,en in ubc2en.items():
        line_out+=write_row([ubc,en])
    
    print(f'{len(ubc2en):,} Unique barcodes outputted... ')
    
    # create fn
    milBarcodes=get_str_million(len(ubc2en))
    outfn=f'{outDictBasename}_ubc2en_filt_{milBarcodes}'
    if minUReadsRequired!=0: 
        outfn+=f'.minUReadsRequired={minUReadsRequired}'
        print(f'\tu  barcodes ignored if reads < {minUReadsRequired} (minUReadsRequired)')
        
    if maxMmReadsAcceptable!=np.inf: 
        outfn+=f'.maxMmReadsAcceptable={maxMmReadsAcceptable}'
        print(f'\tmm barcodes ignored if reads <  {maxMmReadsAcceptable} (maxMmReadsAcceptable)')
        
    if bcTrimLength!=None: 
        outfn+=f'.bcTrimLength={bcTrimLength}'
        print(f'\tbcTrimLength         = {bcTrimLength}')
    outfn+='.tsv'
        
    with open(outfn,'w') as f: f.write(line_out)
    print(f'\n=> Outputted ubc2en to \n\t{outfn}')
    
    return ubc2en,filtbc2en2count



####################################################################################
# Table generation
####################################################################################

def exp_1_count_barcodes(in_fastq,out_path,bcLength,out_basename=None,progressReport=1e6,fastqTotalLines=None,excludeBarcodesWithN=False):

    # Create default basename if one isn't provided
    if not out_basename: out_basename=in_fastq.split('/')[-1]
    
    # Determine size of file
    if not fastqTotalLines:
        print(f'Determining file size of {in_fastq}...')
        fastqTotalLines=int(subprocess.getoutput(f"wc -l < {in_fastq} | bc"))
    else: 
        print(f'File size of {in_fastq} specified by user to be...')
        
    print(f'\t# Lines = {fastqTotalLines:,}')
    print(f'\t# Reads = {int(fastqTotalLines/4):,}')

    
    print(f'\nReading {in_fastq}...\n\t',end='')
    
    out_basename=f'{out_path}/{out_basename}'

    readsToProcess=fastqTotalLines/4
    readsProcessed=0
    readsSkippedWithN=0
    Bc2ReadCount={}
    for seq in get_fastq_reads(in_fastq):
            
            readsProcessed+=1
            
            # report progress
            if readsProcessed%progressReport==0: print(f'{int(readsProcessed/1e6):,}M ({percent(readsProcessed/readsToProcess)}), ',end='')
                
            bc=seq[:bcLength]
            
            if excludeBarcodesWithN:
                if 'N' in bc: 
                    readsSkippedWithN+=1
                    continue
                
            if bc not in Bc2ReadCount:
                Bc2ReadCount[bc]=0
            Bc2ReadCount[bc]+=1

    if excludeBarcodesWithN:
        print(f'\n\n{readsSkippedWithN:,} ({percent(readsSkippedWithN/readsProcessed)}) reads with N skipped in barcode region...')
    print(f'\nWriting...')
    
    # write as pickled python dict for quick loading
    out_fn=f'{out_basename}.Bc2ReadCount.pickle'
    print(f'\t{out_fn}')
    with open(out_fn, 'wb') as handle: pickle.dump(Bc2ReadCount, handle, protocol=pickle.HIGHEST_PROTOCOL)
          
    # write as tsv
    out_fn=f'{out_basename}.Bc2ReadCount.tsv'
    print(f'\t{out_fn}')
    line_out='\n'.join(['\t'.join([bc,str(rc)]) for bc,rc in Bc2ReadCount.items()])
    with open(out_fn, 'w') as f: f.write(line_out)
    
    return Bc2ReadCount
    
def exp_2_create_barcode_table(SampleName2Data, enOrderSheet, ubc2en, outbasename=None):
    
    _,En2EnId=enhancerOrderSheet_from_tsv(enOrderSheet)
        
    if len(En2EnId.values())!=len(set(En2EnId.values())): 
        raise ValueError ('Enhancer Ids are not unique. No duplicate enhancer Ids allowed.')
    
    # Initialize dataframe
    Col2Values={'barcode-seq':[],'enhancer-seq':[],'enhancer-id':[]}
    
    print('Initializing barcode table with dictionary sequences...')
    for bc in ubc2en.keys():
        en=ubc2en[bc]
        Col2Values['barcode-seq'].append(bc)
        Col2Values['enhancer-seq'].append(en)
        Col2Values['enhancer-id'].append(En2EnId[en])
        
    numbctot=len(Col2Values['barcode-seq'])
    print(f'\n\t{numbctot:,} barcodes detected in the dictionary')
        
    sampleNames=list(SampleName2Data.keys())
    
    for sampleName in sampleNames:
        
        Bc2ReadCount=SampleName2Data.pop(sampleName) # remove from memory as you iterate through because you may run into memory issues for larger libs
        
        print(f'\nAdding barcode readcounts from {sampleName}...')
        
        colRc='rc-'+sampleName
        colRpm='rpm-'+sampleName
        Col2Values[colRc]=[]
        Col2Values[colRpm]=[]
        
        readsTotal         = sum(Bc2ReadCount.values())
        readsWithDictBc    = 0
        numBcTotal         = len(Bc2ReadCount)
        numBcInDict        = 0
        numBcNotIndict     = 0

        
        for bc in Col2Values['barcode-seq']:
            
            if bc in Bc2ReadCount: 
                rc=Bc2ReadCount[bc]
                Col2Values[colRc].append(rc)
                Col2Values[colRpm].append(1e6*(rc/readsTotal))
                
                readsTotal+=rc
                readsWithDictBc+=rc
                numBcInDict+=1
                
            else:
                Col2Values[colRc].append(np.NaN)
                Col2Values[colRpm].append(np.NaN)
                numBcNotIndict+=1
                
        readsWithoutDictBc = readsTotal-readsWithDictBc
        
        reportMessage=\
        f'''
        Number Barcodes Analyzed       = {numBcTotal:,}
        Number Barcodes in Dictionary  = {numBcInDict:,} ({percent(numBcInDict/numBcTotal)})
        Number Reads Analyzed          = {readsTotal:,}
        Number Reads in Dict Barcodes  = {readsWithDictBc:,} ({percent(readsWithDictBc/readsTotal)})
        '''
        print(reportMessage)
        
    bcdf=pd.DataFrame(Col2Values)
    
    del Col2Values
    
    if outbasename!=None: 
        print(f'\nWriting out bcdf with basename {outbasename}...')
        bcdf.to_csv(   outbasename+'.bcdf.tsv',sep='\t',index=None)
        # bcdf.to_pickle(outbasename+'.bcdf.pd.df.pickle')
    
    print('\nDone!')
    return bcdf        
        
def exp_3_bcdf_to_endf(bcdf,aggregateKeys=['enhancer-seq','enhancer-id'],outbasename=None):
    
    print('Colappsing on enhancers... This may take some time...')
    AggregateFxns = { col : list for col in bcdf.columns if col not in aggregateKeys}
    endf=bcdf.groupby(aggregateKeys).agg(AggregateFxns).reset_index()
    
    if outbasename!=None: 
        print(f'\nWriting out endf with basename {outbasename}...')
        endf.to_csv(   outbasename+'.endf.tsv',sep='\t',index=None)
        # endf.to_pickle(outbasename+'.endf.pd.df.pickle')
    
    print('\nDone!')
    return endf
    

####################################################################################
# Table Processing
####################################################################################


def mask_na_filter(indf,col,filterType,minVal=None,maxVal=None):
    '''Takes a column of lists and converts real integers or floats to NA fi they are outside the 
    max/min thresholds.
    
    Returns a new column.'''
    if   filterType=='max':
        return indf[col].apply(lambda l: [li if li<=maxVal else np.NaN for li in l])
    
    elif filterType=='min':
        return indf[col].apply(lambda l: [li if li>=minVal else np.NaN for li in l])
    
    elif filterType=='between':
        return indf[col].apply(lambda l: [li if (li<=maxVal) and (li>=minVal) else np.NaN for li in l])
    
    else:
        raise ValueError("filterType must be either 'max','min','between'")


def filter_na_simple(indf,col,treatNa):
    '''Takes a column of lists and filters NA (either by omitting or replacing with 0).
    
    Returns a new column.'''
    if treatNa=='omit':
        return indf[col].apply(lambda l: [li for li in l if pd.notnull(li)])
    
    elif treatNa=='zero':
        return indf[col].apply(lambda l: [li if pd.notnull(li) else 0 for li in l ])
    
    else:
        raise ValueError('treatNa must be omit or zero')
        
def filter_na_matched(indf,col1,col2,treatNaType):
    '''Takes 2 column of lists and filters NA (either by ommitting or replacing with 0).
    If NA are treated as "omit" in col1, then indices with NA values in col1 will also be NA in col2.
    If NA are treated as "zero" in col2, then indices with NA values in col1 will be 0 and col2 will stay the same.
    
    Returns 2 columns (idx 0 = col1-filtered, idx1 = col2-filtered).'''
    
    # make sure treatNa are valid
    for treatNa in treatNaType.values():
        if treatNa not in ['omit','zero']:
            raise ValueError('treatNaType values must be "omit" or "zero"')
            
    # make sure cols are included in treatNaType
    if col1 not in treatNaType or col2 not in treatNaType:
        raise ValueError(f'treatNaType does not contain one of the following columns: {col1} or {col2}.')
        
    new1,new2=[],[]
    
    # for each row
    for l1,l2 in zip(*[indf[c] for c in [col1,col2]]):
        
        # for each item in the row list
        newl1,newl2=[],[]
        for l1i,l2i in zip(l1,l2):
            
            l1i_isnull=pd.isnull(l1i)
            l2i_isnull=pd.isnull(l2i)
            
            omit=False
            if l1i_isnull: 
                if   treatNaType[col1]=='omit': omit=True
                elif treatNaType[col1]=='zero': l1i=0
                
            if l2i_isnull: 
                if   treatNaType[col2]=='omit': omit=True
                elif treatNaType[col2]=='zero': l2i=0
                
            if omit:  
                continue
            else:     
                newl1.append(l1i)
                newl2.append(l2i)
            
        new1.append(newl1)
        new2.append(newl2)
    
    return new1,new2 


def concat_sample_lists(indf,colList):
    '''Takes two columns (each with dtype list) and concatenates the two lists.
    
    Returns a new column.'''
    for col in colList:
        if indf.dtypes[col]!=list: 
            raise ValueError('All columns must be columns of lists')
    
    return indf.loc[:,colList].sum(1)

def concat_sample_values(indf,colList):
    '''Takes two columns (each with dtype string, float, int, etc (any single value)).
    
    Returns a new column.'''
    return indf.apply(lambda row: [row[c] for c in colList],axis=1)


def aggregate_within_samples(indf,col,aggfunc):
    '''Takes a column of lists and aggregates items into a single value using a function (eg mean, median, sd, se, etc).
    
    Returns a new column.'''

    return indf[col].apply(aggfunc)

def ratio(x,y):
    return x/y

def aggregate_matched_across_samples(indf,col1,col2,aggfunc):
    '''Takes two columns (each with dtype list) and in an index-matched fashion performs a user defined function (eg ratio, mean, etc).
    For example: [1,2] [3,4] ==> ratio ==> [1/3, 2/4].
    
    No NA values allowed!!!
    
    Returns a new column.'''
            
    return indf.apply(lambda row: [aggfunc(li1,li2) for li1,li2 in zip(row[col1],row[col2])],axis=1)
        

####################################################################################
# Helper functions to set up comparisons for Reproducibility Plots
####################################################################################

def _get_all_comparison_tuples_2d(comparisonList):
    for i in range(len(comparisonList)):
        for j in range(i+1,len(comparisonList)):
            yield i,j

def _all_combinations_recursive(lists, current_combination=None, index=0):
    if current_combination is None:
        current_combination = []

    if index == len(lists):
        yield tuple(current_combination)
        return

    for item in lists[index]:
        current_combination.append(item)
        yield from _all_combinations_recursive(lists, current_combination, index + 1)
        current_combination.pop()
        
def _get_all_value_tuples(stringScaffold,SampleScaffoldKey2ValueList,comparisonScaffoldKey,outColumnLabelPattern=None):
    
    index2scaffoldKey={}
    indexLists=[]
    index=0
    for scaffoldKey,valueList in SampleScaffoldKey2ValueList.items():
        if comparisonScaffoldKey==scaffoldKey: continue
        index2scaffoldKey[index]=scaffoldKey
        indexLists.append(list(range(len(valueList))))
        index+=1
        
    indexCombinations=list(_all_combinations_recursive(indexLists, current_combination=None, index=0))
    
    outputCombos=[]
    for combo in indexCombinations:
        output=stringScaffold                 
        if outColumnLabelPattern: outputCol=outColumnLabelPattern
        for ci,vi in enumerate(combo):
            scaffoldKey=index2scaffoldKey[ci]
            scaffoldValue=SampleScaffoldKey2ValueList[scaffoldKey][vi]
            
            output=output.replace(scaffoldKey,scaffoldValue)
            if outColumnLabelPattern: outputCol=outputCol.replace(scaffoldKey,scaffoldValue)
            
        if not outColumnLabelPattern: outputCombos.append(output)
        else:                         outputCombos.append((output,outputCol))

    return outputCombos

def label_generator_for_comparing_samples(stringScaffold,SampleScaffoldKey2ValueList,comparisonScaffoldKey,dfColumns):
    comparisonList=SampleScaffoldKey2ValueList[comparisonScaffoldKey]
    
    SampleScaffoldKey2ValueList={k:[str(vi) for vi in v] for k,v in SampleScaffoldKey2ValueList.items()}
    comparisonList=[str(i) for i in comparisonList]
    
    comparisonsMade=set()
    columnsNotInDf=set()
    for valueScaffold in _get_all_value_tuples(stringScaffold,SampleScaffoldKey2ValueList,comparisonScaffoldKey):
        
        for cx,cy in _get_all_comparison_tuples_2d(comparisonList):
            
            
            x=valueScaffold.replace(comparisonScaffoldKey,comparisonList[cx])
            y=valueScaffold.replace(comparisonScaffoldKey,comparisonList[cy])
            
            if x not in dfColumns: 
                columnsNotInDf.add(x)
                continue
            if y not in dfColumns: 
                columnsNotInDf.add(y)
                continue
            comparisonsMade.add((x,y))
            
    if len(columnsNotInDf)>0:
        for i in columnsNotInDf:
            print(f'[WARNING] column {i} not in dataframe')
    return comparisonsMade

def label_generator(stringScaffold,SampleScaffoldKey2ValueList,comparisonScaffoldKey,dfColumns):
    comparisonList=SampleScaffoldKey2ValueList[comparisonScaffoldKey]
    
    SampleScaffoldKey2ValueList={k:[str(vi) for vi in v] for k,v in SampleScaffoldKey2ValueList.items()}
    comparisonList=[str(i) for i in comparisonList]
    
    comparisonsMade=set()
    columnsNotInDf=set()
    outputs=[]
    for valueScaffold in _get_all_value_tuples(stringScaffold,SampleScaffoldKey2ValueList,comparisonScaffoldKey):
        for vi in SampleScaffoldKey2ValueList[comparisonScaffoldKey]:
            oi=valueScaffold.replace(comparisonScaffoldKey,vi)
            outputs.append(oi)
        
    outputs=sorted(set(outputs))    
    finalOutputs=[]
    for oi in outputs:
        if oi not in dfColumns:
            print(f'[WARNING] column {oi} not in dataframe')
            continue
        else:
            finalOutputs.append(oi)
            
    return finalOutputs

# def label_generator_paired(stringScaffoldList,SampleScaffoldKey2ValueList,comparisonScaffoldKey,dfColumns):
    
#     comparisonList=SampleScaffoldKey2ValueList[comparisonScaffoldKey]
    
#     SampleScaffoldKey2ValueList={k:[str(vi) for vi in v] for k,v in SampleScaffoldKey2ValueList.items()}
#     comparisonList=[str(i) for i in comparisonList]
    
#     comparisonsMade=set()
#     columnsNotInDf=set()
#     outputs=[]
#     for valueScaffold in _get_all_value_tuples(stringScaffold,SampleScaffoldKey2ValueList,comparisonScaffoldKey):
#         for vi in SampleScaffoldKey2ValueList[comparisonScaffoldKey]:
#             oi=valueScaffold.replace(comparisonScaffoldKey,vi)
#             if oi not in dfColumns:
#                 print(f'[WARNING] column {oi} not in dataframe')
#                 continue
#             outputs.append(oi)
#     return outputs

def label_generator_for_activity_calculations_paired(stringScaffold,SampleScaffoldKey2ValueList,comparisonScaffoldKey,outColumnLabelPattern,dfColumns):
    
    comparisonList=SampleScaffoldKey2ValueList[comparisonScaffoldKey]
    
    SampleScaffoldKey2ValueList={k:[str(vi) for vi in v] for k,v in SampleScaffoldKey2ValueList.items()}
    comparisonList=[str(i) for i in comparisonList]
    
    comparisonsMade=set()
    for valueScaffold,outcol in _get_all_value_tuples(stringScaffold,SampleScaffoldKey2ValueList,comparisonScaffoldKey,outColumnLabelPattern):
        
        for cx,cy in _get_all_comparison_tuples_2d(comparisonList):
            
            if x not in dfColumns: 
                columnsNotInDf.add(x)
                continue
            if y not in dfColumns: 
                columnsNotInDf.add(y)
                continue
                
            x=valueScaffold.replace(comparisonScaffoldKey,comparisonList[cx])
            y=valueScaffold.replace(comparisonScaffoldKey,comparisonList[cy])
            comparisonsMade.add((outcol,x,y))
            
    return comparisonsMade

# def label_generator_for_activity_calculations(stringScaffold,SampleScaffoldKey2ValueList,comparisonScaffoldKey,outColumnLabelPattern):
    
#     comparisonList=SampleScaffoldKey2ValueList[comparisonScaffoldKey]
    
#     SampleScaffoldKey2ValueList={k:[str(vi) for vi in v] for k,v in SampleScaffoldKey2ValueList.items()}
#     comparisonList=[str(i) for i in comparisonList]
    
#     comparisonsMade=set()
#     for valueScaffold,outcol in _get_all_value_tuples(stringScaffold,SampleScaffoldKey2ValueList,comparisonScaffoldKey,outColumnLabelPattern):
        
#         for cx in label_generator(comparisonList):
            
#             x=valueScaffold.replace(comparisonScaffoldKey,comparisonList[cx])
#             comparisonsMade.add((outcol,x))
            
#     return comparisonsMade
    
####################################################################################
# Reproducibility plots
####################################################################################

def venn(indf,colList,plotDims,plotDpi,vennType,ax=None,**kwargs):
    
    if vennType not in ['enhancer-seq','barcode-seq']: raise ValueError('vennType must be "enhancer-seq" or "barcode-seq".')
    
    if len(colList)==2:
        col1,col2=colList
        bcset1=set(indf.loc[indf[col1].notnull(),vennType])
        bcset2=set(indf.loc[indf[col2].notnull(),vennType])
        
        if len(bcset1&bcset2)==len(bcset1): print(f'[WARNING] Complete overlap between {col1} and {col2} detected. Make sure you replace 0 with np.NaN in the beginning stage (e.g. df=df.replace(0,np.NaN))')
        
        if ax==None: fig,ax=plt.subplots(1,figsize=plotDims,dpi=plotDpi)
        
        venn2([bcset1,bcset2],ax=ax,**kwargs)
        
        if ax!=None: return ax
        
    elif len(colList)==3:
        col1,col2,col3=colList
        bcset1=set(indf.loc[indf[col1].notnull(),vennType])
        bcset2=set(indf.loc[indf[col2].notnull(),vennType])
        bcset3=set(indf.loc[indf[col3].notnull(),vennType])
        if len(bcset1&bcset2)==len(bcset1): print(f'[WARNING] Complete overlap between {col1} and {col2} detected. Make sure you replace 0 with np.NaN in the beginning stage (e.g. df=df.replace(0,np.NaN))')
        
        if ax==None: fig,ax=plt.subplots(1,figsize=plotDims,dpi=plotDpi)
        
        venn3([bcset1,bcset2,bcset3],ax=ax,**kwargs)
        
        if ax!=None: return ax

    else:
        raise ValueError('Venn diagrams only supported for 2 or 3 sample comparisons.')
        

def tile_venn(indf,comparisonList,vennType,plotDims,plotDpi=150,subsample=False,returnCharacter=None,**kwargs):
    
    indf=indf.replace(0,np.NaN)
    if subsample: 
        subsample=min(subsample,len(indf))
        if subsample==len(indf): print('subsample parameter provided is less than total plotted data. Plotting all data instead.')
        indf=indf.sample(int(subsample)).copy(deep=True)
    
    comparisonList=sorted(comparisonList)
    
    # determine tile coords for each subplot
    sample2idx={}
    for x,y in comparisonList:
        if x not in sample2idx: sample2idx[x]=len(sample2idx)
        if y not in sample2idx: sample2idx[y]=len(sample2idx)

    # plot  
    fig,ax=plt.subplots(len(sample2idx),len(sample2idx),figsize=plotDims,dpi=plotDpi)
    
    axUsed=set()
    for x,y in comparisonList:
        position=(sample2idx[x],sample2idx[y])
        axi=ax[position]
        axUsed.add(position)
        
        if returnCharacter: 
            set_labels=[x.replace(returnCharacter,' \n '),
                        y.replace(returnCharacter,' \n ')]
        else:
            set_labels=[x,y]
        venn(indf=indf,
             colList=[x,y],
             plotDims=plotDims,
             plotDpi=plotDpi,
             vennType=vennType,
             set_labels=set_labels,
             ax=axi,
             **kwargs)

        axi.set_title('')
        axi.set_ylabel(x)
        axi.set_xlabel(y)

        axi.spines['top'].set_visible(False)
        axi.spines['right'].set_visible(False)
        
    # turn off other axes
    for x in range(len(sample2idx)):
        for y in range(len(sample2idx)):
            if (x,y) not in axUsed:
                ax[x,y].remove()
        
    plt.tight_layout()
    plt.show()    
    

def scatter(indf,col1,col2,treatNa,plotDims,plotDpi,ax=None,**kwargs):
    
    if treatNa not in ['omit','zero']: raise ValueError('treatNa must be either "omit" or "zero"')
    
    x=[]
    y=[]
    
    if treatNa=='omit':
        for enidx in indf.index:
            xi=indf.at[enidx,col1]
            yi=indf.at[enidx,col2]
            
            if pd.isnull(xi) or pd.isnull(yi):
                continue
            else:
                x.append(xi)
                y.append(yi)
                
            
    if treatNa=='zero':
        for enidx in indf.index:
            xi=indf.at[enidx,col1]
            yi=indf.at[enidx,col2]
            
            if pd.isnull(xi):  x.append(0)
            else:              x.append(xi)

            if pd.isnull(yi):  y.append(0)
            else:              y.append(yi)
                
    if ax==None:    fig,ax=plt.subplots(1,figsize=plotDims,dpi=plotDpi)
    
    ax.scatter(x,y,**kwargs)
    ax.set_ylabel(col2)
    ax.set_xlabel(col1)
    
    if ax!=None: return ax

def tile_scatter(indf,comparisonList,plotDims,plotDpi,subsample=False,**kwargs):
    
    if subsample: 
        subsample=min(subsample,len(indf))
        indf=indf.sample(int(subsample)).copy(deep=True)
        if subsample==len(indf): print('subsample parameter provided is less than total plotted data. Plotting all data instead.')

    comparisonList=sorted(comparisonList)
    
    
    # determine tile coords for each subplot
    sample2idx={}
    for x,y in comparisonList:
        if x not in sample2idx: sample2idx[x]=len(sample2idx)
        if y not in sample2idx: sample2idx[y]=len(sample2idx)

    # plot
    fig,ax=plt.subplots(len(sample2idx),len(sample2idx),figsize=plotDims,dpi=plotDpi)
    
    axUsed=set()
    for x,y in comparisonList:
        position=(sample2idx[x],sample2idx[y])
        axi=ax[position]
        axUsed.add(position)
        
        scatter(indf=indf,
                col1=x,
                col2=y,
                treatNa='zero',
                plotDims=(3,3),
                plotDpi=150,
                ax=axi,
                **kwargs)

        axi.set_title('')
        axi.set_ylabel(x)
        axi.set_xlabel(y)

        axi.spines['top'].set_visible(False)
        axi.spines['right'].set_visible(False)
        
    # turn off other axes
    for x in range(len(sample2idx)):
        for y in range(len(sample2idx)):
            if (x,y) not in axUsed:
                ax[x,y].remove()
        
    plt.tight_layout()
    plt.show()    


    
####################################################################################
# Stats tests
####################################################################################

def stats_test_across_samples(indf,col1,col2,statsFunc,parseValueIdx=None,**kwargs):
    '''Takes two columns (each with dtype list) and perform a stats function (eg stats.mannwhitneyu, or stats.ttest_ind) 
    comparing a single sequences enhancer activiy across multiple samples (eg across tissues, replicates, conditions).
    
    Returns a new column.'''
    
    if parseValueIdx==None:
        return indf.apply(lambda row: statsFunc(row[col1],row[col2],**kwargs) if (len(row[col1])>0) and (len(row[col2])>0) else np.NaN,axis=1)
    
    else:
        return indf.apply(lambda row: statsFunc(row[col1],row[col2],**kwargs)[parseValueIdx] if (len(row[col1])>0) and (len(row[col2])>0) else np.NaN,axis=1)
    
def stats_test_accross_sequences(indf,RefSeqIdx2CompSeqIdList,col,statsFunc,parseValueIdx=None,parseValueColName=None,**kwargs):
    '''Takes a column of lists and performs a stats function comparing two sequences enhancer activities within a single sample.
    
    Returns a pandas dataframe.'''
    if parseValueIdx==None: outcols=['ref-seq-id','comp-seq-id','ref-data','comp-data','result']
        
    else:                   outcols=['ref-seq-id','comp-seq-id','ref-data','comp-data','result',parseValueColName]

    c2v={c:[] for c in outcols}
    
    
    for refseqidx, compseqidxList in RefSeqIdx2CompSeqIdList.items():
        
        refdata=indf.at[refseqidx,col]
        
        for compseqidx in compseqidxList:
            
            compseqdata=indf.at[compseqidx,col]
            
            results=statsFunc(refdata,compseqdata,**kwargs)
            
            c2v['ref-seq-id'] .append(refseqidx)
            c2v['comp-seq-id'].append(compseqidx)
            c2v['ref-data']   .append(refdata)
            c2v['comp-data']  .append(compseqdata)
            c2v['result']     .append(results)
            
            if parseValueIdx!=None: c2v[parseValueColName].append(results[parseValueIdx])
            
    df=pd.DataFrame(c2v)
    return df


#!!!!!!!!!! need to replace this with my own code
def _ecdf(x):
    """No frills empirical cdf used in fdrcorrection."""
    nobs = len(x)
    return np.arange(1, nobs + 1) / float(nobs)

def bh_correction(pvals, alpha=0.05, method='indep'):
    """P-value correction with False Discovery Rate (FDR).
    Correction for multiple comparison using FDR :footcite:`GenoveseEtAl2002`.
    This covers Benjamini/Hochberg for independent or positively correlated and
    Benjamini/Yekutieli for general or negatively correlated tests.
    Parameters
    ----------
    pvals : array_like
        Set of p-values of the individual tests.
    alpha : float
        Error rate.
    method : 'indep' | 'negcorr'
        If 'indep' it implements Benjamini/Hochberg for independent or if
        'negcorr' it corresponds to Benjamini/Yekutieli.
    Returns
    -------
    reject : array, bool
        True if a hypothesis is rejected, False if not.
    pval_corrected : array
        P-values adjusted for multiple hypothesis testing to limit FDR.
    References
    ----------
    .. footbibliography::
    """
    pvals = np.asarray(pvals)
    shape_init = pvals.shape
    pvals = pvals.ravel()

    pvals_sortind = np.argsort(pvals)
    pvals_sorted = pvals[pvals_sortind]
    sortrevind = pvals_sortind.argsort()

    if method in ['i', 'indep', 'p', 'poscorr']:
        ecdffactor = _ecdf(pvals_sorted)
    elif method in ['n', 'negcorr']:
        cm = np.sum(1. / np.arange(1, len(pvals_sorted) + 1))
        ecdffactor = _ecdf(pvals_sorted) / cm
    else:
        raise ValueError("Method should be 'indep' and 'negcorr'")

    reject = pvals_sorted < (ecdffactor * alpha)
    if reject.any():
        rejectmax = max(np.nonzero(reject)[0])
    else:
        rejectmax = 0
    reject[:rejectmax] = True

    pvals_corrected_raw = pvals_sorted / ecdffactor
    pvals_corrected = np.minimum.accumulate(pvals_corrected_raw[::-1])[::-1]
    pvals_corrected[pvals_corrected > 1.0] = 1.0
    pvals_corrected = pvals_corrected[sortrevind].reshape(shape_init)
    reject = reject[sortrevind].reshape(shape_init)
    return [min(1,p) for  p in pvals_corrected]

def bonf_correction(pvalue_list):
    '''a method to take in a list of pvalues and return adjusted pvalues'''
    n=len([pi for pi in pvalue_list if pd.notnull(pi)])
    return [min(1,p*n) for p in pvalue_list]

####################################################################################
# Enhancer Visualizations
####################################################################################

def visualize_different_enhancers_same_sample(indf,plotPointsType,plotDistType,plotDistWidth,enList,actCol,plotDims,plotDpi,plotTitle='',**kwargs):
        
    # get data
    data=[]
    for en in enList:
        dataAdded=indf.at[en,actCol]
        dataAdded=[di for di in dataAdded if pd.notnull(di)]
        data.append(dataAdded)
    
    # initialize figure
    fig,ax=plt.subplots(1,figsize=plotDims,dpi=plotDpi)
    
    # plot points
    if   plotPointsType=='swarm':
        sns.swarmplot(data=data,ax=ax,**kwargs)
        
    elif plotPointsType=='stripplot':
        sns.stripplot(data=data,ax=ax,**kwargs)
        
    elif plotPointsType=='': pass
        
    else: 
        raise ValueError('plotPointsType must be either \'swarm\' or \'stripplot\'')
        
    # initialize boxplots
    if    plotDistType=='box':
        sns.boxplot(data=data,color='white',width=plotDistWidth)
    
    elif  plotDistType=='violin':
        sns.violinplot(data=data,inner=None,color='white',width=plotDistWidth,cut=0)
        
    elif plotDistType=='': pass

    else: 
        raise ValueError('plotPointsType must be either \'violin\' or \'box\'')
        
    ax.set_title(plotTitle)
    ax.set_xticklabels(enList)
        
    return ax

def visualize_same_enhancer_different_samples(indf,plotPointsType,plotDistType,plotDistWidth,en,colList,plotDims,plotDpi,plotTitle='',**kwargs):
        
    # get data
    data=[]
    for actCol in colList:
        dataAdded=indf.at[en,actCol]
        dataAdded=[di for di in dataAdded if pd.notnull(di)]
        data.append(dataAdded)
    
    # initialize figure
    fig,ax=plt.subplots(1,figsize=plotDims,dpi=plotDpi)
    
    # plot points
    if   plotPointsType=='swarm':
        sns.swarmplot(data=data,ax=ax,**kwargs)
        
    elif plotPointsType=='stripplot':
        sns.stripplot(data=data,ax=ax,**kwargs)
        
    elif plotPointsType=='': pass

    else: 
        raise ValueError('plotPointsType must be either \'swarm\' or \'stripplot\'')
        
    # initialize boxplots
    if    plotDistType=='box':
        sns.boxplot(data=data,color='white',width=plotDistWidth)
    
    elif  plotDistType=='violin':
        sns.violinplot(data=data,inner=None,color='white',width=plotDistWidth,cut=0)
        
    elif plotDistType=='': pass


    else: 
        raise ValueError('plotPointsType must be either \'violin\' or \'box\'')
        
    ax.set_title(plotTitle)
    ax.set_xticklabels(colList)
        
    return ax


########################################################################################
# Graveyard
########################################################################################

# def reproducibility_matched_scatter(indf,col1,col2,treatNa,plotDims,plotDpi,**kwargs):
    
#     if treatNa not in ['omit','zero']: raise ValueError('treatNa must be either "omit" or "zero"')
    
#     x=[]
#     y=[]
    
#     if treatNa=='omit':
#         for enidx in indf.index:
#             for xi,yi in zip(indf.at[enidx,col1],indf.at[enidx,col2]):
#                 if pd.isnull(xi) or pd.isnull(yi): 
#                     continue
#                 else:
#                     x.append(xi)
#                     y.append(yi)
                    
#     if treatNa=='zero':
#         for enidx in indf.index:
#             for xi,yi in zip(indf.at[enidx,col1],indf.at[enidx,col2]):
                
#                 if pd.isnull(xi):  x.append(0)
#                 else:              x.append(xi)
                    
#                 if pd.isnull(yi):  y.append(0)
#                 else:              y.append(yi)
                
#     fig,ax=plt.subplots(1,figsize=plotDims,dpi=plotDpi)
    
#     ax.scatter(x,y,**kwargs)
#     ax.set_ylabel(col2)
#     ax.set_xlabel(col1)
    
#     return ax
