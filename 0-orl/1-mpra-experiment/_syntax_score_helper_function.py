#!/usr/bin/env python
# coding: utf-8

# In[2]:


import jsDna as jsd
import jsAff as jsa
import pandas as pd
import js

import numpy as np

# # Define methods

# In[3]:


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


# In[4]:


def get_phrase_str(phrase_key_tuple):
    return '_'.join([str(i) for i in phrase_key_tuple])


# In[5]:


def can_be_converted_to_int(val):
    try: 
        int(val)
        return True
    except ValueError:
        return False

# can_be_converted_to_int('1')


# In[6]:


ori2flip={
    'g':'G',
    'G':'g',
    'e':'E',
    'E':'e'
}

def revcomp(phrase,phraseType):
    
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

# revcomp('g_E','ord_ori')


# ## Get tfbs string from sequence

# In[7]:



def get_tfbs_string(seq,TfbsName2SeqSet,TfbsSeq2NameOri,TfbsName2SpacingStart,TfbsName28mer,TfbsName2Seq2Aff):
    
    ###############################################################
    # Create reference TFBS objects
    ###############################################################
    
    TfbsName2SeqSetRc={tfbs:set() for tfbs in TfbsName2SeqSet}
    TfbsSeq2Name     ={}
    for tfbs,seqSet in TfbsName2SeqSet.items():
        for tfbsSeq in seqSet:
            TfbsName2SeqSetRc[tfbs].add(tfbsSeq)
            TfbsName2SeqSetRc[tfbs].add(jsd.revcomp(tfbsSeq))
            TfbsSeq2Name[tfbsSeq]=tfbs
            TfbsSeq2Name[jsd.revcomp(tfbsSeq)]=tfbs
            
    # Check all seqs within a single tf are same size
    TfbsName2Size={}
    for tfName,tfSeqSet in TfbsName2SeqSetRc.items():
        if len(set([len(seq) for seq in tfSeqSet]))>1:
            raise ValueError("All tf seqs must be same size within a family")
        TfbsName2Size[tfName]=len(list(tfSeqSet)[0])
        
            
    ###############################################################        
    # Get location of all tfbs starts
    ###############################################################
    
    BpIndex2TfbsName={}
    for bpi in range(2,len(seq)-5):
        
        # Check for each sequence
        for tfName,tfSeqSet in TfbsName2SeqSetRc.items():
            kmer=seq[bpi:bpi+TfbsName2Size[tfName]]
            if kmer in tfSeqSet:
                BpIndex2TfbsName[bpi]=TfbsSeq2NameOri[kmer]
                
    ###############################################################        
    # Infer spacing
    ###############################################################
                
    tfbsString = []
    iterList   = sorted([(bpi,tfName) for bpi,tfName in BpIndex2TfbsName.items()])
    
    for i,(bpi,tfName) in enumerate(iterList):
        
        thisSpaceMarker=bpi+TfbsName2SpacingStart[tfName]
        
        if i>0:
            lastSpacing=thisSpaceMarker-lastSpaceMarker
            tfbsString.append(lastSpacing)
            
        startAdj,endAdj=TfbsName28mer[tfName]
        kmer8=seq[bpi+startAdj:bpi+endAdj]
        
        aff=round(TfbsName2Seq2Aff[tfName.upper()][kmer8],2)
        aff=str(aff)[1:]

        tfbsString.append(tfName+aff)
        
        lastSpaceMarker=thisSpaceMarker
        
    return tfbsString
    
seq='AGATAGATCTGAAGCTCGTTATCTCTAACGGAAGTTTTCGAAAAGGAAATTGTTCAATATCTAAGATAGG'

gata=set(['GATA'])
ets=set(['GGAA','GGAT'])
TfbsName2SeqSet={'G':gata,'E':ets}

TfbsSeq2NameOri={'GATA':'G',
                 'TATC':'g',
                 'GGAA':'E','GGAT':'E',
                 'TTCC':'e','ATCC':'e'}

TfbsName2SpacingStart={'G':2,
                       'g':2,
                       'E':2,
                       'e':2}

TfbsName28mer={'G':(-2,6),
               'g':(-2,6),
               'E':(-2,6),
               'e':(-2,6)}

TfbsName2Seq2Aff={'G':jsa.loadAff('supp/parsed_Gata6_3769_contig8mers.txt'),
                  'E':jsa.loadAff('supp/parsed_Ets1_8mers.txt')}

# get_tfbs_string(seq,TfbsName2SeqSet,TfbsSeq2NameOri,TfbsName2SpacingStart,TfbsName28mer,TfbsName2Seq2Aff)


# ## Get grammar features from tfbs lit

# In[8]:


# seq='GATCTGAAGCTCGTTATCTCTAACGGAAGTTTTCGAAAAGGAAATTGTTCAATATCTAAGATAGGA'
# test=get_tfbs_string(seq,TfbsName2SeqSet,TfbsSeq2NameOri,TfbsName2SpacingStart,TfbsName28mer,TfbsName2Seq2Aff)
# test


# ### Ord 

# In[9]:


def parse_ord(sitelist):

    sentence=[]
    for site in sitelist:
        if can_be_converted_to_int(site): continue # Skip spacing
        sentence.append(site.split('.')[0].upper())
    return sentence

# parse_ord(test)


# ### Ord Ori

# In[10]:


# seq='GATCTGAAGCTCGTTATCTCTAACGGAAGTTTTCGAAAAGGAAATTGTTCAATATCTAAGATAGGA'
# test=get_tfbs_string(seq,TfbsName2SeqSet,TfbsSeq2NameOri,TfbsName2SpacingStart,TfbsName28mer,TfbsName2Seq2Aff)
# test


# In[11]:



def parse_ord_ori(sitelist):

    sentence=[]
    for site in sitelist:
        if can_be_converted_to_int(site): continue # Skip spacing
        sentence.append(site.split('.')[0])
    return sentence

# parse_ord_ori(test)


# ### Ord Spc

# In[12]:


def parse_ord_spc(sitelist):

    sentence=[]
    for sitei,site in enumerate(sitelist):
        
        # if spacer
        if can_be_converted_to_int(site):
            sentence.append(site)
            
        # if tfbs
        else:
            sentence.append(site.split('.')[0])
            
    return sentence

# parse_ord_spc(test)


# ### Ord Aff

# In[13]:



def parse_ord_aff(sitelist):

    sentence=[]
    for sitei,site in enumerate(sitelist):
        
        # skip spacing
        if can_be_converted_to_int(site):
            continue
            
        # if tfbs
        else:
            site=site[0].upper()+site[1:]
            sentence.append(site)
            
    return sentence

# parse_ord_aff(test)


# ### Ord Ori Aff

# In[14]:



def parse_ord_ori_aff(sitelist):
    
    sentence=[]
    for sitei,site in enumerate(sitelist):
        
        # skip spacing
        if can_be_converted_to_int(site):
            continue
            
        # if tfbs
        else:
            sentence.append(site)
            
    return sentence

# parse_ord_ori_aff(test)


# ### Ord Spc Aff

# In[15]:



def parse_ord_spc_aff(sitelist):

    sentence=[]
    for sitei,site in enumerate(sitelist):
        
        # skip spacing
        if can_be_converted_to_int(site):
            sentence.append(site)
            
        # if tfbs
        else:
            site=site[0].upper()+site[1:]
            sentence.append(site)
            
    return sentence


# parse_ord_spc_aff(test)


# ### Ord Ori Spc

# In[16]:



def parse_ord_ori_spc(sitelist):
    sentence=[]
    for sitei,site in enumerate(sitelist):
        
        # skip spacing
        if can_be_converted_to_int(site):
            sentence.append(site)
            
        # if tfbs
        else:
            sentence.append(site.split('.')[0])
            
    return sentence


# parse_ord_ori_spc(test)


# ### Ord Ori Spc Aff

# In[17]:



def parse_ord_ori_spc_aff(sitelist):

    sentence=[]
    for sitei,site in enumerate(sitelist):
        
        # skip spacing
        if can_be_converted_to_int(site):
            sentence.append(site)
            
        # if tfbs
        else:
            sentence.append(site)
            
    return sentence


# parse_ord_ori_spc_aff(test)


# ## Wrapper function to get grammar

# In[18]:


parse_ord_ori_spc_aff
Features2Func={
    # tuple(['ord']):parse_ord,
    # ('ord','ori'):parse_ord_ori,
    # ('ord','spc'):parse_ord_spc,
    # ('ord','aff'):parse_ord_aff,
    ('ord','ori','spc'):parse_ord_ori_spc,
    # ('ord','ori','aff'):parse_ord_ori_aff,
    ('ord','ori','spc','aff'):parse_ord_ori_spc_aff
}
Features2Func={tuple(sorted(featureTuple)):func for featureTuple,func in Features2Func.items()}
# Features2Func


# In[19]:


def get_sentence_from_seq(seq,featureTuple=('ord','ori','spc')):
    
    siteList = get_tfbs_string(seq,TfbsName2SeqSet,TfbsSeq2NameOri,TfbsName2SpacingStart,TfbsName28mer,TfbsName2Seq2Aff)
    
    featureTuple = tuple(sorted(featureTuple))
    func=Features2Func[featureTuple]
    
    sentence = func(siteList)
    
    return sentence

# seq='GATCTGAAGCTCGTTATCTCTAACGGAAGTTTTCGAAAAGGAAATTGTTCAATATCTAAGATAGGA'
# get_sentence_from_seq(seq,('ord','spc','ori'))


# ## Wrapper function to get all motifs

# In[20]:


import itertools

def get_all_phrases(sentence,spacing,length):
    
    # If your phrases don't include spacing
    if spacing==False:
        phrase2count={}
        for phrase in itertools.combinations(sentence, length):
            if phrase not in phrase2count:
                phrase2count[phrase]=0
            phrase2count[phrase]+=1
        return phrase2count
    
    # If your phrases include spacing
    if spacing==True:
        phrase2count={}
        sentence_withoutSpacing=[site for site in sentence if type(site)==str]
        for phrase_ilist in itertools.combinations(range(len(sentence_withoutSpacing)), length):
            
            # Get the sentence with spacing
            phrase_ilist=[2*pi for pi in phrase_ilist]
            phrase_withSpacing=[]
            for i,pi in enumerate(phrase_ilist):
                
                # if you've reached the last binding site
                if i==len(phrase_ilist)-1:
                    phrase_withSpacing.append(sentence[pi])
                
                # if there are still other spacings to calculate
                else:
                    
                    # add the first site
                    phrase_withSpacing.append(sentence[pi])
                    
                    # determine the spacing
                    pi_next=phrase_ilist[i+1]
                    
                    # add spacing
                    spacingBetweenFirstLast=sum([s for s in sentence[pi:pi_next] if type(s)==int])
                    phrase_withSpacing.append(spacingBetweenFirstLast)
                
            phrase_withSpacing=tuple(phrase_withSpacing)
            if phrase_withSpacing not in phrase2count: 
                phrase2count[phrase_withSpacing]=0
            phrase2count[phrase_withSpacing]+=1
                
        return phrase2count
        


# In[21]:


# seq='GATCTGAAGCTCGTTATCTCTAACGGAAGTTTTCGAAAAGGAAATTGTTCAATATCTAAGATAGGA'
# sentence=get_sentence_from_seq(seq,tuple(['ord']))
# print(sentence)
# print()
# phrase2count=get_all_phrases(sentence,spacing=False,length=2)
# phrase2count


# # Function to ascribe functional score to enumerated motifs

# ## Load motifs

# In[32]:


# fn='../../20220729-v2-closer-to-publishing/0b-motifs/20220809-grammar-motifs-v4-hypergeom-50k/5-all-phrases-possible-hypergeom-50k-withRealData-stats-sigMotifs-20220810.pandas.df.pickle'
# fmDF=pd.read_pickle(fn)
# fmDF['N']=fmDF['real-data'].apply(lambda t: sum(t))
# fmDF.head(1)


# In[33]:


# sigMotifSet=set(fmDF.index)
# list(sigMotifSet)[:3]


# ## Function

# In[36]:



def get_all_grammar_motifs(seq,sigMotifSet,fmDF,f_effectSizeCol,nf_effectSizeCol='',skip_aff=True,activatorOnly=False):
    
    m_c2v ={'seq':[],'motif':[],'type':[],'fc':[]}#,'n-enhancers':[]}
    
    ########################
    # single motifs
    ########################
    funcPhraseAdded=False
    
    for featureTuple in Features2Func:
        
        # Skip aff if intended
        if skip_aff:
            if 'aff' in featureTuple:
                continue
                
        # Get tfbs list
        sentence=get_sentence_from_seq(seq,featureTuple)
        
        # Format some things
        featureTuple_str='_'.join(featureTuple)
        spacing='spc' in featureTuple
        
        # For all phrase legnths
        sigPhrasesAdded=set()
        
        for phraseLength in [2,3]:
            
            for phrase,count in get_all_phrases(sentence,spacing=spacing,length=phraseLength).items():
                
                for phrase_with_ori in [phrase,revcomp(phrase,featureTuple_str)]:
                    
                    # print(phrase_with_ori)
                    # print(phrase)
                    # print(featureTuple_str)
                    
                    # this accounts for palindromes that are both significant phrases
                    # when phrases are added to this set, we add both fwd and rev so you dont have to try both ways
                    if phrase_with_ori in sigPhrasesAdded: 
                        continue
                    
                    phrase_with_ori=featureTuple_str+'='+get_phrase_str(phrase_with_ori)
                    
                    # only add if functional
                    
                    if phrase_with_ori in sigMotifSet:
                        
                        # add both directions to sigPhrasesAdded
                        sigPhrasesAdded.add(phrase)
                        sigPhrasesAdded.add(revcomp(phrase,featureTuple_str))
                        
                        
                        # if doing func and nonfunc motifs
                        if activatorOnly==False: 
                            
                            f_fc=fmDF.at[phrase_with_ori,f_effectSizeCol]
                            n_fc=fmDF.at[phrase_with_ori,nf_effectSizeCol]
                        
                            if f_fc>1:
                                mType='F'
                                fc=f_fc
                            else:
                                mType='NF'
                                fc=n_fc
                                
                        # if doing func motifs only
                        else:
                            f_fc=fmDF.at[phrase_with_ori,f_effectSizeCol]
                            fc=f_fc
                            mType='F'
                        
                        funcPhraseAdded=True
                        
                        m_c2v['seq'].append(seq)
                        m_c2v['motif'].append(phrase_with_ori)
                        m_c2v['type'].append(mType)
                        m_c2v['fc'].append(fc)
                        # m_c2v['n-enhancers'].append(fmDF.at[phrase_with_ori,'N'])
                        
    if funcPhraseAdded==False:
        if activatorOnly==True: mtype='F'
        else:                   mtype='.'
        m_c2v['seq'].append(seq)
        m_c2v['motif'].append(np.NaN)
        m_c2v['type'].append(mtype)
        m_c2v['fc'].append(np.NaN)
        # m_c2v['n-enhancers'].append(fmDF.at[phrase_with_ori,'N'])
                        
    out_mDF=pd.DataFrame(m_c2v)
    return out_mDF
    
# seq='GATCTGAAGCTCGTTATCTCTAACGGAAGTTTTCGAAAAGGAAATTGTTCAATATCTAAGATAGGA'


# mDF=get_all_grammar_motifs(seq,sigMotifSet,fmDF,'func-enrich','nonfunc-enrich',skip_aff=True)
# mDF


# In[39]:


def get_all_grammar_motifs_BATCH(seqList,sigMotifSet,fmDF,f_effectSizeCol,nf_effectSizeCol,skip_aff=True,activatorOnly=False):
    mDFList=[]
    for seq in seqList:
        mDF=get_all_grammar_motifs(seq,sigMotifSet,fmDF,f_effectSizeCol,nf_effectSizeCol,skip_aff=True,activatorOnly=activatorOnly)
        mDFList.append(mDF)
        
    out_mDF=pd.concat(mDFList)
    
    return out_mDF
     
# seqList=['GATCTGAAGCTCGTTATCTCTAACGGAAGTTTTCGAAAAGGAAATTGTTCAATATCTAAGATAGGA',
#          'ATTTCCTTTGCTCAAGATAGGTAACTTCCGTCATCTGAAGCTCGAGATAACTTTCGAAAATATCTA']
# mDF=get_all_grammar_motifs_BATCH(seqList,sigMotifSet,fmDF,'func-enrich','nonfunc-enrich',skip_aff=True)
# mDF


# # Load motifs

# In[21]:


# fn='20220607-grammar-motifs-v3/5-all-phrases-possible-withRealData-with2PermuteIters-stats-sigMotifs-20220610.pandas.df.pickle'
# fmDF=pd.read_pickle(fn)
# # fmDF.head(1)


# In[22]:


# fmDF['motif-key']=[ti+'='+pi for ti,pi in zip(fmDF['phrase-type'],fmDF['phrase-str'])]
# fmDF=fmDF.set_index('motif-key',drop=True)
# # fmDF.head(1)


# In[23]:


# fn='20220612-grammar-comotifs/sigDF-comotifs.tsv'
# cmDF=pd.read_csv(fn,sep='\t')
# # cmDF.head(1)


# In[24]:


# cmDF['comotif-key']=[(m10,m01) for m10,m01 in zip(cmDF['m10'],cmDF['m01'])]
# cmDF=cmDF.set_index('comotif-key',drop=True)
# # cmDF.head(1)


# In[25]:


# sigMotifSet  =set(fmDF.index)
# # js.sprint(sigMotifSet,3)


# In[26]:


# sigComotifSet=set(cmDF.index)
# # js.sprint(sigComotifSet,3)


# # Annotate Enhancers

# In[27]:


# fn='20220616_Annotate-Enhancer-With-Motifs/all-seqs-to-annotate.tsv'

# skip_aff=True
# line_out=js.write_row(['en-name','en-dataset','en-seq','en-func','motif-type','m10','m10-count','m10-func-fc','m10-nonfunc-fc','m01','m01-count','m01-func-fc','m01-nonfunc-fc'])
# for name,typ,seq,func in js.read_tsv(fn,pc=True,header=True):
    
    ######################################################
    # Get all possible phrases
    ######################################################
    
#     phrase2count={}
#     for featureTuple in Features2Func:
        
#         # Skip aff if intended
#         if skip_aff:
#             if 'aff' in featureTuple:
#                 continue
                
#         # Get tfbs list
#         sentence=get_sentence_from_seq(seq,featureTuple)
        
#         # Format some things
#         featureTuple_str='_'.join(featureTuple)
#         spacing='spc' in featureTuple
        
#         # For all phrase legnths
#         for phraseLength in [2,3,4]:
#             for phrase,count in get_all_phrases(sentence,spacing=spacing,length=phraseLength).items():
                
#                 for phrase_with_ori in [phrase,revcomp(phrase,featureTuple_str)]:
                    
#                     phrase_with_ori=featureTuple_str+'='+get_phrase_str(phrase_with_ori)
                    
#                     # only add if functional
#                     if phrase_with_ori in sigMotifSet:
#                         phrase2count[phrase_with_ori]=count
         
#     ######################################################
#     # Write out single motifs
#     ######################################################
    
#     # Single motifs
#     # phrase2count={}
#     for p,c in phrase2count.items():
#         enrichFunc=fmDF.at[p,'func-enrichment']
#         enrichNonFunc=fmDF.at[p,'nonfunc-enrichment']
#         line_out+=js.write_row( [name,
#                                  typ,
#                                  seq,
#                                  func,
#                                  'single',
#                                  p,c,enrichFunc,enrichNonFunc,
#                                  '','','',''] )
    
#     ######################################################
#     # Write out comotifs
#     ######################################################
    
#     # Comotifs
#     for m10,m01 in sigComotifSet:
#         if m10 in phrase2count and m01 in phrase2count:
            
#             enrich10=[float(i) for i in cmDF.at[(m10,m01),'enrich-11-10'][1:-1].replace(' ','').split(',')]
#             enrich01=[float(i) for i in cmDF.at[(m10,m01),'enrich-11-01'][1:-1].replace(' ','').split(',')]
            
#             funcFC10,_,nonfuncFC10=enrich10
#             funcFC01,_,nonfuncFC01=enrich01
            
#             line_out+=js.write_row( [name,
#                                      typ,
#                                      seq,
#                                      func,
#                                      'comotifs',
#                                      m10,phrase2count[m10],funcFC10,nonfuncFC10,
#                                      m01,phrase2count[m01],funcFC01,nonfuncFC01] )
            
            
            
#     for p,c in phrase2count.items():
#         print(p,c)
#     break
    


# In[28]:


# outfn='20220616_Annotate-Enhancer-With-Motifs/enhancer-motif-annotations-long-table.tsv'
# with open(outfn,'w') as f:
#     f.write(line_out)


# In[ ]:




