"""
accessory tools
"""
import pysam, gzip, os, re
import pandas as pd
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.pylab as pl
import itertools
import seaborn as sns
import cbs
from statsmodels.nonparametric.smoothers_lowess import lowess
from collections import OrderedDict
from general import *
from sklearn import mixture
from scipy import linalg, stats
from numpy import array
from decimal import *

package_directory = os.path.dirname(os.path.abspath(__file__))

def PopFreq_57patient(bg,gene,frequency_table="pop_freq.txt"):
    '''
        rerank the allele combination based on the frequency setup in this fucntion
        the default table is based on the patient in the 56 validation cohort
        new table is based on Stella's paper 
        ### put all the allele frequency here, transform to log-likelihood , add to the calculated freq
        ### remove the allele with frquency of 1
        ### For BG, only if nrow of the bg > 2 or unique index number > 2
        ### if so, start to add the pop freq to each pair bassed on index
        ### keep only the index with the highest freq in the end
    '''
    if 'gz' in frequency_table:
        freq = gzip.open(os.path.join(package_directory, 'database', frequency_table),'rt')
    else:
        freq = open(os.path.join(package_directory, 'database', frequency_table),'rt')
    allele_count={}
    paired_allele_count={}
    total_count={} 
    for line in freq:
        if line.startswith('#'): continue
        d=line.strip().split('\t')
        if d[0]=='Gene': continue ### skip header
        f=int(d[4])
        if d[0] not in allele_count: allele_count[d[0]]={}
        if d[0] not in paired_allele_count: paired_allele_count[d[0]]={}
        if d[0] not in total_count: total_count[d[0]]=0
        if d[1] not in allele_count[d[0]]:
            allele_count[d[0]][d[1]]=f
        else:
            allele_count[d[0]][d[1]]+=f
        if d[1] not in paired_allele_count[d[0]]:
            paired_allele_count[d[0]][d[1]]={}
            paired_allele_count[d[0]][d[1]][d[2]]=f
        else:
            if d[2] not in paired_allele_count[d[0]][d[1]]:
                paired_allele_count[d[0]][d[1]][d[2]]=f
            else:
                paired_allele_count[d[0]][d[1]][d[2]]+=f
        total_count[d[0]]+=f*2
    
    allele_freq={}
    for gene in allele_count:
        if gene not in allele_freq: 
            allele_freq[gene]={}
        for allele in allele_count[gene]:
            allele_freq[gene][allele]={} 
            allele_freq[gene][allele]['count']=allele_count[gene][allele]
            allele_freq[gene][allele]['prop']=allele_count[gene][allele]/float(total_count[gene])
    if bg.shape[0] > 2 or bg['Index'].nunique() > 2:
        for ix in bg['Index'].unique():
            a1=bg.loc[(bg.Index==float(ix)) & (bg.Allele=='A1'),'Bloodtype'].values[0]
            a2=bg.loc[(bg.Index==float(ix)) & (bg.Allele=='A2'),'Bloodtype'].values[0]
            ll=bg.loc[(bg.Index==float(ix)) & (bg.Allele=='A1'),'Likelihood'].values[0]
            if a1 in allele_freq[gene] and a2 in allele_freq[gene]:
                p1=allele_freq[gene][a1]['prop'] if a1 in allele_freq[gene] else 0
                p2=allele_freq[gene][a2]['prop'] if a2 in allele_freq[gene] else 0
                if allele_freq[gene][a1]['count'] > 1 or allele_freq[gene][a2]['count'] > 1:
                    print ("...adjust for population freq")
                    print ("...p1:", p1, " ...p2:", p2, " ...ll:", ll)
                    ll-=np.log(p1)+np.log(p2)
                    print ("...adj.ll:", ll)
            bg.loc[(bg.Index==float(ix)) & (bg.Allele=='A1'),'Likelihood']=ll
            bg.loc[(bg.Index==float(ix)) & (bg.Allele=='A2'),'Likelihood']=ll
           
    return bg          
         

def PopFreq(final,gene,frequency_table="pop_freq.txt"):
   '''
      rerank the allele combination based on the frequency setup in this fucntion
      the default table (see PopFreq_57patient) is based on the patient in the 56 validation cohort
      new table is based on Stella's paper 
   '''
   INDEXls=final.Index.unique().tolist()
   if len(INDEXls) < 1: return final

   fn=os.path.join(package_directory, 'database', frequency_table)
   if 'gz' in frequency_table:
        freq = pd.read_csv(fn, sep='\t', compression='gzip',comment='#',skiprows=1)
   else:
        freq = pd.read_csv(fn, sep='\t', comment='#',skiprows=1)

   pop_freq={}
   INDEXls=final.Index.unique().tolist()
   max_pop_freq=0
   for IDX in INDEXls:
       Alleles=final.loc[final['Index']==IDX,'Bloodtype'].values.tolist()
       status='NA'
       pop_freq_list=[]
       for a in Alleles:
           if a in freq['Allele name'].tolist():
               pop_freq_list.append(freq.loc[ (freq['Allele name']==a), 'PopFreq'].values[0])
           else: pop_freq_list.append(0)
       if IDX not in pop_freq:
           pop_freq[IDX]=pop_freq_list[0] * pop_freq_list[1]
       if pop_freq_list[0] * pop_freq_list[1] > max_pop_freq:
           max_pop_freq=pop_freq_list[0] * pop_freq_list[1]

   index2drop=[]
   for i in pop_freq:
       if pop_freq[i] < max_pop_freq: index2drop.append(i)
   if len(index2drop) > 0: final=final[~final.Index.isin(index2drop)]   
 
   ### keep the pair with one allele with lowest number of mismatch as the candidate
   INDEXls=final.Index.unique().tolist()
   min_variant=final.matchDBn.min()
   index2drop=[]
   if len(INDEXls) > 1:
       for IDX in INDEXls:
          min_var=final.loc[final['Index']==IDX,'matchDBn'].min() 
          if min_var != min_variant:
              index2drop.append(IDX)
   if len(index2drop) > 0: final=final[~final.Index.isin(index2drop)]      

    
   ### keep ce48C and ce733G together, need to use Connie's linkage rules
   #INDEXls=final.Index.unique().tolist()
   #if len(INDEXls)== 1:
   #    idx=INDEXls[0]
   #    alias=final.loc[final['Index']==idx,'Alias'].tolist()
   #    if 'RHCE*ce48C' in alias and 'RHCE*ce733G' in alias:
   #        final.loc[final['Allele']=='A1','Bloodtype']='RHCE*01 or RHCE*ce RHCE*c RHCE*e'
   #        final.loc[final['Allele']=='A1','Alias']='RHCE*ce'
   #        final.loc[final['Allele']=='A1','matchDBn']=0
   #        final.loc[final['Allele']=='A1','matchDB']=''
   #        final.loc[final['Allele']=='A1','matchDBaa']=''
   #        final.loc[final['Allele']=='A2','Bloodtype']='RHCE*01.20.02 RHCE*ce.20.02 RHCE*ceVS.02'
   #        final.loc[final['Allele']=='A2','Alias']='RHCE*ce48C, 733G'  
   #        final.loc[final['Allele']=='A2','matchDBn']=2         
   #        final.loc[final['Allele']=='A2','matchDB']=final.loc[final['Allele']=='A1','AllVariations'].values[0]
   #        final.loc[final['Allele']=='A2','matchDBaa']='48G/C(W/C);733C/G(L/V)'
   return final

    

def filter_annotation(full, filtered, annotation="filtered", LH_adj=-0.00001):
    '''
        Annotate the full ranking with filtered reason in the comment column
    '''
    kept_idx=filtered.Index.unique()
    curr_idx=full.Index.unique()
    #print(annotation)
    #print(kept_idx)
    #print(curr_idx)
    warning=False
    for idx, row in full.iterrows():
        if row['Index'] not in kept_idx and row['Comment']=="":
            full.at[idx, 'Comment']=annotation
            full.at[idx, 'Likelihood']=full.at[idx, 'Likelihood'] + Decimal(LH_adj)
            warning=True
        if row['Index'] in kept_idx:
            if row['Allele']=='A1':
                A1_n=row['Bloodtype']
                A1_fn=filtered.loc[(filtered.Index==row['Index']) & (filtered.Allele=='A1'), 'Bloodtype'].values[0]
                if A1_n != A1_fn:
                    full.loc[(full.Index==row['Index']) & (full.Allele=='A1'), 'Bloodtype']=A1_fn
            if row['Allele']=='A2':
                A2_n=row['Bloodtype']
                A2_fn=filtered.loc[(filtered.Index==row['Index']) & (filtered.Allele=='A2'), 'Bloodtype'].values[0]
                if A2_n != A2_fn:
                    full.loc[(full.Index==row['Index']) & (full.Allele=='A2'), 'Bloodtype']=A2_fn

    #print(full)
    return full, warning



def RHD_RHCE_linked(bg1, bg2fn):
    """
      correct the alleles based on linked information
      INFO based on:
      http://www.isbtweb.org/fileadmin/user_upload/files-2015/red%20cells/blood%20group%20allele%20terminology/allele%20tables/004%20RHCE%20Alleles%20v2.0%20110914.pdf
    """
    links={'RHCE@RHCE_ceTI':'RHD@RHD_IVa.2_(DIV_Type1.0)::RHD_IVa.2_::RHD_IVa_(186T,410T,455C,1048C)',
           'RHCE@RHCE_ceAR':'RHD@RHD_DAR',
           'RHCE@RHCE_ceEK':'RHD@RHD_DAR',
           'RHCE@RHCE_ceMO':'RHD@RHD_DAU-0',
           'RHCE@RHCE_ceBl':'RHD@RHD_DOL',
           'RHCE@RHCE*Ce_733G::RHce_VS::RHCE*Ce_733G_':'RHD@Reference',
           'RHCE@RHce_ces':'RHD_III'               
           }
    resolvable_cases={'RHCE@RHce_ces':'RHCE@ce'}
    solution={'RHCE@RHce_ces':{'a1':'RHCE@RHce_48C::Rhce_e_variant', 'a2':'RHCE@RHCE*Ce_733G::RHce_VS::RHCE*Ce_733G_'}}
    bg1.reset_index(drop=True, inplace=True)
    bg2=pd.read_table(bg2fn, sep="\t", header=0)
    for r in resolvable_cases:
        if r in bg1['Bloodtype'].tolist() and resolvable_cases[r] in bg1['Bloodtype'].tolist():
           if links[r] in bg2['Bloodtype'].tolist():
               print ("No need to correct")
               break
           else:
               print ("Correcting...")
               ID= bg1['ID'].tolist()[0]
               bg1[:]='-'
               bg1['ID']=ID
               bg1.loc[0, 'Bloodtype']=solution[r]['a1'];
               bg1.loc[1, 'Bloodtype']=solution[r]['a2'];
               bg1.loc[0, 'Likelihood']=0
               bg1.loc[1, 'Likelihood']=0
    bg1=bg1.loc[bg1['Likelihood']==bg1['Likelihood'].max()] 
    return bg1
               
def RHD_RHCE_linked_type2(bg1, bg2fn, final_var, link_fn="linked_types.txt"):
    """
      correct the alleles based on linked information
      INFO based on ISBT
    """
    linkf=os.path.join(package_directory, 'database', link_fn)
    lf=pd.read_csv(linkf, sep='\t')
    lf['Allele name'] = lf['Allele name'].str.replace(' ', '_')    
    lf['Allele name'] = lf['Allele name'].astype(str)

    bg1.reset_index(drop=True, inplace=True)
    bg2=pd.read_csv(bg2fn, sep="\t", header=0)
    linked_index=[]
    bk_flag=False    

    

    for id1 in bg1['Index'].tolist():
        if bk_flag: break
        for id2 in bg2['Index'].tolist():
             bg1_a1=bg1.loc[ (bg1['Index']==id1) & (bg1['Allele']=='A1'),'Bloodtype'].values[0]; bg1_a1=re.sub(r'.*@','',bg1_a1)
             bg1_a2=bg1.loc[ (bg1['Index']==id1) & (bg1['Allele']=='A2'),'Bloodtype'].values[0]; bg1_a2=re.sub(r'.*@','',bg1_a2)
             bg2_a1=bg2.loc[ (bg2['Index']==id2) & (bg2['Allele']=='A1'),'Bloodtype'].values[0]; bg2_a1=re.sub(r'.*@','',bg2_a1)
             bg2_a2=bg2.loc[ (bg2['Index']==id2) & (bg2['Allele']=='A2'),'Bloodtype'].values[0]; bg2_a2=re.sub(r'.*@','',bg2_a2)

             #print('bg1', bg1_a1, bg1_a2)
             #print('bg2', bg2_a1, bg2_a2)

             if final_var.shape[0] > 0 and final_var.shape[0]==2 and final_var['gene'].iloc[0]=='RHCE' and 48 in final_var['cpos'].tolist() and 733 in final_var['cpos'].tolist():
                #print("48 733 linkage checking")
                if 'Deletion' in [bg2_a1, bg2_a2] and any(x=='RHD*01' for x in [bg2_a1, bg2_a2]) and \
                any(x=='RHCE*01_or_RHCE*ce_RHCE*c_RHCE*e' for x in [bg1_a1, bg1_a2]) and any(x=='RHCE*01.20.02_RHCE*ce.20.02_RHCE*ceVS.02' for x in [bg1_a1, bg1_a2]):
                      linked_index=[id1]
                      bk_flag=True
                      break
                elif bg2_a1=='RHD*01' and bg2_a2=='RHD*01' and \
                'RHCE*01_or_RHCE*ce_RHCE*c_RHCE*e' in [bg1_a1, bg1_a2] and 'RHCE*01.20.02_RHCE*ce.20.02_RHCE*ceVS.02' in [bg1_a1, bg1_a2]:
                      linked_index=[id1]
                      bk_flag=True
                      break                    

             #print("Other linkage checking")
             if bg1_a1 in lf['Allele name'].tolist() and not all(pd.isna(lf.loc[ (lf['Allele name']==bg1_a1), 'Linked'])):
                # print("check bg1_a1") 
                 bg1_a1_linked=lf.loc[lf['Allele name']==bg1_a1, 'Linked'].dropna().tolist()
                 if any(bg2_a1==x for x in bg1_a1_linked) or any(bg2_a2==x for x in bg1_a1_linked):
                     if id1 not in linked_index:
                        #print("[linked]", id1)
                        #print(bg1_a1_linked)
                        linked_index.append(id1)



             if bg1_a2 in lf['Allele name'].tolist() and not all(pd.isna(lf.loc[ (lf['Allele name']==bg1_a2), 'Linked'])):
                 #print("check bg1_a2") 
                 bg1_a2_linked=lf.loc[lf['Allele name']==bg1_a2, 'Linked'].dropna().tolist()
                 if any(bg2_a1==x for x in bg1_a2_linked) or any(bg2_a2==x for x in bg1_a2_linked):
                     if id1 not in linked_index:
                         #print("[linked]", id1)
                         #print(bg1_a2_linked)
                         linked_index.append(id1) 

    #print("Linked index")
    #print(linked_index)
    #print("Before")
    #print(bg1)

    ### get index with highest score
    if len(linked_index) > 0:
        bg1=bg1.loc[bg1['Index'].isin(linked_index),]
        #print("After")
        #print(bg1)

    return bg1
    



def RHCE_relink(args, final, final_var, cov, Ce_supportN=7, Ce_supportP=0.1, gbuild='hg38'):
    support_n, total_n, support_ratio=RHCe_insertion(args.bam, args.coverage, build=gbuild)
    print ('[109bp] INS support readn: {0}; total_readn: {1}; support_ratio: {2}'.format(support_n, total_n, support_ratio))
    #print ('[exon2] cov', cov.loc[cov['exon']==2, 'RHD_status'].values[0], cov.loc[cov['exon']==2, 'RHCE_status'].values[0])
    ### CE/Ce
    c_types=['RHCE*01','RHCE*03']
    C_types=['RHCE*02','RHCE*04']
    Cstatus=0
    if support_n >= Ce_supportN or support_ratio >= Ce_supportP: ### add coverage criteria here
        ### with one or two Ce allele
        index2drop=[]
        Cstatus=1
        Ccombination='C/c'
        if support_ratio >= Ce_supportP and support_ratio > 0.85: Cstatus=2
        if Cstatus == 1:
                print('[RHCE] C/c identified')
        elif Cstatus == 2:
                Ccombination='C/C'
                print('[RHCE] C/C identified')
        for IDX in final.Index.unique():
            A1=final.loc[ (final['Index']==IDX) & (final['Allele']=='A1'), 'Bloodtype'].values[0]
            A2=final.loc[ (final['Index']==IDX) & (final['Allele']=='A2'), 'Bloodtype'].values[0]
            ### correct for the case that only have 48c variants
            if Cstatus==1 and A1=='RHCE@RHCE*01.01_RHCE*ce.01' and A2=='RHCE@RHCE*01.01_RHCE*ce.01':
                print('correcting ce48C')
                A1='RHCE@RHCE*02_or_RHCE*Ce_RHCE*C_RHCE*e'
                final.loc[ (final['Index']==IDX) & (final['Allele']=='A1'), 'Bloodtype']='RHCE@RHCE*02_or_RHCE*Ce_RHCE*C_RHCE*e'
 
            C_count=0
            c_count=0
            ### Count C and c alleles for each allele pairs
            for c in c_types:
                if c in A1: c_count+=1
                if c in A2: c_count+=1
            for C in C_types:
                if C in A1: C_count+=1 
                if C in A2: C_count+=1

            if Cstatus == 1:
                ### C/c allele
                if C_count!=1:
                    index2drop.append(IDX)

            elif Cstatus == 2:
                ### C/C allele
                if C_count!=2:
                    index2drop.append(IDX)
        #print(index2drop)
        if len(index2drop) > 0 and len(index2drop) != len(final.Index.unique()):
            final=final[~final.Index.isin(index2drop)]
        elif len(index2drop) == len(final.Index.unique()):
            print('[RHCE] all allele pairs do not have', Ccombination, '; no filtering performed ...')

    else:
        ### without Ce allele
        print('[RHCE] c/c identified')
        #print(final)
        index2drop=[]
        for IDX in final.Index.unique():
            if any('RHCE*02' in x for x in final.loc[ (final['Index']==IDX), 'Bloodtype'].tolist()) or \
               any('RHCE*04' in x for x in final.loc[ (final['Index']==IDX), 'Bloodtype'].tolist()):
                if IDX not in index2drop: index2drop.append(IDX)
        if len(index2drop) > 0 and len(index2drop) != len(final.Index.unique()):
            final=final[~final.Index.isin(index2drop)]

    if args.linking_bloodtype:
            #final=RHD_RHCE_linked(final, args.linking_bloodtype) ### for SVM method
            ### for likelihood method
            if args.linking_database:
                print('[Linkage] using user database')
                final=RHD_RHCE_linked_type2(final, args.linking_bloodtype, final_var, link_fn=args.linking_database)
            else:
                print('[Linkage] using ISBT database')
                final=RHD_RHCE_linked_type2(final, args.linking_bloodtype, final_var)
    ### RHCE*cEKK/RHCE*cE
    ### Check the presnece of D1-3 
    final=RHCE_D1_3CE(final, cov) 

    final=final.loc[final['Likelihood']==final['Likelihood'].max()]
    return final


def RHCE_D1_3CE(final, cov):
    '''
        differentiate RHCE*03.02_RHCE*cE.02 and RHCE*03.01_RHCE*cE by checking the coverage of exon1-3
    '''
    if 'RHCE@RHCE*03.02_RHCE*cE.02' not in final.Bloodtype.unique().tolist() and \
       'RHCE@RHCE*03_orRHCE*cERHCE*c_RHCE*E' not in final.Bloodtype.unique().tolist():
        print('No bloodtype need to check for D1_3 exon in RHCE')
        return final

    D1_3=False
    if all( 'loss' in x for x in cov.loc[cov['exon'].isin([1,2,3]), 'RHCE_status'].tolist() ) and \
       all( 'gain' in x for x in cov.loc[cov['exon'].isin([1,2,3]), 'RHD_status'].tolist() ): D1_3=True

    index2drop=[]
    for idx in final.Index.unique():
        A1=final.loc[(final['Index']==idx) & (final['Allele']=='A1'),'Bloodtype'].values[0]
        A2=final.loc[(final['Index']==idx) & (final['Allele']=='A2'),'Bloodtype'].values[0]
        A1=re.sub('.*@','',A1); A2=re.sub('.*@','',A2)
        if A1=='RHCE*03.02_RHCE*cE.02' or A2=='RHCE*03.02_RHCE*cE.02':
            if not D1_3 and idx not in index2drop: index2drop.append(idx)
        else:
            if D1_3 and idx not in index2drop: index2drop.append(idx)
    final=final[~final.Index.isin(index2drop)]
    

    return final


 
def RHD_RHCE_exon_CNV(rhd, rhce, wgs_cov, prefix, ploid=2, verbose=0):
    from math import log
    out={}
    for idx, row in rhd.iterrows():
        acc=idx[0]; gene=idx[1]; exon=idx[2]
        cov=row['avg_QC_readn']
        if exon not in out: out[exon]={}
        if gene not in out[exon]: 
            out[exon][gene]={}
            out[exon][gene]['cov']=cov
            out[exon][gene]['log2R']=round(log(cov/wgs_cov,2),7) if cov/wgs_cov > 0 else 0
    for idx, row in rhce.iterrows():
        acc=idx[0]; gene=idx[1]; exon=idx[2]
        cov=row['avg_QC_readn']
        if exon not in out: out[exon]={}
        if gene not in out[exon]: 
            out[exon][gene]={}
            out[exon][gene]['cov']=cov
            out[exon][gene]['log2R']=round(log(cov/wgs_cov,2),7) if cov/wgs_cov > 0 else 0
    RHDpat={}
    RHCEpat={}
    fn=prefix+'.exonCNV.txt'

    #print out
    out2=[]
    for exon in out:
        RHDstatus=CNV_def(out[exon]['RHD']['log2R'], out[exon]['RHD']['cov'], wgs_cov)    
        RHCEstatus=CNV_def(out[exon]['RHCE']['log2R'], out[exon]['RHCE']['cov'], wgs_cov)
        if RHDstatus == 'normal' or RHCEstatus == 'normal':
            RHDpat[exon]=0
            RHCEpat[exon]=0
            out2.append({'wgs_cov':wgs_cov, 'exon':exon, 'RHD_cov': out[exon]['RHD']['cov'], 'RHD_log2cov':out[exon]['RHD']['log2R'], 'RHD_status':RHDstatus,
                                                     'RHCE_cov':out[exon]['RHCE']['cov'], 'RHCE_log2cov':out[exon]['RHCE']['log2R'],'RHCE_status':RHCEstatus})
            continue
        elif (RHDstatus == 'one-copy gain' and RHCEstatus == 'one-copy loss'):
            RHDpat[exon]=1
            RHCEpat[exon]=-1
        elif (RHDstatus == 'one-copy loss' and RHCEstatus == 'one-copy gain'):
            RHDpat[exon]=-1
            RHCEpat[exon]=1 
        elif (RHDstatus == 'two-copy loss' and RHCEstatus == 'two-copy gain'):
            RHDpat[exon]=-2
            RHCEpat[exon]=2 
        elif (RHDstatus == 'two-copy gain' and RHCEstatus == 'two-copy loss'):
            RHDpat[exon]=2
            RHCEpat[exon]=-2
        out2.append({'wgs_cov':wgs_cov, 'exon':exon, 'RHD_cov': out[exon]['RHD']['cov'], 'RHD_log2cov':out[exon]['RHD']['log2R'], 'RHD_status':RHDstatus,
                                                     'RHCE_cov':out[exon]['RHCE']['cov'], 'RHCE_log2cov':out[exon]['RHCE']['log2R'],'RHCE_status':RHCEstatus})
    out2=pd.DataFrame(out2)
    if verbose>10: print(out2)
    out2.to_csv(fn,sep='\t', index=False)    
    Dpat_blk=[]
    CEpat_blk=[]
    for exon in RHDpat:
        if exon==8:
            ### ignore exon 8 as they are exactly the same between D and CE except for 1136C>T
            Dpat_blk.append('D')
            CEpat_blk.append('CE')
        elif RHDpat[exon]==0:
            Dpat_blk.append('D')
            CEpat_blk.append('CE')
        elif RHDpat[exon]==1:
            Dpat_blk.append('D')
            CEpat_blk.append('D')
        elif RHDpat[exon]==-1:
            Dpat_blk.append('CE')
            CEpat_blk.append('CE') 
        elif RHDpat[exon]==-2:
            Dpat_blk.append('hCE')
            CEpat_blk.append('hCE')
        elif RHDpat[exon]==2:
            Dpat_blk.append('hD')
            CEpat_blk.append('hD')
        else:
            Dpat_blk.append('?D')
            CEpat_blk.append('?CE')
    Dpat=hybrid_call(Dpat_blk)
    CEpat=hybrid_call(CEpat_blk)
    Dpat="" if 'CE' not in Dpat else Dpat
    CEpat="" if 'D' not in CEpat else CEpat

    return Dpat, CEpat, out2

def RHD_RHCE_segment_CNV(rhd_coord, rhce_coord, rhd, rhce, bam, wgs_cov, fasta, prefix, ploid=2, flanking=300, verbose=0, mode='bin_med.wgsratio.log2', chr='chr1'):
    '''
       expand the coverage check region outside of exon
    '''
    from math import log
    print("[RHD_RHCE_segment_CNV]")
    np.seterr(divide='ignore', invalid='ignore')
    RHD_start=rhd_coord.start.min()    
    RHD_end=rhd_coord.end.max()
    RHCE_start=rhce_coord.start.min()
    RHCE_end=rhce_coord.end.max()


    ### flanking region by extending 18700 to start and end positions (ff, revisions)
    ff=0
    RHDbin, RHDgene=perbaseCov(bam, wgs_cov, fasta, prefix=prefix+'.RHD', chr=chr, start=RHD_start-ff, end=RHD_end+ff, coord=rhd_coord, flanking=flanking, save_output=True)
    RHCEbin, RHCEgene=perbaseCov(bam, wgs_cov, fasta, prefix=prefix+'.RHCE', chr=chr, start=RHCE_start-ff, end=RHCE_end+ff, coord=rhce_coord, flanking=flanking, save_output=True)
   
    #RHDbin, RHDgene=perbaseCov(bam, wgs_cov, fasta, prefix=prefix+'.RHD', chr=chr, start=RHD_start, end=RHD_end, coord=rhd_coord, flanking=19000,save_output=True)
    #RHCEbin, RHCEgene=perbaseCov(bam, wgs_cov, fasta, prefix=prefix+'.RHCE', chr=chr, start=RHCE_start, end=RHCE_end, coord=rhce_coord, flanking=19000,save_output=True) 


    RHDbin = gmm_segmentation(RHDbin, mode=mode) 
    RHDbin2= cbs_segmentation(RHDbin, w=10, p=0.001, mode=mode)
   
    RHCEbin = gmm_segmentation(RHCEbin, mode=mode)
    RHCEbin2= cbs_segmentation(RHCEbin, w=10, p=0.001, mode=mode)
 
    #perbinCov(bam, wgs_cov, fasta, bin=100, step=50, prefix=prefix+'RHD', chr=chr, start=RHD_start, end=RHD_end, coord=rhd_coord, flanking=flanking)
    #perbinCov(bam, wgs_cov, fasta, bin=100, step=50, prefix=prefix+'RHCE', chr=chr, start=RHCE_start, end=RHCE_end, coord=rhce_coord, flanking=flanking)

    ### plot
    mode='cbs_segmean'
    RHDexon, RHDperexon=perexon_Cov(RHDbin2, mode=mode)
    RHCEexon, RHCEperexon=perexon_Cov(RHCEbin2, mode=mode)
   
    RHD_RHCE_segment_CNV_plot(RHDbin2, wgs_cov, prefix=prefix+'RHD.bincov')
    RHD_RHCE_segment_CNV_plot(RHCEbin2, wgs_cov, prefix=prefix+'RHCE.bincov')
   
    RHD_s, RHD_e=RHD_RHCE_segment_CNV_plot(RHDbin2, wgs_cov, prefix=prefix+'RHD.bincov', exon_only=False)
    RHCE_s, RHCE_e=RHD_RHCE_segment_CNV_plot(RHCEbin2, wgs_cov, prefix=prefix+'RHCE.bincov', exon_only=False)

    ### define deletion status
    RHDpat={}
    RHCEpat={}
    fn=prefix+'.exonCNV.txt'

    out={}
    for idx, row in rhd.iterrows():
        acc=idx[0]; gene=idx[1]; exon=idx[2]
        cov=row['avg_QC_readn']
        if exon not in out: out[exon]={}
        if gene not in out[exon]:
            out[exon][gene]={}
            out[exon][gene]['cov']=cov
            out[exon][gene]['log2R']=round(log(cov/wgs_cov,2),7) if cov/wgs_cov > 0 else 0
            out[exon][gene][mode]=RHDexon[exon] if exon in RHDexon else 0 
    for idx, row in rhce.iterrows():
        acc=idx[0]; gene=idx[1]; exon=idx[2]
        cov=row['avg_QC_readn']
        if exon not in out: out[exon]={}
        if gene not in out[exon]:
            out[exon][gene]={}
            out[exon][gene]['cov']=cov
            out[exon][gene]['log2R']=round(log(cov/wgs_cov,2),7) if cov/wgs_cov > 0 else 0
            out[exon][gene][mode]=RHCEexon[exon] if exon in RHCEexon else 0


    out2=[]
    
    for exon in out:
	#RHDstatus=CNV_def(out[exon]['RHD'][mode], out[exon]['RHD']['cov'], wgs_cov)
        #RHCEstatus=CNV_def(out[exon]['RHCE'][mode], out[exon]['RHCE']['cov'], wgs_cov)    
        RHDstatus=CNV_def_with_test(exon, RHDperexon, out[exon]['RHD'][mode], out[exon]['RHD']['cov'], wgs_cov)
        RHCEstatus=CNV_def_with_test(exon, RHCEperexon, out[exon]['RHCE'][mode], out[exon]['RHCE']['cov'], wgs_cov)

        if RHDstatus == 'normal' and RHCEstatus == 'normal':
            RHDpat[exon]=0
            RHCEpat[exon]=0
            out2.append({'wgs_cov':wgs_cov, 'exon':exon, 'RHD_cov': out[exon]['RHD']['cov'], 'RHD_log2cov':out[exon]['RHD'][mode], 'RHD_status':RHDstatus,
                                                     'RHCE_cov':out[exon]['RHCE']['cov'], 'RHCE_log2cov':out[exon]['RHCE'][mode],'RHCE_status':RHCEstatus})
            continue
        
        RHDpat[exon]=CNV_def_switch(RHDstatus)
        RHCEpat[exon]=CNV_def_switch(RHCEstatus)
        
        out2.append({'wgs_cov':wgs_cov, 'exon':exon, 'RHD_cov': out[exon]['RHD']['cov'], 'RHD_log2cov':out[exon]['RHD'][mode], 'RHD_status':RHDstatus,
                                                     'RHCE_cov':out[exon]['RHCE']['cov'], 'RHCE_log2cov':out[exon]['RHCE'][mode],'RHCE_status':RHCEstatus})
    out2=pd.DataFrame(out2)
    #if verbose>10: print out2
    out2.to_csv(fn,sep='\t', index=False)
    Dpat_blk=[]
    CEpat_blk=[]
    for exon in RHDpat:
        if exon==8:
           ### ignore exon 8 as they are exactly the same between D and CE except for 1136C>T
           Dpat_blk.append('D')
           CEpat_blk.append('CE')
        elif RHDpat[exon]==0:
            Dpat_blk.append('D')
            CEpat_blk.append('CE')
        elif RHDpat[exon]==1:
            Dpat_blk.append('D')
            CEpat_blk.append('D')
        elif RHDpat[exon]==-1:
            Dpat_blk.append('CE')
            CEpat_blk.append('CE')
        elif RHDpat[exon]==-2:
            Dpat_blk.append('hCE')
            CEpat_blk.append('hCE')
        elif RHDpat[exon]==2:
            Dpat_blk.append('hD')
            CEpat_blk.append('hD')
        else:
            Dpat_blk.append('?D')
            CEpat_blk.append('?CE')
    Dpat=hybrid_call(Dpat_blk)
    CEpat=hybrid_call(CEpat_blk)
    Dpat="" if 'CE' not in Dpat else Dpat
    CEpat="" if 'D' not in CEpat else CEpat

    return Dpat, CEpat, out2, RHD_s, RHD_e, RHCE_s, RHCE_e


def RHD_RHCE_segment_CNV_plot(bin, wgs_cov, prefix='test', mode='bin_med.wgsratio.log2', exon_only=True):
    '''
        Take bin for plotting
    '''
    n = 12
    colors = pl.cm.jet(np.linspace(0,1,n))
    pd_lab=OrderedDict()
    max_wgsRatio=-1
    pos=0
    pd=OrderedDict()
    for p_s in bin: 
        ### limit run to exon region only
        if exon_only and bin[p_s]['bin_exon'] == -1: continue
        pos+=1
        pd_lab[pos]='Exon'+str(bin[p_s]['bin_exon'])+'.'+str(p_s) if bin[p_s]['bin_exon'] != -1 else str(p_s)
        ratio_wgs=bin[p_s][mode]
        if ratio_wgs > max_wgsRatio: max_wgsRatio=ratio_wgs
        if bin[p_s]['bin_exon'] not in pd:
            pd[bin[p_s]['bin_exon']]={}
            pd[bin[p_s]['bin_exon']]['xpos']=[]
            pd[bin[p_s]['bin_exon']]['xlab']=[]
            pd[bin[p_s]['bin_exon']]['y']=[]
            pd[bin[p_s]['bin_exon']]['norm_y_gmm']=[]
            pd[bin[p_s]['bin_exon']]['norm_y_cbs']=[]
            pd[bin[p_s]['bin_exon']]['col']=[]
        pd[bin[p_s]['bin_exon']]['xpos'].append(pos)
        pd[bin[p_s]['bin_exon']]['y'].append(ratio_wgs)
        pd[bin[p_s]['bin_exon']]['col'].append(colors[bin[p_s]['bin_exon']+1])
        pd[bin[p_s]['bin_exon']]['norm_y_gmm'].append(bin[p_s]['gmm_segcenter'])
        pd[bin[p_s]['bin_exon']]['norm_y_cbs'].append(bin[p_s]['cbs_segmean'])

    s_ix, e_ix, break_s, break_e=change_point(bin, wgs_cov, mode='cbs_segmean')

    ### plot normalized coverage
    if exon_only:
        plt.figure(figsize=(20, 7))
    else:
        plt.figure(figsize=(30, 7))
    plt.ticklabel_format(useOffset=False)
    plt.xlim(-0.5,pos+0.5)
    #plt.ylim(-0.1, max_wgsRatio+0.1) if max_wgsRatio+0.1 > 2 else plt.ylim(-0.1, 2)
    plt.ylim(-2, 2)
    for e in pd:
        plt.scatter(pd[e]['xpos'], pd[e]['y'], c=pd[e]['col'], label=e)
        #plt.plot(pd[e]['xpos'], pd[e]['norm_y_gmm'], c='red')
        plt.plot(pd[e]['xpos'], pd[e]['norm_y_cbs'], c='blue')
    if not exon_only:
        if s_ix is not None: plt.axvline(x=s_ix, c='grey')
        if e_ix is not None: plt.axvline(x=e_ix, c='grey')
    if len(list(pd_lab.keys()))> 50:
        plt.xticks(list(pd_lab.keys())[::10],list(pd_lab.values())[::10], rotation=90)
    else:
        plt.xticks(list(pd_lab.keys()),list(pd_lab.values()), rotation=90)
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    outfn=prefix+'.exonpos.png' if exon_only else prefix+'.allpos.png'
    plt.savefig(outfn, bbox_inches = 'tight')
    return break_s, break_e


def perbinCov(bam, wgs_cov, fasta, bin=None, step=None, prefix=None, chr=None, start=None, end=None, coord=None, minimum_base_quality = 15, minimum_mapq = 10, minimum_read_quality = 15, flanking=1000):
        #calculate per-bin coverage
    np.seterr(divide='ignore', invalid='ignore')
    samfile = pysam.AlignmentFile(bam, mode="rb")
    trimchr=False
    if not all('chr' in s for s in samfile.references):trimchr=True
    tarchr=chr.lstrip('chr') if trimchr else chr
    
    out=OrderedDict()
    fa=pysam.FastaFile(fasta)

    test_allbin=[]
    test_exonbin=[]

    for b in range(start-flanking, end+flanking, step):
        bin_start=b
        bin_end=b+bin if (b+bin) < (end+flanking) else end+flanking
        bin_GC=GCpercent(fa.fetch(chr, bin_start, bin_end))
        bin_readn=0
        for read in samfile.fetch(chr, bin_start, bin_end):
            if read.mapping_quality < minimum_mapq:
                continue
            elif read.mapping_quality == 255:
                continue
            elif mean(read.query_qualities) < minimum_read_quality:
                continue
            bin_readn+=1
        if bin_start not in out:
            out[bin_start]={}
            out[bin_start]['bin_end']=bin_end
            out[bin_start]['bin_readN']=bin_readn
            out[bin_start]['bin_GC']=bin_GC
            test_allbin.append(bin_readn) 
    for s in out:
       e=out[s]['bin_end']
       out[s]['bin_normalized_readN']=out[s]['bin_readN']/mean(test_allbin)
       ### add exon info
       if coord is not None:
           for ix, row in coord.iterrows():
               exon=int(row['exon_n'])+1
               if ( row['start']-flanking <= s and s <= row['end']+flanking ) or ( row['start']-flanking <= e and e <= row['end']+flanking ):
                  out[s]['bin_exon']=exon
                  #test_exonbin.append(np.log2(out[s]['bin_normalized_readN']))
                  test_exonbin.append(out[s]['bin_normalized_readN'])
               else:
                  out[s]['bin_exon']=-1
       
    #N, out=Cov_correctForGC(out, prefix, mode='perbin')
       
    #N0 = np.asarray(test_exonbin)
    #L= cbs.segment(N0,p=0.01)
    #S= cbs.validate(N0, L, p=0.01)
    #ax=cbs.draw_segmented_data(N0,  S, title='Circular Binary Segmentation of Data')
    #ax.get_figure().savefig(prefix+'.png')
    return out, test_allbin, test_exonbin

def perbaseCov(bam, wgs_cov, fasta, prefix=None, chr=None, start=None, end=None, coord=None, bin=100, step=50, minimum_base_quality = 15, minimum_mapq = 10, minimum_read_quality = 15, flanking=1000, includeZero=True, exonOnly=True, save_output=False):
    '''
    #   calculate perbase coverage 
    '''
    np.seterr(divide='ignore', invalid='ignore')
    samfile = pysam.AlignmentFile(bam, mode="rb") #, ignore_truncation=False)
    trimchr=False
    if not all('chr' in s for s in samfile.references):trimchr=True
    tarchr=chr.lstrip('chr') if trimchr else chr
    pos_n=0
    total_good_n=0
    RHD_perbase_cov={}
    out={}
    for column in samfile.pileup(chr, start-flanking, end+flanking):
       pos0based=column.pos
       pos1based=column.pos+1
       good_count={'A':0, 'T':0, 'C':0, 'G':0, 'N':0}
       good_n=0
       clipped_n=0
       discordant_n=0
       pos_n+=1
       for pileupread in column.pileups:
          read     = pileupread.alignment
          qpos     = pileupread.query_position
          qname    = pileupread.alignment.query_name
          # skip reads with indels
          if pileupread.is_del: continue
          # skip reads with mapq below threshold
          if pileupread.alignment.mapping_quality < minimum_mapq:
              continue
          elif read.mapping_quality == 255:
              continue
          elif mean(read.query_qualities) < minimum_read_quality:
              continue
          else:
              if read.query_qualities[qpos] < minimum_base_quality:
                  continue
          good_n+=+1
          total_good_n+=1
          if 'S' in read.cigarstring:
              clipped_n+=1
          if not read.is_proper_pair:
              discordant_n+=1
       if pos1based not in out:
           out[pos1based]={}
           out[pos1based]['totalReadN']=column.n
           out[pos1based]['QCtotalReadN']=good_n
           out[pos1based]['QCclippedReadN']=clipped_n
           out[pos1based]['QCdiscordantReadN']=discordant_n
           out[pos1based]['QCtotalRatio']=column.n/wgs_cov
           out[pos1based]['QCclippedRatio']=good_n/wgs_cov
           out[pos1based]['QCdiscordantRatio']=(clipped_n+discordant_n)/wgs_cov
           out[pos1based]['exon']=-1
           out[pos1based]['real_exon']=-1

    if includeZero:
        #### add the zero-covered bases
        for i in range(start-flanking, end+flanking):
            if i not in out:
                out[i]={}
                out[i]['totalReadN']=0
                out[i]['QCtotalReadN']=0
                out[i]['QCclippedReadN']=0
                out[i]['QCdiscordantReadN']=0
                out[i]['QCtotalRatio']=0
                out[i]['QCclippedRatio']=0
                out[i]['QCdiscordantRatio']=0
                out[i]['exon']=-1
                out[i]['real_exon']=-1
             

    ### add exon number
    #exon2=[]
    #exon2_pos=[]

    if coord is not None:
        for ix, row in coord.iterrows():
            exon=int(row['exon_n'])+1
            for p in range(row['start']-flanking, row['end']+flanking):
                if p in out: 
                    out[p]['exon']=exon
                    if row['start'] <= p <= row['end']:
                        out[p]['real_exon']=exon 

                #if exon==2 and row['start'] <= p <= row['end']:
                #    exon2.append(out[p]['QCtotalReadN']) 
                #    exon2_pos.append(p)
    #print(min(exon2_pos), max(exon2_pos), mean(exon2))
            
    if save_output:
       fh=open(prefix+'.binZ'+str(bin)+'.coverage.txt', 'w')
       headers=['prefix', 'position', 'totalReadN','QCtotalReadN','QCclippedReadN','QCdiscordantReadN','QCtotalRatio','QCclippedRatio','QCdiscordantRatio','exon_n_flanking', 'exon', 'wgs_cov']
       fh.write('\t'.join(headers)+'\n')
       for i in range(start-flanking, end+flanking):
           outstr=map(str, [prefix, i, out[i]['totalReadN'], out[i]['QCtotalReadN'], out[i]['QCclippedReadN'], out[i]['QCdiscordantReadN'], out[i]['QCtotalRatio'], out[i]['QCclippedRatio'], out[i]['QCdiscordantRatio'], out[i]['exon'], out[i]['real_exon'], wgs_cov])
           fh.write('\t'.join(outstr)+'\n')
       fh.close()

 
    ### use bin as one datapoint
    out2=OrderedDict()
    fa=pysam.FastaFile(fasta)
    all_bins=[]

    for b in range(start-flanking, end+flanking, step):
        bin_start=b
        if (b+bin) < (end+flanking):
            bin_end=b+bin
        else: 
            bin_end=end+flanking
        bins=[]
        exon=[]
             
        for p in range(bin_start, bin_end):
            if p not in out:
                continue
            bins.append(out[p]['totalReadN'])
            all_bins.append(out[p]['totalReadN'])
            if out[p]['exon'] not in exon: exon.append(out[p]['exon'])

        bins=np.array(bins)
        bin_avg=bins.mean()
        bin_median=np.median(bins)
        bin_GC=GCpercent(fa.fetch(chr, bin_start, bin_end))
        if len(exon) == 1:
            #print bin_start, bin_end, bin_avg, bin_median, exon
            if bin_start not in out2:
                out2[bin_start]={}
                out2[bin_start]['bin_end']=bin_end
                out2[bin_start]['bin_avg']=bin_avg
                out2[bin_start]['bin_avg.wgsratio']=bin_avg/wgs_cov
                out2[bin_start]['bin_avg.wgsratio.log2']=np.log2(bin_avg/wgs_cov) if bin_avg != 0 else np.log2(0.01/wgs_cov)
                out2[bin_start]['bin_med']=bin_median
                out2[bin_start]['bin_med.wgsratio']=bin_median/wgs_cov
                out2[bin_start]['bin_med.wgsratio.log2']=np.log2(bin_median/wgs_cov) if bin_median != 0 else np.log2(0.01/wgs_cov)
                out2[bin_start]['bin_exon']=exon[0]
                out2[bin_start]['bin_GC']=bin_GC
            if -1 in exon and exonOnly: 
                continue
            #if 8 in exon: continue

    ### correct for GC bias, not needed
    #out2=Cov_correctForGC(out2, prefix)
    all_bins=np.array(all_bins)
    all_bins_avg=all_bins.mean()
    all_bins_med=np.median(all_bins)
    for s in out2:
        out2[s]['bin_avg.regionratio']=out2[s]['bin_avg']/all_bins_avg
        out2[s]['bin_med.regionratio']=out2[s]['bin_med']/all_bins_avg
        out2[s]['bin_avg.regionratio.log2']=np.log2(out2[s]['bin_avg']/all_bins_avg) if out2[s]['bin_avg'] != 0 else np.log2(0.01/all_bins_avg) 
        out2[s]['bin_med.regionratio.log2']=np.log2(out2[s]['bin_med']/all_bins_avg) if out2[s]['bin_med'] != 0 else np.log2(0.01/all_bins_avg) 
    return out2, all_bins_med


def change_point(input, mean_cov, mode='gmm_segcenter'):
    '''Given x, Compute the subinterval x[i0:i1] with the maximal segmentation statistic t. 
    Returns t, i0, i1'''
    np.seterr(divide='ignore', invalid='ignore')
    test=[]
    cov=[]
    for s in input:
        test.append(input[s][mode])
        cov.append(input[s]['bin_med'])
    x=np.array(test).reshape(-1,1)
    cov=np.array(cov).reshape(-1,1)
    ### determine the difference from the median 
    L=refine_change_point(x,0,len(x)-1)
    if len(L)==0: return None, None, None, None
    max_t=0
    for vals in L:
        if vals[0] > max_t:
            i=vals[1]
            j=vals[2]

    
    ### mean value between i, j and outside of i,j
    Y=np.median(x[i+1:j])
    Ycov=np.median(cov[i+1:j])
    Z=np.median(np.delete(x, np.arange(i+1,j)))
    Zcov=np.median(np.delete(cov, np.arange(i+1,j)))
    Ystatus=CNV_def(Y, Ycov, mean_cov) 
    Zstatus=CNV_def(Z, Zcov, mean_cov)
 
    ###
    for s, val in enumerate(input):
        if s==i: bk_s=val
        if s==j: bk_e=val
    ### no breakpoint report if no change observed in the maximum interval
    if Ystatus == 'normal' and Zstatus == 'normal':
        i=None; j=None; bk_s=None; bk_e=None 
    return i, j, bk_s, bk_e

def refine_change_point(x,start,end,L=[], w=10, init_t=0, prev_start=None, prev_end=None):
    t0,i,j=cbs.max_tstat(x[start:end])
    if end > len(x)-1: end=len(x)-1
    t=cbs.tstat2(x, start, end)
    if i < w:
        i = 0
    if len(x)-j< w:
        j = len(x)

    if t > init_t:
        L.append([t, start, end])
    if (j-i) > w and (j-i) != (end-start):
        if j-i>0 :
            init_t=t if t > init_t else init_t
            if i != prev_start and j != prev_end:
                refine_change_point(x, start+i, start+j-1, L=L, init_t=init_t, w=w, prev_start=i, prev_end=j)
    return L
     


def gmm_segmentation(x, mode='bin_med.wgsratio.log2'):
    '''
        #use GMM for splitting CNV segment
        #mode: bin_avg.wgsratio, bin_med.wgsratio, bin_avg.regionratio, bin_med.regionratio 
        #      and .log2 version
    '''
    test=[]
    test_exon=[]
          
    for s in x:
        test.append(x[s][mode])
        if x[s]['bin_exon'] != -1:
            test_exon.append(x[s][mode])

    test = np.array(test).reshape(-1,1)
    test_exon=np.array(test_exon).reshape(-1,1)

    lowest_bic=np.infty
    best_k=0
    best_CV_type=None
    bic=[]
    cv_types = ['spherical', 'tied', 'diag', 'full']
    #cv_types = ['tied']
    n_components_range = range(2, 5)

    train=test
    for k in n_components_range:
        for cv_type in cv_types:
            g=mixture.GaussianMixture(n_components=k, covariance_type=cv_type)
            g.fit(train)
            bic.append(g.bic(train))
            #print k, cv_type, g.bic(train)
            if bic[-1] < lowest_bic:
                lowest_bic=bic[-1]
                best_gmm=g
                best_k=k
                best_CV_type=cv_type
    pred = best_gmm.predict(test)
    ix=0

    out2=x
    segn_d={}
    ### match up prediction and coverage
    for p_s in out2:
        out2[p_s]['gmm_predicted']=pred[ix]
        if pred[ix] not in segn_d:
            segn_d[pred[ix]]=[]
        segn_d[pred[ix]].append(out2[p_s][mode])
        ix+=1

    ### print result
    header=['bin.s','bin.e','bin.median','bin.mean','bin.exon','predicted.status','center','ratio.wgs','ratio.region']
    #print '\t'.join(map(str, header))

    for p_s in out2:
        ### limit run to exon region only
        #if out2[p_s]['bin_exon'] != -1:
             center=best_gmm.means_[out2[p_s]['gmm_predicted']][0] if out2[p_s]['gmm_predicted'] is not None else None
             out2[p_s]['gmm_segcenter']=center
             out2[p_s]['gmm_segmean']=np.mean(segn_d[out2[p_s]['gmm_predicted']])
    return out2

def cbs_segmentation(x, mode='bin_med.wgsratio.log2', w=3, p=0.05):
    '''
        #use circular binary segmentation to split CNV segment
        #mode: bin_avg.wgsratio, bin_med.wgsratio, bin_avg.regionratio, bin_med.regionratio 
        #      and log2 version
    '''
    test=[]
    test_exon=[]

    for s in x:
        test.append(x[s][mode])
        if x[s]['bin_exon'] != -1:
            test_exon.append(x[s][mode])

    test = np.array(test).reshape(-1,1)
    test_exon=np.array(test_exon).reshape(-1,1)

    N0 = np.asarray(test)
    L= cbs.segment(N0,p=p, w=w)
    S= cbs.validate(N0, L, p=p)
    segn=0
    segn_d={}
    for ix, breakpoint in enumerate(S[:-1]):
        #if ix==len(S): continue ### skip the end point
        start=S[ix]
        end=S[ix+1]
        for pos, val in enumerate(x):
            if pos >= start and pos <= end:
                x[val]['cbs_predicted']=segn
                if segn not in segn_d:
                    segn_d[segn]=[]
                segn_d[segn].append(x[val][mode])
        segn+=1
    for pos in x:
        x[pos]['cbs_segmean']=np.mean(segn_d[x[pos]['cbs_predicted']])

    #ax=cbs.draw_segmented_data(N0,  S, title='Circular Binary Segmentation of Data')
    return x


def Cov_correctForGC(df, prefix, mode='perbase'):
    '''
         correct coverage for GC bias
    '''
    print('[correct GC bias for coverage]')
    sns.set_style('darkgrid')
    GC=[]
    cov=[]
    for pos in df:
        #print pos, df[pos]['bin_med'], df[pos]['bin_GC'] 
        if mode=='perbase':
            cov.append(df[pos]['bin_avg.wgsratio']) 
        else:
            cov.append(df[pos]['bin_normalized_readN'])
        GC.append(df[pos]['bin_GC'])    

 
    j0=sns.scatterplot(GC,cov)
    t=j0.set_ylabel('Normalized readN')
    t=j0.set_xlabel('GC Proportion')
    t=j0.set_title('Normalized readN vs GC Content')
    fig=j0.get_figure()
    fig.set_size_inches(15, 7)
    plt.savefig(prefix+'.lowess.png')
    plt.clf()

    fig, (j,j1) = plt.subplots(2,1)
    np.seterr(divide='ignore', invalid='ignore')
    print ('[fit lowess]')
    B1=lowess(np.log(cov),GC,frac=0.3)   
    B2=lowess(np.log(cov),GC,frac=0.5)
    np.seterr(divide='warn', invalid='warn')

    sns.scatterplot(GC,cov,ax=j)
    sns.lineplot(B1[:,0],B1[:,1],color='red',label='smoothing parameter f=.3',ax=j)
    sns.lineplot(B2[:,0],B2[:,1],color='green',label='smoothing parameter f=.5',ax=j)
    h,l=j.get_legend_handles_labels()
    j.set_ylabel('Normalized readN')
    j.set_xlabel('GC proportion')
    j.set_title('LOWESS Smoothing for Fragment Counts')
    j.legend(h,l)
    #j.set_ylim([-2.5,2.5])
    N = Cov_lowess(GC,cov,f=0.5)
    sns.scatterplot(GC,N,ax=j1)
    #r=j1.set_ylim([-2.5,2.5])
    j1.set_ylabel('Bias Corrected Counts')
    j1.set_xlabel('GC Proportion')
    s=j1.set_title('Counts corrected by LOWESS with f=0.5')
    fig=j1.get_figure()
    fig.set_size_inches(15, 15)
    plt.savefig(prefix+'.lowess.corrected.png')
    plt.clf()

    Nix=0
    for pos in df:
        df[pos]['bin_GCcorrected_readN']=N[Nix]
        Nix+=1
    return N, df

def Cov_lowess(x, y, f=.5):
    '''
        fit the bin wise cov and GC ratio to a lowess curve
        for GC correction
    '''
    np.seterr(divide='ignore', invalid='ignore')
    jlow=lowess(np.log(y), x, frac=f)
    jz=np.interp(x, jlow[:,0], jlow[:,1])
    return np.log(y)-jz

def GCpercent(seq):
    '''
        calculate the GC percentage of the seq
    '''
    GCnum=0
    ATnum=0
    ALLnum=0
    for i in seq:
        ALLnum+=1
        if i.upper() in ['C','G']:
            GCnum+=1
        elif i.upper() in ['A','T']:
            ATnum+=1
    return GCnum/float(GCnum+ATnum)


def perexon_Cov(df, mode='bin_med.wgsratio.log2'):
    '''
        use per base cov to obtain the exon median coverage
    '''
    exon=OrderedDict()
    for pos in df:
        if df[pos]['bin_exon'] not in exon:
            exon[df[pos]['bin_exon']]=[]
        if mode in  df[pos]:
            exon[df[pos]['bin_exon']].append(df[pos][mode])
    out=OrderedDict()
    for e in exon:
        if e != -1:
            if e not in out:
                out[e]=median(exon[e])
    return out, exon

    

 
def hybrid_call(blk):
    '''
       return hybrid allele if any found
    '''
    Dstr=''
    prevS=''
    begin=''
    for i, val in enumerate(blk):
        exon=i+1
        #print i, exon, val
        if i == 0:
            begin=exon
        else:
            if blk[i] != blk[i-1]:
                Dstr+=blk[i-1] + '(' + str(begin) + '-' + str(exon-1) + ')' if str(begin) != str(exon-1) else blk[i-1] + '(' + str(begin) + ')'
                begin=exon
            elif exon==len(blk):
                Dstr+=blk[i] + '(' + str(begin) + '-' + str(exon) + ')' if str(begin) != str(exon) else blk[i] + '(' + str(begin) + ')'
        prevS=blk[i]
    return Dstr
        

def RHhybrid_rename(final):
    '''
       rename alleles based on the known hybrid alleles (only if no SNP at all)
    '''
    #hybrid_homozy_alleles={'':'RHD*01N.07'}    

    hybrid_hetero_alleles={'CE(1-9)D(10)':'RHD*01N.02',
                    'D(1)CE(2-9)D(10)':'RHD*01N.03',
                    'D(1-2)CE(3-9)D(10)':'RHD*01N.04',
                    'D(1)CE(2-7)D(8-10)':'RHD*01N.05',
                    'D(1-3)CE(4-7)D(8-10)':'RHD*01N.07',
                    'CE(1)D(2-6)CE(7-10)':'RHD*01N.42',
                    'CE(1-3)D(4-10)':'RHD*01N.43',
                    'CE(1-3)D(4-9)':'RHD*01N.43',
                    'CE(1)D(2)CE(3)D(4-10)':'RHD*01N.43',
                    'D(1-4)CE(5-7)D(8-10)':'RHD*01EL.23 RHD*DEL23',
                    'D(1-3)CE(4-9)D(10)':'RHD*01EL.44 RHD*DEL44'
                    }

    hybrid_hetero_alleles_alias={'RHD*01N.02':'RHD*CE(1-9)-D',
                    'RHD*01N.03':'RHD*D-CE(2-9)-D',
                    'RHD*01N.04':'RHD*D-CE(3-9)-D',
                    'RHD*01N.05':'RHD*D-CE(2-7)-D',
                    'RHD*01N.07':'RHD*D-CE(4-7)-D',
                    'RHD*01N.42':'RHD*CE(1)-D(6)-CE(7-10)',
                    'RHD*01N.43':'RHD*CE(1-3)-D',
                    'RHD*01EL.23 RHD*DEL23':'RHD*-CE(5-7)-D',
                    'RHD*01EL.44 RHD*DEL44':'RHD*-CE(4-9)-D'
                    }
   
    if len(final.Index.unique())==0:
        print("[ERROR] [RHhybrid_rename] empty dataframe")

 
    for idx in final.Index.unique():
        A1=final.loc[(final['Index']==idx) & (final['Allele']=='A1'), 'Bloodtype'].values[0]
        A2=final.loc[(final['Index']==idx) & (final['Allele']=='A2'), 'Bloodtype'].values[0]

        ### first correct hybrid for deletion ... should not call hybrid at all
        if 'deletion' in A1.lower() or 'deletion' in A2.lower():
            final.loc[(final['Index']==idx) & (final['Allele']=='A1'), 'CoverageProfile']=""
            final.loc[(final['Index']==idx) & (final['Allele']=='A1'), 'Alias']=""
            final.loc[(final['Index']==idx) & (final['Allele']=='A2'), 'CoverageProfile']=""
            final.loc[(final['Index']==idx) & (final['Allele']=='A2'), 'Alias']=""
    
        C1=final.loc[(final['Index']==idx) & (final['Allele']=='A1'), 'CoverageProfile'].values[0]
        C2=final.loc[(final['Index']==idx) & (final['Allele']=='A2'), 'CoverageProfile'].values[0]
        C1=re.sub(r'\s+','',C1)
        C2=re.sub(r'\s+','',C2)

         
        ### correct the hybrid calling for the known hybrid alleles
        if A1=='RHD*03N.01' and C1 in ['D(1-3)CE(4-7)D(8-10)','D(1-3)CE(4-5)D(6)CE(7)D(8-10)']:
            if A2 != 'RHD*03N.01': 
                final.loc[(final['Index']==idx) & (final['Allele']=='A2'), 'CoverageProfile']=""
                final.loc[(final['Index']==idx) & (final['Allele']=='A2'), 'Alias']=""
                C2=""
        if A2=='RHD*03N.01' and C2 in ['D(1-3)CE(4-7)D(8-10)','D(1-3)CE(4-5)D(6)CE(7)D(8-10)']:
            if A1 != 'RHD*03N.01': 
                final.loc[(final['Index']==idx) & (final['Allele']=='A1'), 'CoverageProfile']=""
                final.loc[(final['Index']==idx) & (final['Allele']=='A1'), 'Alias']=""
                C1=""

        ALLvar1=final.loc[(final['Index']==idx) & (final['Allele']=='A1'), 'AllVariations'].values[0]
        ALLvar2=final.loc[(final['Index']==idx) & (final['Allele']=='A2'), 'AllVariations'].values[0]
        
        ALLvar1=ALLvar1.split(';');
        for x in ALLvar1:
            if 'ins' or 'del' in x:
                ALLvar1.remove(x)
        ALLvar2=ALLvar2.split(';');
        for x in ALLvar2:
            if 'ins' or 'del' in x:
                ALLvar2.remove(x)



        ### check homo alleles first
        ### no example to continue
        ### then check hete alleles
        if 'deletion' not in A1.lower() and 'deletion' not in A2.lower(): 
            if len(ALLvar1)==0 and C2 in hybrid_hetero_alleles:
                final.loc[(final['Index']==idx) & (final['Allele']=='A2'), 'Bloodtype']=hybrid_hetero_alleles[C2]
                final.loc[(final['Index']==idx) & (final['Allele']=='A2'), 'Alias']=hybrid_hetero_alleles_alias[hybrid_hetero_alleles[C2]]
                continue
            if len(ALLvar2)==0 and C1 in hybrid_hetero_alleles:
                final.loc[(final['Index']==idx) & (final['Allele']=='A1'), 'Bloodtype']=hybrid_hetero_alleles[C1]
                final.loc[(final['Index']==idx) & (final['Allele']=='A1'), 'Alias']=hybrid_hetero_alleles_alias[hybrid_hetero_alleles[C1]]
    return final


def CNV_def(log2r, cov, wgs_cov):
    
    '''
        calculate the CNV of each exon
        One copy gain = log2(3/2) = 0.57 (3 copies vs. 2 copies in reference)
        One-copy loss = log2(1/2) = -1
        Two-copy gain = log2(4/2) = 1
    '''

    status='normal'
    if cov/wgs_cov < 0.15:
        status='two-copy loss'
    elif log2r < -2.0 and cov <= 3:
        status='two-copy loss'  
    elif log2r < -2.0 and cov > 3:
        status='one-copy loss'
    elif -2.0 <= log2r < -0.6:
        status='one-copy loss'
    elif -0.6 <= log2r < 0.4:
        status='normal'
    elif 0.4 <= log2r < 0.8:
        status='one-copy gain'
    elif 0.8 <= log2r <= 1.2:
        status='two-copy gain'
    elif log2r < -2 or log2r > 1.5:
        status='abnormal'
    return status

def CNV_def_with_test(exon, perexon, log2r, cov, wgs_cov, alpha=0.05):
    '''
        calculate the CNV of each exon
        One copy gain = log2(3/2) = 0.57 (3 copies vs. 2 copies in reference)
        One-copy loss = log2(1/2) = -1
        Two-copy gain = log2(4/2) = 1
    '''
    target=[]
    non_target=[]
    for e in perexon:
        if e==exon:
            target.extend(perexon[e])
        else:
            non_target.extend(perexon[e])

    tstat, pval=stats.ttest_ind(np.asarray(target), np.asarray(non_target))
    status='normal'
    if cov/wgs_cov < 0.15 and pval < alpha:
        status='two-copy loss'
    elif log2r < -2.0 and cov <= 3 and pval < alpha:
        status='two-copy loss'
    elif log2r < -2.0 and cov > 3 and pval < alpha:
        status='one-copy loss' 
    elif -2.0 <= log2r < -0.6 and pval < alpha:
        status='one-copy loss'
    elif -0.6 <= log2r < 0.4 and pval < alpha:
        status='normal'
    elif 0.4 <= log2r < 0.8 and pval < alpha:
        status='one-copy gain'
    elif 0.8 <= log2r <= 1.2 and pval < alpha:
        status='two-copy gain'
    elif log2r < -2 or log2r > 1.5:
        status='abnormal'
    return status

     

def CNV_def_switch(val):
    switcher = {
        'normal': 0,
        'abnormal': 'abnormal',
        'two-copy loss': -2,
        'one-copy loss': -1,
        'one-copy gain': 1,
        'two-copy gain': 2
    }
    return switcher.get(val, "Invalid CNV status switch")

 
def RHD1136check(df, db, mode):
    '''
        remove the allele combiation with both allele containing the '1136C>T' var when the status is identifed as het 1136 C>T
    '''
    #tdb=db.transpose()
    print(df)
    print(mode)
    df.sort_values(by=['Likelihood', 'Index', 'Allele'], ascending=[0,1,1], inplace=True) 
    indices=df.Index.unique()
    if len(indices) == 1:
         print("[RHD1136] Only a pair of alleles present, not check RHD 1136")
         return df

    drop_indicies=[]


    for i in indices:
        b1=df.loc[(df['Index']==i) & (df['Allele']=='A1'),'Bloodtype']
        b2=df.loc[(df['Index']==i) & (df['Allele']=='A2'),'Bloodtype']
        b1s=db.loc[b1,'1136_C_T'].values[0] 
        b2s=db.loc[b2,'1136_C_T'].values[0]
        #print(mode, b1s, b2s)
        if mode == 'het RHD 1136_C_T':
            if b1s > 0 and b2s > 0:
                print('[Heterozygous RHD1136] drop alleles pairs with non-heterozygous 1136', b1s, b2s, '; index:', i)
                drop_indicies.append(i)
        elif mode == 'hom RHD 1136_C_T':
            if (b1s > 0 and b2s==0) or (b1s == 0 and b2s > 0):
                print('[Homozygous RHD1136] drop alleles pairs with non-homozygous 1136')
                drop_indicies.append(i)
    if len(drop_indicies) < len(df.Index.unique().tolist()):
        df=df[~df.Index.isin(drop_indicies)]
    else:
        print('[RHD1136 check]', 'all allele pairs did not meet filtering criteria, skipping...')
    return df



def RHD_rhesus_box(bam, cov=30, verbose=0):
    """
    Calculate rhesus box coverage
    """
    HG38_box={'upstreamRHbox':{'chr':'chr1', 'start':25258851, 'end':25268086, 'strand':'+'},
         'downstreamRHbox':{'chr':'chr1', 'start':25329004, 'end': 25338415, 'strand':'+'},
         }
    samfile = pysam.AlignmentFile(bam, "rb" )
    type=None
    for box in HG38_box:
        chr=HG38_box[box]['chr']
        start=HG38_box[box]['start']
        end=HG38_box[box]['end']
        basen=0
        totaln=0
        for column in samfile.pileup(chr, start-1, end):
            pos0based=column.pos
            pos1based=column.pos+1
            if start <= pos1based <= end:
                basen+=1
                totaln+=column.n
                #print ("\nstart %d, end %d, coverage at base %s = %s" % (start, end, pos1based, column.n))
        coverage=totaln/basen
        if verbose > 20: print (box, chr, start, end, basen, totaln, coverage)
        if cov*0.25 < coverage < cov*0.75:
            if type:
                type+=';' + box+ 'HeteDel'
            else:
                type=box+ 'HeteDel'
        elif coverage < cov*0.25:
            if type:
                type+=';' + box+ 'HomoDel'
            else:
                type=box+'HomoDel'
    samfile.close()
    return type

def RHCe_insertion(bam, cov, build='hg38', clippedN=5, verbose=0):
    '''
       Check whether the 109bp insertion occured in the 2nd intron
       primer:chr1:25406158-25406177
       primer:chr1:25405515-25405538
       full seg: chr1:25405538-25406158
       narrow range to ensure calling of heterozygous or homozygous CE
    '''
    samfile = pysam.AlignmentFile(bam, mode="rb") #, ignore_truncation=False)
    total_n=0
    sc_r=[]
    sc_n=0
    #chr='chr1'; s=25405500; e=25406100 ### based on primer blast
    #chr='chr1'; s=25405500; e=25405650 ### based on visual observation on WGS data
    chr='chr1'; s=None; e=None ### based on visual observation on WGS data for homozygous CE

    trimchr=False
    if not all('chr' in s for s in samfile.references):trimchr=True
    chr=chr.lstrip('chr') if trimchr else chr

    if build.upper()=='HG38':
        s=25405596; e=25405618
    elif build.upper()=='HG19':
        s=25732087; e=25732109
    for read in samfile.fetch(chr, s-1, e):
        total_n+=1
        if read.cigarstring and 'S' in read.cigarstring:
            for f in read.cigartuples:
                ### f[0] ==4: softcliiped
                if f[0]==4 and f[1] >= clippedN:
                    if verbose >9: print("[109bp INS]", read.reference_start, read.query_name, read.cigarstring)
                    if read.query_name not in sc_r: sc_r.append(read.query_name)
                    sc_n+=1
                    continue
    if verbose>9: print("[109bp INS] TotalN: %d\tclippedN: %d\tratio: %s" % ( total_n, sc_n, str(sc_n/cov) ))
    if total_n > 0:
        return sc_n, total_n, sc_n/float(total_n)
    return sc_n, 0, 0




def gene_exon_del(exon_cov, bam, fasta, coord, exon_cut=0.6, coverage=30, gene='RHD', seqtype='WGS'):
    """
        Calculate gene and exon coverage
        [Note] Cannot use exon average as indicator as it will falsely make hybrid to deletion
    """
    from math import log

    ### get a list of coverage values from perbaseCov
    chr=coord.chr[0]
    start=coord.start.min()
    end=coord.end.max()
    binCov, allbinmed=perbaseCov(bam,coverage, fasta, step=150, bin=300, prefix=None, chr=chr, start=start, end=end, coord=coord, minimum_base_quality = 15, minimum_mapq = 10, minimum_read_quality = 15, flanking=150, includeZero=True)
    allbinCov=[]
    exonbinCov=[]          
    for x in binCov:
        allbinCov.append(binCov[x]['bin_avg'])
        if binCov[x]['bin_exon'] != -1: exonbinCov.append(binCov[x]['bin_avg'])
    

    exon=exon_cov.iloc[exon_cov.index.get_level_values('gene')==gene]
    exon_sd =exon.avg_QC_readn.std()
    gene_avg=exon.avg_QC_readn.mean()
    #gene_avg=gene_cov.iloc[gene_cov.index.get_level_values('gene')==gene, gene_cov.columns.get_loc('avg_QC_readn')].values[0]
    typing=None
    log2r=round(log(gene_avg/coverage,2),7)
   
    homodel_cut=-3
    hetedel_cut=-0.6
    exon_nf=exon_cut
    if seqtype=='WES':
        #cov_list=exon.avg_QC_readn.tolist()
        cov_list=exonbinCov
    elif seqtype=='WGS':
        cov_list=allbinCov
    exon_n=len(cov_list)
    binCov_sd=np.std(cov_list)
    binCov_avg=np.mean(cov_list)
    true_n=exon_n
    #print('gene_avg:', gene_avg, 'exon_sd:', exon_sd)
    #print('bin_n:', exon_n,'bin_avg:', binCov_avg, 'bin_sd:', binCov_sd)
    ###if log2r < homodel_cut:
    for val in cov_list:
            #print('homodel', val)
            if val==0:
                continue
            else: 
               e_log2r=round(log(val/coverage,2),7)
               if e_log2r > homodel_cut: 
                    true_n-=1
    print('Homo-Deleted binN: {0}; ratio={1}'.format(true_n, true_n/exon_n))
    if ( true_n > exon_n * exon_nf ): 
             print('[Detected] Homozygous deletion')
             return 'RHD negative (Homozygous deletion)', 'both_allele_deleted'
    ###elif homodel_cut  <= log2r < hetedel_cut:
    true_n=exon_n
    for val in cov_list:
            #print('hetedel', val)
            if val==0: 
                continue
            else:
                e_log2r=round(log(val/coverage,2),7)
                if homodel_cut > e_log2r or e_log2r >= hetedel_cut: 
                    true_n-=1
                #print(val, e_log2r, true_n)
    print("Hete-Deleted binN: {0}; ratio={1}".format(true_n, true_n/exon_n))
    if ( true_n > exon_n * exon_nf ):
            #print(gene, coverage, gene_avg/float(coverage), log2r, true_n, exon_n)
            print('[Detected] Heterozygous deletion')
            typing="RHD negative (Heterozygous deletion) or low coverage"
    
    abnormal_exon=zip(exon[(exon['avg_QC_readn'] < gene_avg-exon_sd) | (exon['avg_QC_readn'] > gene_avg+exon_sd) ].index, exon[(exon['avg_QC_readn'] < gene_avg-exon_sd) | (exon['avg_QC_readn'] > gene_avg+exon_sd)]['avg_QC_readn'])
    abnormal_exon=[[acc,gene,exon,cov] for (acc, gene, exon), cov in abnormal_exon]
    ab_exon=None
    for e in abnormal_exon:
        if ab_exon:
            ab_exon+='; '+ ':'.join(map(str,e))
        else:
            ab_exon=':'.join(map(str,e))
    return typing, ab_exon

def coverage(df=None, file=None):
    '''
        get average coverage depth for each exon
    '''
    stats = {
        'gpos':['min', 'max'],
        'cpos':['min', 'max'],
        'QC_n':['min', 'max', 'mean', 'median']
    }
    indata=None
    if file is not None:
        indata=pd.read_table(file, sep='\t', header=0)
    elif df is not None:
        indata=df
    else:
        sys.exit('No variant information found ... quit')

    colnames=['genome_start','genome_end','cds_start','cds_end','min_QC_readn','max_QC_readn','avg_QC_readn','med_QC_readn']

    exon_grouped = indata.groupby(['acc','gene','exon']).agg(stats)
    gene_grouped = indata.groupby(['acc','gene']).agg(stats)
                        
    ### rename columns
    exon_grouped.columns = exon_grouped.columns.droplevel(level=0)
    gene_grouped.columns = gene_grouped.columns.droplevel(level=0)
    exon_grouped.columns = colnames
    gene_grouped.columns =colnames
    #grouped.columns = ["_".join(x) for x in grouped.columns.ravel()]
    return exon_grouped, gene_grouped

    


def genome_cov(bam,chr=None, mode='wgs', gbuild='hg38'):
    ''' 
        estimate genome coverage based on chr1 
    '''
    print("[Estimate coverage] based on chr1")
    annotation={'hg38':{#'file': "hg38.refFlat.nonoverlapping.txt.gz",
                        'file': 'hg38.exon.gz',
                        'RHCE':{'start':25362249, 'end':25420914 },
                        'RHD':{'start':25272393, 'end':25330445 }
                       },
                'hg19':{#'file': "hg19.refFlat.nonoverlapping.txt.gz",
                        'file': 'hg19.exon.gz',
                        'RHCE':{'start':25688740, 'end':25747363 },
                        'RHD':{'start':25598981, 'end':25656936}
                       }
               }
    
    ann = gzip.open(os.path.join(package_directory, 'database', annotation[gbuild.lower()]['file']),'rt')
    exons={}
    for line in ann:
        if line.startswith('#'): continue
        flds=line.strip().split('\t')
        flds[0] = re.sub("\.fa","", flds[0])
        ###    chr,     start,   end,    gene
        if chr is not None:
            if flds[0].upper() != chr.upper() and 'CHR'+flds[0].upper() != chr.upper(): continue
        if flds[0] not in exons: exons[flds[0].lower()]=[]
        s=int(flds[2]) if int(flds[1]) > int(flds[2]) else  int(flds[1])
        e=int(flds[1]) if int(flds[1]) > int(flds[2]) else  int(flds[2])
        if mode=='wes':
            if ( s < annotation[gbuild.lower()]['RHCE']['start'] and e < annotation[gbuild.lower()]['RHCE']['end'] ) or ( s > annotation[gbuild.lower()]['RHCE']['start'] and e > annotation[gbuild.lower()]['RHCE']['end'] ):
                continue

        if (s,e) not in exons[flds[0]]:
            exons[flds[0]].append((int(s),int(e)))


    bamfile = pysam.AlignmentFile(bam, "rb")

    chrs = bamfile.header['SQ']
    contigs = {}
    for c in chrs:
        if chr is not None:
            if c['SN'].upper() != chr.upper() and 'CHR'+c['SN'].upper() != chr.upper(): continue
        if 'GL' in c['SN'].upper() and 'GL' not in chr.upper(): continue        
        contigs[c['SN']]={}
        contigs[c['SN']]['Coverage']=0
        contigs[c['SN']]['Length']=c['LN']

    totalSum=0
    totalLen=0

    for contig in contigs:
        covArray = []
        posArray = set()
        if contig in exons:
            ctg_exon=exons[contig]
        elif 'chr' + contig in exons:
            ctg_exon=exons['chr'+contig]
        ex_n=0
        for ex in ctg_exon:
            ex_n+=1
            #if ex_n % 1000 == 0: print ex_n, " reads processed"
            for pile in bamfile.pileup(contig, ex[0], ex[1], truncate=True):
                if pile.reference_pos not in posArray:
                    #print pile.reference_pos, pile.nsegments
                    posArray.add(pile.reference_pos)
                    covArray.append(pile.nsegments)
                #else:
                #    print pile.reference_pos, 'dup'
        try:
            contigs[contig]["Coverage"]=round(sum(covArray) / len(covArray),2)
            contigs[contig]["CoverageSum"]=sum(covArray)
            contigs[contig]["CoverageLen"]=len(covArray)
            totalSum += sum(covArray)
            totalLen += len(covArray)
        except ZeroDivisionError: # should only occur if 0 coverage recorded
            contigs[contig]["Coverage"] = 0
            contigs[contig]["CoverageSum"]=0
            contigs[contig]["CoverageLen"]=0
    
    #print contigs    
    mean_coverage=round((totalSum/totalLen),2)
    median_coverage=median(covArray)
    return mean_coverage, median_coverage

def median(lst):
    n = len(lst)
    if n < 1:
            return None
    if n % 2 == 1:
            return sorted(lst)[n//2]
    else:
            return sum(sorted(lst)[n//2-1:n//2+1])/2.0

def RHD_RHCE_var_compare(RHDvar,RHCEvar, RHcov, altn=7):
    '''
        Compare all the var to determine whether the gene conversion (hybrid) indeed occur
    '''
    for idx, row in RHDvar.iterrows():
        RHDref=row['REF']
        RHDexon=row['exon']
        RHDrefCount=int(row[RHDref])
        RHDalt=row['ALT']
        RHDaltCount=0 if RHDalt in ['.','-'] else int(row[RHDalt])

        RHCEref=RHCEvar.loc[RHCEvar['cpos']==row['cpos'],'REF'].values[0]
        RHCEexon=RHCEvar.loc[RHCEvar['cpos']==row['cpos'],'exon'].values[0]
        RHCErefCount=RHCEvar.loc[RHCEvar['cpos']==row['cpos'], RHCEref].values[0]
        RHCEalt=RHCEvar.loc[RHCEvar['cpos']==row['cpos'],'ALT'].values[0]
        RHCEaltCount=0 if RHCEalt in ['.', '-'] else int(RHCEvar.loc[RHCEvar['cpos']==row['cpos'], RHCEalt].values[0])

        ### adjust call for heterozygous
        RHDaltCount=0  if RHDaltCount < altn else RHDaltCount
        RHCEaltCount=0 if RHCEaltCount < altn else RHCEaltCount
         
        if RHDref == RHCEref:
            ### skip the var that cannot distinguish 
            continue 

        if RHDaltCount != 0 or RHCEaltCount !=0:
            print ("['RHD'] ", RHDexon, row['cpos'], RHDref, RHDrefCount, RHDalt, RHDaltCount)
            print ("['RHCE']", RHCEexon, row['cpos'], RHCEref, RHCErefCount, RHCEalt, RHCEaltCount)
     



def var_report(BG, db, sdb, gene):
    '''
        report the consistent var and inconsistent var for the predicted allele combination
    '''
    BGls=BG['Bloodtype'].tolist()


    tdb=db.transpose()
    sdb2=sdb.copy()
    all_var=[]
    all_aa=[]
    report={}
    for bg in BGls:
        if bg not in tdb.columns:
            ### skip bloodgroup that are not repsent in the database, e.g. deletion ones
            continue
        var=tdb.loc[tdb[bg]>0, bg]
        match_var=[]
        nonmatch_var=[]
        match_aa=[]
        nonmatch_aa=[]
        if var.shape[0] > 0:
            match_var=[]
            nonmatch_var=[]
            for i in var.index:
                cpos, ref, alt=i.split('_')
                ovar=cpos+ref+'/'+alt
                aaovar=ovar
                match=False
                refc=0
                altc=0
                refa=''
                alta=''
                for idx, row in sdb.iterrows():
                    #print row['cpos'], row['REF'], row['ALT'], row['exon'], cpos, ref, alt
                    if int(row['cpos'])==int(cpos) and row['REF']==ref and row['ALT']==alt:
                        if ref in ['A', 'T', 'C', 'G']:
                            refc=row[ref]
                        elif ref =='-':
                            refc=row['INS_n']
                        if alt in ['A', 'T', 'C', 'G']:
                            altc=row[alt]
                        elif alt == '-':
                            altc=row['DEL_n']
                        match=True
                        refa=row['REFaa']
                        alta=row['ALTaa']
                        if idx in sdb2.index:
                            sdb2.drop([idx], inplace=True)
                        break
                if match:
                    ovar=ovar+'('+str(refc)+'/'+str(altc)+')'
                    aaovar=aaovar+'('+str(refa)+'/'+str(alta)+')'
                    match_var.append(ovar)
                    match_aa.append(aaovar)
                    if ovar not in all_var: 
                        all_var.append(ovar)
                        all_aa.append(aaovar)
                else:
                    if refc !=0 and altc !=0:
                        ovar=ovar+'('+str(refc)+'/'+str(altc)+')'
                    if refa !="" and alta !="":
                        aaovar=aaovar+'('+str(refa)+'/'+str(alta)+')'
                    else:
                        aaovar=aaovar+'(./.)'
                    nonmatch_var.append(ovar)
                    nonmatch_aa.append(aaovar)
        report[bg]={}
        report[bg]['match']=';'.join(match_var)
        report[bg]['non-match']=';'.join(nonmatch_var)
        report[bg]['match_aa']=';'.join(match_aa)
        report[bg]['non-match_aa']=';'.join(nonmatch_aa)
    sample_specific=[]
    sample_specific_aa=[]
    exon8_dels=False
    prevalt=""
    prevcpos=-1
    del_start=-1
    del_end=-1
    del_ref=''
    for idx, row in sdb2.iterrows():
        refc=0
        altc=0
        ovar=str(row['cpos'])+row['REF']+'/'+row['ALT']
        aaovar=ovar
        refa=row['REFaa']
        alta=row['ALTaa']
        if row['REF'] in ['A', 'T', 'C', 'G']:
            refc=row[row['REF']]
        elif row['REF'] == '-':
            refc=row['INS_n']
        if row['ALT'] in ['A', 'T', 'C', 'G']:
            altc=row[row['ALT']]
        elif row['ALT'] == '-':
            altc=row['DEL_n']
            del_ref+=row['REF']
            if del_start==-1:
                del_start=row['cpos']
            
            if (prevalt=='-') and (row['ALT']=='-') and (int(row['cpos'])-int(prevcpos)==1):
                prevalt=row['ALT']
                prevcpos=row['cpos']
                del_end=row['cpos']
                if idx != sdb2.shape[0]-1:
                    continue        
            elif (prevalt=='-') and ((row['ALT'] != '-') or (int(row['cpos'])-int(prevcpos) > 1)):
                if del_start == del_end or del_end==-1:
                    delovar=str(del_start) + 'del'
                else:    
                    delovar=str(del_start) + '-' + str(del_end) + 'del'
                sample_specific.append(delovar)
                if delovar not in all_var: all_var.append(delovar)
                del_ref=''
                del_start=-1
                del_end=-1
        ### check whether the last entries are deletion block
        if (del_start != -1) and (del_end != -1) and (idx == sdb2.shape[0]-1):
            if del_start == del_end:
                delovar=str(del_start) + 'del'
            else:    
                delovar=str(del_start) + '-' + str(del_end) + 'del'
            sample_specific.append(delovar)
            if delovar not in all_var: all_var.append(delovar)

        prevalt=row['ALT']
        prevcpos=row['cpos']
        ovar=ovar+'('+str(refc)+'/'+str(altc)+')' 
        aaovar=aaovar+'('+str(refa)+'/'+str(alta)+')' 
        if row['exon']==8 and (gene=='RHD' or gene=='RHCE'):
            # skip exon 8 as it's nearly identical between RHD/RHCE
            exon8_del=True
            continue
        elif row['ALT']=='-':
            # skip the last pos if it's del
            continue
        if row['ALT'] != '-':
            sample_specific.append(ovar)
            sample_specific_aa.append(aaovar)
            if ovar not in all_var: 
                all_var.append(ovar)
                all_aa.append(aaovar)
    #if exon8_dels=True:
    #    sample_specific.append('exon8_dels_not_shown')
       
    report['sample_specific']=sample_specific
    report['sample_specific_aa']=sample_specific_aa
    report['all_variation']=all_var
    report['all_variation_aa']=all_aa
    return report


def pos_coverage(bam, chr, pos, strand, ref, alt, fa, minimum_base_quality = 15, minimum_mapq = 10, minimum_read_quality = 15, ALT_n_LB=None):
    '''
        perform an ad-hoc variant check for intron vars
    '''
    samfile = pysam.AlignmentFile(bam, mode="rb")
    reffasta=pysam.FastaFile(fa)
    trimchr=False
    if not all('chr' in s for s in samfile.references):trimchr=True
    tarchr=chr.lstrip('chr') if trimchr else chr

 
    for column in samfile.pileup(tarchr, pos-1, pos):
        pos0based=column.pos
        pos1based=column.pos+1
        if int(pos) != int(pos1based): continue
        good_count    ={'A':0, 'T':0, 'C':0, 'G':0, 'N':0}
        ins_count     ={'A':0, 'T':0, 'C':0, 'G':0, 'N':0}
        del_count     ={'A':0, 'T':0, 'C':0, 'G':0, 'N':0}
        lowmapq_count ={'A':0, 'T':0, 'C':0, 'G':0, 'N':0}
        nomapq_count  ={'A':0, 'T':0, 'C':0, 'G':0, 'N':0}
        lowreadq_count={'A':0, 'T':0, 'C':0, 'G':0, 'N':0}
        lowbaseq_count={'A':0, 'T':0, 'C':0, 'G':0, 'N':0} 
   
        refbase =reffasta.fetch(chr, pos0based, pos0based+1)
        ### complement base if minus strand
        refbase=complement(refbase) if strand =='-' else refbase
     
        ### process each read in each exon to get the count of bases
        good_n=0
        del_n =0
        ins_n =0
        ins_cpos=pos+1 if strand =='+' else pos-1
        ins_gpos=pos1based+1
        for pileupread in column.pileups:
            read     = pileupread.alignment
            size     = pileupread.indel
            qpos     = pileupread.query_position
            qname    = pileupread.alignment.query_name

            if pileupread.is_del: #or size < 0:
                readbase=''
            else:
                readbase = read.query_sequence[qpos]

            ### complement base if on the minus stand
            readbase=complement(readbase) if strand =='-' else readbase
            readbase='N' if readbase=='.' else readbase

            # skip reads with indels
            if pileupread.is_del:# or pileupread.indel < 0:
                #print skip_message.format(pos1based, qname, 'is del', readbase, pileupread.indel)
                del_count[refbase]+=1
                del_n+=1
                continue

            # skip reads with mapq below threshold
            if pileupread.alignment.mapping_quality < minimum_mapq:
                #print skip_message.format(position, qname, 'low mapq', pileupread.alignment.mapping_quality)
                lowmapq_count[readbase]+=1
                continue
            # skip reads with no mapq specified
            elif read.mapping_quality == 255:
                #print skip_message.format(position, qname, 'no mapq')
                nomapq_count[readbase]+=1
                continue
            # skip mean qscore of the read below threshold
            elif mean(read.query_qualities) < minimum_read_quality:
                #print skip_message.format(position, qname, 'low read quality', mean(read.query_qualities))
                lowreadq_count[readbase]+=1
                continue
            else:
                # check for insertion
                if pileupread.is_refskip: #or pileupread.indel > 0:
                    if read.query_qualities[qpos+1] >= minimum_base_quality:
                        insbase=read.query_sequence[qpos+1]
                        insbase=complement(insbase) if strand =='-' else insbase
                        #print skip_message.format(pos1based, qname, 'is ins', readbase)
                        ins_count[insbase]+=1
                        ins_n+=1
                # skip reads with a base quality below threshold
                if read.query_qualities[qpos] < minimum_base_quality:
                    #print skip_message.format(position, qname, 'low base quality', read.query_qualities[pileupread.query_position])
                    lowbaseq_count[readbase]+=1
                    continue

            ### print details of each passed read
            #print pass_message.format(position, qname, read.query_sequence[qpos])
                
            good_n+=1
            good_count[readbase]+=1
        if alt in ['A','T','C','G','N']:
            if refbase == ref:
                #print("GPOS:", pos, "REF:", ref, ' ', good_count[ref], " ALT:", alt, ' ', good_count[alt])
                return good_count[ref], good_count[alt]
            else: print("GPOS:", gpos, ", REF:", ref, " is not equal to HG ref:", refbase)
        elif alt=='del':
            if refbase in ref:
                return good_count[refbase], del_n
            else: print("GPOS:", pos, "REF:", ref, " is not equal to HG ref:", refbase)
        elif alt in ['ins', 'dup']:
            if refbase in ref:
                return good_count[refbase], ins_n
            else: print("GPOS:", pos, "REF:", ref, " is not equal to HG ref:", refbase)


def table_coverage(sdbC, pos, ref, alt):
    '''
        extract coverage for a position directly from a pd table
    '''
    refC=sdbC.loc[sdbC['cpos']==pos ,ref].values[0]
    altC=sdbC.loc[sdbC['cpos']==pos ,alt].values[0]
    refa=sdbC.loc[sdbC['cpos']==pos ,'REFaa'].values[0]
    alta=sdbC.loc[sdbC['cpos']==pos ,'ALTaa'].values[0]
    return refC, altC, refa, alta

def companyGenecheck(sdbC, gene, bg_nogene, pos, ref, alt, match_types, unmatch_types, match_var, match_aa, nonmatch_var, nonmatch_aa, all_var, all_aa):
    '''
       check var from company gene tables
    '''
    all_types=[]
    all_types.extend(match_types)
    all_types.extend(unmatch_types)
    if bg_nogene in all_types:
        if gene=='RHCE':
            ovar='CEc.'+ str(pos) + ref +'/'+ alt
        elif gene=='RHD':
            ovar='Dc.'+ str(pos) + ref +'/'+ alt
        aaovar=ovar
        if pos in sdbC.cpos.values:
            refC, altC, refa, alta=table_coverage(sdbC, pos, ref, alt)
            ovar=ovar+'('+str(refC)+'/'+str(altC)+')'
            aaovar=aaovar+'('+str(refa)+'/'+str(alta)+')'
            if bg_nogene in match_types:
                match_var.append(ovar)
                match_aa.append(aaovar)
            else:
                nonmatch_var.append(ovar)
                nonmatch_aa.append(aaovar)
            if ovar not in all_var:
                all_var.append(ovar)
                all_aa.append(aaovar)
        else:
            ovar=ovar+'(./.)'
            aaovar=aaovar+'(./.)'
            if bg_nogene in match_types:
                nonmatch_var.append(ovar)
                nonmatch_aa.append(aaovar)
    return match_var, match_aa, nonmatch_var, nonmatch_aa
            
def var_compare(BG, db, sdb, sdbC, adb, gene, bam, fasta, altn_cut=3):
    '''
        report the consistent var and inconsistent var for the predicted allele combination
    '''
    BGls=BG['Bloodtype'].tolist()

    tdb=db
    sdb2=sdb.copy()
    all_var=[]
    all_aa=[]
    report={}
   

    nt_pattern=re.compile("([c,C].\s?[0-9]+\s?[A,T,C,G]\s?>\s?[A,T,C,G])")
    deletion_pattern=re.compile("([c,C].[0-9]*_?[0-9]+\s?(del|DEL)\s?[A,T,C,G]*)")
    insertion_pattern=re.compile("([c,C].[0-9]*_?[0-9]+\s?(ins|INS|dup|DUP)\s?[A,T,C,G]*)")

    for bg in BGls:
        bg_nogene=re.sub(".*@","",bg)
        bg_nogene_ws=re.sub("_"," ",bg_nogene)
        if bg_nogene_ws not in tdb['Allele name'].tolist():
            ### skip bloodgroup that are not repsent in the database, e.g. deletion ones
            continue
        vars=tdb.loc[tdb['Allele name']==bg_nogene_ws, 'Nucleotide'].values[0]
        nts=re.split('; |,|;', vars.strip())
        match_var=[]
        nonmatch_var=[]
        match_aa=[]
        nonmatch_aa=[]
        for nt in nts:
            mut=nt.strip()
            if mut.startswith(('-', '~', '*', '+')): continue ### unknown variant format... skip
            if not mut.startswith('c.'): continue ### non-coding vars ... skip
            if nt_pattern.match(mut.upper()):
                mut=mut.upper() 
                b0, t, cpos, ref, alt, b1=re.split(r'([c,C].)\s?(\d+)\s?([A-Z,a-z])\s?>\s?([A-Z,a-z])\s?', mut)
                ovar=cpos+ref+'/'+alt
                aaovar=ovar 
                match=False
                refc=0
                altc=0
                refa='.'
                alta='.'
                for idx, row in sdb.iterrows():
                    #print(row['cpos'], row['REF'], row['ALT'], row['exon'], cpos, ref, alt)
                    if int(row['cpos'])==int(cpos) and row['REF']==ref and row['ALT']==alt:
                        if ref in ['A', 'T', 'C', 'G']:
                            refc=row[ref]
                        elif ref =='-':
                            refc=row['INS_n']
                        if alt in ['A', 'T', 'C', 'G']:
                            altc=row[alt]
                        elif alt == '-':
                            altc=row['DEL_n']
                        match=True
                        refa=row['REFaa']
                        alta=row['ALTaa']
                        if idx in sdb2.index:
                            sdb2.drop([idx], inplace=True)
                        break
                aaovar=aaovar+'('+str(refa)+'/'+str(alta)+')'
                if match:
                    ovar=ovar+'('+str(refc)+'/'+str(altc)+')'
                    match_var.append(ovar)
                    match_aa.append(aaovar)
                    if ovar not in all_var:
                        all_var.append(ovar)
                        all_aa.append(aaovar)
                else:
                    if ref in ['A', 'T', 'C', 'G']:
                        refc=adb.loc[adb['cpos']==int(cpos),ref].values[0]
                    elif ref == '-':
                        refc=adb.loc[adb['cpos']==int(cpos),'INS_n'].values[0]
                    if alt in ['A', 'T', 'C', 'G']:
                        altc=adb.loc[adb['cpos']==int(cpos),alt].values[0]
                    elif alt == '-':
                        altc=adb.loc[adb['cpos']==int(cpos),'DEL_n'].values[0]
                    ovar=ovar+'('+str(refc)+'/'+str(altc)+')'
                    nonmatch_var.append(ovar)
                    nonmatch_aa.append(aaovar)
            else:
               if deletion_pattern.match(mut.upper()):
                   print("[coverage check] del skipped", mut)
               elif insertion_pattern.match(mut.upper()):
                   print('[coverage check] insertion', mut) 
                   #c.216_217dupCA
                   b0, t, start, end, ins, bases, b1=re.split(r'([c,C].)\s?(\d+)\s?_\s?(\d+)\s?([I,D][N,U][S,P])\s?([A-Z,a-z]+)\s?', mut.upper())
                   print(start, end, bases)
                   gstart=adb.loc[adb['cpos']==int(start), 'gpos'].values[0] 
                   gend  =adb.loc[adb['cpos']==int(end)  , 'gpos'].values[0]
                   ovar=mut
                   aaovar=ovar
                   chr=adb.chr.unique()
                   strand=adb.strand.unique()
                   if len(chr) == 1:
                       all_altc=[]
                       all_refc=[]
                       for gpos in range(gstart,gend):
                           refc, altc=pos_coverage(bam, chr[0], gpos, strand[0], bases, 'ins', fasta, minimum_base_quality = 15, minimum_mapq = 10, minimum_read_quality = 15, ALT_n_LB=None)
                           all_altc.append(altc)
                           all_refc.append(refc)

                       ovar=ovar+'('+','.join(map(str,all_refc))+'/'+ ','.join(map(str,all_altc)) +')'
                       aaovar=aaovar+'(./.)'
                       if any(c >= altn_cut for c in all_altc):
                           match_var.append(ovar)
                           match_aa.append(aaovar)
                           if ovar not in all_var:
                                all_var.append(ovar)
                                all_aa.append(aaovar)
                       else:
                           nonmatch_var.append(ovar)
                           nonmatch_aa.append(aaovar)
                   else:
                       print('[Warning] cannot check insertion coverage due to multiple chr specified')
                       continue 
               else:
                   print('[coverage check] non-indel', mut)
                   b0, t, cpos, shift_direction, shift_basen, ref, alt, b1=re.split(r'([c,C].)\s?(\d+)\s?([+,-])\s?(\d+)\s?([A-Z,a-z]+)\s?>\s?([A-Z,a-z]+)\s?', mut)
                   gpos=adb.loc[adb['cpos']==int(cpos), 'gpos'].values[0]
                   ovar='g.'+str(gpos)+ref+'/'+alt
                   aaovar=ovar
                   #print('\t'.join(map(str,[cpos,gpos,shift_direction,shift_basen,ref, alt])))
                   chr=adb.chr.unique()
                   strand=adb.strand.unique()
                   
                   #print('strand:', strand[0])

                   if strand[0]=='+':
                       if shift_direction == '+':
                           gpos=int(gpos)+int(shift_basen)
                       elif shift_direction == '-':
                           gpos=int(gpos)-int(shift_basen)
                   elif strand[0]=='-':
                       if shift_direction == '+':
                           gpos=int(gpos)-int(shift_basen) 
                       elif shift_direction == '-':
                           gpos=int(gpos)+int(shift_basen)

 
                   if len(chr) == 1:
                       refc, altc=pos_coverage(bam, chr[0], gpos, strand[0], ref, alt, fasta, minimum_base_quality = 15, minimum_mapq = 10, minimum_read_quality = 15, ALT_n_LB=None)
                       ovar=ovar+'('+str(refc)+'/'+str(altc)+')'
                       aaovar=aaovar+'(./.)'
                       if altc >= altn_cut:
                           match_var.append(ovar)
                           match_aa.append(aaovar)
                           if ovar not in all_var:
                                all_var.append(ovar)
                                all_aa.append(aaovar)
                       else:
                           nonmatch_var.append(ovar)
                           nonmatch_aa.append(aaovar) 
                   else:
                       print('[Warning] cannot check intronic var coverage due to multiple chr specified')
                       continue   
        ### extra check
        match_var, match_aa, nomatch_var, nomatch_aa=companyGenecheck(sdbC, 'RHCE', bg_nogene, 733, 'C', 'G', ['RHD*01N.06','RHD*03N.01'], ['RHD*03.04_RHD*DIII.04','RHD*03.01_RHD*DIIIa'], \
                                                                      match_var, match_aa, nonmatch_var, nonmatch_aa, all_var, all_aa)
        match_var, match_aa, nomatch_var, nomatch_aa=companyGenecheck(sdbC, 'RHCE', bg_nogene, 1006, 'G', 'T', ['RHD*01N.06','RHD*03N.01'], ['RHD*03.04_RHD*DIII.04','RHD*03.01_RHD*DIIIa'], \
                                                                      match_var, match_aa, nonmatch_var, nonmatch_aa, all_var, all_aa)

        report[bg]={}
        report[bg]['match']=';'.join(match_var)
        report[bg]['non-match']=';'.join(nonmatch_var)
        report[bg]['match_aa']=';'.join(match_aa)
        report[bg]['non-match_aa']=';'.join(nonmatch_aa)
   

    sample_specific=[]
    sample_specific_aa=[]
    exon8_dels=False
    prevalt=""
    prevcpos=-1
    del_start=-1
    del_end=-1
    del_ref=''
    #print(sdb2)
    for idx, row in sdb2.iterrows():
        refc=0
        altc=0
        ovar=str(row['cpos'])+row['REF']+'/'+row['ALT']
        aaovar=ovar
        refa=row['REFaa']
        alta=row['ALTaa']
        if row['REF'] in ['A', 'T', 'C', 'G']:
            refc=row[row['REF']]
        elif row['REF'] == '-':
            refc=row['INS_n']
        if row['ALT'] in ['A', 'T', 'C', 'G']:
            altc=row[row['ALT']]
        elif row['ALT'] == '-':
            altc=row['DEL_n']
            del_ref+=row['REF']
            if del_start==-1:
                del_start=row['cpos']

            if (prevalt=='-') and (row['ALT']=='-') and (int(row['cpos'])-int(prevcpos)==1):
                prevalt=row['ALT']
                prevcpos=row['cpos']
                del_end=row['cpos']
                if idx != sdb2.shape[0]-1:
                    continue
            elif (prevalt=='-') and ((row['ALT'] != '-') or (int(row['cpos'])-int(prevcpos) > 1)):
                if del_start == del_end or del_end==-1:
                    delovar=str(del_start) + 'del'
                else:
                    delovar=str(del_start) + '-' + str(del_end) + 'del'
                sample_specific.append(delovar)
                if delovar not in all_var: all_var.append(delovar)
                del_ref=''
                del_start=-1
                del_end=-1
        ### check whether the last entries are deletion block
        if (del_start != -1) and (del_end != -1) and (idx == sdb2.shape[0]-1):
            if del_start == del_end:
                delovar=str(del_start) + 'del'
            else:
                delovar=str(del_start) + '-' + str(del_end) + 'del'
            sample_specific.append(delovar)
            if delovar not in all_var: all_var.append(delovar)

        prevalt=row['ALT']
        prevcpos=row['cpos']
        ovar=ovar+'('+str(refc)+'/'+str(altc)+')'
        aaovar=aaovar+'('+str(refa)+'/'+str(alta)+')'
        if row['exon']==8 and (gene=='RHD' or gene=='RHCE') and row['ALT'] == '-':
            # skip exon 8 as it's nearly identical between RHD/RHCE
            exon8_del=True
            continue
        elif row['ALT']=='-':
            # skip the last pos if it's del
            continue

        if row['ALT'] != '-':
            sample_specific.append(ovar)
            sample_specific_aa.append(aaovar)
            if ovar not in all_var:
                all_var.append(ovar)
                all_aa.append(aaovar)
    #if exon8_dels=True:
    #    sample_specific.append('exon8_dels_not_shown')


    report['sample_specific']=sample_specific
    report['sample_specific_aa']=sample_specific_aa
    report['all_variation']=all_var
    report['all_variation_aa']=all_aa


    return report
            

def max_match(final):
    xmDBn={}
    maxDBn=-1
    mmDBn={}
    minDBn=-1
    for idx in final.Index.unique():
        max_not_matchDBn=final.loc[final['Index']==idx,'not_matchDBn'].sum()
        max_matchDBn    =final.loc[final['Index']==idx,'matchDBn'].sum()
        if max_not_matchDBn not in mmDBn: mmDBn[max_not_matchDBn]=[]
        if max_matchDBn     not in xmDBn: xmDBn[max_matchDBn]    =[]
        mmDBn[max_not_matchDBn].append(idx)
        xmDBn[max_matchDBn].append(idx)
        if minDBn==-1: minDBn=max_not_matchDBn
        if maxDBn==-1: maxDBn=max_matchDBn
        if max_not_matchDBn < minDBn: minDBn=max_not_matchDBn
        if max_matchDBn > maxDBn: maxDBn=max_matchDBn
    return xmDBn, maxDBn, mmDBn, minDBn 

