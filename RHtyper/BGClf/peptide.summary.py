#!/usr/bin/env python
"""
    modules to produce epitope matrix 

"""
import pandas as pd
import numpy as np
import re, collections, pysam, sys, os, argparse, glob

### python2
try:
    import commands
except ImportError:
    from subprocess import Popen, PIPE
### python2
try:
    from StringIO import StringIO ## for Python 2
except ImportError:
    from io import StringIO ## for Python 3

package_directory = os.path.dirname(os.path.abspath(__file__))


def process_affinity(file_list, rank_cut=2):
 
   out=pd.DataFrame()
 
   for f in file_list:
      df=pd.read_csv(f, header=None, sep='\s+', names=['Pos','MHC','Peptide','Of','Core','Core_Rel','Identity','Score_EL','Rank_EL_pct','Exp_Bind','BindLevel'])
      ori_size=len(df.index)
      df=df.loc[df.Rank_EL_pct <= rank_cut,]
      if len(out.index)==0:
          out=df
      else:
          out=pd.concat([out,df])
      print(f, ' total:', ori_size, ' filtered:', len(df.index), ' combined:', len(out.index))
   return out

def fa2pd(f):
   out=list()
   with pysam.FastxFile(f) as fh:
      for rec in fh:
         out.append([rec.name, rec.sequence])
   return pd.DataFrame(out, columns=['Identity', 'peptide'])
   
def chech_mut_pos(df,size=15):
   pd.set_option('display.max_columns', 1000)
   df.mutAA_in_core=np.nan
   for idx, row in df.iterrows():
      mutpos_in_full_peptide=int(row['ALT_AApos'])-int(row['AAstart'])
      peptide=np.nan   
      mut_in_peptide=np.nan      
      core_peptide=np.nan
      mut_in_core=np.nan

      if not pd.isnull(row['Pos']):
          peptide=row['Full_peptide'][int(row['Pos'])-1:int(row['Pos'])+size]
          mutpos_in_peptide=mutpos_in_full_peptide-int(row['Pos'])+1
          mut_in_peptide=peptide[mutpos_in_peptide]
 
          core_peptide=peptide[int(row['Of']):int(row['Of'])+9]
          mutpos_in_core=mutpos_in_peptide-int(row['Of'])
          try:
             if mutpos_in_core >=0:
                 mut_in_core=core_peptide[mutpos_in_core]
          except:
             mut_in_core=np.nan
             #print("Mut not found in core peptide")
          if not pd.isnull(mut_in_core) and ((mut_in_core not in mut_in_peptide.split(',')) or (mut_in_peptide not in row['ALTaa'].split(','))):
              print(' '.join(['idx','FullPeptide','Full_peptide_len','ALT_AApos','ALTaa','MutPos_in_Full_peptide','Peptide','Peptide2','Mut_in_peptide','Core','Core_peptide','Mutpos_in_core','Mut_in_core']))
              print(idx, row['Full_peptide'], row['Full_peptide_len'], row['ALT_AApos'], row['ALTaa'], row['Full_peptide'][mutpos_in_full_peptide], row['Peptide'], peptide, mut_in_peptide, row['Core'], core_peptide, mutpos_in_core, mut_in_core)
              sys.exit('[ERROR] aa not matched')
      if not pd.isnull(mut_in_core):
          df.loc[df.index==idx, 'mutAA_in_core']=mut_in_core
      if int(idx)%5000==0:
          print(idx, round(int(idx)*100/len(df.index),2), "%")
      #print(idx, row['Full_peptide'], row['Full_peptide_len'], row['ALT_AApos'], row['ALTaa'], row['Full_peptide'][mutpos_in_full_peptide], row['Peptide'], peptide, mut_in_peptide, row['Core'], core_peptide, mut_in_core) 
   return df

def main():
   parser = argparse.ArgumentParser(description='''Summarize RH epitopes. 
      This script is used for the data generated from peptide.py dry-run mode to interrogate all allele\'s affinity\n. 
      During dry-run, commands for generating the affinity in a nubmer of alleles will be generated that can be used for batch submitting the task for clsuter.''' )
   parser.add_argument('-pre' , '--prefix', help='Output prefix', default='Summary')
   parser.add_argument('-aff_RHD' , '--aff_RHD_path', help='Path to RHD affinity file.', default='RHD')
   parser.add_argument('-aff_RHCE' , '--aff_RHCE_path', help='Path to RHCE affinity file', default='RHCE')
   parser.add_argument('-aff_pattern' , '--affinity_file_pattern', help='pattern of affinity file', default='epitope.affinity')
   parser.add_argument('-info_RHD' , '--peptide_info_RHD', help='File with RHD peptide info generated from peptide.py', default='RHD.epitope.info.txt')
   parser.add_argument('-info_RHCE' , '--peptide_info_RHCE', help='File with RHCE peptide info generated from peptide.py', default='RHCE.epitope.info.txt')
   parser.add_argument('-fa_RHD' , '--peptide_fa_RHD', help='Fasta file with RHD peptide generated from peptide.py', default='RHD.epitope.fa')
   parser.add_argument('-fa_RHCE' , '--peptide_fa_RHCE', help='Fasta file with RHCE peptide generated from peptide.py', default='RHCE.epitope.fa')
   parser.add_argument('-mode' , '--mode', help='Running mode: all, agg_affinity, epitope_merge, summary', default='summary')   

   args=parser.parse_args()

   if args.mode not in ['agg_affinity','epitope_merge','summary','all']:
       sys.exit('[ERROR] unknown mode %s.\nMode available: %s,%s,%s and %s (run all steps)' %(args.mode, 'agg_affinity','epitope_merge','summary','all'))

   RHD_out=args.prefix+'.RHD.agg.affinity.txt'
   RHCE_out=args.prefix+'.RHCE.agg.affinity.txt'
   if args.mode.lower() in ['agg_affinity','all']:
     print('[Affinity] processing')
     RHD_aff=[os.path.join(args.aff_RHD_path,f) for f in os.listdir(args.aff_RHD_path) if re.findall("%s" % args.affinity_file_pattern, f)]
     RHCE_aff=[os.path.join(args.aff_RHCE_path, f) for f in os.listdir(args.aff_RHCE_path) if re.findall("%s" % args.affinity_file_pattern, f)]
    
     RHD_affinity=process_affinity(RHD_aff)
     RHCE_affinity=process_affinity(RHCE_aff)

     RHD_affinity.to_csv(RHD_out, index=False, sep='\t')
     RHCE_affinity.to_csv(RHCE_out, index=False, sep='\t')
     print('[Affinity] processing finished')  
 
   elif args.mode.lower() in ['all', 'epitope_merge']:
     print('[Affinity] loading')
     RHD_affinity=pd.read_csv(RHD_out, sep="\t")
     RHCE_affinity=pd.read_csv(RHCE_out, sep="\t")
     print('[Affinity] loaded')

   RHD_out=args.prefix+'.RHD.Epitope.candidate.txt'
   RHCE_out=args.prefix+'.RHCE.Epitope.candidate.txt'
   if args.mode.lower() in ['epitope_merge','all']:
     print('[FA] loading') 
     RHD_fa=fa2pd(args.peptide_fa_RHD)
     RHCE_fa=fa2pd(args.peptide_fa_RHCE)  
     print('[FA] loaded')
   
     print('[Peptide info] loading')
     RHD_info=pd.read_csv(args.peptide_info_RHD, header=0, sep="\t")
     RHCE_info=pd.read_csv(args.peptide_info_RHCE, header=0, sep="\t")
     print('[Peptide info] loaded')


     print('[Epitope merge] processing')
     RHD_info=RHD_info.merge(RHD_fa, on='peptide', how='outer')
     RHCE_info=RHCE_info.merge(RHCE_fa, on='peptide', how='outer')

     RHD_affinity.Identity=RHD_affinity.Identity.astype('str')
     RHCE_affinity.Identity=RHCE_affinity.Identity.astype('str') 

     RHD_epitope=RHD_info.merge(RHD_affinity, on='Identity', how='outer')
     RHCE_epitope=RHCE_info.merge(RHCE_affinity, on='Identity', how='outer')
   
     RHD_epitope=RHD_epitope.rename(columns=lambda c: re.sub('peptide','Full_peptide', c))
     RHCE_epitope=RHCE_epitope.rename(columns=lambda c: re.sub('peptide','Full_peptide', c))

     RHD_epitope=chech_mut_pos(RHD_epitope)
     RHCE_epitope=chech_mut_pos(RHCE_epitope)

     RHD_epitope.to_csv(RHD_out, index=False, sep='\t') 
     RHCE_epitope.to_csv(RHCE_out, index=False, sep='\t')
     lst2del=[RHD_fa, RHCE_fa, RHD_info, RHCE_info, RHD_affinity, RHCE_affinity]
     del lst2del
     print('[Epitope merge] processing finished')
   elif args.mode.lower() in ['all','summary']:
     print('[Epitope merge] loading')
     RHD_epitope=pd.read_csv(RHD_out, sep="\t")
     RHCE_epitope=pd.read_csv(RHCE_out, sep="\t")
     print('[Epitope merge] loaded')

   RHD_out=args.prefix+'.RHD.Epitope.candidate.summary.txt'
   RHCE_out=args.prefix+'.RHCE.Epitope.candidate.summary.txt' 
   if args.mode.lower() in ['summary','all']: 
     print('[Summary] processing')
    
     RHD_epitope_summary=RHD_epitope.groupby(['Identity','Full_peptide','Full_peptide_len','AAstart','AAend','ALT_AApos','REFaa','ALTaa',
                                              'Pos', 'MHC', 'Peptide','Of','Core','Core_Rel','Score_EL','Rank_EL_pct','BindLevel','mutAA_in_core'])['group'].apply(lambda x: ','.join(set(x.dropna().unique()))).reset_index()
     RHD_epitope_summary.drop_duplicates(inplace=True)
     #print(RHD_epitope_summary)
     RHD_epitope_summary.to_csv(RHD_out, index=False, sep='\t')  
     
     RHCE_epitope_summary=RHCE_epitope.groupby(['Identity','Full_peptide','Full_peptide_len','AAstart','AAend','ALT_AApos','REFaa','ALTaa',
                                              'Pos', 'MHC', 'Peptide','Of','Core','Core_Rel','Score_EL','Rank_EL_pct','BindLevel','mutAA_in_core'])['group'].apply(lambda x: ','.join(set(x.dropna().unique()))).reset_index()
     RHCE_epitope_summary.drop_duplicates(inplace=True)
     #print(RHCE_epitope_summary)
     RHCE_epitope_summary.to_csv(RHCE_out, index=False, sep='\t')
  
     print('[Summary] processing finshed')     

if __name__ == "__main__":
    main()





'''
Pos Residue number (starting from 0)
MHC MHC molecule name
Peptide Amino acid sequence
Of Starting position offset of the optimal binding core (starting from 0)
Core Binding core register
Core_Rel Reliability of the binding core, expressed as the fraction of networks in the ensemble selecting the optimal core
Identity Annotation of the input sequence, if specified
Score_EL Eluted ligand prediction score
%Rank_EL Percentile rank of eluted ligand prediction score
Exp_bind If the input was given in PEPTIDE format with an annotated affinity value (mainly for benchmarking purposes).
Score_BA Predicted binding affinity in log-scale (printed only if binding affinity predictions were selected)
Affinity(nM) Predicted binding affinity in nanomolar IC50 (printed only if binding affinity predictions were selected)
%Rank_BA % Rank of predicted affinity compared to a set of 100.000 random natural peptides. This measure is not affected by inherent bias of certain molecules towards higher or lower mean predicted affinities (printed only if binding affinity predictions were selected)
BindLevel (SB: strong binder, WB: weak binder). The peptide will be identified as a strong binder if the % Rank is below the specified threshold for the strong binders. The peptide will be identified as a weak binder if the % Rank is above the threshold of the strong binders but below the specified threshold for the weak binders.
'''

