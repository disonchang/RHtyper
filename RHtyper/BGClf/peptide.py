#!/usr/bin/env python
"""
    modules to produce epitope matrix 

"""
import pandas as pd
import numpy as np
import re, collections, pysam, sys, os, argparse
from collections import OrderedDict
from general import *
from coordinates import *
from database import *

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

def peptide(cds, cnt, output):
    if (len(cds) != 1):
       sys.exit('[ERROR] mulitple gene in the matrix ... exit')
    
    cds_df=cds[list(cds.keys())[0]]

    #print(cds_df.columns)

    bg=['REF']

    for bloodgroup, row in cnt.iterrows():
        bg.append(str(bloodgroup))
        var=list(row[row>0].index)
        cds_df[bloodgroup+'__NUC']=cds_df['REF']
        for v in var:
           pos, ref, alt=v.split('_')
           cds_df.loc[cds_df['cpos']==int(pos), bloodgroup+'__NUC']=alt

    ### translation
    for bgrp in bg:
        tri_pos=list()
        tri_nuc=list()
        NUC=bgrp+'__NUC'
        PT=bgrp+'__PT'
        AApos=bgrp+'__AApos'
        aapos=0
        premature=False
        if bgrp=='REF': 
            NUC=bgrp
            PT='REFaa'
        cds_df.insert(cds_df.columns.get_loc(NUC)+1, column=PT, value='.')
        if bgrp != 'REF':
             cds_df.insert(cds_df.columns.get_loc(PT), column=AApos, value=0)
        print('[Peptide][bloodgroup] ', bgrp) 
        for idx,row in cds_df.iterrows():
            #print(bgrp, idx, row[NUC])
            if row[NUC] in ['A','T','C','G']:
                if idx not in tri_pos:
                    tri_pos.append(idx)
                    tri_nuc.append(row[NUC])
            if len(tri_pos) == 3:
                aapos+=1
                tri_nuc=''.join(tri_nuc)
                aa=translate(tri_nuc)
                if aa=='*':premature=True
                if premature:aa='*'
                for p in tri_pos:
                    cds_df.loc[cds_df.index==p,PT]=aa
                    if bgrp != 'REF': 
                        cds_df.loc[cds_df.index==p,AApos]=aapos
                #print(bgrp, idx, row[NUC], 'translation', tri_nuc, translate(tri_nuc))
                tri_pos=list()
                tri_nuc=list()
            
    cds_df.to_csv(output, sep="\t", index=False)

    return bg, cds_df


def epitope(bg, cds_df, span=14, outfa="temp.fasta", outtb='temp.txt'):
   epi_id=0
   out=list()
   epitope_fa=list()
   ref_epitope_fa=list()
   for bgrp in bg:
        NUC=bgrp+'__NUC'
        PT=bgrp+'__PT'
        AApos=bgrp+'__AApos'
        if bgrp=='REF': continue
        print('[Epitope][bloodgroup] ', bgrp)

        region2exc=dict()
        
        for baapos in range(1, cds_df[AApos].max()):
             baa=cds_df.loc[cds_df[AApos]==baapos, PT].unique()
             raa=cds_df.loc[cds_df[AApos]==baapos, 'REFaa'].unique()
             if len(baa)==1:
                baa=baa[0]
             else:
                sys.exit('[ERROR] mulitple BG aa annotated at the same position ... exit')
             
             if len(raa)==1:
                raa=raa[0]
             else:
                print('[WARN] mulitple REF aa annotated at the same position')
                raa=','.join(raa)                

             if baa != raa and baa != '*':
                 epi_id+=1
                 s=baapos-span if baapos-span > 0 else 1
                 e=baapos+span if baapos+span <= cds_df[AApos].max() else cds_df[AApos].max()
              
                 baa_all=cds_df.loc[ (cds_df[AApos] >= s) & (cds_df[AApos] <= e) , [AApos,PT,'aapos', 'REFaa']]
 
                 ref_baa_all=baa_all.groupby(['aapos','REFaa']).count().reset_index().rename(columns={0:'count'})
                 baa_all=baa_all.groupby([AApos,PT]).count().reset_index().rename(columns={0:'count'})                    
 

                 epitope=''.join(baa_all[PT].to_list())
                 epitope=re.sub(r'\*.*','',epitope)
 
                 refepitope=''.join(ref_baa_all['REFaa'].to_list()) 
                 refepitope=re.sub(r'\*.*','',refepitope)               

                 outrow=[bgrp, s, e, baapos, raa, baa, len(epitope)]
                 refoutrow=[bgrp+'_REF', s, e, baapos, raa, raa, len(refepitope)]

                 ### NO need to ouput fasta here as many are duplicated peptide sequences, instead output to a unqiue list, then merge later
                 #fa_header='_'.join(map(str,outrow))
                 if epitope not in epitope_fa:
                     epitope_fa.append(epitope)
                 if refepitope not in epitope_fa:
                     epitope_fa.append(refepitope)
                 outrow.append(epitope)
                 refoutrow.append(refepitope)

                 out.append(outrow)
                 out.append(refoutrow) 

   out=pd.DataFrame(out, columns=['group','AAstart','AAend','ALT_AApos','REFaa','ALTaa', 'peptide_len','peptide'])                   
   out.to_csv(outtb, sep="\t") 

   fh=open(outfa,'w')
   epi_n=0
   for fa in epitope_fa:
      epi_n+=1
      fh.write('>'+str(epi_n)+'\n'+fa+'\n')
   fh.close() 

def run_netMHCIIpan(table, fa, allele='DRB1_0101,DRB1_0102', outtb='temp.txt', dry_run=True):
   if not dry_run: 
       table=pd.read_csv(table, sep="\t")
       aff=netMHCIIpan(fa, allele=allele, dry_run=dry_run)
       table=table.merge(aff,how='outer',left_on='epi_id',right_on='Identity')
       table.to_csv(outtb, sep="\t", index=False)
       return outtb
   else:
       print("[DRY_RUN]")
       cmd=netMHCIIpan(fa, allele=allele, dry_run=dry_run) + ' > ' + outtb 
       return cmd

def netMHCIIpan(peptide_fn, allele="DRB1_0101,DRB1_0102", MHC="/home/tchang1/software/netMHCIIpan-4.0/netMHCIIpan", filter_res=1, size=15, rank=20, dry_run=True):
    if dry_run:
        cmd=[MHC,"-f", peptide_fn, "-a", allele, "-filter", str(1), "-rankF", str(100), " | grep -E -v \'^#|^$|^---|Peptide|^\s+Pos|Number\'"]
        return(' '.join(cmd))

    else:
        if os.path.isfile(MHC):
                #### py2.7
                #cmd=MHC + " -f " + peptide_fn + ' -a ' + str(allele) + ' -inptype ' + str(0) + ' | grep -E -v \'^#|^$|^---|^\s+pos|Number\''
                #status, output = commands.getstatusoutput(cmd)
                
                #### py3.7
                cmd=[MHC,"-f", peptide_fn, "-a", allele, "-filter", str(1), "-rankF", str(100)]
                
                print("[cmd] %s" % ' '.join(map(str,cmd)))
                mhc_p= Popen(cmd, stdin=PIPE, stdout=PIPE, stderr=PIPE, universal_newlines=True)
                raw, mhc_p_err=mhc_p.communicate()
                
                grep_p=Popen(["grep", "-E", "-v", "^#|^$|^---|^\s+Pos|Number"],stdin=PIPE,stdout=PIPE,universal_newlines=True)
                output, grep_p_err=grep_p.communicate(raw)
                
                MHC_out=pd.read_table(StringIO(output), header=None, sep='\s+', names=['Pos','MHC','Peptide','Of','Core','Core_Rel','Identity','Score_EL','Rank_EL_pct','Exp_Bind','BindLevel'])
                MHC_out_filtered=MHC_out
                if filter_res == 1:
                        MHC_out_filtered = MHC_out[MHC_out['Rank_EL_pct'] <= int(rank)]
                #pd_print_full(MHC_out_filtered)
                return MHC_out_filtered
        else:
                print("[Error] netMHCIIpan not found ... quit\n")
                quit()

def MHCII_alleles(Locus='all', fn="/home/tchang1/software/RHtyper_dnanexus/resources/RHtyper/RHtyper/BGClf/MHCII_allele_names.txt"):
   df=pd.read_csv(fn, sep="\t") 

   final=list()
   ### DRB
   if Locus in ['all','DRB']:
      DRB=df.loc[df['MHCII']=='DR',]
      DRB_alleles=list()
      for x in DRB.Allele.to_list():
         temp=re.sub(r'\*','_',x)
         temp=re.sub(':','',temp)
         DRB_alleles.append(temp)
      HighQ_DRB='DRB1_0101,DRB1_0102,DRB1_0103,DRB1_0301,DRB1_0305,DRB1_0401,DRB1_0402,DRB1_0403,DRB1_0404,DRB1_0405,DRB1_0408,DRB1_0701,DRB1_0801,DRB1_0803,DRB1_0901,DRB1_1001,DRB1_1101,DRB1_1104,DRB1_1201,DRB1_1301,DRB1_1302,DRB1_1303,DRB1_1401,DRB1_1402,DRB1_1454,DRB1_1501,DRB1_1503,DRB1_1601,DRB3_0101,DRB3_0202,DRB4_0101,DRB4_0103,DRB5_0101,DRB5_0202'
      HighQ_DRB=HighQ_DRB.split(',')
      for x in HighQ_DRB:
         if x not in DRB_alleles:
             print("[WARN] ",x, " not found in DRB_alleles")
      if Locus=='DRB': return DRB_alleles
      final.extend(DRB_alleles)
      

   ### DP
   if Locus in ['all','DP']:
      DPa=df.loc[df['MHCII']=='DP_alpha',]
      DPb=df.loc[df['MHCII']=='DP_beta',]
      DP_alleles=list()
      for a in DPa.Allele.to_list():
         a='HLA-'+re.sub(r'\*|:','',a)
         for b in DPb.Allele.to_list():
            allele=a+'-'+re.sub(r'\*|:','',b)
            DP_alleles.append(allele)
      HighQ_DP="HLA-DPA10103-DPB10201,HLA-DPA10103-DPB10301,HLA-DPA10103-DPB10401,HLA-DPA10103-DPB10402,HLA-DPA10103-DPB10601,HLA-DPA10103-DPB11101,HLA-DPA10103-DPB11701,HLA-DPA10103-DPB12001,HLA-DPA10103-DPB12301"
      HighQ_DP=HighQ_DP.split(',')
      for x in HighQ_DP:
         if x not in DP_alleles:
             print("[WARN] ",x, " not found in DP_alleles")
      if Locus=='DP': return DP_alleles
      final.extend(DP_alleles)

   ### DQ
   if Locus in ['all','DQ']:
      DQa=df.loc[df['MHCII']=='DQ_alpha',]
      DQb=df.loc[df['MHCII']=='DQ_beta',]
      DQ_alleles=list()
      for a in DQa.Allele.to_list():
         a='HLA-'+re.sub(r'\*|:','',a)
         for b in DQb.Allele.to_list():
            allele=a+'-'+re.sub(r'\*|:','',b)
            DQ_alleles.append(allele)
      HighQ_DQ="HLA-DQA10102-DQB10501,HLA-DQA10102-DQB10602,HLA-DQA10102-DQB10604,HLA-DQA10103-DQB10501,HLA-DQA10103-DQB10603,HLA-DQA10201-DQB10201,HLA-DQA10201-DQB10202,HLA-DQA10401-DQB10301,HLA-DQA10501-DQB10201,HLA-DQA10501-DQB10301,HLA-DQA10505-DQB10301"
      HighQ_DQ=HighQ_DQ.split(',')
      for x in HighQ_DQ:
         if x not in DQ_alleles:
            print("[WARN] ",x, " not found in DQ_alleles")
      if Locus=='DQ': return DQ_alleles
      final.extend(DQ_alleles)

   ### Mouse
   if Locus in ['all','Mouse']:
      Mouse=df.loc[df['MHCII']=='Mouse',] 
      Mouse_alleles=Mouse.Allele.to_list()
      if Locus=='Mouse': return Mouse_alleles
      final.extend(Mouse_alleles)
   
   ### 
   return final 
 
def main():
   parser = argparse.ArgumentParser(description='RH peptide to epitopes')
   parser.add_argument('-pre' , '--prefix', help='Output prefix', required=True)
   parser.add_argument('-gene' , '--gene', help='Gene, RHD or RHCE', default='RHD')
   parser.add_argument('-ref' , '--reference', help='Reference fasta file', default='/home/tchang1/software/RHtyper_dnanexus/RHtyper_asset/resources/reference/hg38/GRCh38_no_alt.fa')  
   parser.add_argument('-gbuild','--genome_build', help='Genome build version', default='hg38') 
   parser.add_argument('-mode','--mode', help='mode,[all(epitope + affinity), epitope, affinity]', default='all')
   parser.add_argument('-dry', '--dry_run', help='Not run epitope prediction, but generate commands file only', action="store_true")

   args=parser.parse_args()
   print('[Gene]'+args.gene)
   if args.mode not in ['all','epitope','affinity']:
      sys.exit("[ERROR] unknown mode %s" % (args.mode))
    
   print('[Mode]'+args.mode) 
   g, cds, rh=table_create(args.gene, args.reference, args.genome_build)
   pcnt_wide, cnt_wide, RhBdb=bloodgroupDB(gene=args.gene).isbt_db()
 
   fa_out=args.prefix +'.epitope.fa'
   fa_info_out=args.prefix +'.epitope.info.txt'
 
   if args.mode in ['all', 'epitope']:
       bg, cds=peptide(cds, cnt_wide, args.prefix+'.peptide.txt')
       epitope(bg, cds, outfa=fa_out, outtb=fa_info_out)

   step_n=20
   types=['DRB','DP','DQ'] ## 'DRB','DP','DQ'

   if args.mode in ['all', 'affinity']:
     if args.dry_run:
        cmds_outfn=args.prefix +'.epitope.cmds.txt'
        cmds=list()
        fh=open(cmds_outfn, 'w')
     for MHCII_type in types: 
      alleles=MHCII_alleles(Locus=MHCII_type)
      #print(alleles)
      for allele_n in range(0,len(alleles),step_n):
          print(allele_n, allele_n+step_n)
          affinity_out=args.prefix+'.epitope.affinity.' + MHCII_type + '.' + str(allele_n) + '_' + str(allele_n+step_n) + '.txt'
          cmd=run_netMHCIIpan(fa_info_out, fa_out, allele=','.join(alleles[allele_n:allele_n+step_n]), outtb=affinity_out, dry_run=args.dry_run)
          if args.dry_run:
              fh.write(cmd+'\n')
     if args.dry_run: fh.close()
              

if __name__ == "__main__":
    main()



