#!/usr/bin/env python

'''
   Author: Ti-Cheng Chang
'''

#import re, argparse, sys, collections
#import pysam
import argparse, sys, os, re
import numpy as np
import pandas as pd
#import RHtyper
#from . import *
if sys.version_info >= (3,0):
  import RHtyper.BGClf as BGClf ### py3
  from RHtyper.BGClf import *
else:
  ### python 2
  import BGClf
  from BGClf import *



#import matplotlib.pyplot as plt
#import BGClf
#from collections import OrderedDict
##from matplotlib import style
##style.use("ggplot")

def runRHtyper(args):
    
    if 'est' in str(args.coverage).lower():
       print("\n{0}".format("[Coverage estimation] ..."))
       if args.wes:
           args.coverage, medianCov=BGClf.accessory.genome_cov(args.bam, chr='chr1', gbuild=args.genome_build, mode='wes')
       else:
           args.coverage, medianCov=BGClf.accessory.genome_cov(args.bam, chr='chr1', gbuild=args.genome_build)
       args.alt_read_n=args.coverage*0.1
       print("[Estimated average coverage] {0}; [Median coverage] {1}; [Alternative read N] changed to {2}".format(str(args.coverage), str(medianCov), str(args.alt_read_n)))

    prefix=args.prefix +'.'+ args.gene
    gene, cds, rh=BGClf.coordinates.table_create(args.gene, args.reference, args.genome_build) ### new
    
    
    ### determine whether to use company gene information, useful for RHD/RHCE
    company=None
    if args.gene=='RHD':
        company='RHCE'
        gene_c, cds_c, rh_c=BGClf.coordinates.table_create(company, args.reference, args.genome_build) ### new
    elif args.gene=='RHCE':
        company='RHD' 
        gene_c, cds_c, rh_c=BGClf.coordinates.table_create(company, args.reference, args.genome_build) ### new

    ### call variant
    if args.call_variant:
        print("\n[Variant calling] ... {0}".format(args.gene))
        allvar, filteredvar, finalvar=BGClf.variants.call(rh, cds, args.bam, prefix, gbuild=args.genome_build,ALT_n_LB=args.alt_read_n, verbose=int(args.verbose))
    else:
        print("\n[Variant loading] ... {0}".format(args.gene))
        allvar=pd.read_csv(prefix+'.full.variant.txt', sep="\t", header=0)
        #filteredvar=pd.read_table(prefix+'.filtered.variant.txt', sep="\t", header=0)
        if os.path.isfile(prefix+'.final.variant.txt') and os.stat(prefix+'.final.variant.txt').st_size != 0:
            finalvar=pd.read_csv(prefix+'.final.variant.txt', sep="\t", header=0)
        else:
            finalvar=pd.DataFrame()
    
    ### call variant for company gene
    if company is not None:
        company_prefix=args.prefix+'.'+company
        if args.call_variant:
            print("\n[Variant calling, for paired gene] ... {0}".format(company))
            allvarC, filteredvarC, finalvarC=BGClf.variants.call(rh_c, cds_c, args.bam, company_prefix, gbuild=args.genome_build, ALT_n_LB=args.alt_read_n, verbose=int(args.verbose), write_output=True)
        else:
            print("\n[Variant loading, for paired gene] ... {0}".format(company))
            allvarC=pd.read_csv(company_prefix+'.full.variant.txt', sep="\t", header=0)
            if os.path.isfile(company_prefix+'.final.variant.txt') and os.stat(company_prefix+'.final.variant.txt').st_size != 0:
                finalvarC=pd.read_csv(company_prefix+'.final.variant.txt', sep="\t", header=0)
            else:
                finalvarC=pd.DataFrame()
    
    ### produce bloodtype database
    #bgmut, bgmut_cnt, DB=BGClf.database.bloodgroupDB(gene=args.gene, verbose=int(args.verbose)).bgmut_db() 
    bgmut, bgmut_cnt, rawDB=BGClf.database.bloodgroupDB(gene=args.gene, verbose=int(args.verbose)).isbt_db()     

    ### prepare data and keep only compatible varians
    testdata, untrimmed_test, traindata, untrimmed_traindata, target, weight=BGClf.matrixprep.create(bgmut.transpose(), allvar, args.coverage, gene=args.gene, ALT_n_LB=args.alt_read_n, verbose=args.verbose)
 
    ### classification ... 'likelihood based' method
    print("\n{0}".format("[Typing]..."))
    SNPs, SNPm, SNP_GT, GT, GTm=BGClf.matrixprep.SNP4likelihood(testdata, traindata, args.gene, ALT_n_LB=args.alt_read_n)
    df_merged=None
    a_warn=False; b_warn=False; c_warn=False; d_warn=False; e_warn=False; f_warn=False
    exonCov, geneCov=BGClf.accessory.coverage(df=allvar)    

    print("\n{0}".format("[S1. Detecting large deletion]"))
    if args.wes:
        #exon_cov, bam, fasta, coord, exon_cut=0.6, coverag=30, gene='RHD', seqtype='WGS'
        large_del, abnomral_exon_cov=BGClf.accessory.gene_exon_del(exonCov, args.bam, args.reference, gene, gene=args.gene, coverage=args.coverage, seqtype='WES')
        if company is not None:
            exonCov_C, geneCov_C=BGClf.accessory.coverage(df=allvarC)
            large_del_C, abnomral_exon_cov_C=BGClf.accessory.gene_exon_del(exonCov_C, args.bam, args.reference, gene_c, gene=company, coverage=args.coverage, seqtype='WES')

    else:
        large_del, abnomral_exon_cov=BGClf.accessory.gene_exon_del(exonCov, args.bam, args.reference, gene, gene=args.gene, coverage=args.coverage)
        if company is not None:
            exonCov_C, geneCov_C=BGClf.accessory.coverage(df=allvarC)
            large_del_C, abnomral_exon_cov_C=BGClf.accessory.gene_exon_del(exonCov_C, args.bam, args.reference, gene_c, gene=company, coverage=args.coverage)
 
    ##############################################
     ### test block: speed up the debug process ###
    ##############################################

    if large_del:
        new_data=[]
        if 'Homozygous deletion' in large_del:
            print('hom del')
            new_data.append({'Allele':'A1', 'Bloodtype':large_del, 'Likelihood':0, 'Index':0, "matchDBaa":"","not_matchDBaa":"","SampleExtraVar_aa":"","AllVariations_aa":""})
            new_data.append({'Allele':'A2', 'Bloodtype':large_del, 'Likelihood':0, 'Index':0, "matchDBaa":"","not_matchDBaa":"","SampleExtraVar_aa":"","AllVariations_aa":""})
            df_merged=pd.concat([pd.DataFrame(new_data)], ignore_index=True, sort=True)
        elif 'Heterozygous deletion' in large_del:
            print('het del')
            allele_classified=BGClf.classification.LH(args.bam, args.gene, allvar, GTm, gbuild=args.genome_build, allele_type='hom',ALT_n_LB=args.alt_read_n).calculateLL()
            allele_classified=allele_classified.loc[(allele_classified['Allele']=='A2') & (allele_classified['Likelihood']==allele_classified['Likelihood'].max())]
            allele_idicies=allele_classified['Index'].unique()
            for allele_idx in allele_idicies:
                new_data.append({'Allele':'A1', 'Bloodtype':large_del, 'Likelihood':0, 'Index':allele_idx, "matchDBaa":"","not_matchDBaa":"","SampleExtraVar_aa":"","AllVariations_aa":""})
                df_merged=pd.concat([pd.DataFrame(new_data), allele_classified], ignore_index=True, sort=True)
        df_merged_raw=df_merged.copy()
        df_merged_raw['Comment']=""
    else:
        SNPm_test=SNPm.copy()
        llmode=None
        Var1136=""
        if SNPm_test.shape[0] > 0:
             ### check RHD 1136 first
             if args.gene=='RHD' and '1136_C_T' in SNPm_test.index:
                 print("[RHD1136] checking")
                 #   The RHD 1136 has relatveily low coverage at this position due to exact duplication of exon 8.  It was adjusted to account for that by considering only heterozygrous combination of alleles
                 if (SNPm_test.loc['1136_C_T','alt_cnt']/args.coverage > args.heterozyoteCov) and (SNPm_test.loc['1136_C_T','alt_cnt']/args.coverage < (1-args.heterozyoteCov)):
                     print('[RHD1136] mode changed to hetero RHD 1136_C_T, 1136_cov:', SNPm_test.loc['1136_C_T','alt_cnt']/args.coverage, '; WGS coverage:',args.coverage, '; cutoff: ', args.heterozyoteCov)
                     llmode="het RHD 1136_C_T"
                     Var1136="het RHD 1136_C_T"
                 else:
                     print('[RHD1136] mode changed to hom RHD 1136_C_T, 1136_cov:', SNPm_test.loc['1136_C_T','alt_cnt']/args.coverage, '; WGS coverage:',args.coverage, '; cutoff: ', args.heterozyoteCov) 
                     llmode="hom RHD 1136_C_T"
                     Var1136="hom RHD 1136_C_T"
                 SNPm_test=SNPm_test.drop('1136_C_T')

        #print("TEST")
        #print(GTm, SNPm_test)
        if SNPm_test.shape[0] > 0:
             if GTm.shape[1] < SNPm_test.shape[0]:
                 ### this criteria is used to solve the var not in known database combination
                 llmode='hom/het'
             elif any(SNPm_test['allele_frac'] > args.heterozyoteCov) & any(SNPm_test['allele_frac'] < (1 - args.heterozyoteCov)):
                 llmode='het'
             elif args.gene=='RHD' and (any(((SNPm_test['alt_cnt']/args.coverage) > args.heterozyoteCov)) or any(((SNPm_test['alt_cnt']/args.coverage) < (1-args.heterozyoteCov)))):
                 ### this criteria is used to solve the problematic region on RHD
                 llmode='het'
             else:
                 llmode='hom' if (SNPm_test.shape[0]==1) & (SNPm_test['allele_frac'].max()==1) else 'hom/het'
        else:
             if llmode is None: llmode='hom'

        print ('\n[Zygosity] {0}'.format(llmode))
        if args.gene=='RHD': print ('[VAR1136] {0}'.format(Var1136))
        #print ('[SNP]', SNPm, '\n')



        print ("\n{0}".format("[S2. Calculating likelihood]"))
        if llmode in ['het RHD 1136_C_T','het']:
            df_merged=BGClf.classification.LH(args.bam, args.gene, allvar, GTm, gbuild=args.genome_build, allele_type='het', ALT_n_LB=args.alt_read_n).calculateLL()
            df_merged_raw=df_merged.copy()
            #print(df_merged)
            df_merged=df_merged.loc[df_merged['Likelihood']==df_merged['Likelihood'].max()]
        elif llmode in ['hom', 'hom RHD 1136_C_T']:
            df_merged=BGClf.classification.LH(args.bam, args.gene, allvar, GTm, gbuild=args.genome_build, allele_type='hom', ALT_n_LB=args.alt_read_n).calculateLL()
            df_merged_raw=df_merged.copy()
            #print(df_merged)
            df_merged=df_merged.loc[df_merged['Likelihood']==df_merged['Likelihood'].max()]
        elif llmode in ['hom/het']:
            df_merged=BGClf.classification.LH(args.bam, args.gene, allvar, GTm, gbuild=args.genome_build, ALT_n_LB=args.alt_read_n).calculateLL()
            df_merged_raw=df_merged.copy()
            #print(df_merged)
            df_merged=df_merged.loc[df_merged['Likelihood']==df_merged['Likelihood'].max()]

        df_merged_raw['Comment']=""
        df_merged_raw, a_warn=BGClf.accessory.filter_annotation(df_merged_raw, df_merged, annotation='a', LH_adj=0)


        print("[S2.1 RHD 1136] checking 1136 C/T var")
        #print(df_merged)
        ### remove the combination with both alleles contating 1136C_T
        if Var1136 in ['het RHD 1136_C_T']: 
            df_merged=BGClf.accessory.RHD1136check(df_merged, untrimmed_traindata, Var1136)
            df_merged_raw, b_warn=BGClf.accessory.filter_annotation(df_merged_raw, df_merged, annotation='b')
        #print("after 1136 check")
        #print(df_merged)

    df_merged['ID'] = args.prefix
    df_merged_raw['ID']=args.prefix


    #if args.verbose >9: print "Merged_df\n", df_merged
    #if args.verbose >9: print "Abnormal coverage\n", abnomral_exon_cov
    df_merged.sort_values(inplace=True, axis=0, by=['Likelihood','Index'], ascending=[0,1])
    
    #final=df_merged.head(2)  
    final=df_merged
    final.reset_index(drop=True, inplace=True)
    
    ### output for test purpose  
    #final.to_csv(prefix+'.bloodtyping.test.txt', sep='\t', index=False, mode='w')

    ######################################################
    ### [1] test block end: speed up the debug process ###
    ######################################################
    
    #final=pd.read_csv(prefix+'.bloodtyping.test.txt', sep="\t", header=0)    


    ### other rescue calling for specific type
    print("\n{0}".format("[S3. Exon CNV segmentaion]..."))
    if args.gene=='RHD' or args.gene=='RHCE':
        #rhd, rhce, rhCov=BGClf.accessory.RHD_RHCE_exon_CNV(exonCov, exonCov_C, args.coverage, prefix, verbose=args.verbose)
        if args.gene=='RHD':
            rhd, rhce, rhCov, rhdBreakS, rhdBreakE, rhceBreakS, rhceBreakE=BGClf.accessory.RHD_RHCE_segment_CNV(gene, gene_c, exonCov, exonCov_C, args.bam, args.coverage, args.reference, prefix, verbose=args.verbose)
        elif args.gene=='RHCE':
            rhd, rhce, rhCov, rhdBreakS, rhdBreakE, rhceBreakS, rhceBreakE=BGClf.accessory.RHD_RHCE_segment_CNV(gene_c, gene, exonCov_C, exonCov, args.bam, args.coverage, args.reference, prefix, verbose=args.verbose)
        if args.gene=='RHCE':
            final=BGClf.accessory.RHCE_relink(args, final, finalvar, rhCov, gbuild=args.genome_build)
            df_merged_raw, c_warn=BGClf.accessory.filter_annotation(df_merged_raw, final, annotation='c', LH_adj=0)

    ### output matched and non-matched mutation for prediced alleles
    print("\n{0}".format("[S4. Identifying matched/non-matched vars in database]..."))
    #consistency=BGClf.accessory.var_report(final, untrimmed_traindata, finalvar, args.gene)
    #print(finalvar)
    consistency=BGClf.accessory.var_compare(final, rawDB, finalvar, finalvarC, allvar, args.gene, args.bam, args.reference)
    final.loc[:,'matchDB']=""
    final.loc[:,'matchDBn']=0
    final.loc[:,'not_matchDB']=""
    final.loc[:,'not_matchDBn']=0
    final.loc[:,'matchDBaa']=""
    final.loc[:,'not_matchDBaa']=""
    final.loc[:,'Alias']=""


    for idx, row in final.iterrows():
        BG=re.sub(".*@","",row['Bloodtype'])
        BG=re.sub("_"," ", BG)
        if 'Alias' in rawDB.columns:
            if not rawDB.loc[rawDB['Allele name']==BG,].empty:
                BGalias=rawDB.loc[rawDB['Allele name']==BG,'Alias'].values[0]
                if BG != BGalias: final.loc[idx, 'Alias']=BGalias      
        if row['Bloodtype'] in consistency:
            final.loc[idx,'matchDB']=consistency[row['Bloodtype']]['match']
            final.loc[idx,'not_matchDB']=consistency[row['Bloodtype']]['non-match']
            if len(consistency[row['Bloodtype']]['non-match']) >0:
                final.loc[idx,'not_matchDBn']=len(consistency[row['Bloodtype']]['non-match'].split(";"))
            if len(consistency[row['Bloodtype']]['match']) >0:
                final.loc[idx,'matchDBn']=len(consistency[row['Bloodtype']]['match'].split(";"))
            final.loc[idx,'matchDBaa']=consistency[row['Bloodtype']]['match_aa']
            final.loc[idx,'not_matchDBaa']=consistency[row['Bloodtype']]['non-match_aa']
        ### re-assign the BG name
        final.loc[idx,'Bloodtype']=BG
 
    final.loc[:,'SampleExtraVar']=';'.join(consistency['sample_specific'])
    final.loc[:,'SampleExtraVar_aa']=';'.join(consistency['sample_specific_aa'])
    final.loc[:,'AllVariations']=';'.join(consistency['all_variation'])
    final.loc[:,'AllVariations_aa']=';'.join(consistency['all_variation_aa'])
    outflds=['ID','Allele','Index','Bloodtype','Alias', 'Likelihood','matchDB','matchDBn','matchDBaa', 'not_matchDB','not_matchDBn','not_matchDBaa','SampleExtraVar','SampleExtraVar_aa','AllVariations','AllVariations_aa']
    final=final[outflds]
  
    #print('check1') 
    #print(final)

 
    ### keep only the type with minimum mismatch DBn
    if len(final.Index.unique())> 1: 
        xmDBn, maxDBn, mmDBn, minDBn=BGClf.accessory.max_match(final)
        final=final.loc[final.Index.isin(mmDBn[minDBn])]
    if len(final.Index.unique())> 1: 
        xmDBn, maxDBn, mmDBn, minDBn=BGClf.accessory.max_match(final)
        final=final.loc[final.Index.isin(xmDBn[maxDBn])]

    df_merged_raw, d_warn=BGClf.accessory.filter_annotation(df_merged_raw, final, annotation='d')

    ### identify hybrid alleles
    print("\n{0}".format("[S5. Identifying hybrid alleles]..."))
    if args.gene=='RHD' or args.gene=='RHCE':
        final=final.assign(CoverageProfile="")
        final=final.assign(Break="") 
        if args.gene=='RHD':
             if 'h' not in rhd:
                 linked_BG=['RHD_III_type_4','RHD*03.01_RHD*DIIIa']
                 select_indices=final['Bloodtype'].isin(linked_BG)
                 ### het hybrid
                 if final[final['Bloodtype'].str.contains('RHD_III_type_4')].shape[0] > 0:
                     ### BGMUT
                     ### link hybrid to RHD III (no strict rule here)
                     index=final[final['Bloodtype'].str.contains('RHD_III_type_4')].index
                     final.loc[index,'CoverageProfile']=rhd
                 elif final[final['Bloodtype'].str.contains('RHD*03.01_RHD*DIIIa')].shape[0] > 0:
                     index=final[final['Bloodtype'].str.contains('RHD*03.01_RHD*DIIIa')].index
                     final.loc[index,'CoverageProfile']=rhd
                 else:
                     final=final.assign(CoverageProfile=rhd) 
             else:
                 final=final.assign(CoverageProfile=rhd)
 
             if rhdBreakS is not None and rhdBreakE is not None:
                 final.Break=str(rhdBreakS) + ':' + str(rhdBreakE)
                 
        else:
             if 'h' not in rhce:
                 ### het hybrid
                 final=final.assign(CoverageProfile=rhce)
                 ### not output hybrid ce allele
                 #final=final.assign(CoverageProfile='') 
             else:
                 final=final.assign(CoverageProfile=rhce) 
                 ### not output hybrid ce allele
                 #final=final.assign(CoverageProfile='')
             if rhceBreakS is not None and rhceBreakE is not None:
                 final.Break=str(rhceBreakS) + ':' + str(rhceBreakE)

        ### re-assign allele name for hybrid 
        if args.gene=='RHD': final=BGClf.accessory.RHhybrid_rename(final) 
 
    #print(final) 
 
    ### correct type based on validataion cohort
    if args.popfreq:
        print("\n{0}".format("[S6. Adjust for population frequency..."))
        final=BGClf.accessory.PopFreq(final, args.gene, frequency_table=args.popfreqdb) if args.popfreqdb else BGClf.accessory.PopFreq(final, args.gene)
         
    else:
        if len(final.Index.unique()) > 1:
            print("\n{0}".format("[S6. Adjust for population frequency]..."))
            final=BGClf.accessory.PopFreq(final, args.gene, frequency_table=args.popfreqdb) if args.popfreqdb else BGClf.accessory.PopFreq(final, args.gene)
    #print(final)


    df_merged_raw, e_warn=BGClf.accessory.filter_annotation(df_merged_raw, final, annotation='e')

    ### adjust the deletion allele name
    for idx, row in final.iterrows(): 
        if 'deletion' in row['Bloodtype']: final.loc[idx,'Bloodtype']='Deletion'

    df_merged_raw, f_warn=BGClf.accessory.filter_annotation(df_merged_raw, final, annotation='f')

 
    
    cols=['ID','Index','Allele','Bloodtype','Likelihood','Comment']
    df_merged_raw=df_merged_raw[cols]
    df_merged_raw.sort_values(inplace=True, axis=0, by=['Likelihood','Index','Allele'], ascending=[0,1,1]) 


    if args.verbose >3: 
        print("[Typed result]")
        print(df_merged_raw)
        print(final)

    df_merged_raw.to_csv(prefix+'.allele_ranking.txt', sep='\t', index=False, mode='w')
    final.to_csv(prefix+'.bloodtyping.txt', sep='\t', index=False, mode='w')
   
    xlsxe="xlsxwriter"

    with pd.ExcelWriter(prefix+'.bloodtyping.xlsx', engine=xlsxe) as writer: 
        final.to_excel(writer,sheet_name='bloodtyping', index=False) 
        df_merged_raw.to_excel(writer, sheet_name='allele_ranking', index=False)

        workbook=writer.book

        if e_warn:
            worksheet=writer.sheets['bloodtyping']
            leg="NOTE: Ambiguous/Incomplete SNP phasing. RNA testing would be needed to confirm the SNPs are on the same or separate alleles."
            worksheet.write(4,0,leg)

        worksheet=writer.sheets['allele_ranking']
        
        legends=['Comment legend:', \
                 'a: Allele pair with suboptimal likelihood', \
                 'b: Allele pair with non-heterozygous RHD 1136 C>T', \
                 'c: The RHD-RHCE allele linkage is not preferred', \
                 'd: Allele pair with higher number nucleotide variations mismatched to the allele database compared to the best pair', \
                 'e: Allele pair with lower population frequency']
        li=0
        for l in legends:
            li+=1
            if xlsxe=='xlsxwriter':
              ### xlsxwriter
              worksheet.write(len(df_merged_raw) + 1 + li, 0, l)
            elif xlsxe=='openpyxl':
              ### openpyxl
              worksheet.cell(row=len(df_merged_raw) + 1 + li, column=1).value =l
            else:
              print("[Error] unknown xlsx engine, no xlsx output will be generated")
        writer.save()

    #''' ### [2] test block end: speed up the debug process

    ### generate PDF report
    warning_sec=None
    warning_msg=None
    if e_warn:
        leg="NOTE: Ambiguous/Incomplete SNP phasing. RNA testing would be needed to confirm the SNPs are on the same or separate alleles."
        warning_sec='Genotype'
        warning_msg=leg
    BGClf.report.RHD_RHCE(prefix, args.gene, prefix+'.bloodtyping.txt', prefix+'.final.variant.txt', prefix+'.exonCNV.txt', prefix+'RHD.bincov.allpos.png', prefix+'RHD.bincov.exonpos.png', \
                          prefix=prefix, warning_sec=warning_sec, warning_msg=warning_msg)
    #BGClf.report.HTML(prefix, args.gene, prefix+'.bloodtyping.txt', finalvar, allvar, prefix+'.exonCNV.txt', prefix+'RHD.bincov.allpos.png', prefix+'RHD.bincov.exonpos.png', prefix=prefix).generate()
    
 
    print("\n{0}".format("[Bloodtyping finished]"))

    



def main():
    parser = argparse.ArgumentParser(description='Bloodtyping')
    parser.add_argument('-bam' , '--bam', help='BAM file', required=True)
#    parser.add_argument('-coord', '--gene_table', help='Coordinates of coding exons in BED format\n(can be obtained from UCSC table browser)', required=True)
    parser.add_argument('-gene' , '--gene', help='Name of gene for typing', required=True)
    parser.add_argument('-ref' , '--reference', help='Reference fasta file', required=True)
    parser.add_argument('-gbuild' , '--genome_build', help='Genome build version', required=True, default='hg38') ### new
    parser.add_argument('-pre' , '--prefix', help='Output filename prefix', required=True)
    parser.add_argument('-call' , '--call_variant', help='Call variant. If not call variants, the variant file with name of <output>.full.variant.txt will be loaded for bloodtyping', action='store_true')
    parser.add_argument('-cov' , '--coverage', help='Genome coverage depth. The coverage value is important for calling gene deletion. Use "est" or "estimate" for automatical coverage estimation', required=True)
    parser.add_argument('-altN' , '--alt_read_n', help='Alternative read support N', type=int, required=True, default=3)
    parser.add_argument('-hetCov' , '--heterozyoteCov', help='Cutoff to call as heterozygotes', type=int, required=False, default=0.1)
    parser.add_argument('-wes' , '--wes', help='Input bam are WES', action='store_true')
    parser.add_argument('-phase2' , '--rescue_phase', help='Strategy applied when phasing using SNPs failed (options: ratio/reference). Ratio: The SNPs with alternative allele ratio < 0.4 or > 0.6 will be adjusted to 1. Reference: one allele will be considered with the reference allele always ', required=False, default=None)
    parser.add_argument('-link' , '--linking_bloodtype', help='Link bloodgroup by providing typing result of another bloodtype (useful for RHD/RHCE)', required=False)
    parser.add_argument('-linkdb' , '--linking_database', help='Databse used for bloodtype linkage. NOTE: if not specified, the linkage will be based on ISBT information ', required=False)
    parser.add_argument('-popfreq', '--popfreq', help='adjust prediction based on the population frequency of the allele type', action='store_true')
    parser.add_argument('-popfreqdb', '--popfreqdb', help='Database for allele frequencies. NOTE: if not specified, the default database will be based on African American population', required=False)
    parser.add_argument('-v' , '--verbose', help='Output verbose level', required=False, type=int, default=0)
    args=parser.parse_args()
    if 'est' not in str(args.coverage).lower():
         args.coverage=float(args.coverage)
    runRHtyper(args)


if __name__ == "__main__":
    main()

