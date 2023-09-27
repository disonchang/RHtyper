#!/usr/bin/env python

'''
   Author: Ti-Cheng Chang
'''

#import re, argparse, sys, collections
#import pysam
import argparse, sys, os, re
import numpy as np
import pandas as pd
import os.path

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
   
    args.dtype='WES' if args.wes else 'WGS'
    VAF=None if args.dtype=='WGS' else args.heterozyoteCov
    ALT_p_LB=None if args.dtype=='WGS' else args.heterozyoteCov
 
    if 'est' in str(args.coverage).lower():
       print("\n{0}".format("[Coverage estimation] ..."))
       if args.wes:
           meanCov, medianCov=BGClf.accessory.genome_cov(args.bam, chr='chr1', gbuild=args.genome_build, mode='wes',wes_expand=0, drop_low_coverage=False)
           args.coverage=round(float(meanCov),5)
       else:
           meanCov, medianCov=BGClf.accessory.genome_cov(args.bam, chr='chr1', gbuild=args.genome_build)
           args.coverage=round(float(meanCov),5)
       args.alt_read_n=args.coverage*0.1
       print("[Estimated average coverage] {0}; [Median coverage] {1}; [Alternative read N] changed to {2}".format(meanCov, medianCov, args.alt_read_n))

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
        allvar, filteredvar, finalvar=BGClf.variants.call(rh, cds, args.bam, prefix, gbuild=args.genome_build,ALT_n_LB=args.alt_read_n, verbose=int(args.verbose), VAF=VAF)
    else:
        print("\n[Variant loading] ... {0}".format(args.gene))
        allvar=pd.read_csv(prefix+'.full.variant.txt', sep="\t", header=0)
        #filteredvar=pd.read_table(prefix+'.filtered.variant.txt', sep="\t", header=0)
        if os.path.isfile(prefix+'.final.variant.txt') and os.stat(prefix+'.final.variant.txt').st_size != 0:
            finalvar=pd.read_csv(prefix+'.final.variant.txt', sep="\t", header=0)
        else:
            finalvar=pd.DataFrame()
    
    print("[RHtyper] finalvar:", finalvar)
    
    ### call variant for company gene
    if company is not None:
        company_prefix=args.prefix+'.'+company
        if args.call_variant:
            print("\n[Variant calling, for paired gene] ... {0}".format(company))
            allvarC, filteredvarC, finalvarC=BGClf.variants.call(rh_c, cds_c, args.bam, company_prefix, gbuild=args.genome_build, ALT_n_LB=args.alt_read_n, verbose=int(args.verbose), write_output=True, VAF=VAF)
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
    testdata, untrimmed_test, traindata, untrimmed_traindata, target, weight=BGClf.matrixprep.create(bgmut.transpose(), allvar, args.coverage, gene=args.gene, ALT_n_LB=args.alt_read_n, ALT_p_LB=ALT_p_LB, verbose=args.verbose)
 
    ### classification ... 'likelihood based' method
    print("\n{0}".format("[Typing]..."))
    reduceGT=True if args.dtype=='WES' else True
    fuzzy_dist=2 if args.dtype=='WES' and len(finalvar.index) > 3 else 0
    keep_SNV_GT=True if args.dtype=='WES' else False
    print("[GT reduce strategy] reduce GT: {0}; fuzzy_dist:{1}; keep GT with SNV:{2}".format(reduceGT, fuzzy_dist, keep_SNV_GT))
    SNPs, SNPm, SNP_GT, GT, GTm=BGClf.matrixprep.SNP4likelihood(testdata, traindata, args.gene, ALT_n_LB=args.alt_read_n, reduce_target=reduceGT, fuzzy_dist=fuzzy_dist, keep_SNV_GT=keep_SNV_GT)
  
    df_merged=None
    a_warn=False; b_warn=False; c_warn=False; d_warn=False; e_warn=False; f_warn=False
    exonCov, geneCov=BGClf.accessory.coverage(df=allvar)    

    print("\n{0}".format("[S1. Detecting large deletion]"))

    ### Normalization values option:
    ### use estimated entire coverage as denominator 
    #  bin_avg.wgsratio.log2 (WGS/WES)
    #  bin_med.wgsratio.log2 (WGS/WES)
    ### use all bins as denominator (including introns)
    #  bin_avg.regionratio.log2 (WGS)
    #  bin_med.regionratio.log2 (WGS)
    ### use exon bins as denominator
    #  exon_bin_avg.regionratio.log2 (WES/WGS)
    #  exon_bin_med.regionratio.log2 (WES/WGS)
    ############################3 


    normalization_mode='bin_med.wgsratio.log2' if args.dtype=='WGS' else 'bin_med.wgsratio.log2' 
    cov_mode='avg_QC_readn' if args.dtype=='WGS' else 'avg_QC_readn'
    binz=300 if args.dtype=='WGS' else 20
    step=150 if args.dtype=='WGS'else 10
    flanking=300 if args.dtype=='WGS' else 0
    exonOnly=False if args.dtype=='WGS' else True
    cbs_seg_winZ=10 if args.dtype=='WGS' else 3
    s1_perbaseCov_save=False if args.dtype=='WGS' else True

    print("[S1][normalization mode] {0} [bin size] {1} [step size] {2} [flanking] {3} [exon only] {4} [COV mode] {5}".format(normalization_mode, binz, step, flanking, exonOnly, cov_mode))
    print("[S1][gene1] {0}".format(args.gene))
    large_del, abnomral_exon_cov=BGClf.accessory.gene_exon_del(exonCov, args.bam, args.reference, gene, gene=args.gene, coverage=args.coverage, seqtype=args.dtype, bin=binz, step=step, exonOnly=exonOnly, flanking=flanking, cov_mode=cov_mode, save_output=s1_perbaseCov_save, prefix=args.prefix+'.S1.'+args.gene)
    if company is not None:
        print("[S1][gene2] {0}".format(company))
        exonCov_C, geneCov_C=BGClf.accessory.coverage(df=allvarC)
        large_del_C, abnomral_exon_cov_C=BGClf.accessory.gene_exon_del(exonCov_C, args.bam, args.reference, gene_c, gene=company, coverage=args.coverage, seqtype=args.dtype, bin=binz, step=step, exonOnly=exonOnly, flanking=flanking, cov_mode=cov_mode, save_output=s1_perbaseCov_save, prefix=args.prefix+'.S1.'+company)

    ### extra check for WES coverage:
    if args.dtype == 'WES':
        ### need to add another layer of Cov check for WES if hetero detection failed
        #BGClf.accessory.exon_cov_comparison(exonCov, exonCov_C)
        RHD_ML_fn=args.prefix+'.S1.'+'RHD'+'.binZ'+str(binz)+'.coverage.txt'        
        RHCE_ML_fn=args.prefix+'.S1.'+'RHCE'+'.binZ'+str(binz)+'.coverage.txt'         
        if os.path.exists(RHD_ML_fn) and os.path.exists(RHCE_ML_fn):
            if args.gene=='RHD':
                print("[S1][WES][ML][RHD][DEL] ... predicting RHD deletion status...")
                large_del, hybrid_status, res, prob, prob_class=BGClf.machcls.WES_RHD_DEL_CLF(args.prefix, RHD_ML_fn, RHCE_ML_fn)
                print("[S1][WES][ML][RHD][DEL] {0}; Prob {1}".format(res, max(prob)))    
            elif args.gene=='RHCE':
                print("[S1][WES][ML][RHD][DEL] ... predicting RHD deletion status...")
                large_del=None    
                large_del_C, hybrid_status_C, res_C, prob_C, prob_class_C=BGClf.machcls.WES_RHD_DEL_CLF(args.prefix, RHD_ML_fn, RHCE_ML_fn)   
                print("[S1][WES][ML][RHD][DEL] {0}; Prob {1}".format(res_C, max(prob_C))) 
        else:
            if not os.path.exists(RHD_ML_fn):
               print("[S1][WES][ML][DEL] file not exist:", RHD_ML_fn)
            if not os.path.exists(RHCE_ML_fn):
               print("[S1][WES][ML][DEL] file not exist:", RHCE_ML_fn)
            print("[S1][WES][ML] not run")
 
    print("[S1][large del] {0}".format(large_del))
 
    ##############################################
     ### test block: speed up the debug process ###
    ##############################################


    if large_del:
        new_data=[]
        if 'Homozygous deletion' in large_del:
            print('[S2] hom del')
            new_data.append({'Allele':'A1', 'Bloodtype':large_del, 'Likelihood':0, 'Index':0, "matchDBaa":"","not_matchDBaa":"","SampleExtraVar_aa":"","AllVariations_aa":""})
            new_data.append({'Allele':'A2', 'Bloodtype':large_del, 'Likelihood':0, 'Index':0, "matchDBaa":"","not_matchDBaa":"","SampleExtraVar_aa":"","AllVariations_aa":""})
            df_merged=pd.concat([pd.DataFrame(new_data)], ignore_index=True, sort=True)
        elif 'Heterozygous deletion' in large_del:
            print('[S2] het del')
            if args.dtype=='WES' and len(GT)>15:
                 ### reduce the # of GT for testing
                 GT, GTm=BGClf.matrixprep.reduce_targetGT(SNPm, GTm, fuzzy_dist=0, keep_SNV_GT=keep_SNV_GT)

            allele_classified=BGClf.classification.LH(args.bam, args.gene, allvar, GTm, dtype=args.dtype, gbuild=args.genome_build, allele_type='hom',ALT_n_LB=args.alt_read_n).calculateLL()
            #print(allele_classified)
            allele_classified=allele_classified.loc[(allele_classified['Allele']=='A2') & (allele_classified['Likelihood']==allele_classified['Likelihood'].max())]
            allele_indicies=allele_classified['Index'].unique()
            for allele_idx in allele_indicies:
                new_data.append({'Allele':'A1', 'Bloodtype':large_del, 'Likelihood':0, 'Index':allele_idx, "matchDBaa":"","not_matchDBaa":"","SampleExtraVar_aa":"","AllVariations_aa":""})
                df_merged=pd.concat([pd.DataFrame(new_data), allele_classified], ignore_index=True, sort=True)
        df_merged_raw=df_merged.copy()
        df_merged_raw['Comment']=""
    else:
        SNPm_test=SNPm.copy()
        llmode=None
        Var1136=""
        Var48=""
        if SNPm_test.shape[0] > 0:
             ### check RHD 1136 first
             if args.gene=='RHD' and '1136_C_T' in SNPm_test.index:
                 print("[RHD1136] checking")
                 #   The RHD 1136 has relatveily low coverage at this position due to exact duplication of exon 8.  It was adjusted to account for that by considering only heterozygrous combination of alleles
                 if (SNPm_test.loc['1136_C_T','alt_cnt']/args.coverage > args.heterozyoteCov) and (SNPm_test.loc['1136_C_T','alt_cnt']/args.coverage < (1-args.heterozyoteCov)):
                     print('[RHD1136] mode changed to hetero RHD 1136_C_T. 1136C/T count: {0}; 1136_cov ratio: {1}; {2} coverage: {3}; VAF cutoff: {4}'.format(SNPm_test.loc['1136_C_T','alt_cnt'], SNPm_test.loc['1136_C_T','alt_cnt']/args.coverage, args.dtype, args.coverage, args.heterozyoteCov))
                     llmode="het RHD 1136_C_T"
                     Var1136="het RHD 1136_C_T"
                 else:
                     print('[RHD1136] mode changed to hom RHD 1136_C_T. 1136C/T count: {0}; 1136_cov ratio:{1}; {2} coverage: {3}; VAF cutoff: {4}'.format(SNPm_test.loc['1136_C_T','alt_cnt'], SNPm_test.loc['1136_C_T','alt_cnt']/args.coverage, args.dtype,args.coverage, args.heterozyoteCov) )
                     llmode="hom RHD 1136_C_T"
                     Var1136="hom RHD 1136_C_T"
                 SNPm_test=SNPm_test.drop('1136_C_T')
             ### check RHCE 48 in WES
             if args.dtype in ['WES'] and args.gene=='RHCE' and '48_G_C' in SNPm_test.index: 
                 print("[RHCE48] checking")
                 cutoff=0.2
                 #   The RHCE 48 has relatveily 50% coverage at this position due to exact duplication of exon 1 except for 48C.  It was adjusted to account for that by considering only heterozygrous combination of alleles if 48G/C is present
                 if ((SNPm_test.loc['48_G_C','alt_cnt']/args.coverage > cutoff) and (SNPm_test.loc['48_G_C','alt_cnt']/args.coverage < (1-cutoff))): #and SNPm_test.loc['48_G_C','ref_cnt'] > 0:
                     #((SNPm_test.loc['48_G_C','alt_cnt']/args.coverage > exonCov_C[cov_mode].mean()) and (SNPm_test.loc['48_G_C','alt_cnt']/args.coverage < (1-exonCov_C[cov_mode].mean()))):
                     print('[RHCE48] mode changed to hetero RHCE 48_G_C. 48G/C count: {0}; 48_cov ratio:{1}; {2} coverage: {3}; VAF cutoff: {4}'.format(SNPm_test.loc['48_G_C','alt_cnt'], round(SNPm_test.loc['48_G_C','alt_cnt']/args.coverage,3), args.dtype,args.coverage, cutoff) )
                     llmode="het RHCE 48_G_C"
                     Var48="het RHCE 48_G_C"
                 else:
                     print('[RHCE48] mode changed to hom RHCE 48_G_C. 48G/C count: {0}; 48_cov ratio:{1}; {2} coverage: {3}; VAF cutoff: {4}'.format(SNPm_test.loc['48_G_C','alt_cnt'], round(SNPm_test.loc['48_G_C','alt_cnt']/args.coverage,3), args.dtype,args.coverage, cutoff) )
                     llmode="hom RHCE 48_G_C"
                     Var48="hom RHCE 48_G_C"
                 SNPm_test=SNPm_test.drop('48_G_C')

        ### extra check for WES ML var1136 and var48:
        if args.dtype == 'WES':
            ### need to add another layer of Cov check for WES if hetero detection failed
            #BGClf.accessory.exon_cov_comparison(exonCov, exonCov_C)
            RHD_ML_fn  =args.prefix+'.S1.'+'RHD'+'.binZ'+str(binz)+'.coverage.txt'
            RHCE_ML_fn =args.prefix+'.S1.'+'RHCE'+'.binZ'+str(binz)+'.coverage.txt'
            RHD_VAR_fn =args.prefix+'.RHD.full.variant.txt'
            RHCE_VAR_fn=args.prefix+'.RHCE.full.variant.txt' 
            if os.path.exists(RHD_ML_fn) and os.path.exists(RHCE_ML_fn):
                if args.gene=='RHD':
                    if os.path.exists(RHD_VAR_fn):
                        print("[S1][WES][ML][RHD][1136] ... predicting RHD 1136 status...")
                        Var1136, res, prob, prob_class=BGClf.machcls.WES_RHD_1136_CLF(args.prefix, RHD_ML_fn, RHCE_ML_fn, VAR_fn=RHD_VAR_fn)
                        llmode=Var1136
                        print("[S1][WES][ML][RHD][1136] {0}; Prob {1}".format(Var1136, max(prob)))
                    else:
                        print("[S1][WES][ML][RHD][1136] file not exist:", RHD_VAR_fn)
                if args.gene=='RHCE':
                    if os.path.exists(RHCE_VAR_fn):
                        print("[S1][WES][ML][RHCE][48] ... predicting RHCE 48 status...")
                        Var48, res_C, prob_C, prob_class_C=BGClf.machcls.WES_RHCE_48_CLF(args.prefix, RHD_ML_fn, RHCE_ML_fn, VAR_fn=RHCE_VAR_fn)
                        llmode=Var48
                        print("[S1][WES][ML][RHCE][48] {0}; Prob {1}".format(res_C, max(prob_C)))
                    else:
                        print("[S1][WES][ML][RHCE][48] file not exist:", RHCE_VAR_fn)
            else:
                if not os.path.exists(RHD_ML_fn):
                   print("[S1][WES][ML][VAR] file not exist:", RHD_ML_fn)
                if not os.path.exists(RHCE_ML_fn):
                   print("[S1][WES][ML][VAR] file not exist:", RHCE_ML_fn)
                print("[S1][WES][ML][VAR] not run")

        #print("TEST")
        #print(GTm, SNPm_test)
        llmode2='hom/het'
        if SNPm_test.shape[0] > 0:
             if GTm.shape[1] < SNPm_test.shape[0]:
                 ### this criteria is used to solve the var not in known database combination
                 llmode2='hom/het'
             elif any(SNPm_test['allele_frac'] > args.heterozyoteCov) & any(SNPm_test['allele_frac'] < (1 - args.heterozyoteCov)):
                 llmode2='het'
             elif args.gene=='RHD' and (any(((SNPm_test['alt_cnt']/args.coverage) > args.heterozyoteCov)) or any(((SNPm_test['alt_cnt']/args.coverage) < (1-args.heterozyoteCov)))):
                 ### this criteria is used to solve the problematic region on RHD
                 llmode2='het'
             else:
                 llmode2='hom' if (SNPm_test.shape[0]==1) & (SNPm_test['allele_frac'].max()==1) else 'hom/het'
        #else:
        #     if llmode is None: llmode='hom'
        if llmode is None or llmode=='':
            llmode=llmode2

        print('[TEST] llmode:',llmode, '; llmode2:', llmode2)
        if ('het' in llmode and llmode != 'hom/het' ) or ('het' in llmode2 and llmode2!='hom/het'):
            llmode='het'


        print ('\n[Zygosity] "{0}"'.format(llmode))

        if args.gene=='RHD': print ('[VAR1136] "{0}"'.format(Var1136))
        #print ('[SNP]', SNPm, '\n')



        print ("\n{0}".format("[S2. Calculating likelihood]"))

        problem_snv_type=Var1136 if args.gene=='RHD' else Var48

        if llmode in ['het RHD 1136_C_T','het', 'het RHCE 48_G_C']:
            if args.dtype=='WES':
                if len(GT)>15: fuzzy_dist=0
                GT, GTm=BGClf.matrixprep.reduce_targetGT(SNPm, GTm, fuzzy_dist=fuzzy_dist, keep_SNV_GT=keep_SNV_GT, varmode=Var1136)

            df_merged=BGClf.classification.LH(args.bam, args.gene, allvar, GTm, dtype=args.dtype, gbuild=args.genome_build, allele_type='het', problem_snv_type=problem_snv_type, ALT_n_LB=args.alt_read_n).calculateLL()
            df_merged_raw=df_merged.copy()
            print("het")
            #print(df_merged)
            df_merged=df_merged.loc[df_merged['Likelihood']==df_merged['Likelihood'].max()]
            #print(df_merged)
        elif llmode in ['hom', 'hom RHD 1136_C_T','hom RHCE 48_G_C']:
            if args.dtype=='WES':
                if len(GT)>15: fuzzy_dist=0
                GT, GTm=BGClf.matrixprep.reduce_targetGT(SNPm, GTm, fuzzy_dist=fuzzy_dist, keep_SNV_GT=keep_SNV_GT, varmode=Var48)
            df_merged=BGClf.classification.LH(args.bam, args.gene, allvar, GTm, dtype=args.dtype, gbuild=args.genome_build, allele_type='hom', problem_snv_type=problem_snv_type, ALT_n_LB=args.alt_read_n).calculateLL()
            df_merged_raw=df_merged.copy()
            print("hom")
            #print(df_merged)
            df_merged=df_merged.loc[df_merged['Likelihood']==df_merged['Likelihood'].max()]
            #print(df_merged)
        elif llmode in ['hom/het']:
            if args.dtype=='WES':
                 if len(GT)>15 :
                     ### reduce the # of GT for testing
                     GT, GTm=BGClf.matrixprep.reduce_targetGT(SNPm, GTm, fuzzy_dist=0, keep_SNV_GT=keep_SNV_GT)
            df_merged=BGClf.classification.LH(args.bam, args.gene, allvar, GTm, dtype=args.dtype, gbuild=args.genome_build, ALT_n_LB=args.alt_read_n, problem_snv_type=problem_snv_type).calculateLL()
            df_merged_raw=df_merged.copy()
            print("hom/het")
            #print(df_merged)

        df_merged_raw['Comment']=""
        df_merged_raw, a_warn=BGClf.accessory.filter_annotation(df_merged_raw, df_merged, annotation='a', LH_adj=0)

        #if args.dtype=='WES':
        #    args.alt_read_n=args.wes_alt_read_n



        print("[S2.1 RHD 1136] checking 1136 C/T var")
        print(df_merged)
        ### remove the combination with both alleles contating 1136C_T
        if args.gene=='RHD': 
            df_merged=BGClf.accessory.RHD1136check(df_merged, untrimmed_traindata, Var1136)
            df_merged_raw, b_warn=BGClf.accessory.filter_annotation(df_merged_raw, df_merged, annotation='b')
        #print("after 1136 check")
        #print(df_merged)

        if args.dtype in ['WES'] and args.gene=='RHCE':
            print("[S2.1 RHCE 48] checking 48 G/C var")
            ### remove the combination with both alleles contating 48G_C
            df_merged=BGClf.accessory.RHCE48check(df_merged, untrimmed_traindata, Var48)
            df_merged_raw, b_warn=BGClf.accessory.filter_annotation(df_merged_raw, df_merged, annotation='b')


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

    #print('[S2] debugging .....')
    #print(final)

    ### other rescue calling for specific type
    print("\n{0}".format("[S3. Exon CNV segmentaion]..."))
    print("[S3][normalization mode] {0} [bin size] {1} [step size] {2} [flanking] {3} [exon only] {4}".format(normalization_mode, binz, step, flanking, exonOnly)) 
    
    if args.gene=='RHD' or args.gene=='RHCE':
        if args.gene=='RHD':
            RHDhybrid=False
            if args.dtype=='WES' and hybrid_status and hybrid_status=='Hybrid':
                 RHDhybrid=True
            rhd, rhce, rhCov, rhdBreakS, rhdBreakE, rhceBreakS, rhceBreakE=BGClf.accessory.RHD_RHCE_segment_CNV(gene, gene_c, exonCov, exonCov_C, args.bam, args.coverage, args.reference, prefix, \
                                                                                                                verbose=args.verbose, seg_mode=normalization_mode, bin=binz, step=step, flanking=flanking, exonOnly=exonOnly, cbs_seg_winZ=cbs_seg_winZ, \
                                                                                                                RHDhybrid=RHDhybrid)
        elif args.gene=='RHCE':
            rhd, rhce, rhCov, rhdBreakS, rhdBreakE, rhceBreakS, rhceBreakE=BGClf.accessory.RHD_RHCE_segment_CNV(gene_c, gene, exonCov_C, exonCov, args.bam, args.coverage, args.reference, prefix, \
                                                                                                                verbose=args.verbose, seg_mode=normalization_mode, bin=binz, step=step, flanking=flanking, exonOnly=exonOnly, cbs_seg_winZ=cbs_seg_winZ)
        if args.gene=='RHCE':
            print("{0}".format("[S3. RHCE Relinking]..."))
            Ce_supportP=0.2 if args.dtype=='WES' else 0.1
            Ce_supportN=3
            final=BGClf.accessory.RHCE_relink(args, final, finalvar, rhCov, gbuild=args.genome_build, Ce_supportP=Ce_supportP, binz=binz, Ce_supportN=Ce_supportN)
            df_merged_raw, c_warn=BGClf.accessory.filter_annotation(df_merged_raw, final, annotation='c', LH_adj=0)


    #print(final)
    #print("[RHtyper] finalvar:", finalvar)


    ### output matched and non-matched mutation for prediced alleles
    print("\n{0}".format("[S4. Identifying matched/non-matched vars in database]..."))
    final=BGClf.accessory.var_compare(final, rawDB, finalvar, finalvarC, allvar, args.gene, args.bam, args.reference, rawDB)   

    outflds=['ID','Allele','Index','Bloodtype','Alias', 'Likelihood','matchDB','matchDBn','matchDBaa', 'not_matchDB','not_matchDBn','not_matchDBaa','SampleExtraVar','SampleExtraVar_aa','AllVariations','AllVariations_aa']
    final=final[outflds]

    print(final)  

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
    



    #print(final)                

    #'D(1)CE(2-9)D(10)':'RHD*01N.03',
    #'D(1-2)CE(3-9)D(10)':'RHD*01N.04',
    #'D(1)CE(2-7)D(8-10)':'RHD*01N.05',
    #'CE(1)D(2-6)CE(7-10)':'RHD*01N.42',
    #'CE(1-3)D(4-10)':'RHD*01N.43',
    #'CE(1-3)D(4-9)':'RHD*01N.43',
    #'CE(1)D(2)CE(3)D(4-10)':'RHD*01N.43',
    #'D(1-4)CE(5-7)D(8-10)':'RHD*01EL.23 RHD*DEL23',
    #'D(1-3)CE(4-9)D(10)':'RHD*01EL.44 RHD*DEL44'

    ### keep only the type with minimum mismatch DBn
    if len(final.Index.unique())> 1:
        xmDBn, maxDBn, mmDBn, minDBn, EV, minEVn, hybrid=BGClf.accessory.max_match(final, check_hybrid=True)
        is_hybrid=False
        print('Debugging')
        print(final)
        if rhd != '':
            ### need better rules to control the def of hybrid here
            if 'CE(4-7)' in rhd or 'CE(5-7)' in rhd or 'CE(2-9)' in rhd or 'CE(3-9)' in rhd or 'CE(2-7)' in rhd or 'CE(7-10)' in rhd or 'CE(1-3)' in rhd or 'D(4-10)' in rhd or 'CE(3-9)' in rhd:
                ### force all 4-7 to be hybrid when the string has 4-7
                is_hybrid=True
            elif not rhd.startswith('CE'):
                is_hybrid=True
        final_s5_tmp=final.copy()

        print('[S5][info] RHD coverage profile:{}; is hybrid?:{}'.format(rhd, is_hybrid))
        print('[S5][info] RHCE coverage profile:{}; is hybrid?:{}'.format(rhce, is_hybrid))

        if is_hybrid: 
            final=final.loc[final.Index.isin(hybrid)]
        else:
            final=final.loc[~final.Index.isin(hybrid)]

        print('Debugging ................... 2')
        print(final)

        if len(final.Index.unique()) == 0:
            print('[S5][warning] no allele pairs left after Hybrid correction ... using unfiltred allele pairs')
            final=final_s5_tmp.copy()

        if len(final.Index.unique()) > 1:
            print("{0}".format("[S5] keep best match allele pairs]...init allele_pair#:{0}".format(len(final.Index.unique()))))

        if len(final.Index.unique()) > 1:
            xmDBn, maxDBn, mmDBn, minDBn, EV, minEVn, hybrid=BGClf.accessory.max_match(final)
            final=final.loc[final.Index.isin(mmDBn[minDBn])]
            print("{0}".format("... [keep best match allele pairs]...by minimize variant mismatch ... allele_pair#:{0}".format(len(final.Index.unique()))))
        if len(final.Index.unique()) > 1:
            xmDBn, maxDBn, mmDBn, minDBn, EV, minEVn, hybrid=BGClf.accessory.max_match(final)
            final=final.loc[final.Index.isin(xmDBn[maxDBn])]
            print("{0}".format("... [keep best match allele pairs]...by maximum variant match ... allele_pair#:{0}".format(len(final.Index.unique()))))
        if len(final.Index.unique()) > 1:
            xmDBn, maxDBn, mmDBn, minDBn, EV, minEVn, hybrid=BGClf.accessory.max_match(final)
            final=final.loc[final.Index.isin(EV[minEVn])]
            print("{0}".format("... [keep best match allele pairs]...by minimize extra variants ... allele_pair#:{0}".format(len(final.Index.unique()))))
    #print(final)
    df_merged_raw, d_warn=BGClf.accessory.filter_annotation(df_merged_raw, final, annotation='d')

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


    ### remove underscore
    df_merged_raw.Bloodtype=df_merged_raw.Bloodtype.str.replace('_',' ')
    final.Bloodtype=final.Bloodtype.str.replace('_',' ')


    if args.verbose >3: 
        print("[Typed result]")
        #print(df_merged_raw)
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

