"""
   Module to for future implementation of classification method
"""

def _svaCls():
    weight=None
    df_merged=None
    exonCov, geneCov=BGClf.accessory.coverage(df=allvar)
    large_del, abnomral_exon_cov=BGClf.accessory.gene_exon_del(exonCov, gene=args.gene, coverage=args.coverage)
    if large_del:
        new_data=[]
        new_data.append({'Bloodtype':large_del, 'Rank':0})
        if 'Homozygous deletion' in large_del: 
            new_data.append({'Bloodtype':large_del, 'Rank':0})
        allele_classified=BGClf.classification.CLF(traindata, target, testdata, weight=weight).summary()
        df_merged=pd.concat([pd.DataFrame(new_data), allele_classified], ignore_index=True)
    else:
        allele1, allele2=BGClf.haplotype.build_hap(testdata, untrimmed_test, args.bam, allvar, gene=args.gene, verbose=args.verbose, ALT_n_LB=args.alt_read_n, phase_rescue=args.rescue_phase).build()
        allele1_classified=BGClf.classification.CLF(untrimmed_traindata, target, allele1, weight=weight).summary()
        allele2_classified=BGClf.classification.CLF(untrimmed_traindata, target, allele2, weight=weight).summary()
        df_merged=pd.concat([allele1_classified, allele2_classified])
       
    df_merged['ID'] = args.prefix
   
    if args.verbose >9: print("Merged_df\n", df_merged)
    if args.verbose >9: print("Abnormal coverage\n", abnomral_exon_cov)
    df_merged.sort_values(inplace=True, axis=0, by=['Rank'], ascending=True)

    final=df_merged.head(2)  
    final.reset_index(drop=True, inplace=True)
    

    ### other rescue calling for specific type
    if args.gene=='RHCE':
           support_n, support_ratio=BGClf.accessory.RHCe_insertion(args.bam, args.coverage)
           if support_n >= args.alt_read_n or support_ratio > 0.1:
               ### with one Ce allele
               for idx, row in final.iterrows():
                   if 'RHce_48C' in row['Bloodtype']:
                       ID= final.ix[idx, 'ID']
                       final.iloc[idx]='-'
                       final.ix[idx, 'Bloodtype']='RHCE@RHCe'
                       final.ix[idx, 'Rank']='0'
                       final.ix[idx, 'ID']=ID
                       if support_ratio < 0.85:
                           break
                       else:
                           continue
                   if ((row['Bloodtype']=='RHCE@RHce_ces') and (final.ix[idx+1,'Bloodtype']=='RHCE@ce')) or ((row['Bloodtype']=='RHCE@RHce') and (final.ix[idx+1,'Bloodtype']=='RHCE@ce_ces')):
                       if support_ratio < 0.85:
                           ID= final.ix[idx, 'ID']
                           final.iloc[idx]='-'
                           final.iloc[idx+1]='-'
                           final.ix[idx, 'Bloodtype']='RHCE@RHCe'; 
                           final.ix[idx+1, 'Bloodtype']='RHCE@RHce733G';
                           final['Rank']='0'
                           final['ID']=ID
                       else:
                           print("Both Ce?")
                       break
           else:
               ### without Ce allele
               if args.linking_bloodtype:
                   final=BGClf.accessory.RHD_RHCE_linked(final, args.linking_bloodtype)

    if args.verbose >3: print("Final\n", final)
    final.to_csv(args.prefix+'.bloodtyping.txt', sep='\t', index=False, mode='w', \
                     columns=['ID', 'Bloodtype', 'DecisionLinearSVC', 'DecisionSVClinear', 'SVC_Rank', 'ProbabilityMultinominalNaiveBayes', 'MNB_Rank', 'ProbabilityRandomForest', 'RF_Rank', 'Rank'])


