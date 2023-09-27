"""
   Module for ML classification method
"""

import pandas as pd
import re, os, pickle

package_directory = os.path.dirname(os.path.abspath(__file__))


def WES_RHD_DEL_CLF(sample_id, RHDcov_fn, RHCEcov_fn):
    '''
        Model was trained by XGboost using WGS prediction as ground truth
        Hete	Hete,Hybrid	Homo	Hybrid	Normal
    '''

    model_fn=os.path.join(package_directory, 'database', 'RHDDEL.BorutaFeature.xgb.model') 
    model = pickle.load(open(model_fn, 'rb'))
    res, prob, prob_class=WES_RH_prediction(sample_id, RHDcov_fn, RHCEcov_fn, model)
    hybrid='NotHybrid'
    typing=None
    if res == 'Homo':
        typing='RHD negative (Homozygous deletion)' 
    elif res=='Hete, Hybrid':
        typing='RHD negative (Heterozygous deletion) or low coverage'
        hybrid='Hybrid'
    elif res=='Hete':
        typing='RHD negative (Heterozygous deletion) or low coverage'
    elif res=='Hybrid':
        typing=None
        hybrid='Hybrid'
    elif res=='Normal':
        typing=None
    return typing, hybrid, res, prob, prob_class 


def WES_RHCE_UCLC_CLF(sample_id, RHDcov_fn, RHCEcov_fn, VAR_fn):
    '''
        Model was trained by XGboost using WGS prediction as ground truth
    '''

    model_fn=os.path.join(package_directory, 'database', 'RHCEUCLC.BorutaFeature.xgb.model')
    model = pickle.load(open(model_fn, 'rb'))
    res, prob, prob_class=WES_RH_prediction(sample_id, RHDcov_fn, RHCEcov_fn, model, VAR_fn=VAR_fn)
    typing=res   

    return typing, res, prob, prob_class

def WES_RHD_1136_CLF(sample_id, RHDcov_fn, RHCEcov_fn, VAR_fn):
    '''
        Model was trained by XGboost using WGS prediction as ground truth
    '''

    model_fn=os.path.join(package_directory, 'database', 'RHD1136.BorutaFeature.xgb.model')
    model = pickle.load(open(model_fn, 'rb'))
    res, prob, prob_class=WES_RH_prediction(sample_id, RHDcov_fn, RHCEcov_fn, model, VAR_fn=VAR_fn, cpos2check=1136, ref='C', alt='T')
    typing=""
    if res == 'hom':
        typing="hom RHD 1136_C_T"
    elif res == 'het':
        typing="het RHD 1136_C_T"
    return typing, res, prob, prob_class

def WES_RHCE_48_CLF(sample_id, RHDcov_fn, RHCEcov_fn, VAR_fn):
    '''
        Model was trained by XGboost using WGS prediction as ground truth
    '''

    model_fn=os.path.join(package_directory, 'database', 'RHCE48.BorutaFeature.xgb.model')
    model = pickle.load(open(model_fn, 'rb'))
    res, prob, prob_class=WES_RH_prediction(sample_id, RHDcov_fn, RHCEcov_fn, model, VAR_fn=VAR_fn, cpos2check=48, ref='G', alt='C')
    typing=""
    if res == 'hom':
        typing='hom RHCE 48_G_C'
    elif res == 'het':
        typing='het RHCE 48_G_C'
    return typing, res, prob, prob_class



def WES_RH_prediction(sample_id, RHDcov_fn, RHCEcov_fn, rfmodel,feature_fn="position", VAR_fn=None, cpos2check=48, ref='G', alt='C'):
    '''
        predict sample type by coverage
        rfmodel: RF model loaded from pickle
    '''

    #logging.debug('S1 load test dataset')
    RHCEcov=pd.read_csv(RHCEcov_fn, header=0, sep="\t")
    RHDcov=pd.read_csv(RHDcov_fn, header=0, sep="\t")
    RHcov=pd.concat([RHDcov, RHCEcov])
    #RHcov.prefix=RHcov.prefix.apply(clean_id)

    RHcov['sample']=sample_id
    RHcov['ratio']=RHcov['QCtotalReadN']/RHcov['all_cov']


    #logging.debug('S3 getting selected feature')
    #logging.info('Booster feature names')
    #logging.info(rfmodel.get_booster().feature_names)
    pos=list(rfmodel.get_booster().feature_names)

    if 'ref' in pos: pos.remove('ref')
    if 'alt' in pos: pos.remove('alt')

    pos=[ int(x) for x in pos ]

    keep=RHcov.position.isin(pos)
    selectedRHcov=RHcov[keep]

    #logging.info('Selected feature names')
    #logging.info(selectedRHcov.columns)
    df2=selectedRHcov.pivot_table(index=['sample'], columns='position', values='ratio')

    if VAR_fn is not None and 'ref' in list(rfmodel.get_booster().feature_names) and 'alt' in list(rfmodel.get_booster().feature_names):
        '''
            add var info
        '''
        VAR=pd.read_csv(VAR_fn, header=0, sep="\t")
        ref_c=VAR.loc[VAR.cpos==cpos2check,ref].values[0]
        alt_c=VAR.loc[VAR.cpos==cpos2check,alt].values[0]
        #logging.debug("cpos2check:{0} ... REF:{1}, REFcount:{2}, ALT:{3}, ALTcount:{4}".format(cpos2check, ref, ref_c, alt, alt_c))
        df2['ref']=ref_c
        df2['alt']=alt_c
        #print(df2.columns)

    

    #logging.debug('S4 prediction')
    res =  rfmodel.predict(df2)
    res_prob=rfmodel.predict_proba(df2)


    ### single sample returns
    res=res[0]
    res_prob=res_prob[0]
    res_class=rfmodel.classes_


    return res, res_prob, res_class





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


