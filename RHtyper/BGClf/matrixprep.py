"""
   Module to produce matrix for classification
"""
import pandas as pd

def create(training, input, cov, gene='RHD', QC_n_LB=None, ALT_n_LB=None, ALT_p_LB=None, fillin=True, verbose=0):
    trimmed_training=training.copy()
    trimmed_test    =[]
    untrimmed_test  =[]
    n=0
    n2=0
    for i in training.index.values:
        pos, o_ref, alt=i.split('_')
        n+=1
        ### temporary solution for insertion
        ### use the reference back to refbase
        ref=o_ref
        if o_ref == '-' :
            if int(pos) in input['cpos']:
                ref=input.loc[((input['gene']==gene) & (input['cpos']==int(pos)))]['REF'].values[0]


        if ((input['gene']==gene) & (input['cpos']==int(pos)) & (input['REF']==ref)).any():
            if input.loc[((input['gene']==gene) & (input['REF']==ref) & (input['cpos']==int(pos)))].shape[0] > 0:
                TOTAL_cnt     =input.loc[((input['gene']==gene) & (input['REF']==ref) & (input['cpos']==int(pos)))]['TOTAL_n'].values[0]
                QC_cnt        =input.loc[((input['gene']==gene) & (input['REF']==ref) & (input['cpos']==int(pos)))]['QC_n'].values[0]
                DEL_cnt       =input.loc[((input['gene']==gene) & (input['REF']==ref) & (input['cpos']==int(pos)))]['DEL_n'].values[0]
                INS_cnt       =input.loc[((input['gene']==gene) & (input['REF']==ref) & (input['cpos']==int(pos)))]['INS_n'].values[0]
                ref_cnt       =input.loc[((input['gene']==gene) & (input['REF']==ref) & (input['cpos']==int(pos)))][ref].values[0]
                alt_cnt       =-1
                alt_proportion=-1

                if alt=='-':
                    ### deletion
                    if input.loc[((input['gene']==gene) & (input['REF']==ref) & (input['cpos']==int(pos)) & (input['ALT']=='-'))].shape[0] > 0:
                        if DEL_cnt > 0:
                            ### point, short deletion
                            alt_cnt       =DEL_cnt
                            alt_proportion=input.loc[((input['gene']==gene) & (input['REF']==ref) & (input['cpos']==int(pos)))]['d' + ref +'p'].values[0]
                        else:
                            ### large deletion
                            alt_cnt       =(cov-QC_cnt)
                            alt_proportion=float(alt_cnt/cov)
                    else:
                        alt_cnt       =0
                        alt_proportion=0
                elif o_ref=='-':
                    ### insertion
                    if input.loc[((input['gene']==gene) & (input['REF']==ref) & (input['cpos']==int(pos)))].shape[0] > 0:
                        if 'i' + alt not in input:
                            alt_cnt       =0
                            alt_proportion=0
                        else:
                            alt_cnt       =input.loc[((input['gene']==gene) & (input['REF']==ref) & (input['cpos']==int(pos)))]['i' + alt ].values[0]
                            alt_proportion=input.loc[((input['gene']==gene) & (input['REF']==ref) & (input['cpos']==int(pos)))]['i' + alt +'p'].values[0]
                    else:
                        alt_cnt       =0
                        alt_proportion=0
                else:
                    alt_cnt       =input.loc[((input['gene']==gene) & (input['REF']==ref) & (input['cpos']==int(pos)))][alt].values[0]
                    alt_proportion=input.loc[((input['gene']==gene) & (input['REF']==ref) & (input['cpos']==int(pos)))][alt+'p'].values[0]


                ### asjust probability based on QC read n cutoff for SNV
                if ref != '-' and alt != '-':
                    if QC_n_LB and QC_cnt < QC_n_LB:
                        alt_proportion=0.0
                    if ALT_n_LB and alt_cnt < ALT_n_LB:
                        alt_proportion=0.0
                    if ALT_p_LB and alt_proportion < ALT_p_LB:
                        alt_proportion=0.0
                    if ALT_n_LB and alt_cnt >= ALT_n_LB and ref_cnt < ALT_n_LB:
                        if QC_cnt >= cov * 0.75:
                            ### adjust alt only if coverage is even in the base
                            alt_proportion=1.0
                        else:
                            if verbose > 3: print("Uneven coverage...not adjust ref allele ratio")
                #else:
                    ### put indel filtering here

                trimmed_test.append({'idx':i, 'allele_frac':alt_proportion, 'ref_cnt':ref_cnt, 'alt_cnt':alt_cnt})
            else:
                trimmed_test.append({'idx':i, 'allele_frac':0.0, 'ref_cnt':0, 'alt_cnt':0})
        else:
            if verbose > 5: print(n, gene, pos, ref, alt, "not exist")
            ### replenish trimmed_test for the missing ones, i.e. supposed to be reference ones
            if fillin:
                #print "add missing SNPs into testdata"
                trimmed_test.append({'idx':i, 'allele_frac':0.0, 'ref_cnt':0, 'alt_cnt':0})
            else:
                trimmed_training.drop(i, inplace=True)
    trimmed_test=pd.DataFrame(trimmed_test)
    trimmed_test.set_index('idx', inplace=True)
    ### create untrimmed filtered test data to rescue the SNPs that are not covered in the database
    for idx, row in input.iterrows():
        genen         =row['gene']
        gpos          =row['gpos']
        cpos          =row['cpos']
        TOTAL_cnt     =row['TOTAL_n']
        QC_cnt        =int(row['QC_n'])
        DEL_cnt       =row['DEL_n']
        INS_cnt       =row['INS_n']
        ref           =row['REF']
        alt           =row['ALT']
        ref_cnt       =row[ref]
        alt_cnt       =0
        alt_proportion=0
        if alt=='-':
            ### deletion
            if DEL_cnt > 0:
                ### point, short deletion
                alt_cnt       =DEL_cnt
                alt_proportion=row['d' + ref +'p']
            else:
                ### large deletion
                #print type(cov), type(QC_cnt)
                alt_cnt       =(cov-QC_cnt)
                alt_proportion=float(alt_cnt/cov)
        elif ref=='-':
            ### insertion
            alt_cnt       =row['i' + alt ]
            alt_proportion=row['i' + alt +'p']
        else:
            if alt != '.':
                alt_cnt       =row[alt]
                alt_proportion=row[alt+'p']

        ### asjust probability based on QC read n cutoff for SNV
        if ref != '-' and alt != '-':
            if QC_n_LB and QC_cnt < QC_n_LB:
                alt_proportion=0.0
            if ALT_n_LB and alt_cnt < ALT_n_LB:
                alt_proportion=0.0
            if ALT_p_LB and alt_proportion < ALT_p_LB:
                alt_proportion=0.0
            if ALT_n_LB and alt_cnt >= ALT_n_LB and ref_cnt < ALT_n_LB:
                if QC_cnt >= cov * 0.75:
                    ### adjust alt only if coverage is even in the base
                    alt_proportion=1.0
                #else:
                #    if verbose > 3: print "Uneven coverage...not adjust ref allele ratio"
        idx='_'.join(map(str,[cpos, ref, alt]))
        untrimmed_test.append({'gene':genen, 'gpos':gpos, 'cpos':cpos, 'idx':idx, 'allele_frac':alt_proportion, 'ref':ref, 'alt':alt, 'ref_cnt':ref_cnt, 'alt_cnt':alt_cnt})
    untrimmed_test=pd.DataFrame(untrimmed_test)
    untrimmed_test=untrimmed_test[untrimmed_test['allele_frac'] > 0 ]
    ### keep only columns (variants) with non-zero values and the reference column
    if not fillin:
        cols=list(trimmed_training.columns[trimmed_training.sum()>0].values)
    else:
        cols=list(trimmed_training.columns.values)

    if gene=='RHCE':
       #if gene+'@ce' not in cols: cols.append(gene+'@ce')
       if gene+'@RHCE*01_or_RHCE*ce_RHCE*c_RHCE*e' not in cols: cols.append(gene+'@RHCE*01_or_RHCE*ce_RHCE*c_RHCE*e')
    elif gene=='ABO':
       if gene+'@ABO*A1.01' not in cols: cols.append(gene+'@ABO*A1.01')
    elif gene=='RHD':
       #if gene+'@Reference' not in cols: cols.append(gene+'@Reference')
       if gene+'@RHD*01' not in cols: cols.append(gene+'@RHD*01')
    trimmed_training=trimmed_training[cols] ### keep only columns (variants) with non-zero values


    training_target=trimmed_training.columns.values
    weight=1/trimmed_training.gt(0).sum(axis=1)
    weight=weight.to_dict()
    #if verbose > 3:
    #    training.to_csv('training.txt', sep='\t', index=True, mode='w')
    #    trimmed_training.to_csv('trimmed_training.txt', sep='\t', index=True, mode='w')
    #    trimmed_test.to_csv('trimmed_test.txt', sep='\t', index=True, mode='w')
    #    untrimmed_test.to_csv('untrimmed_test.txt',sep='\t', index=True, mode='w')
    return trimmed_test.transpose(), untrimmed_test, trimmed_training.transpose(), training.transpose(), training_target, weight




def SNP4likelihood(testdata,traindata,gene,QC_n_LB=None, ALT_n_LB=None, ALT_p_LB=None, reduce_target=True, fuzzy_dist=0, keep_SNV_GT=False, skip_indel=False):
    '''
        provide the SNP candidates for calculating likelihood
    '''
    ttestdata=testdata.transpose()
    target_a=list(ttestdata.loc[ttestdata['allele_frac']>0].index)
    if ALT_n_LB is not None:
       target_a=list(ttestdata.loc[ttestdata['alt_cnt'] >= ALT_n_LB].index)

    
    target_am=ttestdata.loc[ttestdata['allele_frac']>0]
    if ALT_p_LB is not None:
        target_am=ttestdata.loc[ttestdata['allele_frac']>=ALT_p_LB]
    if ALT_n_LB is not None:
        target_am=target_am.loc[target_am['alt_cnt'] >= ALT_n_LB]

    target_train=[]
    
    for t in target_a:
        if gene=='RHCE':
            ### make sure the 676 variant is present for RHCE
            if '676_G_C' not in target_a:
                l=list(traindata.loc[(traindata[t]>0) & (traindata['676_G_C']==0)].index)
                target_train.extend(l)
            else:
                l=list(traindata.loc[traindata[t]>0].index)
                target_train.extend(l)
        else:
            l=list(traindata.loc[traindata[t]>0].index)
            target_train.extend(l)


    target_train=list(set(target_train))
    if gene=='RHCE':
       #if gene+'@ce' not in target_train: target_train.append(gene+'@ce')
       if gene+'@RHCE*01_or_RHCE*ce_RHCE*c_RHCE*e' not in target_train: target_train.append(gene+'@RHCE*01_or_RHCE*ce_RHCE*c_RHCE*e')
    elif gene=='ABO':
       if gene+'@ABO*A1.01' not in target_train: target_train.append(gene+'@ABO*A1.01')
    elif gene=='RHD':
       if gene+'@RHD*01' not in target_train: target_train.append(gene+'@RHD*01')

    target_train_mx=traindata.loc[target_train]
    target_train_mx=target_train_mx.loc[:, (target_train_mx != 0).any()] 

    if reduce_target:
        target_train, target_train_mx=reduce_targetGT(target_am, target_train_mx, fuzzy_dist=fuzzy_dist, keep_SNV_GT=keep_SNV_GT)


    target_t=list(target_train_mx)
    return target_a, target_am, target_t, target_train, target_train_mx

def reduce_targetGT(target_am, target_train_mx, fuzzy_dist=0, keep_SNV_GT=False, varmode='all'):
   '''
       reduce the target_train_mx to only the pairs with mismatch number < minimal number of mismatch with the target SNPs + fuzzy_dist
       if keep_SNV_GT: the GT contained SNVs will be kept, set for WES mainly
   '''
   ttarget_train_mx=target_train_mx.transpose()
   target_var=target_am.index.tolist()
   target_var_parameters=target_am.columns.tolist()
   train_var =target_train_mx.columns.tolist()
   train_var_genotypes=target_train_mx.index.tolist()

   dist_mx={}
   dist_min=9999999
   dist_min_withSNV=9999999
   type_with_SNV_ls=list()
   for i in range(0,(len(train_var_genotypes))):
       for j in range(i,len(train_var_genotypes)):
           matched_var=[]
           extra_var=[]
           type1=train_var_genotypes[i]
           type2=train_var_genotypes[j]
           type1_vars=list(ttarget_train_mx.loc[ttarget_train_mx[type1]>0, type1].index)
           type2_vars=list(ttarget_train_mx.loc[ttarget_train_mx[type2]>0, type2].index)
           type_with_SNV=False
           for var in target_var:
               if var in type1_vars or var in type2_vars:
                   if var not in matched_var:
                       matched_var.append(var)
                       if '-' not in var: type_with_SNV=True
           var1=False
           var2=False       
           for var in type1_vars:
               if var not in target_var:
                   if var not in extra_var: extra_var.append(var)
               if '1136_C_T' in varmode:
                   if var=='1136_C_T':
                       var1=True
                       break
               if '48_G_C' in varmode:
                   if var=='48_G_C':
                       var1=True
                       break
           for var in type2_vars:
               if var not in target_var:
                   if var not in extra_var: extra_var.append(var)
               if '1136_C_T' in varmode:   
                   if var=='1136_C_T':
                       var2=True
                       break
               if '48_G_C' in varmode:
                   if var=='48_G_C':
                       var2=True
                       break

           ### check var_mode het or hom
           if 'het' in varmode:
               if var1 and var2:
                   print("[reduce_targetGT] Var mode: {0}; GT {1} and {2} have homozyous zygosity of {0}... skipping".format(varmode, type1, type2))    
                   continue 
           if 'hom' in varmode:
               if not var1 or not var2:
                   print("[reduce_targetGT] Var mode: {0}; GT {1} and {2} have unmatched zygosity of {0}... skipping".format(varmode, type1, type2))
                   continue
             


           dist=len(target_var)-len(matched_var)+len(extra_var)
           problematic=['RHCE@RHCE*01_or_RHCE*ce_RHCE*c_RHCE*e', 'RHCE@RHCE*01.05.01_RHCE*ce.05.01','RHCE@RHCE*01.08_RHCE*ce.08']
           #if type1 in problematic and type2 in problematic:
           #    print("[GT1] {0}; [GT2] {1}; [DIST] {2}; [TARGET_VAR] {3}; [MATCHED_VAR] {4}; [EXTRA_VAR] {5}".format(type1, type2, dist, len(target_var), len(matched_var), len(extra_var)))
           #    print(type1_vars)
           #    print(type2_vars)
           #    print(target_var)
           if dist < dist_min: dist_min=dist
           if dist not in dist_mx: 
               dist_mx[dist]={}
           if type1 not in dist_mx[dist]: 
               dist_mx[dist][type1]=[]
           if type2 not in dist_mx[dist][type1]: 
               dist_mx[dist][type1].append(type2)
           if type_with_SNV:
               if type1 not in type_with_SNV_ls:
                   type_with_SNV_ls.append(type1) 
               if type2 not in type_with_SNV_ls:
                   type_with_SNV_ls.append(type2)  
               if dist < dist_min_withSNV: dist_min_withSNV=dist        


   target_gt=[]
   for d in dist_mx:
      for t1 in dist_mx[d]:
          for t2 in dist_mx[d][t1]:
              #print("[Reduce]", t1,t2,d)
              if d <= dist_min + fuzzy_dist:
                  #print("[reduce_targetGT][Kept]", t1, t2, d)
                  if t1 not in target_gt: target_gt.append(t1)
                  if t2 not in target_gt: target_gt.append(t2)
              if keep_SNV_GT and d <= dist_min_withSNV:
                  if t1 in type_with_SNV_ls or t2 in type_with_SNV_ls:
                      #print("[reduce_targetGT][KeptSNV_GT]", t1, t2, d)
                      if t1 not in target_gt: target_gt.append(t1)
                      if t2 not in target_gt: target_gt.append(t2)

   if len(target_gt)>0:
       target_train_mx=target_train_mx.loc[target_gt]
       target_train_mx=target_train_mx.loc[:, (target_train_mx != 0).any()]
   print("[reduce_targetGT] # of target_gt for likelihood calculation: {0}".format(len(target_gt)))
   return target_gt, target_train_mx 
