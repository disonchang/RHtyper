'''
Classification module
'''
import pandas as pd
import pysam, sys
#from sklearn import preprocessing, linear_model, svm, tree
#from sklearn.naive_bayes import GaussianNB, MultinomialNB
#from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from general import *
from variants import *
from decimal import *

class MLClassifier(object):
    """
    A simple ML classifier
    """
    def __init__ ( self, x_train):
        """Takes in the training dataset, a n_features * n_samples
        (e.g. wavebands * samples) array. Pre calculates a lot of terms
        """
        self.x_train = x_train
        self.n_bloodtypes = x_train.shape[1]
        self.n_features = x_train.shape[0]
        cov_m = np.cov ( x_train )
        self.i_cov_m = np.linalg.inv ( cov_m )
        self.det_cov = np.linalg.det ( cov_m )
        self.mu = x_train.mean ( axis=1 )
        self.constant_term = -0.5*np.log( self.det_cov ) - self.n_features*0.5*np.log(2.*np.pi)
        

    def classify ( self, x_test, threshold = None ):
        """Classifies (ie gives the probability of belonging to a 
        class defined by the `__init__` training set) for a number
        of test data vectors. These vectors are n_features*n_samples.
        If `threshold` is specified, it selects samples with a probability
        larger or equal to it's value
        """
        s = x_test.T - self.mu
        log_prob =  -  0.5*  (np.dot(s, self.i_cov_m) *  s).sum(1)
        prob = np.exp ( log_prob )

        if threshold is not None:
            return np.where ( prob >= threshold, True, False )
        else:
            return prob

class LH(object):
    def __init__(self, bam, gene, SNPm, GTm, gbuild='hg38', dtype='WGS', allele_type='all', problem_snv_type='NA', minimum_base_quality = 15, minimum_mapq = 10, minimum_read_quality = 15, ALT_n_LB=None, verbose=0, het_cov=0.1):
        '''
        Takes bam, target SNP and SNPmatrix for reverting back to genome pos
        '''
        self.verbose=verbose
        self.bam=bam
        self.GTm=GTm
        self.gene=gene
        self.snp=list(GTm)
        self.snp_matrix=SNPm
        self.selected_snp_matrix=pd.DataFrame()
        self.genotype_matrix=GTm.transpose()
        self.allele_type=allele_type ### hete, homo, all
        self.problem_snv_type=problem_snv_type
        self.het_cov=het_cov
        self.build=gbuild
        self.ALT_n_LB=ALT_n_LB
        self.dtype=dtype

        for s in self.snp:
            a=s.split('_')
            t=self.snp_matrix.loc[self.snp_matrix['cpos']==int(a[0])]
            self.selected_snp_matrix=self.selected_snp_matrix.append(t)
        if self.selected_snp_matrix.shape[0] > 0:
            self.strand        = self.selected_snp_matrix['strand'].iloc[0]
            self.snp_d         = self.selected_snp_matrix.loc[:,['gpos','cpos','REF','ALT']]
            self.span_chr      = self.selected_snp_matrix['chr'].iloc[0]
            self.span_snp_gpos = list(self.selected_snp_matrix['gpos'])
            self.span_snp_gpos.sort()
            self.span_start    = min(self.span_snp_gpos) - 1
            self.span_end      = max(self.span_snp_gpos)
            self.minimum_base_quality  = minimum_base_quality 
            self.minimum_mapq =minimum_mapq
            self.minimum_read_quality=minimum_read_quality
            self.ALT_n_LB = ALT_n_LB



    def geno_ll(self):
        '''
            Calculate geno likelihood of each base
            ref for error probability
            https://www.drive5.com/usearch/manual/quality_score.html
        '''
        if self.selected_snp_matrix.shape[0] == 0:
             return "REF"
        ll_out={}
        samfile = pysam.AlignmentFile(self.bam, mode="rb") #, ignore_truncation=False)
        ### option1: per read
        #for read in samfile.fetch(self.span_chr, self.span_start, self.span_end):
        #    qname=read.query_name
        #    qaln =read.query_alignment_sequence
        #    qqual=read.query_alignment_qualities
        #    qpos=read.get_reference_positions()
        #    qmap=zip(qpos, qaln, qqual)
        #    
        #    for idx in range(0,len(qmap)):
        #         gpos=qmap[idx][0]
        #         b  =qmap[idx][1]
        #         Q   =qmap[idx][2]
        #         e   = 10 ** (-Q/float(10))
        #         if gpos in self.span_snp_gpos:
        #             cpos=self.snp_d.loc[self.snp_d['gpos']==gpos,['cpos']].values[0][0]
        #             for i in range(0, self.genotype_matrix.shape[1]):
        #                for j in range(i, self.genotype_matrix.shape[1]):
        #                    type1=list(self.genotype_matrix)[i]
        #                    type2=list(self.genotype_matrix)[j]
        #                    for idx, row in self.genotype_matrix.iterrows():
        #                        GTcpos, GTref, GTalt=idx.split('_')
        #                        if int(cpos)==int(GTcpos):
        #                            GTf1=row[type1]
        #                            GTf2=row[type2]
        #                            GTa1=GTref if GTf1==0 else GTalt
        #                            GTa2=GTref if GTf2==0 else GTalt
        #                            ll=None
        #                            if GTa1 == GTa2 and GTa1 == b:
        #                                print 'case1'
        #                                ll=1-e
        #                            elif GTa1 != GTa2 and (GTa1 == b or GTa2 == b):
        #                                print 'case2'
        #                                ll=((1-e)/float(2))+(((e/float(3)))/float(2))
        #                            elif b != GTa1 and b != GTa2:
        #                                print 'case3'
        #                                ll=e/float(3)
        #                            print qname, gpos, cpos, b, Q, e, type1, GTf1, GTa1, type2, GTf2, GTa2, ll
                                         
                            
                            
                            
                     #ref= self.selected_snp_matrix.loc[self.selected_snp_matrix['gpos']==gpos,['REF']].values[0]
                     #if self.verbose > 0: 
                     #    if len(ref) > 1: print "[Warning geno_all] multiple reference base at the same pos %d" % gpos
                     #for refb in ref:
                     #    for idx2 in range(0, len(self.snp_d)):
                     #        if gpos==self.snp_d[idx2][3]:
                     #            print qname, gpos, refb, b, Q, E
        ### option2: per position            
        for column in samfile.pileup(self.span_chr, self.span_start, self.span_end):
            pos0based=column.pos
            pos1based=column.pos+1
            if pos1based in self.span_snp_gpos:
                #print "calculating genoLL for %d" % pos1based
                for pileupread in column.pileups:
                    read     = pileupread.alignment
                    gpos     = pos1based
                    qname    = pileupread.alignment.query_name
                    qpos     = pileupread.query_position

                    # skip reads with indels
                    #if pileupread.is_del: #or pileupread.indel < 0:
                    #    continue

                    if pileupread.alignment.mapping_quality < self.minimum_mapq:
                        continue
                    # skip reads with no mapq specified
                    elif read.mapping_quality == 255:
                        continue
                    # skip mean qscore of the read below threshold
                    elif mean(read.query_qualities) < self.minimum_read_quality:
                        continue
                    else:
                        # check for insertion
                        if pileupread.is_refskip:# or pileupread.indel > 0:
                            # skip reads with a base quality below threshold
                            if read.query_qualities[qpos] < self.minimum_base_quality:
                                continue
                        if qpos and read.query_qualities[qpos] < self.minimum_base_quality:
                            continue

                    ### skip misaligned reads
                    if pos1based in [25321796,25321822,25321905,25321928] or \
                       pos1based in [25648287,25648313,25648396,25648419]:
                        if reads_misaligned(pileupread, self.build):
                            continue


                    if pileupread.is_del:
                        b='-'
                        Q=40
                    else:
                        b        = read.query_sequence[qpos]
                        Q        = read.query_qualities[qpos]
                    e        = 10 ** (-Q/float(10))
                    cpos     = self.snp_d.loc[self.snp_d['gpos']==gpos,['cpos']].values[0][0]

                    if self.strand == '-':
                        b=complement(b)
                    if b in ['A','T','C','G']:
                         bp=self.snp_matrix.loc[self.snp_matrix['gpos']==gpos,b+'p'].values[0]
                         bn=self.snp_matrix.loc[self.snp_matrix['gpos']==gpos,b].values[0]
                         if self.het_cov is not None and (0 < bp < self.het_cov) and \
                            self.ALT_n_LB is not None and (0 < bn < self.ALT_n_LB):
                             print("[genoLL] Skip low freq reads:", qname, 'at cpos:', cpos, '; bp:', bp, '; bn:', bn)
                             continue
                    for i in range(0, self.genotype_matrix.shape[1]):
                        for j in range(i, self.genotype_matrix.shape[1]):
                            type1=list(self.genotype_matrix)[i]
                            type2=list(self.genotype_matrix)[j]
                            for idx, row in self.genotype_matrix.iterrows():
                                GTcpos, GTref, GTalt=idx.split('_')
                                
                                if int(cpos)==int(GTcpos):
                                    GTf1=row[type1]
                                    GTf2=row[type2]
                                    GTa1=GTref if GTf1==0 else GTalt
                                    GTa2=GTref if GTf2==0 else GTalt
                                    ll=None
                                     
                                    if GTa1 == GTa2 and GTa1 == b:
                                        #print 'case1'
                                        ll=1-e
                                    elif GTa1 != GTa2 and (GTa1 == b or GTa2 == b):
                                        #print 'case2'
                                        ll=((1-e)/float(2))+(((e/float(3)))/float(2))
                                    elif b != GTa1 and b != GTa2:
                                        #print 'case3'
                                        ll=e/float(3)
                                    else:
                                        print('ll not determined')
                                    if ll is not None and ll != 0:
                                        ### resolve the issue for 25317062 (HG38) for RHD that always have only 50% coverage
                                        llori=ll
                                        checkpos=25317062
                                        if self.build.upper()=='HG19':
                                            checkpos=25643553

                                        #if gpos==checkpos and self.allele_type in ['het','het RHD 1136_C_T']: ### major changes
                                        if gpos==checkpos and self.problem_snv_type in ['het RHD 1136_C_T']:
                                            #if type1 == "RHD@RHD_DAU-0" or type2 == "RHD@RHD_DAU-0":
                                               #ll=0.5 * Decimal(ll)
                                               if (GTa1 == GTa2 and GTa1 == b) or (GTa1 != GTa2 and (GTa1 == b or GTa2 == b)):
                                                   ll=((1-e)/float(2))+(((e/float(3)))/float(2))

                                        ### resolve the issue for 25420739 (HG38) for RHCE that have only 50% coverage in WES
                                        if self.dtype in ['WES']:
                                            checkpos=25420739
                                            if self.build.upper()=='HG19':
                                                checkpos=25747230
                                            if gpos==checkpos and self.problem_snv_type in ['het RHCE 48_G_C']:
                                               if (GTa1 == GTa2 and GTa1 == b) or (GTa1 != GTa2 and (GTa1 == b or GTa2 == b)):
                                                   ll=((1-e)/float(2))+(((e/float(3)))/float(2))
 

                                        ### the print section below is for debugging
                                        #if gpos==25317062: 
                                        #    if b != GTa1 and b != GTa2:
                                        #        print(qname, gpos, cpos, b, Q, e, type1, GTa1, type2, GTa2, llori, ll)
                                        #if cpos==835:
                                        #    if b != GTa1 and b != GTa2:
                                        #        print qname, gpos, cpos, b, Q, e, type1, GTa1, type2, GTa2, ll
                                        #if type1 == 'RHCE@RHCE*Ce_733G::RHce_VS::RHCE*Ce_733G_' or type2 == 'RHCE@RHCE*Ce_733G::RHce_VS::RHCE*Ce_733G_':
                                        #    if b != GTa1 and b != GTa2:
                                        #        print qname, gpos, cpos, b, Q, e, type1, GTa1, type2, GTa2, ll
                                        #if type1 == "RHD@RHD_DAU-3" and type2=="RHD@RHD_DAU-0":
                                        #    print qname, gpos, cpos, b, Q, e, type1, GTa1, type2, GTa2, ll
                                        #if type1 == "RHD@RHD_DAU-3" and type2=="RHD@Reference":
                                        #    print qname, gpos, cpos, b, Q, e, type1, GTa1, type2, GTa2, ll
                                        #if type2=='RHCE@RHCe_polypeptide;_Rh30':
                                        #    print qname, gpos, cpos, b, Q, e, type1, GTa1, type2, GTa2, ll
                                        #if type1=='RHD@RHD*XX.XXRHD(L390L)' and type2=='RHD@Reference':
                                        #    print(qname, gpos, cpos, b, Q, e, type1, GTa1, type2, GTa2, ll)
                                        #if type1=='RHD@RHD*01W.41.0.1RHD*weak_D_type_41.0.1' and type2=='RHD@Reference':
                                        #    print(qname, gpos, cpos, b, Q, e, type1, GTa1, type2, GTa2, ll)

                                        if type1 not in ll_out:
                                            ll_out[type1]={}
                                            ll_out[type1][type2]={}
                                            ll_out[type1][type2][gpos]={}
                                            ll_out[type1][type2][gpos]['rawLL']=[]
                                            ll_out[type1][type2][gpos]['rawLL'].append(ll)
                                            ll_out[type1][type2][gpos]['prodLL']=ll
                                        else:
                                            if type2 not in ll_out[type1]:
                                                ll_out[type1][type2]={}
                                                ll_out[type1][type2][gpos]={}
                                                ll_out[type1][type2][gpos]['rawLL']=[]
                                                ll_out[type1][type2][gpos]['rawLL'].append(ll)
                                                ll_out[type1][type2][gpos]['prodLL']=ll
                                            else:
                                                if gpos not in ll_out[type1][type2]:
                                                    ll_out[type1][type2][gpos]={}
                                                    ll_out[type1][type2][gpos]['rawLL']=[]
                                                    ll_out[type1][type2][gpos]['rawLL'].append(ll)
                                                    ll_out[type1][type2][gpos]['prodLL']=ll
                                                else:
                                                    ll_out[type1][type2][gpos]['rawLL'].append(ll)
                                                    ori=ll_out[type1][type2][gpos]['prodLL']
                                                    ll_out[type1][type2][gpos]['prodLL']=Decimal(ll_out[type1][type2][gpos]['prodLL']) * Decimal(ll)
                                                    if ll_out[type1][type2][gpos]['prodLL'] == 0:
                                                        print(qname, gpos, cpos, b, Q, e, type1, GTf1, GTa1, type2, GTf2, GTa2, ori, ll)
                                                        sys.exit("become 0 here")
                                    else:
                                        print("ll == 0 for %s at pos %d" %(read, gpos))
                                   
        '''
           calculate log geno likelihood
        '''
        geno_LL={}
        for BG1 in ll_out:
            geno_LL[BG1]={}
            for BG2 in ll_out[BG1]:
                LLgeno=0
                for pos in ll_out[BG1][BG2]:
                    ll=ll_out[BG1][BG2][pos]['prodLL']
                    if ll != 0:
                        logll=Decimal(ll).log10()
                        LLgeno+=logll
                    else:
                        print('%s\t%s\t%d\tll == 0, skip' %(BG1, BG2, pos))
                geno_LL[BG1][BG2]=LLgeno
                #print('genoLL', BG1, BG2, LLgeno)
        if not geno_LL: return "REF"
        return geno_LL
        
         
    def phase_ll(self):
        '''
            Calculate phase likelihood of each base
        '''
        if self.selected_snp_matrix.shape[0] == 0:
             return 'REF'
        ll_out={}
        samfile = pysam.AlignmentFile(self.bam, mode="rb")
        #print self.span_snp_gpos
        for column in samfile.pileup(self.span_chr, self.span_start, self.span_end):
            pos0based=column.pos
            pos1based=column.pos+1
            gpos     =pos1based
            if pos1based in self.span_snp_gpos: 
                #print "calculating genoLL for %d" % pos1based
                readn=0
                for pileupread in column.pileups:
                    read     = pileupread.alignment
                    # skip reads with indels
                    #if pileupread.is_del: #or pileupread.indel < 0:
                    #    continue
                    if pileupread.alignment.mapping_quality < self.minimum_mapq:
                        continue
                    # skip reads with no mapq specified
                    elif read.mapping_quality == 255:
                        continue
                    # skip mean qscore of the read below threshold
                    elif mean(read.query_qualities) < self.minimum_read_quality:
                        continue
                    else:
                        # check for insertion
                        if pileupread.is_refskip: #or pileupread.indel > 0:
                            # skip reads with a base quality below threshold
                            if read.query_qualities[qpos] < self.minimum_base_quality:
                                continue

                    ### skip misaligned reads
                    if pos1based in [25321796,25321822,25321905,25321928] or \
                       pos1based in [25648287,25648313,25648396,25648419]:
                        if reads_misaligned(pileupread, self.build):
                            continue


                    if self.span_snp_gpos.index(gpos) + 1 < len(self.span_snp_gpos):
                        gpos2    = self.span_snp_gpos[self.span_snp_gpos.index(gpos) + 1]
                        qname    = pileupread.alignment.query_name
                        qpos     = pileupread.query_position
                        qaln     = pileupread.alignment.query_sequence
                        qqual    = pileupread.alignment.query_qualities
                        qrpos_0  = read.get_reference_positions()
                        qrpos_1  = [ x+1 for x in qrpos_0 ]
                        qmap     = zip(qrpos_1, qaln, qqual)
                        
                        if pileupread.is_del:
                            b1='-'
                            Q1=42
                        else:
                            b1 = read.query_sequence[qpos]
                            Q1 = read.query_qualities[qpos]
                        b11= None; b2=None
                        Q11= None; Q2=None
                        for p in qmap:
                            if p[0]==gpos:
                                b11=p[1]
                                Q11=p[2]
                            if p[0]==gpos2:
                                b2=p[1]
                                Q2=p[2]
                        if b1=='-' and b11 is None:
                            b11='-'
                            Q11=42 

                        if b2 is not None:
                            if self.strand == '-':
                                b1=complement(b1)
                                b11=complement(b11)
                                b2=complement(b2)
                            #print gpos, b1, Q1, b11, Q11, gpos2, b2, Q2
                            if b11 == b1:
                                Qerr  = 0.01 ### out-of-phase error rate
                                cpos1 = int(self.snp_d.loc[self.snp_d['gpos']==gpos,['cpos']].values[0][0])
                                cpos2 = int(self.snp_d.loc[self.snp_d['gpos']==gpos2,['cpos']].values[0][0])

                                if b1 in ['A','T','C','G'] and b2 in ['A','T','C','G']:
                                    ### AF based
                                    if self.het_cov is not None:
                                        b1p=self.snp_matrix.loc[self.snp_matrix['gpos']==gpos,b1+'p'].values[0]
                                        b2p=self.snp_matrix.loc[self.snp_matrix['gpos']==gpos2,b2+'p'].values[0]
                                        if (0 < b2p < self.het_cov) and (0 < b1p < self.het_cov):
                                            print("[PhaseLL] Skip low freq reads:", qname, ' at cpos(1,2):', cpos1, cpos2)
                                            continue
                                    ### read N based filtering
                                    if self.ALT_n_LB is not None:
                                        b1n=self.snp_matrix.loc[self.snp_matrix['gpos']==gpos,b1].values[0]
                                        b2n=self.snp_matrix.loc[self.snp_matrix['gpos']==gpos2,b2].values[0]
                                        if (0 < b2n < self.ALT_n_LB) and (0 < b1n < self.ALT_n_LB):
                                            print("[PhaseLL] Skip low freq reads:", qname, ' at cpos(1,2):', cpos1, cpos2)
                                            continue


                                readn+=1
                                
                                for i in range(0, self.genotype_matrix.shape[1]):
                                    for j in range(i, self.genotype_matrix.shape[1]):
                                        type1=list(self.genotype_matrix)[i]
                                        type2=list(self.genotype_matrix)[j]
                                        GTa1={}
                                        GTa2={}
                                        for idx, row in self.genotype_matrix.iterrows():
                                            GTcpos, GTref, GTalt=idx.split('_')
                                            # pos1
                                            if int(cpos1)==int(GTcpos):
                                                GT1f1=float(row[type1])
                                                GT1f2=float(row[type2])
                                                GT1a1=GTref if GT1f1==0 else GTalt
                                                GT1a2=GTref if GT1f2==0 else GTalt
                                                if gpos not in GTa1:
                                                    GTa1['g1']=[]
                                                    GTa1['g1'].append(GT1a1)
                                                else:
                                                    GTa1['g1'].append(GT1a1)
                                                if gpos not in GTa2:
                                                    GTa2['g2']=[]
                                                    GTa2['g2'].append(GT1a2)
                                                else:
                                                    GTa2[gpos].append(GT1a2)
                                            # pos2
                                            if int(cpos2)==int(GTcpos):
                                                GT2f1=float(row[type1])
                                                GT2f2=float(row[type2])
                                                GT2a1=GTref if GT2f1==0 else GTalt
                                                GT2a2=GTref if GT2f2==0 else GTalt
                                                if gpos2 not in GTa1:
                                                    GTa1['g1_next']=[]
                                                    GTa1['g1_next'].append(GT2a1)
                                                else:
                                                    GTa1.append(GT1a1)
                                                if gpos2 not in GTa2:
                                                    GTa2['g2_next']=[]
                                                    GTa2['g2_next'].append(GT2a2)
                                                else:
                                                    GTa2[gpos2].append(GT2a2)
                                        if len(GTa1['g1']) > 1 or len(GTa1['g1_next']) > 1 or len(GTa2['g2']) > 1 or len(GTa2['g2_next']) > 1:
                                            print('allele1', type1, 'g1=', gpos, GTa1['g1'], 'g1_next=', gpos2, GTa1['g1_next'])
                                            print('allele2', type2, 'g2=', gpos, GTa2['g2'], 'g2_next=', gpos2, GTa2['g2_next'])
                                        ### phase likelihood
                                        phase_ll=None
                                        phase_type=None
                                        if (GTa1['g1'][0] == b1 and GTa1['g1_next'][0]==b2) and (GTa2['g2'][0]==b1 and GTa2['g2_next'][0]==b2):
                                            phase_ll=1-Qerr
                                            phase_type=1
                                        elif (GTa1['g1'][0] != GTa2['g2'][0] or GTa1['g1_next'][0] != GTa2['g2_next'][0]) and ((GTa1['g1'][0] == b1 and GTa1['g1_next'][0]==b2) or (GTa2['g2'][0]==b1 and GTa2['g2_next'][0]==b2)):
                                            phase_ll=((1-Qerr)/2)+((Qerr/15)/2)
                                            phase_type=2
                                            
                                        elif ((GTa1['g1'][0] != b1 or GTa1['g1_next'][0] != b2) and (GTa2['g2'][0] != b1 or GTa2['g2_next'][0] != b2)):
                                            phase_ll=Qerr/15
                                            phase_type=3
                                        else:
                                            print('phase ll not determined')
                                            print(gpos, b1, Q1, b11, Q11, gpos2, b2, Q2)
                                            print('allele1', type1, 'g1=', gpos, GTa1['g1'], 'g1_next=', gpos2, GTa1['g1_next'])
                                            print('allele2', type2, 'g2=', gpos, GTa2['g2'], 'g2_next=', gpos2, GTa2['g2_next'])
 
                                        ll=phase_ll
                                        
                                        if ll is not None and ll != 0:
                                            ### resolve the issue for 25317062 for RHD that always have only 50% coverage
                                            problems=['RHD@Reference','RHD@RHD_DAU-0','RHD@RHD*10.00_RHD*DAU0','RHD@RHD*01']
                                            llori=ll
                                            checkpos=25317062

                                            if self.build.upper()=='HG19':
                                                checkpos=25643553

                                            if ( gpos==checkpos or gpos2==checkpos ) and ( type1 in problems or type2 in problems or self.allele_type in ['het','het RHD 1136_C_T'] ):
                                                if (GTa1['g1'][0] != GTa2['g2'][0])  or (GTa1['g1_next'][0] != GTa2['g2_next'][0]):
                                                    ### adjust ll only if heterozygous
                                                    #ll=2 * Decimal(ll)
                                                    ll=1-Qerr
                                            
                                            ### resolve the issue for 25420739 (HG38) for RHCE that have only 50% coverage in WES
                                            if self.dtype in ['WES']:
                                                checkpos=25420739
                                                if self.build.upper()=='HG19':
                                                    checkpos=25747230
                                                if ( gpos==checkpos or gpos2==checkpos ) and ( self.problem_snv_type in ['het RHCE 48_G_C'] ):
                                                   if (GTa1['g1'][0] != GTa2['g2'][0])  or (GTa1['g1_next'][0] != GTa2['g2_next'][0]):
                                                       ### adjust ll only if heterozygous
                                                       ll=1-Qerr    

 
                                            ### test print
                                            #if type1 in ['RHD@RHD*04.04_RHD*DIV.4','RHD@RHD*49_RHD*DWN'] and type2=='RHD@RHD*01':
                                            #     print(type1, cpos1, type2, cpos2, ll, phase_type)
                                            #     print(gpos, b1, Q1, b11, Q11, gpos2, b2, Q2)
                                            #     print('allele1', type1, 'g1=', gpos, GTa1['g1'], 'g1_next=', gpos2, GTa1['g1_next'])
                                            #     print('allele2', type2, 'g2=', gpos, GTa2['g2'], 'g2_next=', gpos2, GTa2['g2_next'])
                                            #if type1 == 'RHD@RHD(L390L)' and type2 =='RHD@Reference':
                                            #         print type1, cpos1, type2, cpos2, ll, phase_type
                                            #         print gpos, b1, Q1, b11, Q11, gpos2, b2, Q2
                                            #         print 'allele1', type1, 'g1=', gpos, GTa1['g1'], 'g1_next=', gpos2, GTa1['g1_next']
                                            #         print 'allele2', type2, 'g2=', gpos, GTa2['g2'], 'g2_next=', gpos2, GTa2['g2_next']
                                            #if type1 == 'RHD@RHD_weak_D_330GTdel' and type2 =='RHD@Reference':
                                            #         print type1, cpos1, type2, cpos2, ll, phase_type
                                            #         print gpos, b1, Q1, b11, Q11, gpos2, b2, Q2
                                            #         print 'allele1', type1, 'g1=', gpos, GTa1['g1'], 'g1_next=', gpos2, GTa1['g1_next']
                                            #         print 'allele2', type2, 'g2=', gpos, GTa2['g2'], 'g2_next=', gpos2, GTa2['g2_next']
                                            #if gpos==25420681:
                                            #  if ll < 0.01:
                                            #    if type1 == 'RHCE@RHcE_polypeptide;_Rh30' and type2 =='RHCE@RHce_ces':
                                            #        print qname, type1, cpos1, type2, cpos2, ll, phase_type
                                            #        print gpos, b1, Q1, b11, Q11, gpos2, b2, Q2
                                            #        print 'allele1', type1, 'g1=', gpos, GTa1['g1'], 'g1_next=', gpos2, GTa1['g1_next']
                                            #        print 'allele2', type2, 'g2=', gpos, GTa2['g2'], 'g2_next=', gpos2, GTa2['g2_next']
                                            #    if type1 == 'RHCE@RHcE_polypeptide;_Rh30' and type2 =='RHCE@RHCE*Ce_733G::RHce_VS::RHCE*Ce_733G_':
                                            #        print qname, type1, cpos1, type2, cpos2, ll, phase_type
                                            #        print gpos, b1, Q1, b11, Q11, gpos2, b2, Q2
                                            #        print 'allele1', type1, 'g1=', gpos, GTa1['g1'], 'g1_next=', gpos2, GTa1['g1_next']
                                            #        print 'allele2', type2, 'g2=', gpos, GTa2['g2'], 'g2_next=', gpos2, GTa2['g2_next']
                                            #        print '\n'
                                            #if type1 == 'RHD@RHD_DAU-3' and type2 =='RHD@RHD_DAU-0':
                                            #        print qname, type1, cpos1, type2, cpos2, 'originalLL:', llori, 'adjusatedLL:', ll, phase_type
                                            #        print gpos, b1, Q1, b11, Q11, gpos2, b2, Q2
                                            #        print 'allele1', type1, 'g1=', gpos, GTa1['g1'], 'g1_next=', gpos2, GTa1['g1_next']
                                            #        print 'allele2', type2, 'g2=', gpos, GTa2['g2'], 'g2_next=', gpos2, GTa2['g2_next']
                                            #if type1 == 'RHD@RHD_DAU-0' and type2 =='RHD@Reference':
                                            #        print qname, type1, cpos1, type2, cpos2, 'originalLL:', llori, 'adjusatedLL:', ll, phase_type
                                            #        print gpos, b1, Q1, b11, Q11, gpos2, b2, Q2
                                            #        print 'allele1', type1, 'g1=', gpos, GTa1['g1'], 'g1_next=', gpos2, GTa1['g1_next']
                                            #        print 'allele2', type2, 'g2=', gpos, GTa2['g2'], 'g2_next=', gpos2, GTa2['g2_next']
                                            #        print '\n'



                                            if type1 not in ll_out:
                                                ll_out[type1]={}
                                                ll_out[type1][type2]={}
                                                ll_out[type1][type2][gpos]={}
                                                ll_out[type1][type2][gpos]['rawLL']=[]
                                                ll_out[type1][type2][gpos]['rawLL'].append(ll)
                                                ll_out[type1][type2][gpos]['prodLL']=ll
                                            else:
                                                if type2 not in ll_out[type1]:
                                                    ll_out[type1][type2]={}
                                                    ll_out[type1][type2][gpos]={}
                                                    ll_out[type1][type2][gpos]['rawLL']=[]
                                                    ll_out[type1][type2][gpos]['rawLL'].append(ll)
                                                    ll_out[type1][type2][gpos]['prodLL']=ll
                                                else:
                                                    if gpos not in ll_out[type1][type2]:
                                                        ll_out[type1][type2][gpos]={}
                                                        ll_out[type1][type2][gpos]['rawLL']=[]
                                                        ll_out[type1][type2][gpos]['rawLL'].append(ll)
                                                        ll_out[type1][type2][gpos]['prodLL']=ll
                                                    else:
                                                        ll_out[type1][type2][gpos]['rawLL'].append(ll)
                                                        ori=ll_out[type1][type2][gpos]['prodLL']
                                                        ll_out[type1][type2][gpos]['prodLL']=Decimal(ll_out[type1][type2][gpos]['prodLL']) * Decimal(ll)
                                                        if ll_out[type1][type2][gpos]['prodLL'] == 0:
                                                            print(qname, gpos, cpos, b, Q, e, type1, GTf1, GTa1, type2, GTf2, GTa2, ori, ll)
                                                            sys.exit("become 0 here")
                                        else:
                                            print("ll == 0 for %s at pos %d" %(qname, gpos))
                                            print(qname, type1, gpos, GTa1, type2, gpos2, GTa2, ll)
                                            sys.exit("become 0 here")
       
                            else:
                                print("Allele base is discordant when calling haplotype: read> %s, gpos> %d, b1> %s, b11> %s" % (qname, gpos, b1, b11))
                #print pos1based, readn                
                #print self.snp_d
        '''
             calculate log phase likelihood
        '''
        #type1 = 'RHCE@RHcE_polypeptide;_Rh30'
        #type2 = 'RHCE@RHCE*Ce_733G::RHce_VS::RHCE*Ce_733G_'
        #type3 = 'RHCE@RHce_ces'
        #gpos  = 25420681
        #print ll_out[type1][type2][gpos]['rawLL']
        #print ll_out[type1][type3][gpos]['rawLL']
        #print ll_out[type1][type2][gpos]['prodLL']
        #print ll_out[type1][type3][gpos]['prodLL']

        phase_LL={}
        for BG1 in ll_out:
            phase_LL[BG1]={}
            for BG2 in ll_out[BG1]:
                LLphase=0
                for pos in ll_out[BG1][BG2]:
                    ll=ll_out[BG1][BG2][pos]['prodLL']
                    if ll != 0:
                        logll=Decimal(ll).log10()
                        LLphase+=logll
                    else:
                        print('%s\t%s\t%d\tll == 0, skip' %(BG1, BG2, pos))
                phase_LL[BG1][BG2]=LLphase
                #print('phaseLL', BG1, BG2, LLphase)
        return phase_LL
                  
    def calculateLL(self):
        genoLL=self.geno_ll()
        if genoLL=='REF':
            bt='REF'
            if self.gene=='RHCE':
                #bt='RHCE@RHce'
                bt='RHCE@RHCE*01_or_RHCE*ce_RHCE*c_RHCE*e'
            elif self.gene=='ABO':
                bt='ABO@ABO*A1.01'
            else:
                #bt='RHD@Reference'
                bt='RHD@RHD*01'
             
            res=pd.DataFrame([[1, 'A1', bt, 0,'','','',''],[1, 'A2', bt, 0,'','','','']], columns=['Index', 'Allele', 'Bloodtype', 'Likelihood', 'matchDBaa', 'not_matchDBaa','SampleExtraVar_aa','AllVariations_aa'])
            return res
        else:
            phaseLL=self.phase_ll()
            totalLL=pd.DataFrame(columns=['Allele', 'Bloodtype', 'Likelihood'])
            index=0
            for BG1 in genoLL:
                for BG2 in genoLL[BG1]:
                    if 'het' in self.allele_type:
                        if BG1 == BG2:
                            continue
                    elif 'hom' in self.allele_type:
                        if BG1 != BG2:
                            continue
                        
                    index+=1
                    gLL=genoLL[BG1][BG2]
                    pLL=phaseLL[BG1][BG2] if BG1 in phaseLL and BG2 in phaseLL[BG1] else 0
                    getcontext().prec = 7
                    tLL=Decimal(gLL)+Decimal(pLL)
                    res=pd.DataFrame([[index, 'A1', BG1, tLL], [index, 'A2',BG2, tLL]], columns=['Index', 'Allele', 'Bloodtype', 'Likelihood'])
                    totalLL=totalLL.append(res, ignore_index=True, sort=True)
            totalLL.sort_values(by=['Likelihood', 'Index', 'Allele'], ascending=[0,1,1], inplace=True)
            totalLL.reset_index(inplace=True, drop=True)
            return totalLL

         


class CLF(object):
    def __init__(self, train, target, test, weight=None):
        """
        Takes training, target, test
        """
        self.train =train.values
        self.target=target
        self.test  =test.values
        self.train_bloodtype_n=train.shape[1]
        self.train_feature_n  =train.shape[0]
        self.weight           =weight
        self.random_state     =1
    
    
    def SVM(self, probout=False, maxiter=1000):
        #print test
        clf=svm.LinearSVC(penalty='l2', loss='squared_hinge', dual=True, tol=0.0001, C=1.0, multi_class='ovr',fit_intercept=True, intercept_scaling=1, class_weight=self.weight, verbose=0, random_state=self.random_state, max_iter=maxiter)
        clf.fit(self.train, self.target)
        deci=pd.DataFrame(zip(clf.classes_, clf.decision_function(self.test)[0]))
        #print len(clf.decision_function(test)[0])

        clf2=svm.SVC(kernel='linear', probability=False, decision_function_shape='ovr', class_weight=self.weight)
        clf2.fit(self.train, self.target)
        deci2=pd.DataFrame(zip(clf2.classes_, clf2.decision_function(self.test)[0]))
        #print len(clf2.decision_function(test)[0])

        deci.columns=['Bloodtype','DecisionLinearSVC']
        deci2.columns=['Bloodtype','DecisionSVClinear']
        res=deci.merge(deci2, on='Bloodtype', how='outer').fillna(0)

        best_n=''
        n=2
        if probout:
            clf3=svm.SVC(kernel='linear', probability=True, class_weight=self.weight)
            clf3.fit(self.train, self.target)
            probs=clf3.predict_proba(self.test)
            best_n=sorted( zip(clf3.classes_, probs[0] ), key=lambda x:x[1] )[-n:]
            prob=pd.DataFrame(zip(clf3.classes_, clf3.predict_proba(self.test)[0]))
            #print len(clf3.predict_proba(test)[0])
            prob.columns=['Bloodtype','ProbabilitySVClinear']
            res=res.merge(prob, on='Bloodtype', how='outer').fillna(0)
        res=res.sort_values('DecisionLinearSVC', ascending=False)
        res['SVC_Rank'] = res['DecisionLinearSVC'].rank(ascending=False)
        return res


    def decisionTree(self, method='gini'):
        # model = tree.DecisionTreeRegressor() for regression
        clf = tree.DecisionTreeClassifier(criterion=method, random_state=self.random_state) # for classification, gini or entropy
        # Train the model using the training sets and check score
        clf.fit(self.train, self.target)
        #Predict Output
        prob=pd.DataFrame(zip(clf.classes_, clf.predict_proba(self.test)[0]))
        return prob


    def MNBS(self):
        #clf = GaussianNB()
        clf = MultinomialNB(class_prior=self.weight)
        clf.fit(self.train, self.target)
        prob=pd.DataFrame(zip(clf.classes_, clf.predict_proba(self.test)[0]))
        prob.columns=['Bloodtype','ProbabilityMultinominalNaiveBayes']
        res=prob.sort_values('ProbabilityMultinominalNaiveBayes', ascending=False)
        res['MNB_Rank'] = res['ProbabilityMultinominalNaiveBayes'].rank(ascending=False)
        return res

    def RF(self):
        clf=RandomForestClassifier(class_weight=self.weight, random_state=self.random_state)
        clf.fit(self.train, self.target)
        prob=pd.DataFrame(zip(clf.classes_, clf.predict_proba(self.test)[0]))
        prob.columns=['Bloodtype','ProbabilityRandomForest']
        res=prob.sort_values('ProbabilityRandomForest', ascending=False)
        res['RF_Rank'] = res['ProbabilityRandomForest'].rank(ascending=False)
        return res

    def GBM(self):
        clf=GradientBoostingClassifier(n_estimators=100, learning_rate=1.0, max_depth=1, random_state=self.random_state)
        clf.fit(self.train, self.target)
        prob=pd.DataFrame(zip(clf.classes_, clf.predict_proba(self.test)[0]))
        prob.columns=['Bloodtype','ProbabilityGBM']
        res=prob.sort_values('ProbabilityGBM', ascending=False)
        res['GBM_Rank'] = res['ProbabilityGBM'].rank(ascending=False)
        return res

    def summary(self):
        svm_res =self.SVM()
        mnbs_res=self.MNBS()
        rf_res  =self.RF()
        #gbm_res =self.GBM(traindata.values, target, test.values)
        dfs = [svm_res, mnbs_res, rf_res]
        df_merged=reduce(lambda  left,right: pd.merge(left,right,on=['Bloodtype'], how='outer'), dfs).fillna('void')
        df_merged['Rank']=df_merged['RF_Rank'] + df_merged['MNB_Rank'] + df_merged['SVC_Rank']
        df_merged.sort_values('Rank', ascending=True)
        df_merged=df_merged.drop_duplicates()
        summary=df_merged.nsmallest(1, 'Rank')
        return summary


