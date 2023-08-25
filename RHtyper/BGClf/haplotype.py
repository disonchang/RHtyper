"""
Haplotype modules
    
"""
import pandas as pd
import pysam
from general import complement
 
class build_hap(object):
    def __init__(self, df, untrimmed_df, bam, vars, gene='RHD',  verbose=0, minimum_base_quality = 15, ALT_n_LB=None, phase_rescue=None):
        self.df=df
        self.untrimmed_df=untrimmed_df
        self.bam=bam
        self.vars=vars
        self.gene=gene
        self.verbose=verbose
        self.minimum_base_quality=minimum_base_quality
        self.ALT_n_LB=ALT_n_LB
        self.phase_rescue=phase_rescue
        
    
    def hapblock_smooth(self, hapblock):
        '''
           check whether in allele called are all ref on one allele and alt on the allele.
           If yes, always put ref on the allele1 and alt on the allele2
    
        '''
        df=hapblock
        if self.verbose >0: print("Checking allele strands...")
        max_blk=[]
        max_blk_n=None
        for blk in df:
            block_a1_strand=None
            consistent_orientation=True
            varn=0
            for pos in df[blk]:
                varn+=1
                a1=df[blk][pos]['a1'] if 'a1' in df[blk][pos] else '.'
                a2=df[blk][pos]['a2'] if 'a2' in df[blk][pos] else '.'
                for ref in self.vars.loc[(self.vars['gene']==self.gene) & (self.vars['gpos']==pos)]['REF'].values:
                    if not block_a1_strand:
                        block_a1_strand='ref' if a1==ref else 'alt'
                    else:
                        if a1==ref:
                             if block_a1_strand !='ref':
                                 consistent_orientation=False
                        else:
                             if block_a1_strand !='alt':
                                 consistent_orientation=False
            if consistent_orientation:
                print("Consistent strand for block %d" % blk)
                for pos in df[blk]:
                     a1=df[blk][pos]['a1'] if 'a1' in df[blk][pos] else '.'
                     a2=df[blk][pos]['a2'] if 'a2' in df[blk][pos] else '.'
                     for ref in self.vars.loc[(self.vars['gene']==self.gene) & (self.vars['gpos']==pos)]['REF'].values:
                         if a1 != ref:
                             #print "Switching...alleles"
                             #print "Original", df[blk][pos]
                             a_temp=df[blk][pos]['a1']
                             df[blk][pos]['a1']=df[blk][pos]['a2']
                             df[blk][pos]['a2']=a_temp
                             #print "After", df[blk][pos]
            else:
                print ("Inconsistent strand for block %d" % blk)
            if self.verbose > 3: print('Block %d contatins %d Vars' % (blk, varn))
            if not max_blk_n:
                max_blk=[blk]
                max_blk_n=varn
            else:
                if varn > max_blk_n:
                    max_blk_n=varn
                    max_blk=[blk]
                elif varn==max_blk_n:
                    max_blk.append(blk)
        if self.verbose >0: print('Max haploblock:', ','.join(map(str,max_blk)), '; variation N:', max_blk_n)
        return df, max_blk, consistent_orientation
    
    
    def haplotype_complement_base(self, gpos, base):
        '''
        Return the complement alternative base given a variation positoin in a gene
        Helpful for defining the alternative allelic base
        '''
        cBase=None
        REFbase=self.vars.loc[(self.vars['gene']==self.gene) & (self.vars['gpos']==gpos)]['REF'].values
        ALTbase=self.vars.loc[(self.vars['gene']==self.gene) & (self.vars['gpos']==gpos)]['ALT'].values
        #if len(REFbase) != 1: return None
        #if len(ALTbase) != 1: return None
        if len(REFbase) != 1:
            print("[Warn Multiple REF]", gpos, REFbase)
            print("[Warn Multiple ALT]", gpos, ALTbase)
        for i in range(0,len(REFbase)):
           if ALTbase[i] == '.': continue ### skip the ones that have no varitaion, i.e. ABO 260 case
           if base == REFbase[i]:
               cBase=ALTbase[i]
           elif base == ALTbase[i]:
               cBase=REFbase[i]
        return cBase

    def BAF_adjust(self, allele1):
        '''
          adjust the bases that cannot be haplotyped to complement proportion
          since the BAF was calculated in the main matrix, should use allele1 for adjust here
        '''
        for idx, row in allele1.iterrows():
            if 0 < row['allele_frac'] < 1:
                BAF=allele1.loc[idx, 'allele_frac']
                allele1.loc[idx, 'allele_frac']=1-BAF
        return allele1

    
    def ABO_261Gdel_block(self, hap_block_pos, allele1, allele2):
        ''' 
           adjust ratio for ABO 261Gdel 
        '''
        if self.verbose >0: print("Adjust ratio for alleles in the 261G del block")
        targetblock=-1
        for blk in hap_block_pos:
            for pos in hap_block_pos[blk]:
                if pos==133257522: 
                    targetblock=blk
        if targetblock>=0:
            blk=targetblock
            for pos in hap_block_pos[blk]:
                a1=hap_block_pos[blk][pos]['a1'] if 'a1' in hap_block_pos[blk][pos] else '.'
                a2=hap_block_pos[blk][pos]['a2'] if 'a2' in hap_block_pos[blk][pos] else '.'
                cposs=self.vars.loc[(self.vars['gene']==self.gene) & (self.vars['gpos']==pos)]['cpos'].values
                for cpos in cposs:
                    ref=self.vars.loc[(self.vars['gene']==self.gene) & (self.vars['cpos']==cpos)]['REF'].values[0]
                    alt=self.vars.loc[(self.vars['gene']==self.gene) & (self.vars['cpos']==cpos)]['ALT'].values[0]
                    if alt=='.': continue ### avoid insertion that have uplicated gpos
                    ref_n=99; alt_n=99
                    if ref != '-': ref_n=self.vars.loc[(self.vars['gene']==self.gene) & (self.vars['cpos']==cpos)][ref].values[0]  
                    if alt != '-': alt_n=self.vars.loc[(self.vars['gene']==self.gene) & (self.vars['cpos']==cpos)][alt].values[0]
                    if ref_n < self.ALT_n_LB or alt_n < self.ALT_n_LB: continue ### skip low confidence variant          

                    idx=str(cpos)+'_'+ref+'_'+alt
                    if a1=='.': a1=alt if a2==ref else ref
                    if a2=='.': a2=alt if a1==ref else ref
                    if self.verbose > 3: print(blk, pos, cpos, idx, 'a1:'+a1, 'a2:'+a2, 'ref:'+ref, 'alt:'+alt)

                    if self.verbose > 3:
                        print('original allele1:', allele1[allele1['idx']==idx]['idx'].values[0], allele1[allele1['idx']==idx]['allele_frac'].values[0])
                        print('original allele2:', allele2[allele2['idx']==idx]['idx'].values[0], allele2[allele2['idx']==idx]['allele_frac'].values[0])
                    if a1==ref and a2==alt:
                        allele1.loc[(allele1['idx']==idx),['allele_frac']]=1.00
                        allele2.loc[(allele2['idx']==idx),['allele_frac']]=0
                    elif a1==alt and a2==ref:
                        allele1.loc[(allele1['idx']==idx),['allele_frac']]=0.00
                        allele2.loc[(allele2['idx']==idx),['allele_frac']]=1
                    else:
                        print('[Warning] REF/ALT not match phased nucleotides')

                    if self.verbose > 3:
                        print('adjusted allele1', allele1[allele1['idx']==idx]['idx'].values[0], allele1[allele1['idx']==idx]['allele_frac'].values[0])
                        print('adjusted allele2', allele2[allele2['idx']==idx]['idx'].values[0], allele2[allele2['idx']==idx]['allele_frac'].values[0])
        return allele1, allele2

    
    def build(self):
        dt=self.df.transpose()
        dt.reset_index(inplace=True)
        dt[['cpos','ref','alt']]=dt['idx'].str.split('_', expand=True)
        dt['cpos']=dt['cpos'].astype('int')
        dt=dt.sort_values(by=['cpos'])
        dt.reset_index(drop=True, inplace=True)
        dt_f=dt[(dt['ref'] != dt['alt']) & (0 < dt['allele_frac']) & (dt['allele_frac'] < 1) & (self.ALT_n_LB<=dt['alt_cnt'])]
        dt_f.reset_index(drop=True, inplace=True)
        dt_f=dt_f.assign(InDB='Y')
        cpos_max=dt_f['cpos'].max()
        cpos_min=dt_f['cpos'].min()
    
        untrimmed_df=self.untrimmed_df[self.untrimmed_df['gene']==self.gene]
        untrimmed_df=untrimmed_df.sort_values(by=['cpos'])
        untrimmed_df.reset_index(drop=True, inplace=True)
        untrimmed_df_f=untrimmed_df[(untrimmed_df['allele_frac'] < 1) & (self.ALT_n_LB<=untrimmed_df['alt_cnt']) & (untrimmed_df['cpos']<=cpos_max) & (untrimmed_df['cpos']>=cpos_min)]
        untrimmed_df_f.reset_index(drop=True, inplace=True)
    
        if dt_f.shape[0] != untrimmed_df_f.shape[0]:
            if self.verbose > 3: print("Adding more SNPs for phasing")
            dt_f=dt_f.merge(untrimmed_df_f, how='outer')
            dt_f=dt_f.sort_values(by=['cpos'])
            dt_f.reset_index(drop=True, inplace=True)
            dt_f['InDB'].fillna(value='N', inplace=True)
            if self.verbose > 3: print(dt_f)
    
        allele1=dt.copy()
        allele2=dt.copy()
        if dt_f.shape[0] == 0:
            if self.verbose > 0: print("Haplotype, no segregating variations available")
        if dt_f.shape[0] == 1:
            if self.verbose > 0: print("Haplotype will be determined by a single variant")
            ### if dt_f.shape[0] == 0, both alleles are the same 
            ### if dt_f.shape[0] == 1, this single variation will be able to distinguish allele1 and allele2
            ### There are supposed one freq for ref and one freq for alt here ... check
            q='cpos==' + str(dt_f['cpos'].values[0])
            dt_alleles=dt.query(q)
            if self.verbose > 0:
                print('Two haplotypes found' if dt_alleles.shape[0] == 1 else 'More than two haplotypes found ... use the top 2 as alleles')
            ### get two most freqent alleles
            dt_alleles=dt_alleles.nlargest(1, 'allele_frac')
            idxes=dt_alleles['idx'].values
    
            ### adjust ratios
            if self.verbose > 3:
                print('original allele1:', allele1[allele1['idx']==idxes[0]]['idx'].values[0], allele1[allele1['idx']==idxes[0]]['allele_frac'].values[0])
                print('original allele2:', allele2[allele2['idx']==idxes[0]]['idx'].values[0], allele2[allele2['idx']==idxes[0]]['allele_frac'].values[0])
    
    
            allele1.loc[(allele1['idx']==idxes[0]),['allele_frac']]=1.00
            allele2.loc[(allele2['idx']==idxes[0]),['allele_frac']]=0
    
            if self.verbose > 3:
                print('adjusted allele1', allele1[allele1['idx']==idxes[0]]['idx'].values[0], allele1[allele1['idx']==idxes[0]]['allele_frac'].values[0])
                print('adjusted allele2', allele2[allele2['idx']==idxes[0]]['idx'].values[0], allele2[allele2['idx']==idxes[0]]['allele_frac'].values[0])
    
        elif dt_f.shape[0] > 1:
            ### if dt_f.shape[0] > 1, need to check whether varations can be combined to allele haplotypes based on 
            ### the distance between the alleles
            if self.verbose >0: print("Haplotype will be determined by multiple variants")
            if self.verbose >0: print(dt_f)
            hap_block=-1
            hap_block_pos={}
            hap_block_reads={}
            for idx, row in dt_f.iterrows():
                if idx+1 < dt_f.shape[0]:
                    pos     =row['cpos']
                    pos_next=dt_f.iloc[idx+1]['cpos']
                    strand=self.vars.loc[(self.vars['gene']==self.gene) & (self.vars['cpos']==pos)]['strand'].values[0]
                    chr   =self.vars.loc[(self.vars['gene']==self.gene) & (self.vars['cpos']==pos)]['chr'].values[0]
                    if strand == '+':
                        gpos=int(self.vars.loc[(self.vars['gene']==self.gene) & (self.vars['cpos']==pos)]['gpos'].values[0])
                        gpos_next=int(self.vars.loc[(self.vars['gene']==self.gene) & (self.vars['cpos']==pos_next)]['gpos'].values[0])
                    else:
                        gpos=int(self.vars.loc[(self.vars['gene']==self.gene) & (self.vars['cpos']==pos_next)]['gpos'].values[0])
                        gpos_next=int(self.vars.loc[(self.vars['gene']==self.gene) & (self.vars['cpos']==pos)]['gpos'].values[0])
    
                    REF=self.vars.loc[(self.vars['gene']==self.gene) & (self.vars['cpos']==pos)]['REF'].values[0]
                    ALT=self.vars.loc[(self.vars['gene']==self.gene) & (self.vars['cpos']==pos)]['ALT'].values[0]
                    if self.verbose > 10: print(hap_block, pos, row['idx'], pos_next, dt_f.iloc[idx+1]['idx'])
                    ### check support reads of haplotype
                    samfile = pysam.AlignmentFile(self.bam, mode="rb") #, ignore_truncation=False)
                    block_bases={}
                    for column in samfile.pileup(chr, gpos-1, gpos_next):
                        pos0based=column.pos
                        pos1based=column.pos+1
                        if pos1based == gpos or pos1based == gpos_next:
                            for pileupread in column.pileups:
                                read     = pileupread.alignment
                                qname    = pileupread.alignment.query_name
                                if self.gene=='ABO' and pos1based==133257522:
                                    ''' handle the particular case of 261G del in ABO here'''
                                    qpos     = pileupread.query_position
                                    readbase = read.query_sequence[qpos-1]
                                    readbaseQ= read.query_qualities[qpos-1]
                                    #print "[[Checking]]", pos1based, REF, ALT, qname, readbase, read.cigarstring                                    
                                    if readbase=='T':
                                        ### no insertion on HG38, interpretated as deletion
                                        readbase='-'
                                        readbaseQ=99
                                elif pileupread.is_del or pileupread.indel < 0:
                                    ### deletion
                                    readbase ='-'
                                    readbaseQ=99
                                elif pileupread.is_refskip or pileupread.indel > 0:
                                    ### insertion
                                    qpos     = pileupread.query_position
                                    readbase = read.query_sequence[qpos] ### need check +1 or -1
                                    readbaseQ= read.query_qualities[qpos]
                                else:
                                    qpos     = pileupread.query_position
                                    readbase = read.query_sequence[qpos]
                                    readbaseQ= read.query_qualities[pileupread.query_position]
                                if self.verbose > 99: print(pos1based, REF, ALT, qname, readbase)
                                rb=readbase       if readbaseQ >= self.minimum_base_quality else 'N'
                                rb=complement(readbase) if strand == '-' and rb != '-' else readbase
                                if qname not in block_bases:
                                    block_bases[qname]={}
                                    block_bases[qname][pos1based]=rb
                                else:
                                    block_bases[qname][pos1based]=rb
                    ### construct haplotypes
                    #print block_bases
                    for read in block_bases:
                        poses=list(block_bases[read].keys()); n_pos=len(poses)
                        #if hap_block not in hap_block_pos: hap_block_pos[hap_block]={}
                        if n_pos==2:
                              ### all the comparison should be in pairs
                              pos1=poses[0] if poses[0] < poses[1] else poses[1]
                              pos2=poses[1] if poses[0] < poses[1] else poses[0]
                              base1=block_bases[read][pos1]; base2=block_bases[read][pos2]
                              if self.verbose > 9: print(pos1, pos2, base1, base2, read)
    
                              if base1 == 'N' or base2 =='N':
                                  ### ignore low quality block
                                  if self.verbose > 9: print('skip low Q reads', read)
                                  continue
    
                              ### initiate
                              if hap_block==-1:
                                  hap_block=0
                                  hap_block_pos[hap_block]={}
                                  hap_block_pos[hap_block][pos1]={}
                                  hap_block_pos[hap_block][pos2]={}
                                  hap_block_pos[hap_block][pos1]['a1']=base1
                                  hap_block_pos[hap_block][pos2]['a1']=base2
                                  ### fill in alternative allele
                                  if self.haplotype_complement_base(pos1, base1): hap_block_pos[hap_block][pos1]['a2']=self.haplotype_complement_base(pos1, base1)
                                  if self.haplotype_complement_base(pos2, base2): hap_block_pos[hap_block][pos2]['a2']=self.haplotype_complement_base(pos2, base2)
    
                                  hap_block_reads[hap_block]={}
                                  hap_block_reads[hap_block]['allele1_reads']=[read]
                                  hap_block_reads[hap_block]['allele2_reads']=[]
                                  hap_block_reads[hap_block]['confliced_reads']=[]
                                  if self.verbose > 9: print('0. Initiate Hapblock', hap_block, 'allele1');
                                  if self.verbose > 9: print(hap_block_pos,'\n')
    
                              if pos1 not in hap_block_pos[hap_block] and pos2 not in hap_block_pos[hap_block]:
                                  hap_block+=1
                                  ### RE-initiate broken hap block
                                  hap_block_pos[hap_block]={}
                                  hap_block_pos[hap_block][pos1]={}
                                  hap_block_pos[hap_block][pos2]={}
                                  hap_block_pos[hap_block][pos1]['a1']=base1
                                  hap_block_pos[hap_block][pos2]['a1']=base2
                                  ### fill in alternative allele
                                  if self.haplotype_complement_base(pos1, base1): hap_block_pos[hap_block][pos1]['a2']=self.haplotype_complement_base(pos1, base1)
                                  if self.haplotype_complement_base(pos2, base2): hap_block_pos[hap_block][pos2]['a2']=self.haplotype_complement_base(pos2, base2)
    
                                  hap_block_reads[hap_block]={}
                                  hap_block_reads[hap_block]['allele1_reads']=[read]
                                  hap_block_reads[hap_block]['allele2_reads']=[]
                                  hap_block_reads[hap_block]['confliced_reads']=[]
                                  if self.verbose > 9: print('1. Broken... Re-Initiate Hapblock', hap_block, 'allele1');
                                  if self.verbose > 9: print(hap_block_pos,'\n')
                              elif pos1 in hap_block_pos[hap_block] and pos2 in hap_block_pos[hap_block]:
                                  ### check of support evidence of existing haplotype
                                  ### create 2nd haplotype if not existing
                                  if base1 == hap_block_pos[hap_block][pos1]['a1'] and base2 == hap_block_pos[hap_block][pos2]['a1']:
                                      if read not in hap_block_reads[hap_block]['allele1_reads']:
                                          hap_block_reads[hap_block]['allele1_reads'].append(read)
                                          if self.verbose > 9: print('support hap1', hap_block, 'allele1')
                                          if self.verbose > 9: print(hap_block_pos, '\n')
                                  elif base1 != hap_block_pos[hap_block][pos1]['a1'] and base2 != hap_block_pos[hap_block][pos2]['a1']:
                                      if 'a2' not in hap_block_pos[hap_block][pos1]:
                                          hap_block_pos[hap_block][pos1]['a2']=base1
                                          hap_block_pos[hap_block][pos2]['a2']=base2
                                          ### fill in alternative allele
                                          if 'a1' not in hap_block_pos[hap_block][pos1]:
                                              if self.haplotype_complement_base(pos1, base1): hap_block_pos[hap_block][pos1]['a1']=self.haplotype_complement_base(pos1, base1)
                                          if 'a1' not in hap_block_pos[hap_block][pos2]:
                                              if self.haplotype_complement_base(pos2, base2): hap_block_pos[hap_block][pos2]['a1']=self.haplotype_complement_base(pos2, base2)
    
                                          hap_block_reads[hap_block]['allele2_reads'].append(read)
                                          if self.verbose > 9: print('2. Initiate Hapblock', hap_block, 'allele2')
                                          if self.verbose > 9: print(hap_block_pos, '\n')
                                      else:
                                          #print hap_block_pos
                                          if base1 == hap_block_pos[hap_block][pos1]['a2']:
                                              if 'a2' not in hap_block_pos[hap_block][pos2]:
                                                  hap_block_pos[hap_block][pos2]['a2']=base2
                                                  ### fill in alternative allele
                                                  if 'a1' not in hap_block_pos[hap_block][pos2]:
                                                      if self.haplotype_complement_base(pos2, base2): hap_block_pos[hap_block][pos2]['a1']=self.haplotype_complement_base(pos2, base2)
                                                  if read not in hap_block_reads[hap_block]['allele2_reads']:
                                                      hap_block_reads[hap_block]['allele2_reads'].append(read)
                                                  if self.verbose >9: print('3. Extend Hapblock', hap_block, 'allele2')
                                                  if self.verbose >9: print(hap_block_pos,'\n')
                                              else:
                                                  if base2 == hap_block_pos[hap_block][pos2]['a2']:
                                                      if read not in hap_block_reads[hap_block]['allele2_reads']:
                                                          hap_block_reads[hap_block]['allele2_reads'].append(read)
                                                      if self.verbose > 9: print('support hap2', hap_block, 'allele2')
                                                      if self.verbose > 9: print(hap_block_pos, '\n')
                                                  else:
                                                      if read not in hap_block_reads[hap_block]['confliced_reads']:
                                                          hap_block_reads[hap_block]['confliced_reads'].append(read)
                                                      if self.verbose >3: print('conflict reads')
                              elif pos1 in hap_block_pos[hap_block] and pos2 not in hap_block_pos[hap_block]:
                                  ### extend haplotype
                                  hap_block_pos[hap_block][pos2]={}
                                  if base1 == hap_block_pos[hap_block][pos1]['a1']:
                                      hap_block_pos[hap_block][pos2]['a1']=base2
                                       ### fill in alternative allele
                                      if 'a2' not in hap_block_pos[hap_block][pos1]:
                                              if self.haplotype_complement_base(pos1, base1): hap_block_pos[hap_block][pos1]['a2']=self.haplotype_complement_base(pos1, base1)
                                      if 'a2' not in hap_block_pos[hap_block][pos2]:
                                              if self.haplotype_complement_base(pos2, base2): hap_block_pos[hap_block][pos2]['a2']=self.haplotype_complement_base(pos2, base2)
                                      if read not in hap_block_reads[hap_block]['allele1_reads']:
                                          hap_block_reads[hap_block]['allele1_reads'].append(read)
                                      if self.verbose >9: print('4. Extend Hapblock', hap_block, 'allele1')
                                      if self.verbose >9: print(hap_block_pos,'\n')
                                  else:
                                      if 'a2' not in hap_block_pos[hap_block][pos1]:
                                          hap_block_pos[hap_block][pos1]['a2']=base1
                                          hap_block_pos[hap_block][pos2]['a2']=base2
                                          ### fill in alternative allele
                                          if 'a1' not in hap_block_pos[hap_block][pos1]:
                                              if self.haplotype_complement_base(pos1, base1): hap_block_pos[hap_block][pos1]['a1']=self.haplotype_complement_base(pos1, base1)
                                          if 'a1' not in hap_block_pos[hap_block][pos2]:
                                              if self.haplotype_complement_base(pos2, base2): hap_block_pos[hap_block][pos2]['a1']=self.haplotype_complement_base(pos2, base2)
    
                                          if read not in hap_block_reads[hap_block]['allele2_reads']:
                                              hap_block_reads[hap_block]['allele2_reads'].append(read)
                                          if self.verbose >9:print('4. Extend Hapblock', hap_block, 'allele2')
                                          if self.verbose >9:print(hap_block_pos,'\n')
                                      else:
                                          if base1 == hap_block_pos[hap_block][pos1]['a2']:
                                              hap_block_pos[hap_block][pos2]['a2']=base2
                                              ### fill in alternative allele
                                              if 'a1' not in hap_block_pos[hap_block][pos1]:
                                                  #print '5. Need fill in a1 for', pos1, haplotype_complement_base(gene, pos1, base1, vars)
                                                  if self.haplotype_complement_base(pos1, base1): hap_block_pos[hap_block][pos1]['a1']=self.haplotype_complement_base(pos1, base1)
                                              if 'a1' not in hap_block_pos[hap_block][pos2]:
                                                  #print '5. Need fill in a1 for', pos2, haplotype_complement_base(gene, pos2, base2, vars)
                                                  if self.haplotype_complement_base(pos2, base2): hap_block_pos[hap_block][pos2]['a1']=self.haplotype_complement_base(pos2, base2)
                                              if read not in hap_block_reads[hap_block]['allele2_reads']:
                                                  hap_block_reads[hap_block]['allele2_reads'].append(read)
                                              if self.verbose >9:print('5. Extend Hapblock', hap_block, 'allele2')
                                              if self.verbose >9:print(hap_block_pos,'\n')
                                          else:
                                              print("[error] extra alleles")
                              elif pos1 not in hap_block_pos[hap_block] and pos2 in hap_block_pos[hap_block]:
                                  ### extend haplotype
                                  hap_block_pos[hap_block][pos1]={}
                                  if base2 == hap_block_pos[hap_block][pos2]['a1']:
                                      hap_block_pos[hap_block][pos1]['a1']=base1
                                       ### fill in alternative allele
                                      if 'a2' not in hap_block_pos[hap_block][pos1]:
                                              if self.haplotype_complement_base(pos1, base1): hap_block_pos[hap_block][pos1]['a2']=self.haplotype_complement_base(pos1, base1)
                                      if 'a2' not in hap_block_pos[hap_block][pos2]:
                                              if self.haplotype_complement_base(pos2, base2): hap_block_pos[hap_block][pos2]['a2']=self.haplotype_complement_base(pos2, base2)
                                      if read not in hap_block_reads[hap_block]['allele1_reads']:
                                          hap_block_reads[hap_block]['allele1_reads'].append(read)
                                      if self.verbose >9: print('6. Extend Hapblock', hap_block, 'allele1')
                                      if self.verbose >9: print(hap_block_pos,'\n')
                                  else:
                                      if 'a2' not in hap_block_pos[hap_block][pos2]:
                                          hap_block_pos[hap_block][pos1]['a2']=base1
                                          hap_block_pos[hap_block][pos2]['a2']=base2
                                          ### fill in alternative allele
                                          if 'a1' not in hap_block_pos[hap_block][pos1]:
                                              if self.haplotype_complement_base(pos1, base1): hap_block_pos[hap_block][pos1]['a1']=self.haplotype_complement_base(pos1, base1)
                                          if 'a1' not in hap_block_pos[hap_block][pos2]:
                                              if self.haplotype_complement_base(pos2, base2): hap_block_pos[hap_block][pos2]['a1']=self.haplotype_complement_base(pos2, base2)

                                          if read not in hap_block_reads[hap_block]['allele2_reads']:
                                              hap_block_reads[hap_block]['allele2_reads'].append(read)
                                          if self.verbose >9:print('6. Extend Hapblock', hap_block, 'allele2')
                                          if self.verbose >9:print(hap_block_pos,'\n')
                                      else:
                                          if base2 == hap_block_pos[hap_block][pos2]['a2']:
                                              hap_block_pos[hap_block][pos1]['a2']=base1
                                              ### fill in alternative allele
                                              if 'a1' not in hap_block_pos[hap_block][pos1]:
                                                  #print '5. Need fill in a1 for', pos1, haplotype_complement_base(gene, pos1, base1, vars)
                                                  if self.haplotype_complement_base(pos1, base1): hap_block_pos[hap_block][pos1]['a1']=self.haplotype_complement_base(pos1, base1)
                                              if 'a1' not in hap_block_pos[hap_block][pos2]:
                                                  #print '5. Need fill in a1 for', pos2, haplotype_complement_base(gene, pos2, base2, vars)
                                                  if self.haplotype_complement_base(pos2, base2): hap_block_pos[hap_block][pos2]['a1']=self.haplotype_complement_base(pos2, base2)
                                              if read not in hap_block_reads[hap_block]['allele2_reads']:
                                                  hap_block_reads[hap_block]['allele2_reads'].append(read)
                                              if self.verbose >9:print('7. Extend Hapblock', hap_block, 'allele2')
                                              if self.verbose >9:print(hap_block_pos,'\n')
                                          else:
                                              print("[error] extra alleles")

                              else:
                                  print('unknown condition', pos1, pos2, base1, base2, read)
                #else:
                #    print hap_block, row['cpos'], row['idx'], 'stop'
            ### adjust ratio based on haplotype
            if self.verbose > 3:
                print("Haplotype:")
                if hap_block_pos:
                    print("%d blocks found" % len(hap_block_pos))
                    print(hap_block_pos)
                else:
                    print("Not able to phase by short reads")
    
            if not hap_block_pos:
                if self.phase_rescue and self.phase_rescue=='ratio':
                    if self.verbose > 0: print("Haplotpye will be determined by allelic ratios of the variants")
                    ### if dt_f.shape[0] > 1 and no SNPs were able to phased by short reads, i.e. distacne of SNP too far
                    ### The allelic ratio can be adjusted based on the allelic ratio in theory
                    for idx, row in dt_f.iterrows():
                        q='cpos=='+str(row['cpos'])
                        dt_alleles=dt.query(q)
                        if self.verbose > 5: print(dt_alleles)
                        if self.verbose > 0:
                            print('Two haplotypes found') if dt_alleles.shape[0] == 1 else 'More than two haplotypes found ... use the top 2 as alleles'
                            ### get two most freqent alleles
                            dt_alleles=dt_alleles.nlargest(1, 'allele_frac')
                            idxes=dt_alleles['idx'].values
    
                        ### adjust ratios
                        if self.verbose > 3:
                            print ('original allele1:', allele1[allele1['idx']==idxes[0]]['idx'].values[0], allele1[allele1['idx']==idxes[0]]['allele_frac'].values[0])
                            print ('original allele2:', allele2[allele2['idx']==idxes[0]]['idx'].values[0], allele2[allele2['idx']==idxes[0]]['allele_frac'].values[0])
    
                        if float(row['allele_frac']) > 0.5:
                            allele1.loc[(allele1['idx']==idxes[0]),['allele_frac']]=0
                            allele2.loc[(allele2['idx']==idxes[0]),['allele_frac']]=1
                        else:
                            allele1.loc[(allele1['idx']==idxes[0]),['allele_frac']]=1
                            allele2.loc[(allele2['idx']==idxes[0]),['allele_frac']]=0
    
                        if self.verbose > 3:
                            print ('adjusted allele1', allele1[allele1['idx']==idxes[0]]['idx'].values[0], allele1[allele1['idx']==idxes[0]]['allele_frac'].values[0])
                            print ('adjusted allele2', allele2[allele2['idx']==idxes[0]]['idx'].values[0], allele2[allele2['idx']==idxes[0]]['allele_frac'].values[0])
                elif self.phase_rescue and self.phase_rescue=='ref':
                    if self.verbose > 0: print ("Haplotpye will be determined by reference, i.e. reference one will be reference")
                    for idx, row in dt_f.iterrows():
                        q='cpos=='+str(row['cpos'])
                        dt_alleles=dt.query(q)
                        if self.verbose > 0:
                            print ('Two haplotypes found' if dt_alleles.shape[0] == 1 else 'More than two haplotypes found ... use the top 2 as alleles')
                            ### get two most freqent alleles
                            dt_alleles=dt_alleles.nlargest(1, 'allele_frac')
                            idxes=dt_alleles['idx'].values
    
                        ### adjust ratios
                        if self.verbose > 3:
                            print ('original allele1:', allele1[allele1['idx']==idxes[0]]['idx'].values[0], allele1[allele1['idx']==idxes[0]]['allele_frac'].values[0])
                            print ('original allele2:', allele2[allele2['idx']==idxes[0]]['idx'].values[0], allele2[allele2['idx']==idxes[0]]['allele_frac'].values[0])
    
                        ### reference always reference, alernative always alternative
                        allele1.loc[(allele1['idx']==idxes[0]),['allele_frac']]=0
                        allele2.loc[(allele2['idx']==idxes[0]),['allele_frac']]=1
    
                        if self.verbose > 3:
                            print ('adjusted allele1', allele1[allele1['idx']==idxes[0]]['idx'].values[0], allele1[allele1['idx']==idxes[0]]['allele_frac'].values[0])
                            print ('adjusted allele2', allele2[allele2['idx']==idxes[0]]['idx'].values[0], allele2[allele2['idx']==idxes[0]]['allele_frac'].values[0])
            else:
              hap_block_pos, max_blk, orientation=self.hapblock_smooth(hap_block_pos)
              if self.verbose >0: print ("Adjust ratio for alleles in the largest hap block")
              for blk in max_blk:
                for pos in hap_block_pos[blk]:
                    a1=hap_block_pos[blk][pos]['a1'] if 'a1' in hap_block_pos[blk][pos] else '.'
                    a2=hap_block_pos[blk][pos]['a2'] if 'a2' in hap_block_pos[blk][pos] else '.'
                    a1_n=len(hap_block_reads[blk]['allele1_reads'])
                    a2_n=len(hap_block_reads[blk]['allele2_reads'])
                    ac_n=len(hap_block_reads[blk]['confliced_reads'])
                    cposs=self.vars.loc[(self.vars['gene']==self.gene) & (self.vars['gpos']==pos)]['cpos'].values
                    for cpos in cposs:
                    #for ref in self.vars.loc[(self.vars['gene']==self.gene) & (self.vars['gpos']==pos)]['REF'].values:
                    #    for alt in self.vars.loc[(self.vars['gene']==self.gene) & (self.vars['gpos']==pos)]['ALT'].values:

                            ref=self.vars.loc[(self.vars['gene']==self.gene) & (self.vars['cpos']==cpos)]['REF'].values[0]
                            alt=self.vars.loc[(self.vars['gene']==self.gene) & (self.vars['cpos']==cpos)]['ALT'].values[0]
                            if alt=='.': continue ### avoid insertion that have uplicated gpos
                            idx=str(cpos)+'_'+ref+'_'+alt
                            if a1=='.': a1=alt if a2==ref else ref
                            if a2=='.': a2=alt if a1==ref else ref
                            if self.verbose > 3: print (blk, pos, cpos, idx, 'a1:'+a1, 'a2:'+a2, 'ref:'+ref, 'alt:'+alt, 'SN1:'+str(a1_n), 'SN2:'+str(a2_n), 'CN:'+str(ac_n))
                            ### adjust ratios
                            ref_n=99; alt_n=99
                            if ref != '-': ref_n=self.vars.loc[(self.vars['gene']==self.gene) & (self.vars['cpos']==cpos)][ref].values[0]
                            if alt != '-': alt_n=self.vars.loc[(self.vars['gene']==self.gene) & (self.vars['cpos']==cpos)][alt].values[0]
                            if ref_n < self.ALT_n_LB or alt_n < self.ALT_n_LB: continue ### skip low confidence variant      
    
                            ### untrimmed_df[(untrimmed_df['allele_frac'] < 1) & (ALT_n_LB<=untrimmed_df['alt_cnt']) & (untrimmed_df['cpos']<=cpos_max) & (untrimmed_df['cpos']>=cpos_min)]
                            if dt_f[(dt_f['cpos']==cpos)]['InDB'].values=='N':
                                if self.verbose > 9: print ("Skip SNPs not in DB: %d" % cpos)
                                continue ### ignore the added SNP that help for phasing 
                            if self.verbose > 3:
                                print ('original allele1:', allele1[allele1['idx']==idx]['idx'].values[0], allele1[allele1['idx']==idx]['allele_frac'].values[0])
                                print ('original allele2:', allele2[allele2['idx']==idx]['idx'].values[0], allele2[allele2['idx']==idx]['allele_frac'].values[0])
                            if a1==ref and a2==alt:
                                allele1.loc[(allele1['idx']==idx),['allele_frac']]=1.00
                                allele2.loc[(allele2['idx']==idx),['allele_frac']]=0
                            elif a1==alt and a2==ref:
                                allele1.loc[(allele1['idx']==idx),['allele_frac']]=0.00
                                allele2.loc[(allele2['idx']==idx),['allele_frac']]=1
                            else:
                                print ('[Warning] REF/ALT not match phased nucleotides')
    
                            if self.verbose > 3:
                                print ('adjusted allele1', allele1[allele1['idx']==idx]['idx'].values[0], allele1[allele1['idx']==idx]['allele_frac'].values[0])
                                print ('adjusted allele2', allele2[allele2['idx']==idx]['idx'].values[0], allele2[allele2['idx']==idx]['allele_frac'].values[0])
              if self.gene=="ABO":
                  allele1, allele2=self.ABO_261Gdel_block(hap_block_pos, allele1, allele2)
              ### adjust BAF
              #allele1=self.BAF_adjust(allele1)
    
        ### reformat alleles
        allele1.drop_duplicates(inplace=True)
        allele2.drop_duplicates(inplace=True)
        allele1.drop(['cpos','ref','alt'], axis=1, inplace=True); allele1.set_index('idx', drop=True, inplace=True); allele1=allele1.reindex(self.df.transpose().index)
        allele2.drop(['cpos','ref','alt'], axis=1, inplace=True); allele2.set_index('idx', drop=True, inplace=True); allele2=allele2.reindex(self.df.transpose().index)
    
        if self.verbose > 9:
            allele1.to_csv('allele1.txt', sep='\t', index=True, mode='w')
            allele2.to_csv('allele2.txt', sep='\t', index=True, mode='w')
    
        return allele1.transpose(), allele2.transpose()
    
