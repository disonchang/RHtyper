"""
Modules to generate a table with all variations in a BAM
"""
import pysam, re
import collections
import pandas as pd
from general import *
from collections import OrderedDict

def call(gene, cds, bam, output, gbuild='hg38', minimum_base_quality = 15, minimum_mapq = 10, minimum_read_quality = 15, verbose=5, ALT_n_LB=None, VAF=None, write_output=True):
    skip_message = 'SKIPPED: {} {} {} {} {}'
    pass_message = 'PASSED: {} {} {}'
    out=[]
    out_filtered=[]

    for g in gene:
        CDSseq=''
        genetable=gene[g]
        #print genetable
        cdstable=cds[g]
        cdstable['ALT']='.'
        cdstable['REFaa']='.'; cdstable['ALTaa']='.'
        cdstable['TOTAL_n']='0'; cdstable['QC_n']='0'; cdstable['DEL_n']='0'; cdstable['INS_n']='0'
        cdstable['A']=0; cdstable['T']=0; cdstable['C']=0; cdstable['G']=0
        cdstable['Ap']=0; cdstable['Tp']=0; cdstable['Cp']=0; cdstable['Gp']=0
        cdstable['dA']=0; cdstable['dT']=0; cdstable['dC']=0; cdstable['dG']=0
        cdstable['dAp']=0; cdstable['dTp']=0; cdstable['dCp']=0; cdstable['dGp']=0
        cdstable['iA']=0; cdstable['iT']=0; cdstable['iC']=0; cdstable['iG']=0
        cdstable['iAp']=0; cdstable['iTp']=0; cdstable['iCp']=0; cdstable['iGp']=0
        cdstable['max']=0; cdstable['max_del']=0; cdstable['max_ins']=0
        n=0
        for idx, row in genetable.iterrows():
            #n+=1
            #if n == 1: ### test purpose
              samfile = pysam.AlignmentFile(bam, mode="rb") #, ignore_truncation=False)

              trimchr=False
              if not all('chr' in s for s in samfile.references):trimchr=True
              tarchr=row['chr'].lstrip('chr') if trimchr else row['chr'] 

              for column in samfile.pileup(tarchr, row['start']-1, row['end']):
                  pos0based=column.pos
                  pos1based=column.pos+1
                  #print cdstable.loc[(cdstable['chr']==row['chr']) & (cdstable['gpos']==pos1based) & (cdstable['gene']==row['gene']),['cpos']].values
                  for cposs in cdstable.loc[(cdstable['chr']==tarchr) & (cdstable['gpos']==pos1based) & (cdstable['gene']==row['gene']),['cpos']].values:
                      for cpos in cposs:
                           if int(row['start']) <= int(pos1based) <= int(row['end']):
                              refbase= str(cdstable.loc[(cdstable['chr']==tarchr) & (cdstable['gpos']==pos1based) & (cdstable['cpos']==cpos) & (cdstable['gene']==row['gene']),['REF']].values[0][0])
                              strand = str(cdstable.loc[(cdstable['chr']==tarchr) & (cdstable['gpos']==pos1based) & (cdstable['cpos']==cpos) & (cdstable['gene']==row['gene']),['strand']].values[0][0])
                              CDSseq += refbase
                              #if verbose > 99: print("\ngene:%s, acc:%s, exon:%s, region:%s:%s-%s, coverage at base %s:%s = %s" % (row['gene'], row['acc'], row['exon_n'], row['chr'], row['start'], row['end'], pos1based, refbase, column.n))

                              good_count    ={'A':0, 'T':0, 'C':0, 'G':0, 'N':0}
                              ins_count     ={'A':0, 'T':0, 'C':0, 'G':0, 'N':0}
                              del_count     ={'A':0, 'T':0, 'C':0, 'G':0, 'N':0}
                              lowmapq_count ={'A':0, 'T':0, 'C':0, 'G':0, 'N':0}
                              nomapq_count  ={'A':0, 'T':0, 'C':0, 'G':0, 'N':0}
                              lowreadq_count={'A':0, 'T':0, 'C':0, 'G':0, 'N':0}
                              lowbaseq_count={'A':0, 'T':0, 'C':0, 'G':0, 'N':0}

                              ### process each read in each exon to get the count of bases
                              good_n=0
                              del_n =0
                              ins_n =0
                              ins_cpos=cpos+1 if strand =='+' else cpos-1
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
                                  
                                  #if pos1based==25321928:print(pos1based, qname, size, qpos, read.cigarstring, readbase)

                                  # skip reads with indels
                                  if pileupread.is_del:# or pileupread.indel < 0:
                                      #if pos1based==25321928: print(skip_message.format(pos1based, qname, 'is del', readbase, pileupread.indel))
                                      del_count[refbase]+=1
                                      del_n+=1
                                      continue

                                  # skip reads with mapq below threshold
                                  if pileupread.alignment.mapping_quality < minimum_mapq:
                                      #if pos1based==25321928: print(skip_message.format(position, qname, 'low mapq', pileupread.alignment.mapping_quality))
                                      lowmapq_count[readbase]+=1
                                      continue
                                  # skip reads with no mapq specified
                                  elif read.mapping_quality == 255:
                                      #if pos1based==25321928: print(skip_message.format(position, qname, 'no mapq'))
                                      nomapq_count[readbase]+=1
                                      continue
                                  # skip mean qscore of the read below threshold
                                  elif mean(read.query_qualities) < minimum_read_quality:
                                      #if pos1based==25321928: print(skip_message.format(position, qname, 'low read quality', mean(read.query_qualities)))
                                      lowreadq_count[readbase]+=1
                                      continue
                                  else: 
                                      # check for insertion
                                      if pileupread.is_refskip: #or pileupread.indel > 0:
                                          if read.query_qualities[qpos+1] >= minimum_base_quality:
                                              insbase=read.query_sequence[qpos+1]
                                              insbase=complement(insbase) if strand =='-' else insbase
                                              #if pos1based==25321928: print(skip_message.format(pos1based, qname, 'is ins', readbase))
                                              ins_count[insbase]+=1
                                              ins_n+=1
                                      # skip reads with a base quality below threshold
                                      if read.query_qualities[qpos] < minimum_base_quality:
                                          #if pos1based==25321928: print(skip_message.format(position, qname, 'low base quality', read.query_qualities[pileupread.query_position]))
                                          lowbaseq_count[readbase]+=1
                                          continue

                                  ### skip reads if misaligned
                                  if pos1based in [25321796,25321822,25321905,25321928] or \
                                     pos1based in [25648287,25648313,25648396,25648419]: ### hg19
                                      if reads_misaligned(pileupread, gbuild):
                                          continue
 

                                  ### print details of each passed read
                                  #print pass_message.format(position, qname, read.query_sequence[qpos])
                                  good_n+=1
                                  good_count[readbase]+=1
 
                              ### otuput the summary for each position of the exon
                              cdstable.loc[(cdstable['chr']==tarchr) & (cdstable['gpos']==pos1based) & (cdstable['cpos']==cpos) & (cdstable['gene']==row['gene']),['TOTAL_n', 'QC_n', 'DEL_n']]   = [column.n, good_n, del_n]
                              cdstable.loc[(cdstable['chr']==tarchr) & (cdstable['gpos']==ins_gpos) & (cdstable['cpos']==ins_cpos) & (cdstable['gene']==row['gene']),['INS_n']]   = [ins_n]
                              #cdstable.loc[(cdstable['chr']==tarchr) & (cdstable['gpos']==pos1based) & (cdstable['cpos']==cpos) & (cdstable['gene']==row['gene']),['INS_n']]   = [ins_n]
                              cdstable.loc[(cdstable['chr']==tarchr) & (cdstable['gpos']==pos1based) & (cdstable['cpos']==cpos) & (cdstable['gene']==row['gene']),['A','T','C','G']]     = [good_count['A'], good_count['T'], good_count['C'], good_count['G']]
                              if good_n > 0:
                                  cdstable.loc[(cdstable['chr']==tarchr) & (cdstable['gpos']==pos1based) & (cdstable['cpos']==cpos) & (cdstable['gene']==row['gene']),['Ap','Tp','Cp','Gp']] = [float(good_count['A'])/float(good_n),
                                                                                                                                                            float(good_count['T'])/float(good_n),
                                                                                                                                                            float(good_count['C'])/float(good_n),
                                                                                                                                                            float(good_count['G'])/float(good_n)]
                              cdstable.loc[(cdstable['chr']==tarchr) & (cdstable['gpos']==pos1based) & (cdstable['cpos']==cpos) & (cdstable['gene']==row['gene']),['dA','dT','dC','dG']] = [del_count['A'], del_count['T'], del_count['C'], del_count['G']]
                              if del_n > 0:
                                  cdstable.loc[(cdstable['chr']==tarchr) & (cdstable['gpos']==pos1based) & (cdstable['cpos']==cpos) & (cdstable['gene']==row['gene']),['dAp','dTp','dCp','dGp']] = [float(del_count['A'])/(float(del_n)+float(good_n)),
                                                                                                                                                                    float(del_count['T'])/(float(del_n)+float(good_n)),
                                                                                                                                                                    float(del_count['C'])/(float(del_n)+float(good_n)),
                                                                                                                                                                    float(del_count['G'])/(float(del_n)+float(good_n))]
                              cdstable.loc[(cdstable['chr']==tarchr) & (cdstable['gpos']==ins_gpos) & (cdstable['cpos']==ins_cpos) & (cdstable['gene']==row['gene']),['iA','iT','iC','iG']] = [ins_count['A'], ins_count['T'], ins_count['C'], ins_count['G']]
                              #cdstable.loc[(cdstable['chr']==tarchr) & (cdstable['gpos']==pos1based) & (cdstable['cpos']==cpos) & (cdstable['gene']==row['gene']),['iA','iT','iC','iG']] = [ins_count['A'], ins_count['T'], ins_count['C'], ins_count['G']]
                              if ins_n > 0:
                                  cdstable.loc[(cdstable['chr']==tarchr) & (cdstable['gpos']==ins_gpos) & (cdstable['cpos']==ins_cpos) & (cdstable['gene']==row['gene']),['iAp','iTp','iCp','iGp']] = [float(ins_count['A'])/float(good_n),
                                  #cdstable.loc[(cdstable['chr']==tarchr) & (cdstable['gpos']==pos1based) & (cdstable['cpos']==cpos) & (cdstable['gene']==row['gene']),['iAp','iTp','iCp','iGp']] = [float(ins_count['A'])/float(good_n),
                                                                                                                                                                    float(ins_count['T'])/float(good_n),
                                                                                                                                                                    float(ins_count['C'])/float(good_n),
                                                                                                                                                                    float(ins_count['G'])/float(good_n)]
                              ### assign alt allele
                              ALT='.'
                              if good_n > 0:
                                  cdstable.loc[(cdstable['chr']==tarchr) & (cdstable['gpos']==pos1based) & (cdstable['cpos']==cpos) & (cdstable['gene']==row['gene']),['max']]     = max(float(good_count['A'])/float(good_n),
                                                                                                                                                     float(good_count['T'])/float(good_n),
                                                                                                                                                     float(good_count['C'])/float(good_n),
                                                                                                                                                     float(good_count['G'])/float(good_n))
                                  ALT=rank_allele(cdstable.loc[(cdstable['chr']==tarchr) & (cdstable['gpos']==pos1based) & (cdstable['cpos']==cpos) & (cdstable['gene']==row['gene']),['Ap','Tp','Cp','Gp']], refbase)
                                  if str(ALT) == str(refbase):
                                      ALT='.'

                              if del_n >= ALT_n_LB:
                                  ALT='-'
                                  cdstable.loc[(cdstable['chr']==tarchr) & (cdstable['gpos']==pos1based) & (cdstable['cpos']==cpos) & (cdstable['gene']==row['gene']),['max_del']] = max(float(del_count['A'])/(float(del_n)+float(good_n)),
                                                                                                                                                     float(del_count['T'])/(float(del_n)+float(good_n)),
                                                                                                                                                     float(del_count['C'])/(float(del_n)+float(good_n)),
                                                                                                                                                     float(del_count['G'])/(float(del_n)+float(good_n)))
                              if column.n <= 3:
                                  ALT='-'
                              if ins_n >= ALT_n_LB:
                                  #print "Need to work on inertions variant calling, the REF allele not changed right now"
                                  cdstable.loc[(cdstable['chr']==tarchr) & (cdstable['gpos']==ins_gpos) & (cdstable['cpos']==ins_cpos) & (cdstable['gene']==row['gene']),['max_ins']] = max(float(ins_count['A'])/float(good_n),
                                  #cdstable.loc[(cdstable['chr']==tarchr) & (cdstable['gpos']==pos1based) & (cdstable['cpos']==cpos) & (cdstable['gene']==row['gene']),['max_ins']] = max(float(ins_count['A'])/float(good_n),
                                                                                                                                                     float(ins_count['T'])/float(good_n),
                                                                                                                                                     float(ins_count['C'])/float(good_n),
                                                                                                                                                     float(ins_count['G'])/float(good_n))
                              cdstable.loc[(cdstable['chr']==tarchr) & (cdstable['gpos']==pos1based) & (cdstable['cpos']==cpos) & (cdstable['gene']==row['gene']),['ALT']]     =ALT

                              ### translation
                              REFcodon='.'
                              ALTcodon='.'
                              if cpos % 3 == 0:
                                  #print "3rd base"
                                  REFcodon=''.join(map(str,[cdstable.loc[(cdstable['chr']==tarchr) & (cdstable['cpos']==cpos-2) & (cdstable['gene']==row['gene']),['REF']].values[0][0],
                                                    cdstable.loc[(cdstable['chr']==tarchr) & (cdstable['cpos']==cpos-1) & (cdstable['gene']==row['gene']),['REF']].values[0][0],
                                                    cdstable.loc[(cdstable['chr']==tarchr) & (cdstable['cpos']==cpos) & (cdstable['gene']==row['gene']),['REF']].values[0][0]]))
                                  ALTcodon=''.join(map(str,[cdstable.loc[(cdstable['chr']==tarchr) & (cdstable['cpos']==cpos-2) & (cdstable['gene']==row['gene']),['REF']].values[0][0],
                                                    cdstable.loc[(cdstable['chr']==tarchr) & (cdstable['cpos']==cpos-1) & (cdstable['gene']==row['gene']),['REF']].values[0][0],
                                                    ALT]))
                              elif cpos %3 == 1:
                                  #print "1st base"
                                  REFcodon=''.join(map(str,[cdstable.loc[(cdstable['chr']==tarchr) & (cdstable['cpos']==cpos) & (cdstable['gene']==row['gene']),['REF']].values[0][0],
                                                   cdstable.loc[(cdstable['chr']==tarchr) & (cdstable['cpos']==cpos+1) & (cdstable['gene']==row['gene']),['REF']].values[0][0],
                                                   cdstable.loc[(cdstable['chr']==tarchr) & (cdstable['cpos']==cpos+2) & (cdstable['gene']==row['gene']),['REF']].values[0][0]]))
                                  ALTcodon=''.join(map(str,[ALT,
                                                   cdstable.loc[(cdstable['chr']==tarchr) & (cdstable['cpos']==cpos+1) & (cdstable['gene']==row['gene']),['REF']].values[0][0],
                                                   cdstable.loc[(cdstable['chr']==tarchr) & (cdstable['cpos']==cpos+2) & (cdstable['gene']==row['gene']),['REF']].values[0][0]]))

                              elif cpos %3 == 2:
                                  #print "2nd base"
                                  REFcodon=''.join(map(str,[cdstable.loc[(cdstable['chr']==tarchr) & (cdstable['cpos']==cpos-1) & (cdstable['gene']==row['gene']),['REF']].values[0][0],
                                                   cdstable.loc[(cdstable['chr']==tarchr) & (cdstable['cpos']==cpos) & (cdstable['gene']==row['gene']),['REF']].values[0][0],
                                                   cdstable.loc[(cdstable['chr']==tarchr) & (cdstable['cpos']==cpos+1) & (cdstable['gene']==row['gene']),['REF']].values[0][0]]))
                                  ALTcodon=''.join(map(str,[cdstable.loc[(cdstable['chr']==tarchr) & (cdstable['cpos']==cpos-1) & (cdstable['gene']==row['gene']),['REF']].values[0][0],
                                                   ALT,
                                                   cdstable.loc[(cdstable['chr']==tarchr) & (cdstable['cpos']==cpos+1) & (cdstable['gene']==row['gene']),['REF']].values[0][0]]))
                              ## reverse complement the codon if minus strand
                              #REFcodon=complement(REFcodon) if strand == '-' else REFcodon
                              #ALTcodon=complement(ALTcodon) if strand == '-' else ALTcodon   

                              REFaa='.';ALTaa='.'
                              if '.' in  REFcodon:
                                  REFaa='?'
                              elif '-' in REFcodon:
                                  REFaa='del'
                              else:
                                  REFaa=translate(REFcodon)

                              if '-' in  ALTcodon:
                                  ALTaa='del'
                              elif '.' in ALTcodon:
                                  ALTaa='.'
                              elif 'N' in ALTcodon:
                                  ALTaa='*'
                              else:
                                  ALTaa=translate(ALTcodon)


                              cdstable.loc[(cdstable['chr']==tarchr) & (cdstable['gpos']==pos1based) & (cdstable['cpos']==cpos) & (cdstable['gene']==row['gene']),['REFaa']]     = REFaa
                              cdstable.loc[(cdstable['chr']==tarchr) & (cdstable['gpos']==pos1based) & (cdstable['cpos']==cpos) & (cdstable['gene']==row['gene']),['ALTaa']]     = ALTaa
        
        pd.options.display.max_columns = 50
        #pd.options.display.float_format = '{:,.4f}'.format
        #final = cdstable.loc[(cdstable['max'] < 1) | (cdstable['max_del'] > 0) | (cdstable['max_ins'] > 0) | (cdstable['ALT'] != '.')]
        filtered = cdstable.loc[(cdstable['ALT'] != '.')]
         
        #print final
        if len(cdstable) >0:
            out.append(cdstable)
        if len(filtered) >0 :
            out_filtered.append(filtered)
        CDSseq=CDSseq if strand == '+' else reverse(CDSseq)
        #if verbose > 10: print CDSseq

    outtotal=pd.concat(out, sort=True)
    #outtotal.apply(pd.to_numeric, errors='ignore')
    outtotal[['TOTAL_n', 'QC_n', 'DEL_n', 'INS_n', 'exon']] = outtotal[['TOTAL_n', 'QC_n', 'DEL_n', 'INS_n','exon']].apply(pd.to_numeric)
    out_filtered_total=pd.concat(out_filtered, sort=True)
    if outtotal.loc[outtotal['gene']=='ABO'].shape[0]>0:
        if verbose > 0: "[ABO] adjust 261 del"
        outtotal=ABO_261_del(outtotal) 
        out_filtered_total=ABO_261_del(out_filtered_total)

    outtotal.sort_values(by=['gene', 'cpos'], inplace=True)
    outtotal.reset_index(drop=True, inplace=True)
    out_filtered_total.sort_values(by=['gene', 'cpos'], inplace=True)
    out_filtered_total.reset_index(drop=True, inplace=True)
    column_order=['gene','acc','chr','exon','start','end','strand','gpos','cpos','REF','ALT','aapos','REFaa','ALTaa',\
                  'TOTAL_n','QC_n','DEL_n','INS_n',\
                  'A','T','C','G','Ap','Tp','Cp','Gp','dA','dT','dC','dG','dAp','dTp','dCp','dGp','iA','iT','iC','iG','iAp','iTp','iCp','iGp','max','max_del','max_ins']

    out_filtered_total[['TOTAL_n', 'QC_n', 'DEL_n', 'INS_n', 'exon']] = out_filtered_total[['TOTAL_n', 'QC_n', 'DEL_n', 'INS_n','exon']].apply(pd.to_numeric) 
    outtotal=outtotal[column_order]
    out_filtered_total=out_filtered_total[column_order]
    if int(outtotal['exon'].min())==0: 
        outtotal.exon=outtotal.exon+1
        out_filtered_total.exon=out_filtered_total.exon+1
    if write_output:
        outtotal.to_csv(output+'.full.variant.txt', sep='\t', index=False, mode='w')
        out_filtered_total.to_csv(output+'.filtered.variant.txt', sep='\t', index=False, mode='w')

    ### determine the final set of allele with polymorphsim
    final=[]
    for idx, row in out_filtered_total.iterrows():
       if row['ALTaa']=='del' or row['ALTaa']=='ins':
           final.append(row)
       else:
           alt=row['ALT']
           if out_filtered_total.at[idx, alt] < ALT_n_LB:
               continue
           if VAF is not None:
               if (out_filtered_total.at[idx, alt]/out_filtered_total.at[idx, 'QC_n']) < VAF:
                   continue 
           final.append(row)
    
    final=pd.DataFrame(final)
    if final.shape[0] > 0:
        final.sort_values(by=['gene', 'cpos'], inplace=True)
        final.reset_index(drop=True, inplace=True)
        final=final[column_order]
        final[['TOTAL_n', 'QC_n', 'DEL_n', 'INS_n', 'exon']] = final[['TOTAL_n', 'QC_n', 'DEL_n', 'INS_n','exon']].apply(pd.to_numeric)
        if int(outtotal['exon'].min())==0: final.exon=final.exon+1
        if write_output:
            final.to_csv(output+'.final.variant.txt', sep='\t', index=False, mode='w')

    return outtotal, out_filtered_total, final

def ABO_261_del(df):
    '''
       adjust 261 indel for HG38
    '''
    gene='ABO'
    pos=261
    insGr=float(df.loc[(df['cpos']==pos) & (df['gene']==gene),['iGp']].values[0][0])
    insN=float(df.loc[(df['cpos']==pos) & (df['gene']==gene),['INS_n']].values[0][0])
    good_n=int(df.loc[(df['cpos']==pos) & (df['gene']==gene),['QC_n']].values[0][0])
    if insGr == 0:
        print("261G full deleted")
        #df.loc[(df['cpos']==pos) & (df['gene']==gene),['REFaa']]     = REFaa
        df.loc[(df['cpos']==pos) & (df['gene']==gene),['ALT']]       = '-'
        df.loc[(df['cpos']==pos) & (df['gene']==gene),['ALTaa']]     = 'del'
        ###  change del proportion
        Tn=df.loc[(df['cpos']==pos) & (df['gene']==gene),['T']].values[0][0]
        Gn=df.loc[(df['cpos']==pos) & (df['gene']==gene),['G']].values[0][0]
        df.loc[(df['cpos']==pos) & (df['gene']==gene),['DEL_n']]=Tn-Gn ### this value is important for matrix calculation
        df.loc[(df['cpos']==pos) & (df['gene']==gene),['dG']]=Tn-Gn
        df.loc[(df['cpos']==pos) & (df['gene']==gene),['dGp']]=float((Tn-Gn)/good_n) 
        df.loc[(df['cpos']==pos) & (df['gene']==gene),['T']]=0                                           
        df.loc[(df['cpos']==pos) & (df['gene']==gene),['Tp']]=0 
        df.loc[(df['cpos']==pos) & (df['gene']==gene),['max_del']]=float((Tn-Gn)/good_n)
        ### adjust other codon bases
        df.loc[(df['cpos']==pos-1) & (df['gene']==gene),['ALTaa']]   = 'del' 
        df.loc[(df['cpos']==pos-2) & (df['gene']==gene),['ALTaa']]   = 'del'
    elif insGr > 0:
        print("261G partially deleted")
        df.loc[(df['cpos']==pos) & (df['gene']==gene),['REFaa']]     = 'V'
        df.loc[(df['cpos']==pos) & (df['gene']==gene),['ALT']]       = '-'
        df.loc[(df['cpos']==pos) & (df['gene']==gene),['ALTaa']]     = 'del'
        Tn =int(df.loc[(df['cpos']==pos) & (df['gene']==gene),['T']].values[0][0])
        Gn =int(df.loc[(df['cpos']==pos) & (df['gene']==gene),['iG']].values[0][0])
        df.loc[(df['cpos']==pos) & (df['gene']==gene),['DEL_n']]=good_n-Gn
        df.loc[(df['cpos']==pos) & (df['gene']==gene),['G']]=Gn
        df.loc[(df['cpos']==pos) & (df['gene']==gene),['T', 'Tp']]=0
        df.loc[(df['cpos']==pos) & (df['gene']==gene),['Gp']]=Gn/float(good_n)
        df.loc[(df['cpos']==pos) & (df['gene']==gene),['dG']]=good_n-Gn
        df.loc[(df['cpos']==pos) & (df['gene']==gene),['dGp']]=(good_n-Gn)/float(good_n)
        df.loc[(df['cpos']==pos) & (df['gene']==gene),['iG', 'iGp']]=0
        df.loc[(df['cpos']==pos) & (df['gene']==gene),['max']]=Gn/float(good_n)
        df.loc[(df['cpos']==pos) & (df['gene']==gene),['max_ins']]=0
        df.loc[(df['cpos']==pos) & (df['gene']==gene),['max_del']]=(good_n-Gn)/float(good_n)
        df.loc[(df['cpos']==pos) & (df['gene']==gene),['INS_n']]=0
        ### adjust other codon bases
        df.loc[(df['cpos']==pos-1) & (df['gene']==gene),['ALTaa']]   = 'del' 
        df.loc[(df['cpos']==pos-2) & (df['gene']==gene),['ALTaa']]   = 'del'
    return df

def rank_allele(ls, ref=None):
    '''
        rank the allele frequency to predict ALT allele
    '''

    ### transform pandas column to dictionary
    ls = ls.to_dict('records')
    if len(ls) == 1:
        ls_sorted = OrderedDict(sorted(ls[0].items(), key=lambda value: value[1], reverse=True))
        for k in list(ls_sorted.keys())[:2]:
            base=re.sub('[idp]','',str(k))
            if ref:
                if str(ref) != base:
                    return base if ls_sorted[k] > 0 else '.'
            else:
                return base if ls_sorted[k] > 0 else '.'

def reads_misaligned(pileupread, gbuild):
        '''
            check whether the read is misaligned given a piled up read info
        '''
        check_pos={'hg38':{25321796:{'ref':'T','alt':'C'},
                           25321822:{'ref':'C','alt':'T'},
                           25321905:{'ref':'T','alt':'C'},
                           25321928:{'ref':'A','alt':'T'}},
                   'hg19':{25648287:{'ref':'T','alt':'C'},
                           25648313:{'ref':'C','alt':'T'},
                           25648396:{'ref':'T','alt':'C'},
                           25648419:{'ref':'A','alt':'T'}}}
        qname    = pileupread.alignment.query_name
        qpos     = pileupread.query_position
        qaln     = pileupread.alignment.query_sequence
        qqual    = pileupread.alignment.query_qualities
        qrpos_0  = pileupread.alignment.get_reference_positions()
        qrpos_1  = [ x+1 for x in qrpos_0 ]
        qmap     = zip(qrpos_1, qaln, qqual)
        misalnbN=0
        alnbN=0
        for p in qmap:
            if p[0] in check_pos[gbuild]:
                alnbN+=1
                if p[1]==check_pos[gbuild][p[0]]['alt']:
                    misalnbN+=1
        if misalnbN > 1:
            print(qname, 'misaligned; mis-matched base N:', misalnbN, '/', alnbN)
            return True
        return False


