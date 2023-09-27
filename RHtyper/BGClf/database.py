"""
Database module
"""
import pandas as pd
import re, os

package_directory = os.path.dirname(os.path.abspath(__file__))

class bloodgroupDB(object):
    def __init__(self, gene='RHD', DB="ISBT", verbose=0):
       
        self.BGMUTf=os.path.join(package_directory, 'database', 'bgmut_data_export-2017-10-20-bgmut_alleles.txt')
        self.BGMUTf2=os.path.join(package_directory, 'database', 'bgmut_data_export-2017-10-20-bgmut_alleles_submitted.txt')
        self.ISBT=os.path.join(package_directory, 'database', 'ISBT.tsv')
        #self.DB=DB
        self.gene=gene
        self.verbose=verbose

    def BloodDB(self):
        ''' deprecated'''
        db=None
        if self.DB=='BGMUT':
            db_in1=self.BGMUTf
            db_in2=self.BGMUTf2
            db1   =pd.read_csv(db_in1, sep='|', header=0)
            db2   =pd.read_csv(db_in2, sep='|', header=0)
            ### add more variations from other db
            extra_db=pd.DataFrame([['RHCE', 'ce254G', '254C>G'], \
                                   ['RHD','DUC-3','48G>C'], \
                                   ['RHCE','48C_105T','48G>C;105C>T']], 
                                  columns=['gene', 'name', 'nt_change'])
            frames= [db1, db2, extra_db]
            ### put all together
            db    =pd.concat(frames, sort=True)
            db    =db.drop_duplicates()
            ### remove weird entries
            db.drop(db.loc[(db['gene']=='ABO') & (db['alias']=='O01 intronic')].index, inplace=True)
        elif self.DB=="ISBT":
            db_in=self.ISBT
            db1   =pd.read_csv(db_in, sep='\t', header=0)
            extra_db=pd.DataFrame([['RHD', 'RHD*XX.XX RHD(L390L)', 'c.1170T>C','RHD(L390L)'], \
                                   ['RHD', 'RHD*XX.XX RHD648D', 'c.648G>C','RHD648D'], \
                                   ['RHCE', 'RHCE*01.XX RHce(48C,105T)','c.48G>C;c.105C>T','RHce(48C,105T)'], \
                                   ['RHD', 'RHD*XX.XX RHD674T','c.674C>T','RHD674T'], \
                                   ['RHD', 'RHD*XX.XX RHDIIIa(2)', 'c.186G>T; c.410C>T; c.455A>C; c.602C>G; c.667T>G', 'RHDIIIa(2)'], \
                                   ['RHD', 'RHD*XX.XX RHD841C' , 'c.841G>C','RHD841C']
                                  ],
                                  columns=['Gene','Allele name','Nucleotide','Alias'])
            db    =pd.concat([db1, extra_db], sort=True)
            db    =db.drop_duplicates()
        else:
            print("[Error] Unknown database")

        return db

    def isbt_db(self, include_non_mut=False):
        self.DB='ISBT'
        gene=self.gene
        Bdb=self.BloodDB()
        verbose=self.verbose
        RhBdb=Bdb.loc[ (Bdb['Gene']==gene), ['Gene', 'Allele name', 'Nucleotide','Amino Acid', 'rs number', 'Allele detail', 'Exon', 'Alias']]
        if verbose > 9: RhBdb.to_csv(gene+'raw_database'+'.txt')
        nuclbase=['A', 'T', 'C', 'G']
        nt_pattern=re.compile("(C.\s?[0-9]+\s?[A,T,C,G]\s?>\s?[A,T,C,G])")
        indel_pattern=re.compile("(c.[0-9]*_?[0-9]+\s?(del|ins)\s?[A,T,C,G]*)")
        data = []
        prob = []
        for idx, row in RhBdb.iterrows():
            nts=re.split('; |,|;', row['Nucleotide'].strip())
            for nt in nts:
                mut=nt.strip()
                if len(mut) ==0: continue
                if mut.startswith(('-', '~', '*', '+')):
                    if verbose > 9: print('\t\t'.join([row['Gene'], row['Allele name'], row['Nucleotide'], mut, 'unknown variant format... skip'])) 
                    continue
                if not mut.startswith('c.'):
                    if verbose > 9: print('\t\t'.join([row['Gene'], row['Allele name'], row['Nucleotide'], mut, 'non-coding vars ... skip']))
                    continue

                pattern=-1
                if nt_pattern.match(mut.upper()):
                    pattern=0
                    mut=mut.upper()
                    #if verbose > 9: print '\t\t'.join([row['Gene'], row['Allele name'], mut, 'SNPs ... process'])
                elif indel_pattern.match(mut):
                    pattern=1
                    #if verbose > 9: print '\t\t'.join([row['Gene'], row['Allele name'], row['Nucleotide'], mut, 'Indels ... process'])
                else:
                    if 'Intr' in row['Exon']:
                        pattern=-1
                        #if verbose > 9: print('\t\t'.join([row['Gene'], row['Allele name'], row['Nucleotide'], mut, 'Intron skipped... skip']))
                    else:
                        pattern=-1
                        #if verbose > 9: print('\t\t'.join([row['Gene'], row['Allele name'], row['Nucleotide'], mut, 'Intron or splicing ... skip']))

                if pattern ==-1:
                    #if verbose > 9: print '\t\t'.join([row['gene'], row['name'], row['nt_change'], mut, 'unknown variant format... skip'])
                    continue
                else:
                    if pattern==0:
                        b0, t, pos, ref, alt, b1=re.split(r'([c,C].)\s?(\d+)\s?([A-Z,a-z])\s?>\s?([A-Z,a-z])\s?', mut)
                        #print '\t\t'.join([row['Gene'], row['Allele name'], row['Nucleotide'], mut, pos, ref, alt])
                        
                        data.append({'gene':row['Gene'], 'name':row['Allele name'], 'pos' :pos, 'ref' :ref, 'alt' :alt, 'detail':row['Allele detail'] })
                        ### calculate frequency
                        for n in nuclbase:
                            for n2 in nuclbase:
                                if n != n2:
                                    mut = 1 if ref == n and alt == n2 else 0
                                    #print '\t\t'.join(map(str, [row['gene'], row['name'], pos, n, n2, mut]))
                                    prob.append({'gene' : row['Gene'], 'name' : row['Allele name'], 'pos'  : pos, 'ref'  : n, 'alt'  : n2, 'mut'  : mut,  'detail':row['Allele detail']})
                    elif pattern==1:
                        if 'del' in mut.lower():
                            b0, t, pos, ref, b1=re.split(r'([c,C].)\s?(.*)\s?del\s?([A-Z]*)\s?', mut)
                        
                            if '_' in pos:
                                start, end =map(int, re.split('_', pos))
                                refn=0
                                for p in range(start, end+1):
                                    if refn < len(ref):
                                        #print '\t\t'.join(map(str,[row['Gene'], row['Allele name'], row['Nucleotide'], mut, str(p), ref[refn], '-']))
                                        data.append({'gene':row['Gene'], 'name':row['Allele name'], 'pos' :p, 'ref' :ref[refn], 'alt' :'-' , 'detail':row['Allele detail']})
                                        ### calculate frequency
                                        for n in nuclbase:
                                            mut = 1 if ref[refn] == n else 0
                                            prob.append({'gene' : row['Gene'], 'name' : row['Allele name'], 'pos'  : p, 'ref'  : n, 'alt'  : '-', 'mut'  : mut, 'detail':row['Allele detail']})
                                    refn+=1
                            else:
                                #print '\t\t'.join([row['Gene'], row['Allele name'], row['Nucleotide'], mut, pos, ref, '-'])
                                start, end=int(pos), int(pos)+len(ref)
                                refn=0
                                for p in range(start, end):
                                    data.append({'gene':row['Gene'], 'name':row['Allele name'], 'pos' :p, 'ref' :ref[refn], 'alt' :'-' , 'detail':row['Allele detail']})
                                    ### calculate frequency
                                    for n in nuclbase:
                                        mut = 1 if ref[refn] == n else 0
                                        prob.append({'gene' : row['Gene'], 'name' : row['Allele name'], 'pos'  : p, 'ref'  : n, 'alt'  : '-', 'mut'  : mut, 'detail':row['Allele detail']})
                                    refn+=1
                        
                        elif 'ins' in mut.lower():
                            b0, t, pos, alt, b1=re.split(r'([c,C].)\s?([0-9]*_?[0-9]+)\s?ins\s?([A-Z]*)\s?', mut)
                            
                            if '_' in pos:
                                start, end =map(int, re.split('_', pos))
                                altn=0
                                for p in range(start, end+1):
                                    if altn<len(alt):
                                        #print '\t\t'.join([row['gene'], row['name'], row['nt_change'], mut, str(p), '-', alt[altn]])
                                        data.append({'gene':row['Gene'],'name':row['Allele name'], 'pos' :p, 'ref' :'-', 'alt' :alt[altn] , 'detail':row['Allele detail']})
                                        ### calculate frequency
                                        prob.append({'gene' : row['Gene'], 'name' : row['Allele name'], 'pos'  : p, 'ref'  : '-', 'alt'  : alt[altn], 'mut'  : 1 , 'detail':row['Allele detail']})
                                    altn+=1
                            else:
                                #print '\t\t'.join([row['gene'], row['name'], row['nt_change'], mut, pos, ref, '-'])
                                start, end=int(pos), int(pos)+len(alt)
                                altn=0
                                for p in range(start, end):
                                    data.append({'gene':row['Gene'], 'name':row['Allele name'], 'pos' :p, 'ref' :'-', 'alt' :alt[altn] , 'detail':row['Allele detail']})
                                    ### calculate frequency
                                    for n in nuclbase:
                                        prob.append({'gene' : row['Gene'], 'name' : row['Allele name'], 'pos'  : p, 'ref'  : '-', 'alt'  : alt[altn], 'mut'  : 1, 'detail':row['Allele detail']})
                                    altn+=1
                        else:
                            if verbose > 9: print('\t\t'.join([row['Gene'], row['Allele name'], row['Allele detail'], mut, 'Unknown type']))
        data=pd.DataFrame(data)
        prob=pd.DataFrame(prob)
        prob["index"] = prob["gene"].map(str) + '@' + prob["name"].map(str).str.replace('\s+','_')
        prob=prob.pivot_table(index=['index'], columns=['pos','ref','alt'], values=['mut'], aggfunc=max).fillna(0)
        #prob=prob.pivot_table(index=['gene','name'], columns=['pos','ref','alt'], values=['mut'], aggfunc=max).fillna(0)
        prob.columns = prob.columns.droplevel(level=0)
        prob.columns = ["_".join(map(str,x)) for x in prob.columns.ravel()]
        prob=prob.loc[:, (prob != 0).any(axis=0)] ### remove position with non information
        cnt_wide=prob.copy()
        if not include_non_mut:
            ### add reference back with all zero
            ### [NOTE] matrixprep.py need to change at the same time
            cnt_wide=cnt_wide.transpose()
            if gene=='RHCE':
                #cnt_wide[gene+'@ce']=0
                cnt_wide[gene+'@RHCE*01_or_RHCE*ce_RHCE*c_RHCE*e']=0
            elif gene=='ABO':
                cnt_wide[gene+'@ABO*A1.01']=0
            else:
                #cnt_wide[gene+'@Reference']=0
                cnt_wide[gene+'@RHD*01']=0
            cnt_wide=cnt_wide.transpose()

        pcnt_wide=cnt_wide.apply(lambda x: x/x.sum(), axis=0) ### calculate percentage
        #pcnt_wide=cnt_wide.apply(lambda x: x/sum(x>0), axis=1) ### percentage divided by marker in each BG
        #pcnt_wide=cnt_wide

        if verbose > 20:
            cnt_wide.to_csv('test.database.txt', sep='\t', index=True, mode='w')
            pcnt_wide.to_csv('test.database.pcnt.txt', sep='\t', index=True, mode='w')
            data.to_csv('test.database2.txt', sep='\t', index=False, mode='w')
        return pcnt_wide, cnt_wide, RhBdb#, pcnt_pos

                        
                       

    def bgmut_db(self, include_non_mut=False):
        self.DB='BGMUT'
        gene=self.gene
        verbose=self.verbose
        Bdb=self.BloodDB()
        RhBdb=Bdb.loc[ (Bdb['gene']==gene), ['gene', 'name', 'aa_change','nt_change', 'genbank_id']]
        ### ignore hybrid ones for now
        if verbose > 3: print ('Hybrid nt_change ignore for now')
        RhBdb=RhBdb.loc[~RhBdb['nt_change'].str.contains('hybrid')]
        RhBdb=RhBdb.groupby(['gene','nt_change']).agg(lambda x: '::'.join(list(set(x)))).reset_index()

        if verbose > 9: RhBdb.to_csv(gene+'raw_database'+'.txt')

        nuclbase=['A', 'T', 'C', 'G']
        nt_pattern=re.compile("([0-9]+[A,T,C,G]>[A,T,C,G])")
        nt_pattern2=re.compile("([0-9]+[A,T,C,G]/[A,T,C,G])")
        nt_pattern3=re.compile("(nt[0-9]+[A,T,C,G]>[A,T,C,G])")
        indel_pattern=re.compile("([0-9]*-?[0-9]+(del|ins)[A,T,C,G]*)")
        indel_pattern2=re.compile("([0-9]*-?[0-9]+[A,T,C,G]*?(del|ins))")
        data = []
        prob = []
        for idx, row in RhBdb.iterrows():
            #print "This is me", row['nt_change']
            nts=re.split(';|,', row['nt_change'].strip())
            for nt in nts:
                mut=nt.strip()
                if len(mut) ==0: continue
                if mut.startswith(('-', '~', '*', '+')):
                   #if verbose > 9: print '\t\t'.join([row['gene'], row['name'], row['nt_change'], mut, 'unknown variant format... skip'])
                   continue
                if 'intron' in mut.lower():
                   #if verbose > 9: print '\t\t'.join([row['gene'], row['name'], row['nt_change'], mut, 'non-coding vars ... skip'])
                   continue
                pattern=-1
                mut=re.sub("[Ee]xon [0-9]+", "", mut)
                mut=re.sub("[Ee]xon\s?[0-9]+:", "", mut)
                mut=re.sub(' \(reading frame is restored by the latter deletion\)', "", mut)
                mut=re.sub(' or mutation \+ del',"",mut)
                mut=re.sub('[() ]','',mut)
                if nt_pattern.match(mut.upper()):
                    pattern=0
                    mut=mut.upper()
                    #if verbose > 9: print '\t\t'.join([row['gene'], row['name'], mut, 'SNPs ... process'])
                elif nt_pattern2.match(mut.upper()):
                    mut=mut.replace('/','>').upper()
                    pattern=0
                    #if verbose > 9: print '\t\t'.join([row['gene'], row['name'], mut, 'SNPs ... adjust format ... process'])
                elif nt_pattern3.match(mut):
                    mut=mut[2:]
                    pattern=0
                   #if verbose > 9: print '\t\t'.join([row['gene'], row['name'], mut, 'SNPs ... adjust format ... process'])
                elif indel_pattern.match(mut):
                    pattern=1
                    #if verbose > 9: print '\t\t'.join([row['gene'], row['name'], row['nt_change'], mut, 'Indels ... process'])
                elif indel_pattern2.match(mut):
                    pattern=2
                    #if verbose > 9: print '\t\t'.join([row['gene'], row['name'], mut, 'Indels ... process'])
                if pattern ==-1:
                    #if verbose > 9: print '\t\t'.join([row['gene'], row['name'], row['nt_change'], mut, 'unknown variant format... skip'])
                    continue
                else:
                    if pattern==0:
                        b0, pos, ref, alt, b1=re.split(r'(\d+)([A-Z,a-z])>([A-Z,a-z])', mut)
                        #print '\t\t'.join([row['gene'], row['name'], row['nt_change'], mut, pos, ref, alt])
                        data.append({'gene':row['gene'], 'name':row['name'], 'pos' :pos, 'ref' :ref, 'alt' :alt })
                        ### calculate frequency
                        for n in nuclbase:
                            for n2 in nuclbase:
                                if n != n2:
                                    mut = 1 if ref == n and alt == n2 else 0
                                    #print '\t\t'.join(map(str, [row['gene'], row['name'], pos, n, n2, mut]))
                                    prob.append({'gene' : row['gene'], 'name' : row['name'], 'pos'  : pos, 'ref'  : n, 'alt'  : n2, 'mut'  : mut })
                    elif pattern==1:
                        if 'del' in mut.lower():
                            b0, pos, ref, b1=re.split(r'(.*)del([A-Z]*)', mut)
                            if '-' in pos:
                                start, end =map(int, re.split('-', pos))
                                refn=0
                                for p in range(start, end+1):
                                    if refn < len(ref):   
                                        #print '\t\t'.join(map(str,[row['gene'], row['name'], row['nt_change'], mut, str(p), ref[refn], '-']))
                                        data.append({'gene':row['gene'], 'name':row['name'], 'pos' :p, 'ref' :ref[refn], 'alt' :'-' })
                                        ### calculate frequency
                                        for n in nuclbase:
                                            mut = 1 if ref[refn] == n else 0
                                            prob.append({'gene' : row['gene'], 'name' : row['name'], 'pos'  : p, 'ref'  : n, 'alt'  : '-', 'mut'  : mut})
                                    refn+=1
                            else:
                                #print '\t\t'.join([row['gene'], row['name'], row['nt_change'], mut, pos, ref, '-'])
                                start, end=int(pos), int(pos)+len(ref)
                                refn=0
                                for p in range(start, end):
                                    data.append({'gene':row['gene'], 'name':row['name'], 'pos' :p, 'ref' :ref[refn], 'alt' :'-' })
                                    ### calculate frequency
                                    for n in nuclbase:
                                        mut = 1 if ref[refn] == n else 0
                                        prob.append({'gene' : row['gene'], 'name' : row['name'], 'pos'  : p, 'ref'  : n, 'alt'  : '-', 'mut'  : mut})
                                    refn+=1
                        elif 'ins' in mut.lower():
                            if len(mut)==0: continue
                            b0, pos, alt, b1=re.split(r'([0-9]*-?[0-9]+)ins([A-Z]*)', mut)
                            if '-' in pos:
                                start, end =map(int, re.split('-', pos))
                                altn=0
                                for p in range(start, end+1):
                                    if altn<len(alt):
                                        #print '\t\t'.join([row['gene'], row['name'], row['nt_change'], mut, str(p), '-', alt[altn]])
                                        data.append({'gene':row['gene'],'name':row['name'], 'pos' :p, 'ref' :'-', 'alt' :alt[altn] })
                                        ### calculate frequency
                                        prob.append({'gene' : row['gene'], 'name' : row['name'], 'pos'  : p, 'ref'  : '-', 'alt'  : alt[altn], 'mut'  : 1 })
                                    altn+=1
                            else:
                                #print '\t\t'.join([row['gene'], row['name'], row['nt_change'], mut, pos, ref, '-'])
                                start, end=int(pos), int(pos)+len(alt)
                                altn=0
                                for p in range(start, end):
                                    data.append({'gene':row['gene'], 'name':row['name'], 'pos' :p, 'ref' :'-', 'alt' :alt[altn] })
                                    ### calculate frequency
                                    for n in nuclbase:
                                        prob.append({'gene' : row['gene'], 'name' : row['name'], 'pos'  : p, 'ref'  : '-', 'alt'  : alt[altn], 'mut'  : 1})
                                    altn+=1
                        else:
                            if verbose > 9:print('\t\t'.join([row['gene'], row['name'], mut, 'Unknown type']))
                        
                    elif pattern==2:
                        if 'del' in mut.lower():
                            b0, pos, ref, b1=re.split(r'([0-9]*-?[0-9]+)([A-Z]*)del', mut)                            
                            if '-' in pos:
                                start, end =map(int, re.split('-', pos))
                                refn=0
                                for p in range(start, end+1):
                                    if refn < len(ref):
                                        #print '\t\t'.join(map(str,[row['gene'], row['name'], row['nt_change'], mut, str(p), ref[refn], '-']))
                                        data.append({'gene':row['gene'],'name':row['name'],'pos' :p,'ref' :ref[refn], 'alt' :'-'})
                                        ### calculate frequency
                                        for n in nuclbase:
                                            mut = 1 if ref[refn] == n else 0
                                            prob.append({'gene' : row['gene'], 'name' : row['name'], 'pos'  : p,'ref'  : n, 'alt'  : '-', 'mut'  : mut })
                                    refn+=1
                            else:
                                #print '\t\t'.join([row['gene'], row['name'], row['nt_change'], mut, pos, ref, '-'])
                                start, end=int(pos), int(pos)+len(ref)
                                refn=0
                                for p in range(start, end):
                                    data.append({'gene':row['gene'], 'name':row['name'], 'pos' :p, 'ref' :ref[refn], 'alt' :'-' })
                                    ### calculate frequency
                                    for n in nuclbase:
                                        mut = 1 if ref[refn] == n else 0
                                        prob.append({'gene' : row['gene'], 'name' : row['name'], 'pos'  : p, 'ref'  : n, 'alt'  : '-', 'mut'  : mut})
                                    refn+=1

                        elif 'ins' in mut.lower():
                            b0, pos, alt, b1=re.split(r'([0-9]*-?[0-9]+)ins([A-Z]*)', mut)             
                            if '-' in pos:
                                    #print '\t\t'.join([row['gene'], row['name'], row['nt_change'], mut, pos,'-', alt])
                                    start, end =map(int, re.split('-', pos))
                                    altn=0
                                    for p in range(start, end+1):
                                        if altn<len(alt):
                                            if verbose >0: '\t\t'.join([row['gene'], row['name'], row['nt_change'], mut, str(p), '-', alt[altn]])
                                            data.append({'gene':row['gene'], 'name':row['name'], 'pos' :p, 'ref' :'-', 'alt' :alt[altn]})
                                            ### calculate frequency
                                            prob.append({'gene' : row['gene'], 'name' : row['name'], 'pos'  : p, 'ref'  : '-', 'alt'  : alt[altn], 'mut'  : 1})
                                        altn+=1
                            else:
                                #print '\t\t'.join([row['gene'], row['name'], row['nt_change'], mut, pos, ref, '-'])
                                start, end=int(pos), int(pos)+len(alt)
                                altn=0
                                for p in range(start, end):
                                    data.append({'gene':row['gene'], 'name':row['name'], 'pos' :p, 'ref' :'-', 'alt' :alt[altn] })
                                    ### calculate frequency
                                    for n in nuclbase:
                                        prob.append({'gene' : row['gene'], 'name' : row['name'], 'pos'  : p, 'ref'  : '-', 'alt'  : alt[altn], 'mut'  : 1})
                                    altn+=1
                        else:
                            if verbose > 9: print('\t\t'.join([row['gene'], row['name'], mut, 'Unknown type']))
                            continue 
        data=pd.DataFrame(data)
        prob=pd.DataFrame(prob)
        prob["index"] = prob["gene"].map(str) + '@' + prob["name"].map(str).str.replace('\s+','_')
        prob=prob.pivot_table(index=['index'], columns=['pos','ref','alt'], values=['mut'], aggfunc=max).fillna(0)
        #prob=prob.pivot_table(index=['gene','name'], columns=['pos','ref','alt'], values=['mut'], aggfunc=max).fillna(0)
        prob.columns = prob.columns.droplevel(level=0)
        prob.columns = ["_".join(map(str,x)) for x in prob.columns.ravel()]
        prob=prob.loc[:, (prob != 0).any(axis=0)] ### remove position with non information
        cnt_wide=prob.copy()
        if not include_non_mut:
            ### add reference back with all zero
            ### [NOTE] matrixprep.py need to change at the same time
            cnt_wide=cnt_wide.transpose()
            if gene=='RHCE':
                #cnt_wide[gene+'@ce']=0
                cnt_wide[gene+'@RHCE*01_or_RHCE*ce_RHCE*c_RHCE*e']=0
            elif gene=='ABO':
                cnt_wide[gene+'@ABO*A1.01']=0
            else:
                #cnt_wide[gene+'@Reference']=0
                cnt_wide[gene+'@RHD*01']=0
            cnt_wide=cnt_wide.transpose()

        pcnt_wide=cnt_wide.apply(lambda x: x/x.sum(), axis=0) ### calculate percentage
        #pcnt_wide=cnt_wide.apply(lambda x: x/sum(x>0), axis=1) ### percentage divided by marker in each BG
        #pcnt_wide=cnt_wide

        if verbose > 20:
            cnt_wide.to_csv('test.database.txt', sep='\t', index=True, mode='w')
            pcnt_wide.to_csv('test.database.pcnt.txt', sep='\t', index=True, mode='w')
            data.to_csv('test.database2.txt', sep='\t', index=False, mode='w')
        return pcnt_wide, cnt_wide#, pcnt_pos


