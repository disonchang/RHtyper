"""
    modules to produce a matrix with all nucleotides and translation at each position

"""
import pandas as pd
import re, collections, pysam, sys, os
from collections import OrderedDict
from general import complement

package_directory = os.path.dirname(os.path.abspath(__file__))


def gene_coordinates(gene, build='hg38'):
    """
       Dict for gene positions/coordinates
    """
    Hg38={#'RHD':'BGClf/coordinates/RHD_RHCE_UCSC_gene_table.txt',
          #'RHCE':'BGClf/coordinates/RHD_RHCE_UCSC_gene_table.txt',
          'RHD':os.path.join(package_directory,'coordinates', 'RHD.txt'),
          'RHCE':os.path.join(package_directory,'coordinates','RHCE.txt'),
          'ABO':os.path.join(package_directory,'coordinates','ABO.txt')}
    Hg19={
          'RHD':os.path.join(package_directory,'coordinates', 'RHD.hg19.txt'),
          'RHCE':os.path.join(package_directory,'coordinates','RHCE.hg19.txt'),
          'ABO':os.path.join(package_directory,'coordinates','ABO.hg19.txt')}

    if build.upper()=='HG38':
        return Hg38[gene.upper()]
    elif build.upper()=='HG19':
        return Hg19[gene.upper()]


def genetable(fn):
    df=pd.read_csv(fn, sep="\s+", header=None)
    df.columns = ['chr', 'start', 'end', 'exon', 'frame', 'strand']
    return df

def format_table(df, gene):
    '''
       reformat the table into an appropriate format
    '''
    ### exon
    df['acc']='.'
    df['exon_n']='0'
    df['exon_size']='0'
    df['gene']='.'
    for idx, row in df.iterrows():
        df.loc[idx,'acc']=re.sub(r"(NM_[0-9]+.*?)_.*", r"\1", row['exon'])
        df.loc[idx,'exon_n']=re.sub(r".*_cds_([0-9]+)_.*", r"\1", row['exon'])
        df.loc[idx,'exon_size']=row['end']-row['start']
        df.loc[idx,'gene']=gene.upper()
        df.loc[idx,'start'] = df.loc[idx,'start'] + 1 ### UCSC 0-based
    df=df.sort_values(by=['acc', 'start'])
    intron=[]
    for idx, row in df.iterrows():
        intron_row=row.copy()
        if (idx+1) < df.shape[0]:
            if (df.loc[idx,'acc']==df.loc[idx+1,'acc']):
                intron_row['start']=df.loc[idx,'end']+1
                intron_row['end']=df.loc[idx+1,'start']-1
                intron_row['intron_n']=int(df.loc[idx,'exon_n'])+1
                intron_row['intron_size']=intron_row['end']-intron_row['start']
                intron.append(intron_row)
    intron=pd.DataFrame(intron)
    intron.drop(labels=['exon', 'exon_n','exon_size'], axis=1, inplace=True)

    ### put each gene coordinates to dict
    final={}
    final_intron={}
    d=df.sort_values(by=['start'])
    d.reset_index(drop=True, inplace=True)
    i=intron.sort_values(by=['start'])
    i=i.sort_values(by=['start'])
    i.reset_index(drop=True, inplace=True)
    a=df.loc[0,'acc']
    final[a]=d
    final_intron[a]=i
    return final, final_intron

def CDS_genome_coordinate(df):
    '''
        generate coorespondent CDS position
    '''
    final=collections.OrderedDict()
    for acc in df:
        d=df[acc]
        gene=d.gene.unique()
        strand=d.strand.unique()
        chr=d.chr.unique()
        if len(gene)==1 and len(strand)==1 and len(chr)==1:
             gene=gene[0]; strand=strand[0]; chr=chr[0]
        else:
             sys.exit('[Error] mutlitple genes in one df: %s' % ','.join(map(str,gene)))
        cds=collections.OrderedDict()
        cds_pos=1
        if strand=='-':
            d=d.sort_values('end', ascending=False)
            d=d.reset_index(drop=True)
            for idx, row in d.iterrows():
                cds[row['exon_n']]=collections.OrderedDict()
                cds[row['exon_n']]['start']=int(row['start'])
                cds[row['exon_n']]['end']=int(row['end'])
                cds[row['exon_n']]['cpos']=collections.OrderedDict()
                for genome_pos in range(int(row['end']), int(row['start'])-1, -1):
                    if genome_pos not in cds[row['exon_n']]['cpos']:
                        cds[row['exon_n']]['cpos'][genome_pos]=cds_pos
                        cds_pos+=1
                        if chr=='chr9' or chr==9:
                            if genome_pos==133257522:
                                ### This is the ABO where a Gdel occured in the Reference HG38
                                ### inrease the CDS pos by one for insertion a ref base
                                cds_pos+=1
                    else:
                        sys.exit('[Error] mutlitple position mapped: %d' % genome_pos)
        else:
            d=d.reset_index(drop=True)
            for idx, row in d.iterrows():
                cds[row['exon_n']]=collections.OrderedDict()
                cds[row['exon_n']]['start']=int(row['start'])
                cds[row['exon_n']]['end']=int(row['end'])
                cds[row['exon_n']]['cpos']=collections.OrderedDict()
                for genome_pos in range(int(row['start']), int(row['end'])+1, 1):
                    if genome_pos not in cds[row['exon_n']]['cpos']:
                        cds[row['exon_n']]['cpos'][genome_pos]=cds_pos
                        cds_pos+=1
                    else:
                        sys.exit('[Error] mutlitple position mapped: %d' % genome_pos)
        final[gene]={}
        final[gene]['acc']=acc
        final[gene]['chr']=chr
        final[gene]['cds']=cds
        final[gene]['strand']=strand
    return final

def CDS_genome_coordinate_todf(df, ref, build='hg38'):
    '''
    convert CDS dict to daraframe
    NOTE: Need to correct refbase at here for alleles
    '''
    accs={}
    reffasta=pysam.FastaFile(ref)
    trimchr=False
    if not all('chr' in r for r in reffasta.references):
        trimchr=True
        
    for gene in df:
        out=[]
        acc=df[gene]['acc']
        chr=df[gene]['chr'].lstrip('chr') if trimchr else df[gene]['chr']
        
        strand=df[gene]['strand']
        for exon_n in df[gene]['cds']:
            start=df[gene]['cds'][exon_n]['start']
            end  =df[gene]['cds'][exon_n]['end']
            for gpos in df[gene]['cds'][exon_n]['cpos']:
                cpos=df[gene]['cds'][exon_n]['cpos'][gpos]
                aapos=(cpos/3)+1 if cpos % 3 != 0 else cpos/3
                aapos=int(aapos)
                pos0based=gpos-1
                
                refbase =reffasta.fetch(chr, pos0based, pos0based+1)

                ### complement base if minus strand
                refbase=complement(refbase) if strand =='-' else refbase

                ### correct the reference for genes
                if gene.upper() == 'RHD' and cpos==1136:
                    refbase='C'
                if gene.upper() == 'RHCE' and cpos==48:
                    refbase='G'
                if chr=='chr9' or chr==9:
                        if build.upper()=='HG38' and gpos==133257522:
                            ### This is the ABO where a G del occured in the Reference HG38
                            out.append({'gene':gene, 'acc':acc, 'exon':exon_n, 'chr': chr, 'start':start, 'end': end, 'strand':strand, 'gpos':gpos, 'cpos':261, 'REF': 'G', 'aapos': 87})
                        elif build.upper()=='HG19' and gpos==136132909:
                            ### This is the ABO where a G del occured in the Reference HG38
                            out.append({'gene':gene, 'acc':acc, 'exon':exon_n, 'chr': chr, 'start':start, 'end': end, 'strand':strand, 'gpos':gpos, 'cpos':261, 'REF': 'G', 'aapos': 87}) 
                        
                out.append({'gene':gene, 'acc':acc, 'exon':exon_n, 'chr': chr, 'start':start, 'end': end, 'strand':strand, 'gpos':gpos, 'cpos':cpos, 'REF': refbase, 'aapos': aapos})
        out=pd.DataFrame(out)
        out[['cpos','exon']]=out[['cpos','exon']].apply(pd.to_numeric)
        if out.loc[out['cpos']==1,'exon'].values[0] > 1:
            exon_num=out.exon.max()
            out.exon=(exon_num+1)-out.exon
        accs[acc]=out
    return accs

def table_create(gene, reference, gbuild):
    """
       Create table in one go
    """
    g=genetable(gene_coordinates(gene, build=gbuild))
    rh, rh_intron=format_table(g, gene)
    cds=CDS_genome_coordinate_todf(CDS_genome_coordinate(rh), reference, build=gbuild)
    return g, cds, rh


