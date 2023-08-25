#!/usr/bin/env python

import re, os
import pandas as pd
import matplotlib, argparse

### fpdf
from pylab import title, figure, xlabel, ylabel, xticks, bar, legend, axis, savefig
from fpdf import FPDF

### html
import plotly
import plotly.plotly as py
import plotly.tools as ptls
import plotly.graph_objs as go
import plotly.figure_factory as ff
import pkg_resources, base64
from plotly import tools
from plotly.graph_objs import Scatter, Layout, Data, Figure

class PDF(FPDF):
    def header(self):
        # Logo
        logo_path="SJ_Tag_V_C.jpg"
        img_fn=pkg_resources.resource_filename(__name__, logo_path)
        self.image(img_fn, 5,4,33)
        
        self.set_font('Arial', 'B', 15)
        # Move to the right
        self.cell(80)
        # Title
        self.cell(20, 20, 'RHtyper', 0, 0, 'C')
        # Line break
        self.ln(20)
        # add header boundary
        self.line(5,self.get_y()+3,200,self.get_y()+3)
        self.ln(3)

    def footer(self):
        # Position at 1.5 cm from bottom
        self.set_y(-15)
        # Arial italic 8
        self.set_font('Arial', 'I', 8)
        # Page number
        self.cell(0, 10, 'Page ' + str(self.page_no()) + '/{nb}', 0, 0, 'C')

class HTML():
    def __init__ (self, sample, gene, BGtable, SNPtable, fullSNPtable, CNVtable, allpos_CNVplot, exon_CNVplot, prefix='test'):
        self.sample=re.sub("\."+gene,"",sample)
        self.gene=gene
        self.prefix=prefix
        self.SNPtable=SNPtable
        self.fullSNPtable=fullSNPtable
    
    def generate(self):
        ### plot title
        title=self.report_title()
        
        ### SNP information
        snp_table, snp_pos=self.var_count_list()

        ### SNP plot
        snp_trace=self.var_plot(1)
        snp_trace2=self.var_plot(2, snp_pos)

         

        data=[snp_table, snp_trace, snp_trace2]
        layout = self.report_layout()
        fig = dict(data=data, layout=layout)
       
 
        #fig=tools.make_subplots(rows=3, cols=1, 
        #                        subplot_titles=('title','SNPs','SNP trace'))

        #fig.append_trace(title,1,1)
        #fig.append_trace(snp_table,2,1)
        #fig.append_trace(snp_trace,3,1)

        #fig['layout'].update(showlegend=False, title='Title')

        plotly.offline.plot(fig, filename = self.prefix+'.html', auto_open=False, config=dict(displaylogo=False, modeBarButtonsToRemove=['sendDataToCloud'])) 

    def report_layout(self):
        logo_path="SJ_Tag_V_C.png"
        image=pkg_resources.resource_filename(__name__, logo_path)
        with open(image, "rb") as image_file:
            encoded_string = base64.b64encode(image_file.read()).decode()
        # Add the prefix that plotly will want when using the string as source
        encoded_image = "data:image/png;base64," + encoded_string
        #image="https://raw.githubusercontent.com/cldougl/plot_images/add_r_img/vox.png"

        axis=dict(
            showline=True,
            zeroline=False,
            showgrid=True,
            mirror=True, 
            ticklen=4, 
            gridcolor='#ffffff',
            tickfont=dict(size=10)
        )


        title  = go.Layout(images=[go.layout.Image(dict(source=encoded_image,xref="paper",yref="paper",x=1, y=1, sizex=0.15, sizey=0.15, xanchor="right",yanchor="bottom"))],
                           autosize=True, 
                           #height=800, width=1200,
                           #margin=dict(r=20, l=300, b=75, t=125),
                           title=go.layout.Title(text='Bloodtyping<br>Blood Group: <i>'+self.gene+'</i><br>Sample: <i>'+self.sample+'</i><br><br>',
                                                 font=dict(family='Courier New, monospace', size=18, color='#7f7f7f'),
                                                 xanchor="left",yanchor="bottom",x=0, xref="paper" 
                                                 ),
                           #yaxis=go.layout.YAxis(automargin=True),
                           showlegend = False,
                           xaxis1=dict(axis, **dict(domain=[0, 1], anchor='y1', showticklabels=False)),
                           xaxis2=dict(axis, **dict(domain=[0, 1], anchor='y2', showticklabels=True)),        
                           xaxis3=dict(axis, **dict(domain=[0, 1], anchor='y3', showticklabels=True)), 
                           yaxis1=dict(axis, **dict(domain=[2 * 0.21 + 0.02 + 0.02, 0.68], anchor='x1')),  
                           yaxis2=dict(axis, **dict(domain=[0.21 + 0.02, 2 * 0.21 + 0.02], anchor='x2',hoverformat='.2f')),
                           yaxis3=dict(axis, **dict(domain=[0.0, 0.21], anchor='x3',hoverformat='.2f')),
                           plot_bgcolor='rgba(228, 222, 249, 0.65)'
                           )
        return title

    def report_title(self):
        trace=go.Table(
            header=dict(values=['Sample', self.sample]),
            cells =dict(values=['Gene', self.gene]))
        sep=dict(type='scatter', x=[1,5], y=[0,0], mode='lines', line=dict(color='red'))
        
        return trace
    
    def var_count_list(self):
        snp=self.SNPtable
        ls_ref=[]
        ls_alt=[]
        ls_ref_aa=[]
        ls_alt_aa=[]
        ls_ref_cnt=[]
        ls_alt_cnt=[]
        for i in range(0, len(snp)):
            ref=snp['REF'].ix[i]
            alt=snp['ALT'].ix[i]
            ls_ref.append(ref)
            ls_ref_aa.append(snp['REFaa'].ix[i])
            ls_alt.append(alt) 
            ls_alt_aa.append(snp['ALTaa'].ix[i])
            ls_ref_cnt.append(snp[ref].ix[i])
            altn=snp[alt].ix[i] if alt in ['A','T','C','G'] else '.'
            ls_alt_cnt.append(altn)

        snp_table=go.Table(
            domain=dict(x=[0, 1], y=[0.7, 1.0]), 
            columnwidth=[2,2,2,2,1,1,2,2,1,1],
            columnorder=[0,1,2,3,4,5,6,7,8,9],
            header=dict(
                height=50,
                values=list(['Chromosome','Genome.position','Coding.position', 'AA.position', 'Ref','Alt','RefCount','AltCount','RefAA','AltAA']),
                line=dict(color='rgb(50, 50, 50)'),
                align = ['left'] * 10,
                font=dict(color=['rgb(45, 45, 45)'] * 10, size=14), 
                fill=dict(color='#d562be')),
            cells =dict(
                height=27,
                values=[snp.chr, snp.gpos, snp.cpos, snp.aapos, ls_ref, ls_alt, ls_ref_cnt, ls_alt_cnt, ls_ref_aa, ls_alt_aa], 
                line = dict(color='#506784'),
                align = ['left'] * 10,
                font = dict(color=['rgb(40, 40, 40)'] * 10, size=12),
                fill = dict(color=['rgb(235, 193, 238)', 'rgba(228, 222, 249, 0.65)'])))
        return snp_table, list(snp.gpos)

    def var_plot(self, plot_object=1, snp_pos=None):
        snp=self.fullSNPtable
        mafs=[]
        info=[]
        for i in range(0, len(snp)):
            ref=snp['REF'].ix[i]
            alt=snp['ALT'].ix[i]
            refaa=snp['REFaa'].ix[i]
            altaa=snp['ALTaa'].ix[i]
            refn=int(snp[ref].ix[i])
            altn=int(snp[alt].ix[i]) if alt in ['A','T','C','G'] else '.'
            if snp_pos is not None:
                 if snp['gpos'].ix[i] in snp_pos: 
                     maf=refn/float(refn+altn) if altn != '.' else 0
                 else:
                     maf=None
            else:
                 maf=refn/float(refn+altn) if altn != '.' else None
            mafs.append(maf)
            if maf is not None and maf >= 0:
                infotext='[' + str(snp['gpos'].ix[i]) + '] [' + str(snp['cpos'].ix[i]) + ' ' + ref +'/' + alt + '] [' + str(snp['aapos'].ix[i]) + ' ' + refaa + '/' + altaa + ']'
            else:
                infotext=None 
            info.append(infotext)
        trace=go.Scatter(
            x=snp.gpos,
            y=mafs,
            text=info, 
            mode="markers",
            xaxis='x'+str(plot_object),
            yaxis='y'+str(plot_object)
        )
        return trace
  

def RHD_RHCE (sample, gene, BGtable, SNPtable, CNVtable, allpos_CNVplot, exon_CNVplot, prefix='test', border='TB', warning_sec=None,warning_msg='Warning'):
    '''
       generate report using FPDF
    '''
    bg=pd.read_csv(BGtable, sep="\t")   if os.path.isfile(BGtable) and os.stat(BGtable).st_size != 0 else pd.DataFrame()
    snp=pd.read_csv(SNPtable, sep="\t") if os.path.isfile(SNPtable) and os.stat(SNPtable).st_size != 0 else pd.DataFrame() 
    cnv=pd.read_csv(CNVtable, sep="\t") if os.path.isfile(CNVtable) and os.stat(CNVtable).st_size != 0 else pd.DataFrame()


    ### plot initiate
    pdf = PDF()
    pdf.set_fill_color(89,89,89)
    pdf.alias_nb_pages()
    #pdf.add_font('DejaVu', '', 'DejaVuSansCondensed.ttf', uni=True)
    pdf.set_font('Arial', 'B', 11)
    pdf.add_page()
    ### cell format
    ### width, height, text, border, ln, align, fill, link
    sample=re.sub('\.'+gene,"", sample)
    pdf.set_text_color(0, 0,0)
    pdf.cell(0, 10, 'Sample: ' + str(sample), 0, 1)
    pdf.cell(0, 10, 'Gene: ' + str(gene), 0, 1)
    
    ### bloodtyping
    #print bg
    pdf.cell(0, 10, 'Genotype:', 0, 1)
    pdf.set_font('Arial', '', 10); pdf.set_text_color(255,255,255)
    pdf.cell(10,7,'Allele','TB',0,'C', fill=True)
    pdf.cell(90,7,'Genotype','TB',0,'C', fill=True)
    pdf.cell(60,7,'Alias','TB',1,'C', fill=True)
    #pdf.cell(50,7,'CoverageProfile',1,0,'C')
    #pdf.cell(40,7,'Break' ,1,1,'C')
    pdf.set_text_color(0, 0,0)
    for i in range(0, len(bg)):
        pdf.cell(10,7, '%s' % (bg['Allele'].iloc[i]),'TB',0,'C')
        pdf.cell(90,7, '%s' % (bg['Bloodtype'].iloc[i].encode('ascii','ignore').decode('ascii')),'TB',0,'C')
        if bg['Alias'].iloc[i] != bg['Bloodtype'].iloc[i] and pd.notna(bg['Alias'].iloc[i]):
            pdf.cell(60,7, '%s' % (bg['Alias'].iloc[i]),'TB',1,'C')
        else:
            pdf.cell(60,7, '%s' % (''),'TB',1,'C')

    if warning_sec is not None and warning_sec=="Genotype":
       pdf.set_font('Arial', '', 7) 
       pdf.cell(0, 3, warning_msg, 0, 1)
    if pdf.get_y() > 300: pdf.add_page()

    ### hybrid info
    pdf.set_font('Arial', 'B', 11)
    pdf.cell(0, 10, 'CoverageProfile:', 0, 1)
    pdf.set_font('Arial', '', 10); pdf.set_text_color(255,255,255)
    pdf.cell(10,7,'Allele','TB',0,'C', fill=True)
    pdf.cell(40,7,'CoverageProfile','TB',0,'C', fill=True)
    pdf.cell(40,7,'Break' ,'TB',1,'C', fill=True)
    pdf.set_text_color(0, 0,0)
    for i in range(0, len(bg)):
        allele=bg['Allele'].iloc[i]
        cp=bg['CoverageProfile'].iloc[i] if pd.notna(bg['CoverageProfile'].iloc[i]) else ''
        bk=bg['Break'].iloc[i] if pd.notna(bg['Break'].iloc[i]) else ''
        pdf.cell(10,7, '%s' % (allele),'TB',0,'C')
        pdf.cell(40,7, '%s' % (cp),'TB',0,'C')
        pdf.cell(40,7, '%s' % (bk),'TB',1,'C')
    #if pdf.get_y() > 300: pdf.add_page()
    pdf.add_page()

    ### variations
    #print snp
    pdf.set_font('Arial', 'B', 11)
    pdf.cell(0, 10, 'Variations:', 0, 1)
    pdf.set_font('Arial', '', 8); pdf.set_text_color(255,255,255)
    pdf.cell(20,7,'Exon','TB',0,'C', fill=True) 
    pdf.cell(20,7,'Cpos','TB',0,'C', fill=True)
    pdf.cell(20,7,'AApos','TB',0,'C', fill=True) 
    pdf.cell(20,7,'Ref' ,'TB',0,'C', fill=True)
    pdf.cell(20,7,'Alt' ,'TB',0,'C', fill=True)
    pdf.cell(20,7,'Ref_aa' ,'TB',0,'C', fill=True)
    pdf.cell(20,7,'Alt_aa' ,'TB',0,'C', fill=True)
    pdf.cell(20,7,'RefN' ,'TB',0,'C', fill=True)
    pdf.cell(20,7,'AltN' ,'TB',1,'C', fill=True)
    pdf.set_text_color(0, 0,0)
    for i in range(0, len(snp)):
        ref=snp['REF'].iloc[i]; refaa=snp['REFaa'].iloc[i]
        alt=snp['ALT'].iloc[i]; altaa=snp['ALTaa'].iloc[i]
        refn=snp[ref].iloc[i]
        altn=snp[alt].iloc[i] if alt in ['A','T','C','G'] else '.'
        pdf.cell(20,7, '%s' % (snp['exon'].iloc[i]),'TB',0,'C') 
        pdf.cell(20,7, '%s' % (int(snp['cpos'].iloc[i])),'TB',0,'C')
        pdf.cell(20,7, '%s' % (int(snp['aapos'].iloc[i])),'TB',0,'C')  
        pdf.cell(20,7, ref,'TB',0,'C')
        pdf.cell(20,7, alt,'TB',0,'C')
        pdf.cell(20,7, refaa,'TB',0,'C')
        pdf.cell(20,7, altaa,'TB',0,'C')
        pdf.cell(20,7, str(refn),'TB',0,'C')
        pdf.cell(20,7, str(altn),'TB',1,'C')
    #if pdf.get_y() > 300: pdf.add_page()
    pdf.set_font('Arial', '', 7)
    pdf.cell(0, 3, 'Exon: exon number', 0, 1)
    pdf.cell(0, 3, 'Cpos: coding position of the variant', 0, 1)
    pdf.cell(0, 3, 'AApos: amino acid position of the variant', 0, 1)
    pdf.cell(0, 3, 'Ref: reference nucleotide', 0, 1)
    pdf.cell(0, 3, 'Alt: alternative nucleotide', 0, 1)
    pdf.cell(0, 3, 'Ref_aa: reference amino acid', 0, 1)
    pdf.cell(0, 3, 'Alt_aa: alternative amino acid', 0, 1)
    pdf.cell(0, 3, 'RefN: Number of reads with reference nucelotide', 0, 1)
    pdf.cell(0, 3, 'AltN: Number of reads with alternative nucelotide', 0, 1)

    ### CNV table
    #print(cnv)
    pdf.set_font('Arial', 'B', 11)
    pdf.cell(0, 10, 'Copy number variation of the exon regions:', 0, 1)  
    pdf.set_font('Arial', '', 10) 
    pdf.cell(0, 10, 'Mean coverage:'+ str(cnv['wgs_cov'].iloc[0]),0,1)
    pdf.set_font('Arial', '', 10); pdf.set_text_color(255,255,255)
    pdf.cell(20,7,'Exon','TB',0,'C', fill=True)
    pdf.cell(40,7,'RHD_log2R','TB',0,'C', fill=True)
    pdf.cell(40,7,'RHD_status','TB',0,'C', fill=True)
    pdf.cell(40,7,'RHCE_log2R' ,'TB',0,'C', fill=True)
    pdf.cell(40,7,'RHCE_status' ,'TB',1,'C', fill=True)
    pdf.set_text_color(0, 0,0)
    for i in range(0, len(cnv)):
        pdf.cell(20,7, '%s' % (cnv['exon'].iloc[i]),border,0,'C')
        pdf.cell(40,7, '%s' % (round(cnv['RHD_log2cov'].ix[i],2)),border,0,'C')
        pdf.cell(40,7, '%s' % (cnv['RHD_status'].iloc[i]),border,0,'C')
        pdf.cell(40,7, '%s' % (round(cnv['RHCE_log2cov'].ix[i],2)),border,0,'C')
        pdf.cell(40,7, '%s' % (cnv['RHCE_status'].iloc[i]),border,1,'C')
    pdf.add_page()

        

    ### all pos CNV
    pdf.set_text_color(0, 0,0)
    pdf.set_font('Arial', 'B', 11)
    pdf.cell(0, 10, 'Normalized coverage of the gene region', 0, 1)
    pdf.set_font('Arial', '', 10)
    pdf.image(allpos_CNVplot, x=pdf.get_x()-2, y=pdf.get_y()+3,w=180,h=60)
    pdf.set_y(pdf.get_y()+60)
    pdf.ln(3) 

    if pdf.get_y() + 60 > 300: pdf.add_page()

    ### exon pos CNV
    pdf.set_text_color(0, 0,0)
    pdf.set_font('Arial', 'B', 11)
    pdf.cell(0, 10, 'Normalized coverage of the exon regions', 0, 1)
    pdf.set_font('Arial', '', 10)
    pdf.image(exon_CNVplot, x=pdf.get_x()-2, y=pdf.get_y()+3,w=190,h=60)
    pdf.set_y(pdf.get_y()+60)
    pdf.set_font('Arial', '', 9)
    pdf.cell(0, 10, 'Note: the points are colored by exons.', 0, 1)
    pdf.ln(3)    


    pdf.output(prefix+'.bloodtyping.pdf', 'F')

def main():
    parser = argparse.ArgumentParser(description='Bloodtyping')
    parser.add_argument('-sample' , '--sample', help='sample', required=True)
    parser.add_argument('-g' , '--gene', help='gene name', required=True)
    parser.add_argument('-bg' , '--typing', help='blood typing table', required=True)
    parser.add_argument('-snp' , '--snp', help='snp table', required=True)
    parser.add_argument('-cnv' , '--cnv', help='cnv table', required=True)
    parser.add_argument('-ac' , '--all_cnv', help='CNV plot', required=True)
    parser.add_argument('-ec' , '--exon_cnv', help='exon CNV plot', required=True)
    
    args=parser.parse_args() 
    RHD_RHCE(args.sample, args.gene, args.typing, args.snp, args.cnv, args.all_cnv, args.exon_cnv) 

if __name__ == "__main__":
    main()

