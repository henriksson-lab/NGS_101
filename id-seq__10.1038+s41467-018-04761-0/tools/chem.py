"""Published ID-seq two-stage reporter barcoding; incomplete source-local model."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
import illumina as il
from chemdraw import Construct,Scene,Segment,feature,revcomp,panel,workflow_from_protocol_graph,workflow_panel
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,Terminal,Disposition,SourceItem
from pcr import PCRGeometry
from page import head,info
PAGE_FILE='chem.html'
TITLE='ID-seq — released antibody reporters to plate-index PCR'
SOURCE='Published DOI10.1038/s41467-018-04761-0, PMC6008431 Methods sample barcoding; published SI Note1/Fig1 and SupplementaryData3'
LIMIT='Private partial model: antibody conjugation, modified released-dsDNA ends, complete PCR product strands and sequencing-primer/run geometry remain unassembled. A production template oligo is not silently equated with its manufactured antibody-linked dsDNA product.'
WORKBOOK='s41467-018-04761-0_supp_41467_2018_4761_MOESM6_ESM.xlsx.txt'
FORWARD=il.P5[:19]
ANCHOR='ATCAGTCAACAGATAAGCG'
PRODUCTION_REVERSE='TCGCTTATCTGTTGACTGAT'
WELL_BARCODE='AGAGTCGCGT'
WELL_HANDLE='GATCGGAAGAGCACACGTC'
PLATE_PRIMER_INDEX='ATCAGT'

def production_template():
    return Construct([Segment('source P5 prefix',il.P5[:25]),Segment('Read1 region',il.TRUSEQ_READ1),
        Segment('FLD0001 diversity spacer','A'),Segment('reporter UMI','N'*15,placeholder=True,feature=feature('idseq.umi','umi','random')),
        Segment('FLD0001 antibody barcode','GTATCGTCGT',feature=feature('idseq.antibody','feature_barcode','fixed')),
        Segment('template anchor',ANCHOR)],'Source production template FLD0001; not released dsDNA')

def well_primer():
    return Construct([Segment('second reverse handle',revcomp(WELL_HANDLE)),
        Segment('FLD0151 well barcode reverse strand',revcomp(WELL_BARCODE),feature=feature('idseq.well','sample_index','fixed')),
        Segment('production reverse segment',PRODUCTION_REVERSE)],'Source well reverse primer FLD0151')

def plate_primer():
    return Construct([Segment('P7',il.P7),Segment('source plate25 index in primer orientation',PLATE_PRIMER_INDEX,'idx',feature=feature('idseq.plate','sample_index','fixed')),
        Segment('source truncated Read2 region',il.TRUSEQ_READ2[:-1])],'Source plate reverse primer25')

def unknown_scene():
    sc=Scene();sc.strand('unknown',[Segment('released dsDNA ends and complete product','[unassembled molecular geometry]',placeholder=True)])
    sc.note('unknown','Source-scope assembly boundary; not a claim of missing disclosure.')
    return sc

def cleanup(g,key,x,label):
    sel=g.input(key+'_selection',label,K.SELECTION)
    kept=g.product(key+'_kept','Retained DNA after '+label,K.MOLECULAR)
    waste=g.product(key+'_waste','Physically removed '+key+' fractions',K.WASTE)
    g.react(key,label,S.CLEANUP,((x,R.TEMPLATE),(sel,R.SELECTION)),(kept,waste));g.discard(waste,'Removed only at this source-described separation')
    return kept

def protocol_graph():
    g=ProtocolGraph(TITLE,authority=SOURCE)
    fw=g.input('forward','Source19nt forward primer',K.OLIGO,state=Construct([Segment('forward',FORWARD)]))
    enzyme=g.input('pfu','PFU polymerase, source buffer and dNTPs',K.CATALYST)
    branches=[]
    for key in ('a','b'):
        x=g.input(key+'_released','Released antibody dsDNA from symbolic cell-population well '+key,K.MOLECULAR,primary=True)
        rv=g.input(key+'_well_primer','Distinct well-specific reverse primer; FLD0151 inventory is an example, not assigned to this symbolic well',K.OLIGO)
        spike=g.input(key+'_spike','Source spike-in DNA barcodes; identities not assigned',K.MOLECULAR)
        y=g.product(key+'_barcoded','Well-barcoded reporters and source spike-ins in '+key,K.MOLECULAR)
        g.react(key+'_well_pcr','PCR1 directly adds well barcode; no scID15cycle preamplification',S.PCR,
            ((x,R.TEMPLATE),(spike,R.SPIKE_IN),(fw,R.STRUCTURAL),(rv,R.STRUCTURAL),(enzyme,R.CATALYST)),(y,),
            geometry=PCRGeometry.unavailable(unknown_scene(),reason='Full manufactured reporter and well-specific product strands unassembled; printed primer inventory is separate.'))
        branches.append(y)
    pool=g.product('pool','Pooled well-barcoded products from one plate',K.MOLECULAR)
    g.react('plate_pool','Pool well-specific PCR products per plate',S.POOL,tuple((x,R.TEMPLATE) for x in branches),(pool,))
    clean=cleanup(g,'first_column',pool,'Qiagen PCR purification column')
    exo=g.input('exoi','ExonucleaseI',K.CATALYST)
    dig=g.product('digested','Pooled DNA with primer digestion products still present',K.MOLECULAR)
    g.react('primer_digest','ExonucleaseI37C30min',S.DIGESTION,((clean,R.TEMPLATE),(exo,R.CATALYST)),(dig,))
    heat=g.input('heat','80C20min inactivation',K.CONDITION)
    inactive=g.product('inactive','Heat-inactivated reaction; digestion products not physically discarded',K.MOLECULAR)
    g.react('exo_inactivation','Heat-inactivate ExonucleaseI',S.THERMAL_TREATMENT,((dig,R.TEMPLATE),(heat,R.CONDITION)),(inactive,))
    rv=g.input('plate_reverse','Example source plate25 reverse primer',K.OLIGO,state=plate_primer())
    indexed=g.product('indexed','Plate-index PCR reaction; complete duplex unassembled',K.MOLECULAR)
    g.react('plate_pcr','PCR2 adds plate index and second adapter',S.PCR,((inactive,R.TEMPLATE),(fw,R.STRUCTURAL),(rv,R.STRUCTURAL),(enzyme,R.CATALYST)),(indexed,),
        geometry=PCRGeometry.unavailable(unknown_scene(),reason='Source plate25 primer known; full pooled product/end assembly and exact selected experimental well remain unassigned.'))
    second=cleanup(g,'second_column',indexed,'Qiagen PCR purification after PCR2')
    finaldig=g.product('final_digest','Second ExonucleaseI reaction including digestion products',K.MOLECULAR)
    g.react('final_primer_digest','ExonucleaseI45min',S.DIGESTION,((second,R.TEMPLATE),(exo,R.CATALYST)),(finaldig,))
    selected=cleanup(g,'pkis_size_selection',finaldig,'PKIS branch only: Zymo selection for fragments greater than150bp')
    out=cleanup(g,'minelute',selected,'Final Qiagen MinElute purification')
    g.mark_terminal(Terminal(out.id,Disposition.UNRESOLVED,LIMIT))
    g.declare_reviewed_scope(tuple(g.materials)+tuple(g.steps),SOURCE)
    g.add_source_item(SourceItem('unassembled_source','Remaining source and full reporter/read assembly',LIMIT))
    return g

def render_page():
    g=protocol_graph()
    def inventory(x):
        sc=Scene();sc.strand('oligo',x.segments,label=x.name);return sc.rows()
    return '\n'.join([head(TITLE),'<main class="wrap"><h1>'+TITLE+'</h1>',info(LIMIT),'<p><a href="01_id_seq.html">Research notes</a></p>',
        workflow_panel(workflow_from_protocol_graph(g,'minelute_kept')),panel(inventory(production_template()),caption='Verbatim production-template component example FLD0001; not manufactured dsDNA geometry'),
        panel(inventory(well_primer()),caption='Source well-primer FLD0151 inventory; source-field orientation preserved'),'</main>'])

def selftest():
    from checks import Check
    from storage import sources_root
    c=Check();g=protocol_graph()
    c('Private incomplete model',g.publishable,False)
    c('No false completion marker','data-protocol-complete=' in render_page(),False)
    c('Well barcoding directly uses released reporters',g.steps['a_well_pcr'].inputs[0].material,'a_released')
    c('Perplate pooling follows both independent wellPCRs',[i.material for i in g.steps['plate_pool'].inputs],['a_barcoded','b_barcoded'])
    c('Digest debris remains through heat inactivation',g.steps['exo_inactivation'].inputs[0].material,'digested')
    p=sources_root()/Path(__file__).resolve().parents[1].name/WORKBOOK
    if p.exists():
        source=''.join(p.read_text().split())
        for x in (production_template(),well_primer(),plate_primer()):c('Published Data3 transcription '+x.name,x.top() in source,True)
        c('Published forward transcription',FORWARD in source,True)
    else:c.skip('Published Data3 transcriptions','External original unavailable')
    c.report()
