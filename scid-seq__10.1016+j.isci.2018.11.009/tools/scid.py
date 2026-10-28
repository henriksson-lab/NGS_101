"""scID-seq post-sort protein reporter branch; no bulk RNA/ChIP recipe reuse."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
import illumina as il
from chemdraw import Construct, Scene, Segment, feature, panel, protocol_graph_panels
from protocol_graph import ProtocolGraph, MaterialKind as K, InputRole as R, StepKind as S, SequencingPlan
from pcr import PCRGeometry
from page import head, info
PAGE_FILE='scid.html'
SOURCE='https://doi.org/10.1016/j.isci.2018.11.009; published Transparent Methods antibody/library sections and TableS1'
FORWARD=il.TRUSEQ_READ1[-19:]
REVERSE='TCGCTTATCTGTTGACTGAT'
FW_LONG=il.P5+il.TRUSEQ_READ1[4:]
ACTIN_BARCODE='GTATCGTCGT'
UMI_NT=15

def identity_components():
    sc=Scene()
    sc.strand('antibody',[Segment('Actin barcode1',ACTIN_BARCODE,feature=feature('scid.antibody','feature_barcode','fixed'))],label='source identity fragment')
    sc.strand('umi',[Segment('UMI','N'*UMI_NT,placeholder=True,feature=feature('scid.umi','umi','random'))],label='source UMI width')
    sc.strand('well',[Segment('well identity','[well barcode]',placeholder=True,
        feature=feature('scid.cell.well','cell_barcode','combinatorial',group='scid.cell',part='well'))],label='placement/length unassigned')
    sc.strand('plate',[Segment('plate identity','[plate barcode]',placeholder=True,
        feature=feature('scid.cell.plate','cell_barcode','combinatorial',group='scid.cell',part='plate'))],label='placement/length unassigned')
    sc.note('antibody','Separate identity components; not a complete reporter or final library')
    return sc

def cleanup(g,key,template,label,kind=K.MOLECULAR):
    r=g.input(key+'_column','QIAquick purification and source cleanup; exact recovery not predicted',K.SELECTION)
    out=g.product(key+'_kept',label,kind);w=g.product(key+'_waste','Unretained '+key+' material',K.WASTE)
    g.react(key,key.replace('_',' '),S.CLEANUP,((template,R.TEMPLATE),(r,R.SELECTION)),(out,w))
    g.discard(w,'Removed at column purification/oil separation as applicable');return out

def well_branch(g,key,forward,reverse,enzyme,release):
    # Two symbolic input classes expose the actual plate pooling operation. No
    # experimental well IDs, barcode bases or study-specific cells are invented.
    cell=g.input(key+'_cell','Source-fixed immunostained single cells: symbolic well class '+key,K.BIOLOGICAL,primary=True)
    released=g.product(key+'_released','DTT-released antibody dsDNA in '+key,K.MOLECULAR)
    g.react(key+'_release','Release tags in10mM DTT/15mMTris,pH8.8 under Vapor-lock',S.CHEMICAL_CONVERSION,
        ((cell,R.TEMPLATE),(release,R.CONDITION)),(released,))
    pre=g.product(key+'_preamplified','Un-cell-barcoded preamplified reporters in '+key,K.MOLECULAR)
    geom=PCRGeometry.unavailable(identity_components(),reason='Full antibody dsDNA template/anchor order is unassigned pending cited ID-seq sources; identity components are not a template reconstruction.')
    g.react(key+'_preamp','PCR1:15cycles BEFORE well barcoding',S.PCR,
        ((released,R.TEMPLATE),(forward,R.STRUCTURAL),(reverse,R.STRUCTURAL),(enzyme,R.CATALYST)),(pre,),geometry=geom)
    wellprimer=g.input(key+'_well_primer','Distinct reverse well-barcoding primer for '+key+'; sequence/length unassigned',K.OLIGO)
    tagged=g.product(key+'_barcoded','Well-barcoded antibody reporters in '+key,K.MOLECULAR)
    g.react(key+'_well_pcr','PCR2:5cycles of well-specific barcoding',S.PCR,
        ((pre,R.TEMPLATE),(forward,R.STRUCTURAL),(wellprimer,R.STRUCTURAL),(enzyme,R.CATALYST)),(tagged,),geometry=geom)
    return tagged

def protocol_graph():
    g=ProtocolGraph('scID-seq — post-sort plate barcoding',authority=SOURCE)
    fw=g.input('forward','Source forward amplification primer',K.OLIGO,state=Construct([Segment('forward',FORWARD)]))
    rv=g.input('reverse','Source reverse amplification primer',K.OLIGO,state=Construct([Segment('reverse',REVERSE)]))
    enzyme=g.input('herculase','HerculaseII Fusion polymerase with source buffer/dNTPs',K.CATALYST)
    release=g.input('release_buffer','10mMDTT in15mMTris,pH8.8; source plate release buffer',K.CONDITION)
    a=well_branch(g,'a',fw,rv,enzyme,release);b=well_branch(g,'b',fw,rv,enzyme,release)
    pooled=g.product('pooled','Pooled per-plate well-barcoded reporters',K.MOLECULAR)
    g.react('plate_pool','Pool well-specific products per plate',S.POOL,((a,R.TEMPLATE),(b,R.TEMPLATE)),(pooled,))
    clean=cleanup(g,'plate_cleanup',pooled,'Purified pooled DNA after Vapor-lock removal')
    exo=g.input('exonuclease','ExonucleaseI',K.CATALYST)
    treated=g.product('exo_treated','Pooled DNA after free-primer digestion',K.MOLECULAR)
    g.react('exo','Digest remaining primers with EXO1',S.DIGESTION,((clean,R.TEMPLATE),(exo,R.CATALYST)),(treated,))
    reclean=cleanup(g,'post_exo_cleanup',treated,'Repurified pooled DNA')
    long=g.input('fw_long','Source Fw_long primer',K.OLIGO,state=Construct([Segment('Fw_long',FW_LONG)]))
    plate=g.input('plate_reverse','Plate-specific reverse primer; exact identity/index unassigned',K.OLIGO)
    amplified=g.product('plate_amplified','Plate-barcoded library products',K.MOLECULAR)
    g.react('plate_pcr','PCR3:5cycles with plate-specific reverse primer',S.PCR,
        ((reclean,R.TEMPLATE),(long,R.STRUCTURAL),(plate,R.STRUCTURAL),(enzyme,R.CATALYST)),(amplified,),
        geometry=PCRGeometry.unavailable(identity_components(),reason='Plate-specific reverse oligos and complete intermediate DNA require cited ID-seq source/version review.'))
    lib=cleanup(g,'final_cleanup',amplified,'Purified scID protein-reporting library (source expected~185bp)',K.LIBRARY)
    g.sequence(lib,SequencingPlan(acquisition='unspecified',
        limitation='NextSeq500 and bcl2fastq stated; antibody/UMI/well anchor analysis described. Read lengths, plate index placement and full priming geometry remain unassigned.'))
    g.declare_reviewed_scope(('forward','reverse','herculase','release_buffer',
        'a_cell','a_released','a_release','a_preamplified','a_preamp','a_well_primer','a_barcoded','a_well_pcr',
        'b_cell','b_released','b_release','b_preamplified','b_preamp','b_well_primer','b_barcoded','b_well_pcr',
        'pooled','plate_pool','plate_cleanup_column','plate_cleanup_kept','plate_cleanup_waste','plate_cleanup',
        'exonuclease','exo_treated','exo','post_exo_cleanup_column','post_exo_cleanup_kept','post_exo_cleanup_waste','post_exo_cleanup',
        'fw_long','plate_reverse','plate_amplified','plate_pcr','final_cleanup_column','final_cleanup_kept','final_cleanup_waste','final_cleanup'),SOURCE)
    return g

def render_page():
    return '\n'.join([head('scID-seq — bounded plate model'),'<div class="wrap"><h1>scID-seq</h1>',
        '<p class="research-notes"><a href="01_scid_seq.html">Research notes</a></p>',
        info('Post-sort protein reporter branch. Two symbolic well classes show plate pooling; they do not assign experimental cells or barcode values. Complete reporter/well/plate primer geometry remains unassigned. Bulk RNA/ChIP validation libraries are outside this branch.'),
        protocol_graph_panels(protocol_graph(),cls='long'),panel(identity_components().rows(),caption='Separate source identity components, not an assembled oligo or final library.'),'</div>'])

def selftest():
    from checks import Check
    c=Check();g=protocol_graph();c('Protein-only scID output',len(g.sequencing),1)
    c('Preamplification precedes well barcoding',g.steps['a_well_pcr'].inputs[0].material,'a_preamplified')
    c('Pool occurs only after well barcoding',[p.material for p in g.steps['plate_pool'].inputs],['a_barcoded','b_barcoded'])
    from storage import sources_root
    folder=sources_root()/Path(__file__).resolve().parents[1].name
    p=folder/'j.isci.2018.11.009_supp_mmc1.pdf.txt'
    if p.exists():
        text=''.join(p.read_text().split())
        for seq in (FORWARD,REVERSE,FW_LONG):c('Transparent Methods transcription '+seq,seq in text)
    else:c.skip('Primer transcription','External Transparent Methods unavailable')
    p=folder/'j.isci.2018.11.009_supp_mmc2.xlsx.txt'
    if p.exists():c('TableS1 Actin barcode1 mapping',ACTIN_BARCODE+'\t1\t1\tActin' in p.read_text())
    else:c.skip('Antibody barcode transcription','External TableS1 unavailable')
    c.report(exit_on_fail=True)
