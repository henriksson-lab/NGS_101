"""scID-seq post-sort protein reporter branch; no bulk RNA/ChIP recipe reuse."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
import illumina as il
from chemdraw import Construct, Scene, Segment, feature, panel, protocol_graph_panels
from protocol_graph import ProtocolGraph, MaterialKind as K, InputRole as R, StepKind as S, SequencingPlan, Terminal, Disposition
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

def well_branch(g,key,forward,reverse,enzyme,nts,buffer,release):
    # Two symbolic input classes expose the actual plate pooling operation. No
    # experimental well IDs, barcode bases or study-specific cells are invented.
    cell=g.input(key+'_cell','Source-fixed immunostained single cells: symbolic well class '+key,K.BIOLOGICAL,primary=True)
    release_conditions=g.input(key+'_release_conditions','Source7uLVaporlock over releasebuffer; molecular sideproducts unassembled',K.CONDITION)
    released=g.product(key+'_released','DTT-released antibody DNA with copresent cellular/protein material unassembled in '+key,K.MOLECULAR)
    g.react(key+'_release','Release tags in10mM DTT/15mMTris,pH8.8 under Vapor-lock',S.CHEMICAL_CONVERSION,
        ((cell,R.TEMPLATE),(release,R.STRUCTURAL),(release_conditions,R.CONDITION)),(released,))
    pre=g.product(key+'_preamplified','Un-cell-barcoded preamplified reporters in '+key,K.MOLECULAR)
    geom=PCRGeometry.unavailable(identity_components(),reason='Full antibody dsDNA template/anchor order is unassigned pending cited ID-seq sources; identity components are not a template reconstruction.')
    g.react(key+'_preamp','PCR1:15cycles BEFORE well barcoding',S.PCR,
        ((released,R.TEMPLATE),(forward,R.STRUCTURAL),(reverse,R.STRUCTURAL),(enzyme,R.CATALYST),(nts,R.STRUCTURAL),(buffer,R.CONDITION)),(pre,),geometry=geom)
    wellprimer=g.input(key+'_well_primer','Distinct reverse well-barcoding primer for '+key+'; sequence/length unassigned',K.OLIGO)
    tagged=g.product(key+'_barcoded','Well-barcoded antibody reporters in '+key,K.MOLECULAR)
    well_mix=g.input(key+'_well_mix','ActualsourcePCR2addition:1uL Herculasebuffer/forward/reversewellprimer mix; enzyme/dNTP carryover fromPCR1',K.CONDITION)
    g.react(key+'_well_pcr','PCR2:5cycles of well-specific barcoding',S.PCR,
        ((pre,R.TEMPLATE),(forward,R.STRUCTURAL),(wellprimer,R.STRUCTURAL),(enzyme,R.CATALYST),(nts,R.STRUCTURAL),(well_mix,R.CONDITION)),(tagged,),geometry=geom,reason='Forward/wellprimer ports are componentviews of the actual1uLfreshPCR2buffer/primeraddition; enzyme/dNTPs refer to carriedPCR1inventory, not fresh extra reagent doses.')
    return tagged

def protocol_graph():
    g=ProtocolGraph('scID-seq — post-sort plate barcoding',authority=SOURCE,representation='illustrative_example',example_description='Two invented symbolic well classes illustrate the published per-well procedure and per-plate pooling. They are not experimental cells, barcode values or a reconstruction of occupied wells in a96well plate.')
    fw=g.input('forward','Source forward amplification primer',K.OLIGO,state=Construct([Segment('forward',FORWARD)]))
    rv=g.input('reverse','Source reverse amplification primer',K.OLIGO,state=Construct([Segment('reverse',REVERSE)]))
    enzyme=g.input('herculase','Actual HerculaseII Fusion polymerase supply for sourcePCR1; reused catalytic view inPCR2, not a fresh addition',K.CATALYST)
    nts=g.input('dntp','Actual separately supplied dNTP inventory in4uLperwellPCR1mix; PCR2carryover',K.MOLECULAR)
    buffer=g.input('pcr1_buffer','Source5xHerculasebuffer component of4uLperwellPCR1mix',K.CONDITION)
    release=g.input('release_buffer','Actual source1uL perwell releasebuffer:10mMDTT in15mMTris,pH8.8',K.MOLECULAR)
    a=well_branch(g,'a',fw,rv,enzyme,nts,buffer,release);b=well_branch(g,'b',fw,rv,enzyme,nts,buffer,release)
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
    final_enzyme=g.input('plate_enzyme','FreshPCR3polymerase supply; identity/formulation unassigned in selected20uL source recipe',K.CATALYST)
    final_nt=g.input('plate_nt','RequiredPCR3nucleotide supply; packaging/dose unassigned',K.MOLECULAR)
    final_buffer=g.input('plate_buffer','PCR3buffer supply; formulation/packaging unassigned',K.CONDITION)
    g.react('plate_pcr','PCR3:5cycles with plate-specific reverse primer',S.PCR,
        ((reclean,R.TEMPLATE),(long,R.STRUCTURAL),(plate,R.STRUCTURAL),(final_enzyme,R.CATALYST),(final_nt,R.STRUCTURAL),(final_buffer,R.CONDITION)),(amplified,),
        geometry=PCRGeometry.unavailable(identity_components(),reason='Plate-specific reverse oligos and complete intermediate DNA require cited ID-seq source/version review.'))
    lib=cleanup(g,'final_cleanup',amplified,'Purified scID protein-reporting library (source expected~185bp)',K.LIBRARY)
    qc=g.product('final_qc','Analytical material for sourceagarose/BioAnalyzer; amounts unassigned',K.QC)
    retained=g.product('retained_library','Library retained after analytic selection',K.LIBRARY)
    g.react('select_qc','Reserve source library-characterization material',S.SPLIT,((lib,R.TEMPLATE),),(qc,retained),changes_molecule=False)
    g.mark_terminal(Terminal(qc.id,Disposition.QC,'Expected185bp checked by sourcegel/BioAnalyzer; assaysetup/aliquot amounts outside this preparation model.'))
    loading=g.product('loading_aliquot','Selected sequencing loading aliquot; volume unassigned',K.LIBRARY)
    remainder=g.product('library_remainder','Untransferred library stock; amount unassigned',K.LIBRARY)
    g.react('select_loading','Select library contribution for sequencing',S.SPLIT,((retained,R.TEMPLATE),),(loading,remainder),changes_molecule=False)
    g.mark_terminal(Terminal(remainder.id,Disposition.UNRESOLVED,'No whole-stock loading or returned QC material is assumed.'))
    g.sequence(loading,SequencingPlan(acquisition='unspecified',
        limitation='NextSeq500 and bcl2fastq stated; antibody/UMI/well anchor analysis described. Read lengths, plate index placement and full priming geometry remain unassigned.'))
    g.declare_reviewed_scope(('dntp','pcr1_buffer','a_release_conditions','b_release_conditions','a_well_mix','b_well_mix','plate_enzyme','plate_nt','plate_buffer','select_qc','final_qc','select_loading','loading_aliquot','library_remainder','forward','reverse','herculase','release_buffer',
        'a_cell','a_released','a_release','a_preamplified','a_preamp','a_well_primer','a_barcoded','a_well_pcr',
        'b_cell','b_released','b_release','b_preamplified','b_preamp','b_well_primer','b_barcoded','b_well_pcr',
        'pooled','plate_pool','plate_cleanup_column','plate_cleanup_kept','plate_cleanup_waste','plate_cleanup',
        'exonuclease','exo_treated','exo','post_exo_cleanup_column','post_exo_cleanup_kept','post_exo_cleanup_waste','post_exo_cleanup',
        'fw_long','plate_reverse','plate_amplified','plate_pcr','final_cleanup_column','final_cleanup_kept','final_cleanup_waste','final_cleanup'),SOURCE)
    return g.require_publishable()

def render_page():
    from protocol_views import assemble_page,WorkflowSection,OligoSection,OligoEntry
    return assemble_page(protocol_graph(),(head('scID-seq — single-cell protein reporters'),'<div class="wrap" data-protocol-complete="true"><h1>scID-seq</h1>',
        '<p class="research-notes"><a href="01_scid_seq.html">Research notes</a></p>',info('<a href="https://doi.org/10.1016/j.isci.2018.11.009">Published scID-seq article, TransparentMethods and TableS1</a>'),
        '<h2>Three amplification stages preserve cell and protein identity</h2>',
        '<p>The selected branch starts with fixed, immunostained cells already sorted into wells. DTT releases the antibody reporters. Fifteen amplification cycles precede five cycles that add a well barcode. Products are pooled per plate, oil and free primers are removed, then five further cycles add plate identity.</p>',
        info('ILLUSTRATIVE EXAMPLE — two invented symbolic well classes show this operation. The source sorts into96-well plates; this view does not claim two experimental cells,96 occupied wells, specific well/plate sequences or an actual study plate. BulkRNA/ChIP validation libraries and upstream staining/conjugation are outside this post-sort graph.'),
        '<h2>What the known sequences establish</h2>',
        OligoSection('Source amplification primer orders',(
            OligoEntry('Forward',(Segment('Forward',FORWARD),)),OligoEntry('Reverse',(Segment('Reverse',REVERSE),)),OligoEntry('Fw_long',(Segment('Fw_long',FW_LONG),)))),
        panel(identity_components().rows(),caption='Separate source identity components: Actinbarcode1,15-baseUMI and well/plate identity. Their relative positions, complete reporter polymers and final library are unassembled; no equal placeholder strings prove binding.'),
        '<h2>Stock and read boundaries</h2>',
        '<p>Each source well has1 µL releasebuffer under7 µL Vapor-lock. PCR1 adds4 µL enzyme/dNTP/buffer/primer mix. PCR2 adds1 µL buffer/primer mix and uses carried enzyme/nucleotides; it does not receive an invented fresh fullPCR mix. All illustrated well products enter the plate pool before oil removal and QIAquick/ExoI cleanup. FreshPCR3reaction supplies remain partly unassigned in the inspected20 µL recipe.</p>',
        '<p>The source checks an expected185 bp library by agarose/BioAnalyzer and sequences on NextSeq500. Read lengths, loaded sequencing primers, exact well/plate spans and sequencing-index placement remain unassigned. Final analytical material and a selected loading aliquot have explicit retained/untransferred counterparts.</p>',
        WorkflowSection('Full illustrated post-sort preparation and molecular limits'),'</div>'))

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
    c('symbolic wells are explicitly illustrative',g.representation,'illustrative_example')
    c('source DTT is molecular, not only a condition',any(p.material=='release_buffer' and p.role is R.STRUCTURAL for p in g.steps['a_release'].inputs),True)
    c('PCR2uses PCR1catalyst view, not fresh enzyme stock',any(p.material=='herculase' for p in g.steps['a_well_pcr'].inputs),True)
    c('post-purification PCR3 does not presume enzyme carryover',any(p.material=='herculase' for p in g.steps['plate_pcr'].inputs),False)
    c('final assay material terminates analytically',g.terminals['final_qc'].disposition,Disposition.QC)
    c('untransferred stock survives sequencing',g.terminals['library_remainder'].disposition,Disposition.UNRESOLVED)
    c.report(exit_on_fail=True)
