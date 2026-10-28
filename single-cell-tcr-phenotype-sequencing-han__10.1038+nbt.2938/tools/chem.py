"""Bounded, source-local molecular workflow; not a complete protocol review."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Construct, Scene, Segment, feature, panel, workflow_from_protocol_graph, workflow_panel
from protocol_graph import ProtocolGraph, MaterialKind as K, InputRole as R, StepKind as S, SequencingPlan, Terminal, Disposition, SourceItem
from pcr import PCRGeometry
from page import head, info
PAGE_FILE='chem.html'

def unknown_template(label):
    return Construct([Segment(label,'X'*30,placeholder=True)])

def clean(g,key,material,label,selection,kind=K.MOLECULAR):
    sel=g.input(key+'_selection',selection,K.SELECTION)
    kept=g.product(key+'_kept',label,kind)
    waste=g.product(key+'_waste','Separated material from '+key,K.WASTE)
    g.react(key,label,S.CLEANUP,((material,R.TEMPLATE),(sel,R.SELECTION)),(kept,waste))
    g.discard(waste,'Removed only at this source-described physical separation; recovery is not predicted')
    return kept

def pcr(g,key,template,forward,reverse,label,conditions):
    f=g.input(key+'_forward',forward,K.OLIGO)
    r=g.input(key+'_reverse',reverse,K.OLIGO)
    cat=g.input(key+'_enzyme',conditions,K.CATALYST)
    out=g.product(key+'_product',label,K.MOLECULAR)
    g.react(key,label,S.PCR,((template,R.TEMPLATE),(f,R.STRUCTURAL),(r,R.STRUCTURAL),(cat,R.CATALYST)),(out,),
        geometry=PCRGeometry.unavailable(unknown_template('unassigned target sequence'),
            reason='This capped model does not assign the study-specific target substrate and complete primer-pool landing geometry.'))
    return out

def aliquot(g,key,material,label):
    selected=g.product(key+'_selected',label,K.MOLECULAR)
    remainder=g.product(key+'_remainder','Untransferred '+key+' material',K.MOLECULAR)
    g.react(key,label,S.SPLIT,((material,R.TEMPLATE),),(selected,remainder))
    g.mark_terminal(Terminal(remainder.id,Disposition.UNRESOLVED,'The source specifies an aliquot; remaining reaction material is not claimed discarded.'))
    return selected
SOURCE='https://doi.org/10.1038/nbt.2938; Online Methods, PCR primer design, Supplementary Figure1'
TITLE='Single-cell TCR/phenotype sequencing (Han)'
LIMIT='Cell identity links separate TCR and phenotype molecules. Full primer-pool/target geometry and Illumina endpoint bases remain unassigned; no molecular chimera or UMI is implied.'
HANDLE5='CCAGGGTTTTCCCAGTCACGAC'
HANDLE3='AGCGGATAACAATTTCACACAGGA'
def barcode_components():
    # Components only: this is not a contiguous final sequencing-library strand.
    return [Segment(name,'N'*5,'cbc',placeholder=True,
        feature=feature('han.cell.'+part,'cell_barcode','combinatorial',group='han.cell',part=part))
        for name,part in (('plate barcode','plate'),('row barcode','row'),('column barcode','column'))]

def protocol_graph():
    g=ProtocolGraph(TITLE,authority=SOURCE)
    cell=g.input('cell','Sorted single T cell in RT-PCR buffer',K.BIOLOGICAL,primary=True)
    f=g.input('rt_forward','Multiplex V-region and phenotype forward primers',K.OLIGO)
    r=g.input('rt_reverse','Multiplex C-region and phenotype reverse primers',K.OLIGO)
    rt=g.input('rt_enzyme','Qiagen OneStep RT-PCR mix;50°C30min',K.CATALYST)
    cdna=g.product('cdna','First-strand TCR and phenotype cDNA in original well',K.MOLECULAR)
    g.react('rt','Reverse transcribe sorted-cell transcripts',S.REVERSE_TRANSCRIPTION,
        ((cell,R.TEMPLATE),(r,R.STRUCTURAL),(rt,R.CATALYST)),(cdna,))
    first=g.product('first','Multiplex first-PCR products in original well',K.MOLECULAR)
    g.react('first_pcr','First multiplex amplification:25cycles',S.PCR,
        ((cdna,R.TEMPLATE),(f,R.STRUCTURAL),(r,R.STRUCTURAL),(rt,R.CATALYST)),(first,),
        geometry=PCRGeometry.unavailable(unknown_template('cell transcript targets'),reason='Multiplex target/primer-pool geometry is unassigned.'))
    tcr=g.product('tcr_aliquot','1µl first-PCR aliquot for TCR nestedPCR',K.MOLECULAR)
    pheno=g.product('phenotype_aliquot','Separate1µl first-PCR aliquot for phenotype nestedPCR',K.MOLECULAR)
    remainder=g.product('first_remainder','Untransferred first-PCR material',K.MOLECULAR)
    g.react('branch','Transfer separate aliquots to TCR and phenotype reactions',S.SPLIT,((first,R.TEMPLATE),),(tcr,pheno,remainder))
    g.mark_terminal(Terminal(remainder.id,Disposition.UNRESOLVED,'First reaction is20µl; source transfers1µl per branch, without specifying remainder fate.'))
    t2=pcr(g,'tcr_nested',tcr,'Nested V primers with common5prime handle','Nested Cα/Cβ primers','Nested TCR products','HotStarTaq;25cycles')
    p2=pcr(g,'phenotype_nested',pheno,'Phenotype5prime primers with common handle','Phenotype3prime primers with common handle','Nested phenotype products','HotStarTaq;35cycles')
    t3=pcr(g,'tcr_barcoding',aliquot(g,'tcr_transfer',t2,'1µl nestedTCR transfer'),
        'Plate5+row5 barcode primer pool','Column5+C-region primer pool','Cell-identified TCR products','HotStarTaq;36cycles; combinedα/β tumor/colon branch')
    p3=pcr(g,'phenotype_barcoding',aliquot(g,'phenotype_transfer',p2,'1µl nestedphenotype transfer'),
        'Same plate5+row5 barcode primer pool','Column5+phenotype common-handle primer pool','Cell-identified phenotype products','HotStarTaq;36cycles')
    pool=g.product('pooled','Equal-volume combined TCR and phenotype PCR products',K.MOLECULAR)
    g.react('pool','Combine separate branch products',S.POOL,((t3,R.TEMPLATE),(p3,R.TEMPLATE)),(pool,))
    lib=clean(g,'gel',pool,'Pooled TCR/phenotype MiSeq library','1.2% agarose gel; excise350–380bp; QIAquick gel extraction',K.LIBRARY)
    g.sequence(lib,SequencingPlan(acquired_roles=('Read 1','Read 2'),limitation='Source specifies MiSeq paired-end. Exact run-cycle allocation and final primer coordinates are deferred.'))
    g.declare_reviewed_scope(REVIEWED,SOURCE)
    g.add_source_item(SourceItem('unassembled_source','Full source primer-pool, intermediate and final-strand assembly',LIMIT+' Source material is available for follow-up; not an absence finding.'))
    return g

def detail():
    sc=Scene()
    for i,s in enumerate(barcode_components()):sc.strand(s.name,[s],label=s.name+' (identity component, not a complete strand)')
    return panel(sc.rows(),caption='Three source-defined components jointly identify a cell')

REVIEWED=('cell', 'rt_forward', 'rt_reverse', 'rt_enzyme', 'cdna', 'first', 'tcr_aliquot', 'phenotype_aliquot', 'first_remainder', 'tcr_nested_forward', 'tcr_nested_reverse', 'tcr_nested_enzyme', 'tcr_nested_product', 'phenotype_nested_forward', 'phenotype_nested_reverse', 'phenotype_nested_enzyme', 'phenotype_nested_product', 'tcr_transfer_selected', 'tcr_transfer_remainder', 'tcr_barcoding_forward', 'tcr_barcoding_reverse', 'tcr_barcoding_enzyme', 'tcr_barcoding_product', 'phenotype_transfer_selected', 'phenotype_transfer_remainder', 'phenotype_barcoding_forward', 'phenotype_barcoding_reverse', 'phenotype_barcoding_enzyme', 'phenotype_barcoding_product', 'pooled', 'gel_selection', 'gel_kept', 'gel_waste', 'rt', 'first_pcr', 'branch', 'tcr_nested', 'phenotype_nested', 'tcr_transfer', 'tcr_barcoding', 'phenotype_transfer', 'phenotype_barcoding', 'pool', 'gel')
def extra_checks(c,g):
    from storage import sources_root
    source=sources_root()/Path(__file__).resolve().parents[1].name
    original=source/SOURCE_FILE
    if original.exists():
        text=''.join(original.read_text().split())
        for oligo in TRANSCRIBED:
            c('Verified original transcription '+oligo,oligo in text)
    else:c.skip('Original transcription checks','External original unavailable')
    source_checks(c,g)


SOURCE_FILE='nbt.2938_41587_2014_BFnbt2938_MOESM24_ESM.pdf.txt'
TRANSCRIBED=(HANDLE5,HANDLE3)
def source_checks(c,g):
    c('Three barcode components form one cell identity',tuple(x.feature.part for x in barcode_components()),('plate','row','column'))
    c('Phenotype and TCR aliquots remain separate',tuple(g.steps['branch'].outputs),('tcr_aliquot','phenotype_aliquot','first_remainder'))

def render_page():
    return '\n'.join([head(TITLE),'<div class="wrap"><h1>'+TITLE+'</h1>',
        '<p class="research-notes"><a href="01_han_tcr.html">Research notes</a></p>',
        info(LIMIT),info('Source wording discrepancy: main Methods calls the common5′ handle23bases, while its printed literal and SupplementaryFigure1 give '+str(len(HANDLE5))+'bases. This model retains the literal '+HANDLE5+' and does not infer an omitted base. The combined TCR primer concentration ratio is also left unassigned because the main text gives conflicting prose and numeric values.'),detail(),''.join(workflow_panel(workflow_from_protocol_graph(protocol_graph(),endpoint)) for endpoint in ['gel_kept', 'phenotype_barcoding_product']),'</div>'])

def selftest():
    from checks import Check
    c=Check(); g=protocol_graph()
    issues=g.publication_problems()
    c('Only explicit unassembled-source coverage blocks publication',issues,["authoritative source item 'unassembled_source' is not accounted for"])
    c('Private preview does not claim protocol completion','data-protocol-complete="true"' in render_page(),False)
    extra_checks(c,g)
    c.report(exit_on_fail=True)
