"""Source-local bounded partial model; no inherited kit recipe or finished library claim."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from dataclasses import dataclass,replace
from chemdraw import Construct,Scene,Segment,feature,revcomp,panel,workflow_from_protocol_graph,workflow_panel
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,Terminal,Disposition,SourceItem
from pcr import PCRGeometry
from page import head,info
PAGE_FILE='chem.html'

def unknown(label):return Construct([Segment(label,'[unassembled]',placeholder=True)])
def cleanup(g,key,x,label,procedure,state=None):
    sel=g.input(key+'_selection',procedure,K.SELECTION)
    kept=g.product(key+'_kept',label,K.MOLECULAR,state=state)
    waste=g.product(key+'_waste','Physically separated '+key+' rejects',K.WASTE)
    g.react(key,label,S.CLEANUP,((x,R.TEMPLATE),(sel,R.SELECTION)),(kept,waste))
    g.discard(waste,'Source-described physical separation; purity and quantitative recovery are not assumed')
    return kept

def aliquot(g,key,x,label):
    chosen=g.product(key+'_chosen',label,K.MOLECULAR)
    rest=g.product(key+'_remainder','Untransferred '+key+' material',K.MOLECULAR)
    g.react(key,label,S.SPLIT,((x,R.TEMPLATE),),(chosen,rest))
    g.mark_terminal(Terminal(rest.id,Disposition.UNRESOLVED,'Source aliquot remainder fate is not specified in the inspected section; no disposal inferred'))
    return chosen

def limit(g,x):
    g.mark_terminal(Terminal(x.id,Disposition.UNRESOLVED,LIMIT))
    g.declare_reviewed_scope(REVIEWED,SOURCE)
    g.add_source_item(SourceItem('unassembled_source','Remaining source library/input chemistry and full strand/read assembly',LIMIT))
    return g

def render_page():
    g=protocol_graph()
    return '\n'.join([head(TITLE),'<main class="wrap"><h1>'+TITLE+'</h1>',
        '<p class="research-notes"><a href="'+NOTE+'">Research notes</a></p>',
        info('Private incomplete source model. '+LIMIT),detail(),
        ''.join(workflow_panel(workflow_from_protocol_graph(g,x)) for x in ENDPOINTS),'</main>'])

def selftest():
    from checks import Check
    from storage import sources_root
    c=Check();g=protocol_graph()
    c('Unassembled source and endpoint remain unpublished',g.publishable,False)
    c('Private preview has no completion marker','data-protocol-complete=' in render_page(),False)
    p=sources_root()/Path(__file__).resolve().parents[1].name/SOURCE_FILE
    if p.exists():
        source=''.join(p.read_text().split()).upper()
        for seq in TRANSCRIBED:c('Original transcription '+seq,seq.upper() in source,True)
    else:c.skip('Source transcription','External original unavailable')
    source_checks(c,g)
    c.report()
import illumina as il
TITLE='MULTI-seq — published supplemental protocol variant'
NOTE='01_multi_seq.html'
SOURCE='Published SI DOI10.1038/s41592-019-0433-8 MOESM1 PDFp33–35 protocol; associated accepted manuscript PMC6837808 used only as separately named comparison'
LIMIT='Published SI protocol route only. Historical10xv2 capture/RT and full CBC/UMI cDNA/primer landing remain unassembled. Accepted-manuscript oligo variant is separate; final VOR equivalence and exact study-specific version remain unresolved. No anchor/primer recipe is merged across versions.'
SOURCE_FILE='s41592-019-0433-8_41592_2019_433_MOESM1_ESM.pdf.txt'
CO='AGTGACAGCTGGATCGTTAC'
BARCODE_ARM=il.TRUSEQ_SMALL_RNA_RA3
HANDLE=revcomp(BARCODE_ARM)
ANCHOR_SI=BARCODE_ARM+revcomp(CO)
ANCHOR_MANUSCRIPT=revcomp(CO)+BARCODE_ARM
ADDITIVE_SI='CTTGGCACCCGAGAATTCC'
ADDITIVE_MANUSCRIPT='CCTTGGCACCCGAGAATTCC'
P5=il.NEBNEXT_UNIVERSAL_PRIMER
RPI_TAIL='GTGACTGGAGTTCCTTGGCACCCGAGAATTCCA'
TRANSCRIBED=(ANCHOR_SI,CO,HANDLE+'N'*8+'A30',ADDITIVE_SI,il.P7+'N'*6+RPI_TAIL,P5)
ENDPOINTS=('library_cleanup_kept','assembly')
@dataclass(frozen=True)
class MembraneComplex:
    version:str='published_si'
    def __post_init__(self):
        if self.version not in {'published_si','accepted_manuscript'}:raise ValueError('Version identity must be explicit')
    @property
    def anchor(self):return ANCHOR_SI if self.version=='published_si' else ANCHOR_MANUSCRIPT
    @property
    def additive(self):return ADDITIVE_SI if self.version=='published_si' else ADDITIVE_MANUSCRIPT
    def scene(self):
        sc=Scene()
        arms=[Segment('barcode-binding arm',BARCODE_ARM),Segment('co-anchor-binding arm',revcomp(CO))]
        if self.version=='accepted_manuscript':arms.reverse()
        sc.strand('anchor',arms,mod5='lipid (version-specific source)',label=self.version+' anchor; noncovalent membrane association')
        sc.anneal('co',[Segment('co-anchor',CO)],to='anchor',pair=('co-anchor','co-anchor-binding arm'),mod3='lipid (version-specific source)')
        sc.anneal('barcode',[Segment('sample-tag handle',HANDLE),Segment('sample tag','N'*8,'idx',placeholder=True,
            feature=feature('multi.biological_sample','sample_index','unknown',note='Membrane sample tag; distinct from Illumina library index')),
            Segment('polyA capture tail','A'*30)],to='anchor',pair=('sample-tag handle','barcode-binding arm'),unpaired=('sample tag','polyA capture tail'))
        sc.note('barcode','Sample tag and polyA are unpaired overhangs in this partial binding view, not a covalent continuation of anchor/co-anchor')
        return sc

def rpi_primer():return Construct([Segment('P7',il.P7),Segment('sequencer library index','N'*6,'idx',placeholder=True,feature=feature('multi.library_index','sample_index','unknown')),
    Segment('RPIX constant tail',RPI_TAIL)])

def protocol_graph():
    g=ProtocolGraph(TITLE,authority=SOURCE)
    cells=g.input('cells','Source prepared single-cell suspension; independent sample label per pool',K.BIOLOGICAL,primary=True)
    anchor=g.input('anchor','Published SI5prime-lipid anchor sequence',K.OLIGO,state=Construct([Segment('anchor',ANCHOR_SI)]))
    barcode=g.input('barcode','Published SI sample-tag DNA with8base tag and A30',K.OLIGO,state=Construct([Segment('handle',HANDLE),
        Segment('sample tag','N'*8,'idx',placeholder=True,feature=feature('multi.biological_sample','sample_index','unknown')),
        Segment('polyA','A'*30)]))
    assembly=g.product('assembly','Noncovalent anchor/barcode duplex from1:1 stock; co-anchor supplied later',K.MOLECULAR)
    g.react('assemble','Pre-anneal anchor and sample-barcode oligos in1:1 molar ratio',S.HYBRIDIZATION,((anchor,R.TEMPLATE),(barcode,R.STRUCTURAL)),(assembly,))
    bound=g.product('bound','Anchor/barcode membrane-associated cells plus free reagent fractions',K.MOLECULAR)
    cold=g.input('cold','Published SI5min on ice; no BSA or serum during initial labeling',K.CONDITION)
    g.react('label','Associate anchor/barcode with cell membranes',S.BINDING,((cells,R.TEMPLATE),(assembly,R.STRUCTURAL),(cold,R.CONDITION)),(bound,),changes_molecule=False)
    co=g.input('co','Published SI3prime-lipid co-anchor',K.OLIGO,state=Construct([Segment('co-anchor',CO)]))
    stabilized=g.product('stabilized','Co-anchor-associated membrane sample tags plus free oligos',K.MOLECULAR,state=MembraneComplex().scene(),chemical_state=MembraneComplex())
    g.react('co_anchor','Add co-anchor;5min further on ice',S.HYBRIDIZATION,((bound,R.TEMPLATE),(co,R.STRUCTURAL),(cold,R.CONDITION)),(stabilized,))
    bsa=g.input('bsa','1ml ice-cold1%BSA inPBS per SI protocol',K.CONDITION)
    quenched=g.product('quenched','BSA-quenched sample; free LMO/BSA complexes remain until washes',K.MOLECULAR)
    g.react('quench','Quench free lipid-oligo labeling using BSA',S.QUENCHING,((stabilized,R.TEMPLATE),(bsa,R.CONDITION)),(quenched,))
    washed=cleanup(g,'cell_washes',quenched,'Washed labeled-cell pools','Centrifuge4°C;at least two ice-cold1%BSA/PBS washes;combine/filter pools')
    tenx=g.input('tenx','Published SI referenced10xworkflow through post-RT cDNA retrieval; historical molecular recipe unassembled',K.CONDITION)
    captured=g.product('post_rt','Source mixed endogenous cDNA and synthetic-DNA barcode-derived products; shared source cell identity/UMI',K.MOLECULAR)
    g.react('tenx_boundary','Source10x capture/RT handoff',S.MAGIC,((washed,R.TEMPLATE),(tenx,R.CONDITION)),(captured,),reason='Source explicitly follows10x workflow through this point; complete historical bead/cDNA polymer and copy geometry unreviewed. Synthetic DNA tags are not relabeled as RNA.')
    mix=g.input('amplification_mix','Source10xcDNA amplification master mix',K.CATALYST)
    tenxprimers=g.input('tenx_primers','Historical10xcDNA primer mix; exact members unassembled',K.OLIGO)
    additive=g.input('additive','Published SI19nt MULTI-seq additive primer',K.OLIGO,state=Construct([Segment('additive',ADDITIVE_SI)]))
    amplified=g.product('amplified','Amplified endogenous and barcode-derived populations in same reaction',K.MOLECULAR)
    g.react('cdna_amplification','Amplify with source10xprimer mix plus SI-specific additive primer',S.PCR,((captured,R.TEMPLATE),(tenxprimers,R.STRUCTURAL),(additive,R.STRUCTURAL),(mix,R.CATALYST)),(amplified,),geometry=PCRGeometry.unavailable(unknown('postRTmixedtemplate'),reason='Known additive sequence is not sufficient to assemble original10xcDNA substrate or all primer landings.'))
    spri=g.input('partition_spri','0.6×SPRI size selection',K.SELECTION)
    endogenous=g.product('endogenous','Bead-bound endogenous cDNA fraction retained for separate10x library preparation',K.MOLECULAR)
    tags=g.product('tag_supernatant','Barcode-derived DNA supernatant; RETAINED, not discarded',K.MOLECULAR)
    g.react('size_partition','Separate endogenous bead fraction from barcode supernatant',S.SELECTION,((amplified,R.TEMPLATE),(spri,R.SELECTION)),(endogenous,tags))
    g.mark_terminal(Terminal(endogenous.id,Disposition.UNRESOLVED,'Original10x endogenous library preparation and run geometry remain unassembled; this fraction is physically retained'))
    clean=cleanup(g,'barcode_enrichment',tags,'Recovered barcode-derived cDNA','3.2×finalSPRI plus1.8×isopropanol;ethanol washes;EBelution')
    selected=aliquot(g,'barcode_aliquot',clean,'3.5ng barcode-derived cDNA into libraryPCR')
    p5=g.input('p5','Published SI UniversalI5/P5 primer',K.OLIGO,state=Construct([Segment('P5 primer',P5)]))
    rpi=g.input('rpi','Published SI TruSeqRPIX with unassigned6base library index',K.OLIGO,state=rpi_primer())
    kapa=g.input('kapa','KAPAHiFiHotStart;SI prints repeat98/60/72°C steps8–12times; total-cycle interpretation unassigned',K.CATALYST)
    library=g.product('library_pcr','Barcode library-PCR products with source sample tags and distinct sequencer index',K.MOLECULAR)
    g.react('library_amplify','Amplify barcode library using SI-specific workflow',S.PCR,((selected,R.TEMPLATE),(p5,R.STRUCTURAL),(rpi,R.STRUCTURAL),(kapa,R.CATALYST)),(library,),geometry=PCRGeometry.unavailable(unknown('barcode-derivedcDNA'),reason='Full first-strand/cDNA end geometry and run-primer landing remain unassembled; two source index meanings stay distinct.'))
    final=cleanup(g,'library_cleanup',library,'Purified MULTI-seq barcode-library material','1.6×SPRI cleanup;ethanol washes;EBelution')
    return limit(g,final)

def detail():
    return info('Version boundary: published SI protocol anchor='+ANCHOR_SI+' and additive='+ADDITIVE_SI+'. Associated accepted manuscript anchor='+ANCHOR_MANUSCRIPT+' and additive='+ADDITIVE_MANUSCRIPT+'. The primary workflow follows published SI only; manuscript interface is a separately labelled comparison, not a published-error finding or proof of which reagent was used.')+panel(MembraneComplex().scene().rows(),caption='Published SI molecular-interface variant')+panel(MembraneComplex('accepted_manuscript').scene().rows(),caption='Associated accepted-manuscript interface variant; not merged into SI workflow')
def source_checks(c,g):
    for variant in ('published_si','accepted_manuscript'):
        state=MembraneComplex(variant);sc=state.scene()
        c('Separate physical strands retained for '+variant,len(sc.strands),3)
        c('Version-specific anchor conserved through Scene for '+variant,''.join(s.top for s in sc.strands['anchor'].segs),state.anchor)
    c('Barcode supernatant supplies barcode enrichment',g.steps['barcode_enrichment'].inputs[0].material,'tag_supernatant')
    c('Endogenous cDNA retained independently',g.terminals['endogenous'].disposition,Disposition.UNRESOLVED)
    c('Biological sample tag differs from library index',rpi_primer().get('sequencer library index').feature.id,'multi.library_index')

REVIEWED=('cells', 'anchor', 'barcode', 'assembly', 'cold', 'bound', 'co', 'stabilized', 'bsa', 'quenched', 'cell_washes_selection', 'cell_washes_kept', 'cell_washes_waste', 'tenx', 'post_rt', 'amplification_mix', 'tenx_primers', 'additive', 'amplified', 'partition_spri', 'endogenous', 'tag_supernatant', 'barcode_enrichment_selection', 'barcode_enrichment_kept', 'barcode_enrichment_waste', 'barcode_aliquot_chosen', 'barcode_aliquot_remainder', 'p5', 'rpi', 'kapa', 'library_pcr', 'library_cleanup_selection', 'library_cleanup_kept', 'library_cleanup_waste', 'assemble', 'label', 'co_anchor', 'quench', 'cell_washes', 'tenx_boundary', 'cdna_amplification', 'size_partition', 'barcode_enrichment', 'barcode_aliquot', 'library_amplify', 'library_cleanup')
