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
import rt
TITLE='Gel-seq — pre-separation RT and independent genomic/cDNA branches'
NOTE='01_gel_seq.html'
SOURCE='Published LabChip2017 VOR DOI10.1039/c7lc00430c; institutional qt4120h56x.pdf pp2620–2624'
LIMIT='Low-input100–1000cell route, not a single-cell barcode method. Historical SMART-Seqv4 oligos, full cDNA geometry, kit-specific AMPure workup and two NexteraXT library preparations/read layouts remain unassembled. Gel retention is qualitative source evidence, not a pore-size, diffusion or purity calculation.'
SOURCE_FILE='qt4120h56x.pdf.txt'
RANDOM_PRIMER=rt.SMART_HANDLE+'AC'+'N'*6
TRANSCRIBED=(RANDOM_PRIMER,)
ENDPOINTS=('gdna_concentrated','cdna_purified_kept')
@dataclass(frozen=True)
class GelFraction:
    population:str
    position:str
    def __post_init__(self):
        expected={'genomic DNA':'low-density loading boundary','RNA/cDNA hybrids':'high-density interface'}
        if self.population not in expected or self.position!=expected[self.population]:raise ValueError('Only the two source-described retention populations/positions are represented')
    def scene(self):
        sc=Scene();sc.strand('population',[Segment(self.population,'[variable '+self.population+']',placeholder=True)],label=self.population+' at '+self.position)
        sc.note('population','Source-declared enrichment/retention; no assertion of perfect separation or quantified diffusion')
        return sc
G=GelFraction('genomic DNA','low-density loading boundary')
C=GelFraction('RNA/cDNA hybrids','high-density interface')

def protocol_graph():
    g=ProtocolGraph(TITLE,authority=SOURCE)
    cells=g.input('cells','100–1000 intact cells in source diluted SMART-Seqv4 lysis buffer',K.BIOLOGICAL,primary=True)
    lys=g.input('lysis','SMART-Seqv4 lysis/RNase inhibitor;72°C3min',K.CONDITION)
    lysate=g.product('lysate','Coexisting genomic DNA and cellular RNA after lysis',K.MOLECULAR)
    g.react('lyse','Lyse cells while limiting genomic DNA shearing',S.GENERIC,((cells,R.TEMPLATE),(lys,R.CONDITION)),(lysate,),changes_molecule=False,reason='Physical cell lysis; no purification or DNA amplification is inferred')
    primer=g.input('rt_primers','3prime SMART-Seq CDSPrimerII and source random hexamer+SMART handle',K.OLIGO)
    tso=g.input('tso','Source SMART-Seqv4 oligos; polymers unassembled',K.OLIGO)
    enzyme=g.input('reverse_transcriptase','SMARTScribe;42°C90min then70°C10min',K.CATALYST)
    mixed=g.product('mixed','Genomic DNA and template-switched RNA/cDNA hybrid populations before separation',K.MOLECULAR)
    g.react('rt','Reverse transcribe RNA BEFORE gel separation; preserve genomic material',S.REVERSE_TRANSCRIPTION,((lysate,R.TEMPLATE),(primer,R.STRUCTURAL),(tso,R.STRUCTURAL),(enzyme,R.CATALYST)),(mixed,))
    gel=g.input('gel_device','Source4%T/3%C low-density and30%T/5%C high-density polyacrylamide layers;210V30min',K.SELECTION)
    gdna=g.product('gdna_gel','Excised genomic-DNA-containing gel region',K.MOLECULAR,state=G.scene(),chemical_state=G)
    cdna=g.product('cdna_gel','Separate excised RNA/cDNA-hybrid gel region',K.MOLECULAR,state=C.scene(),chemical_state=C)
    rejected=g.product('gel_reject','Uncollected gel/buffer regions',K.WASTE)
    g.react('separate','Electrophorese then excise two spatially distinct source fractions',S.SELECTION,((mixed,R.TEMPLATE),(gel,R.SELECTION)),(gdna,cdna,rejected));g.discard(rejected,'Physical gel-region collection boundary; retained fractions are not assumed pure')
    recovered={}
    for key,x,volume in (('genomic',gdna,40),('transcript',cdna,80)):
        soak=g.input(key+'_soak',str(volume)+'µl water;crush/soak37°C8–12h',K.CONDITION)
        soluble=g.product(key+'_soaked','Eluted '+key+' nucleic acids mixed with gel fragments',K.MOLECULAR)
        g.react(key+'_elution','Diffuse nucleic acids out of excised/crushed gel',S.EXTRACTION,((x,R.TEMPLATE),(soak,R.CONDITION)),(soluble,))
        recovered[key]=cleanup(g,key+'_filter',soluble,'Gel-filtered '+key+' fraction','8µm mesh filter;2600RCF5min')
    protease=g.input('protease','Qiagen protease0.9AU/ml;50°C15min then70°C15min',K.CATALYST)
    deproteinized=g.product('deproteinized','Genomic DNA plus proteolysis products; nucleosomes depleted, not physically purified by protease',K.MOLECULAR)
    g.react('deproteinize','Proteolyse genomic fraction before NexteraXT',S.DEPROTEINIZATION,((recovered['genomic'],R.TEMPLATE),(protease,R.CATALYST)),(deproteinized,))
    gcon=g.product('gdna_concentrated','5µl concentrated genomic fraction, including nonvolatile proteolysis products',K.MOLECULAR)
    water=g.product('genomic_water','Evaporated water fraction',K.WASTE)
    vacuum=g.input('vacuum','Source vacufuge; do not equate concentration with DNA purification',K.CONDITION)
    g.react('genomic_concentration','Reduce genomic sample volume to5µl',S.CLEANUP,((deproteinized,R.TEMPLATE),(vacuum,R.CONDITION)),(gcon,water));g.discard(water,'Water physically evaporated; proteins/reaction products not assumed removed')
    tcon=g.product('cdna_concentrated','10µl concentrated RNA/cDNA fraction',K.MOLECULAR)
    twater=g.product('transcript_water','Evaporated transcript-fraction water',K.WASTE)
    g.react('transcript_concentration','Reduce RNA/cDNA sample volume to10µl',S.CLEANUP,((recovered['transcript'],R.TEMPLATE),(vacuum,R.CONDITION)),(tcon,twater));g.discard(twater,'Only water removal claimed')
    single=g.input('pcr_primer','Source PCRPrimerIIA12µM, one declared species; full polymer unassembled',K.OLIGO)
    poly=g.input('pcr_mix','KAPA SYBRFast;20–30cycles chosen from qPCR curves',K.CATALYST)
    amplified=g.product('amplified_cdna','Amplified cDNA plus unassigned original RNA and reaction-material fractions',K.MOLECULAR)
    g.react('cdna_pcr','Amplify cDNA before its NexteraXT library preparation',S.PCR,((tcon,R.TEMPLATE),(single,R.STRUCTURAL),(poly,R.CATALYST)),(amplified,),
        geometry=PCRGeometry.unavailable(unknown('source SMART cDNA'),symmetric_species=True,reason='Source specifies one PCRPrimerIIA species; exact full historical cDNA handles/primer polymer are unassembled, not supplied by platform name.'))
    clean=cleanup(g,'cdna_purified',amplified,'AMPure-recovered amplified cDNA; carryover unquantified','Original SMART-Seq manual AMPure procedure; detailed ratios unresolved')
    g.mark_terminal(Terminal(gcon.id,Disposition.UNRESOLVED,'Half-volume original NexteraXT genomic library preparation/read geometry remains unassembled; concentration retained nonvolatile products'))
    return limit(g,clean)

def detail():return panel(G.scene().rows(),caption='Genomic fraction')+panel(C.scene().rows(),caption='Separate transcript-derived fraction; no shared molecule is implied')
def source_checks(c,g):
    c('Reverse transcription precedes physical fractionation',g.steps['separate'].inputs[0].material,'mixed')
    c('Distinct genomic and transcript fractions survive selection',tuple(g.steps['separate'].outputs),('gdna_gel','cdna_gel','gel_reject'))
    c('Proteolysis is restricted to genomic branch',g.steps['deproteinize'].inputs[0].material,'genomic_filter_kept')
    c('Random hexamer is priming sequence, not a declared UMI',len(RANDOM_PRIMER)-len(rt.SMART_HANDLE+'AC'),6)

REVIEWED=('cells', 'lysis', 'lysate', 'rt_primers', 'tso', 'reverse_transcriptase', 'mixed', 'gel_device', 'gdna_gel', 'cdna_gel', 'gel_reject', 'genomic_soak', 'genomic_soaked', 'genomic_filter_selection', 'genomic_filter_kept', 'genomic_filter_waste', 'transcript_soak', 'transcript_soaked', 'transcript_filter_selection', 'transcript_filter_kept', 'transcript_filter_waste', 'protease', 'deproteinized', 'vacuum', 'gdna_concentrated', 'genomic_water', 'cdna_concentrated', 'transcript_water', 'pcr_primer', 'pcr_mix', 'amplified_cdna', 'cdna_purified_selection', 'cdna_purified_kept', 'cdna_purified_waste', 'lyse', 'rt', 'separate', 'genomic_elution', 'genomic_filter', 'transcript_elution', 'transcript_filter', 'deproteinize', 'genomic_concentration', 'transcript_concentration', 'cdna_pcr', 'cdna_purified')
