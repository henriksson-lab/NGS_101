"""Protocol-owned source-limited model; no finished-publication marker."""
from dataclasses import dataclass
from pathlib import Path
import sys
HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE), str(HERE.parents[1] / 'lib')]
from chemdraw import Scene, Segment, Construct, feature, revcomp, oligo, panel, workflow_from_protocol_graph, workflow_panel
from protocol_model import role, molecule
from protocol_graph import ProtocolGraph, InputRole as R, MaterialKind as K, StepKind as S, Terminal, Disposition, SourceItem, SequencingPlan
from adapter_ligation import AdapterLigationGeometry, reported_adapter_end
from pcr import PCRGeometry
from page import head, info

def mixture(*species):
    sc = Scene()
    for n, (label, parts) in enumerate(species):
        sc.strand('s' + str(n), list(parts), label=label)
    return sc

def cleanup(g, ident, source, template, state):
    reagent = g.input(ident + '_condition', source, K.CONDITION)
    out = g.product(ident, 'recovered source preparation', K.MOLECULAR, state=state)
    waste = g.product(ident + '_waste', 'physically separated workup fractions; exact molecular partition unquantified', K.WASTE)
    g.react(ident + '_step', source, S.CLEANUP, ((template, R.TEMPLATE), (reagent, R.CONDITION)), (out, waste))
    g.discard(waste, 'Actual source physical separation; no perfect recovery/depletion asserted.')
    return out

def end_boundary(label):
    sc = molecule('unassembled ligation-ready geometry', (role(label),))
    return AdapterLigationGeometry('unavailable', 'unavailable', (sc,), label + '; no source-absence claim.')

def render_page():
    import html
    g = protocol_graph()
    parts = [head(TITLE), '<main class="wrap"><h1>' + html.escape(TITLE) + '</h1>',
             '<p><a href="' + NOTES + '">Research notes</a></p>', info(SOURCE),
             info('Private incomplete model. ' + CAVEAT)]
    if callable(globals().get('oligos')):
        parts.extend(['<h2>Source ordering oligos</h2>', *oligos()])
    for label, terminal in PREVIEWS:
        parts.extend(['<h2>' + html.escape(label) + '</h2>', workflow_panel(workflow_from_protocol_graph(g, terminal))])
    parts.extend([info(READOUT), '</main>'])
    return '\n'.join(parts)

def selftest():
    from checks import Check
    c = Check(); g = protocol_graph()
    c('Pending source/assembly boundary prevents finished publication', bool(g.publication_problems()))
    c('Private preview has no finished marker', 'data-protocol-complete=' not in render_page())

    x=MethylatedDNA('AGCGCA',frozenset({2}),frozenset({2,4}))
    c('GCG physicalmark has overlapping possible origins',x.physical_5mc,frozenset({2,4}))
    c('GCG provenance remains ambiguous',x.context_class(2),'ambiguous GCG')
    c('GpC outsideCpG supplies distinct accessibility context',x.context_class(4),'GCH accessibility context')
    c('Marking does not convert genomicC intoT',x.sequence,'AGCGCA')
    c('Source acquisition retained despite unknown primer geometry',len(g.sequencing),1)
    c.report()

def pop(label):
    return mixture(('source DNA',(role(label),)),('other physical material',(role('other original DNA/RNA/protein or reaction remnants; exact partition and recovery unquantified'),)))
def op(g,key,label,kind,x,out_state,reagent_label,reagent_kind=K.CONDITION,reagent_role=R.CONDITION):
    reagent=g.input(key+'_stock',reagent_label,reagent_kind)
    product=g.product(key+'_product',label,K.MOLECULAR,state=out_state)
    kw={'reason':'Source operation known; exact molecular substrate/product relation and mixture fractions unassembled.'} if kind in {S.MAGIC,S.GENERIC} else {}
    if kind is S.GENERIC:kw['changes_molecule']=False
    g.react(key,label,kind,((x,R.TEMPLATE),(reagent,reagent_role)),(product,),**kw)
    return product

PAGE_FILE='nano-nome.html'
TITLE='nanoNOMe'
NOTES='01_nano-nome.html'
SOURCE='<a href="https://doi.org/10.1038/s41592-020-01000-7">Lee et al.2020</a>; associated manuscript PMC7704922 main71 and Methods139–145; published SI Figure1/call-interpretation scope.'
CAVEAT='Selected approximately10kb original LSK-SQK108 nanopore route;20kb/Circulomics samples are separate source options. Validation WGBS/PCR and ChIP-qPCR are not native nanoNOMe preparation. Native and exogenous5mC chemical identity is the same, GCG origin is ambiguous and excluded from source computational analysis. No complete methylation/accessibility specificity or accurate-call guarantee; original historical kit/adaptor/read geometry and finalVOR method applicability remain unreviewed.'
READOUT='Source unamplified long nanopore molecules retain methylated cytosines. Nanopolish dual-context calls and accessibility/footprint/phasing analyses are computational; an unmethylated call is not proof of protein occupancy. Source MinION Mk1b/GridION/PromethION, FLO-MIN106 orPRO-002, up to72h.'
PREVIEWS=(('Selected10kb native methylation route','library'),)
@dataclass(frozen=True)
class MethylatedDNA:
    sequence:str
    native_cpg:frozenset[int]=frozenset()
    exogenous_gpc:frozenset[int]=frozenset()
    def __post_init__(self):
        if not self.sequence or set(self.sequence)-set('ACGT'):raise ValueError('explicit source-independent DNA substrate required')
        for values,context in ((self.native_cpg,'cpg'),(self.exogenous_gpc,'gpc')):
            if not isinstance(values,frozenset) or any(type(i)is not int or not 0<=i<len(self.sequence) or self.sequence[i]!='C' or (context=='cpg' and (i+1>=len(self.sequence) or self.sequence[i+1]!='G')) or (context=='gpc' and (i==0 or self.sequence[i-1]!='G')) for i in values):raise ValueError('selected5mC requires an actual cytosine in the declared context')
    @property
    def physical_5mc(self):return self.native_cpg|self.exogenous_gpc
    def context_class(self,i):
        if type(i)is not int or not 0<=i<len(self.sequence) or self.sequence[i]!='C':raise ValueError('actual cytosine position required')
        if i==0 or i+1==len(self.sequence):return 'incomplete context'
        left,right=self.sequence[i-1],self.sequence[i+1]
        if left=='G' and right=='G':return 'ambiguous GCG'
        if left=='G':return 'GCH accessibility context'
        if right=='G':return 'HCG endogenous context'
        return 'other context'

@dataclass(frozen=True)
class NOMePopulation:
    stage:str
    def rows(self):
        if self.stage not in {'nuclei','methylated','proteolyzed','sheared','end_prepared','library'}:raise ValueError('unsupported nanoNOMe state')
        sc=Scene()
        sc.strand('genomic',[role('source leaderDNA-linked genomic library with nativeCpG/GpC5mC; exact strand/junction geometry unreviewed' if self.stage=='library' else 'original genomic DNA with nativeCpG and added accessibleGpC5mC classes; sites/fractions unassigned' if self.stage!='nuclei' else 'original genomic DNA with native endogenous methylation; sites/fractions unassigned')],label='genomic DNA, noC-to-T base conversion')
        sc.strand('other',[role('proteinase peptides and possible other DNA/RNA/carryover' if self.stage in {'proteolyzed','sheared','end_prepared','library'} else 'source nuclear proteins/RNA and other material')],label='other physical material')
        if self.stage=='library':sc.strand('motor',[role('source motor protein associated with the adaptor-linkedDNA; no covalent protein-DNA topology assigned')],label='motor-associated component, not a separate unligatedDNA adaptor')
        return sc.rows()

def protocol_graph():
    g=ProtocolGraph(TITLE,authority=SOURCE,representation='published_instance')
    cells=g.input('cells','source individual cell-line preparation; cell lines not physically pooled',K.MOLECULAR,state=NOMePopulation('nuclei'),primary=True)
    current=op(g,'nuclei_lysis','Source trypsinized cells/0.25%NP40 resuspensionbuffer5minonice',S.MAGIC,cells,cells.state,'source100mMTris/100mMNaCl/30mMMgCl2 andNP40; nuclei isolation, not chromatin fixation')
    current=cleanup(g,'nuclei_pellet','Source500g5min4C nuclear recovery',current,current.state)
    enzyme=g.input('mcvi','source200UM.CviPI per500000nuclei/500µl',K.CATALYST)
    sam=g.input('sam','source96µMSAM, replenished at7.5min, with300mMsucrose/1xbuffer',K.MOLECULAR,state=molecule('methyl donor',(role('sourceSAM methyl donor; product molecules unassembled'),)))
    marked=g.product('marked','source15min37C1000rpm methylated nuclear preparation',K.MOLECULAR,state=NOMePopulation('methylated'))
    g.react('gpc_mark','Source GpC methyltransferase marking; no ideal access/saturation guarantee',S.CHEMICAL_CONVERSION,((current,R.TEMPLATE),(enzyme,R.CATALYST),(sam,R.STRUCTURAL)),(marked,))
    current=op(g,'stop','Source equal-volumeTris/600mMNaCl/1%SDS/10mMEDTA stop solution',S.QUENCHING,marked,marked.state,'source stop solution; no physical removal at quench')
    current=op(g,'proteinase','Source proteinaseK55C over2h; peptide remnants retained until extraction',S.DEPROTEINIZATION,current,NOMePopulation('proteolyzed'),'sourceProteinaseK',K.CATALYST,R.CATALYST)
    current=cleanup(g,'organic','Source phenol:chloroform extraction; wide-orifice tips avoid excessive fragmentation',current,current.state)
    current=cleanup(g,'ethanol','Source ethanol precipitation recovery',current,current.state)
    current=op(g,'shear','Selected source G-tube approximately10kb:2–3µg,5000g1min/invert/repeat',S.FRAGMENTATION,current,NOMePopulation('sheared'),'sourceCovarisG-tube; optional20kb/Circulomics branch unassembled separately')
    current=op(g,'end_prep','Source NEBNextUltraII end repair and dA-tail',S.MAGIC,current,NOMePopulation('end_prepared'),'historicalNEBNextUltraIIend-repair module; exact sequence/end populations unreviewed')
    current=cleanup(g,'end_cleanup','Source1xAMPureXP after end preparation',current,current.state)
    donor_scene=Scene();donor_scene.strand('leader',[role('sourceleaderDNA adaptor; exact polymer/end geometry unreviewed')],label='leaderDNA');donor_scene.strand('motor',[role('source motor-protein component; DNA association topology unreviewed')],label='motor-associated component, not a DNA base')
    donor=g.input('leader_motor','source LSK-SQK108 leaderDNA and motor-protein adaptor',K.ADAPTER,state=donor_scene)
    ligase=g.input('ligase','source NEBBlunt/TALigaseMasterMix',K.CATALYST)
    library=g.product('ligated','source leader/motor adaptor-ligated methylated DNA',K.LIBRARY,state=NOMePopulation('library'))
    geometry=AdapterLigationGeometry('unavailable','unavailable',(donor_scene,),'Source end-repaired/dA-tailed substrate and leader/motor adaptor stated; exact cohesive/blunt strand assembly unreviewed, no source-absence claim.')
    g.react('ligate','Source leaderDNA/motor adaptor ligation',S.LIGATION,((current,R.TEMPLATE),(donor,R.STRUCTURAL),(ligase,R.CATALYST)),(library,),geometry=geometry)
    final=cleanup(g,'library','Source0.4xAMPureXP/sequencingkit-reagent cleanup; greater400ngloaded per flowcell',library,library.state)
    g.sequence(final,SequencingPlan(limitation='Source native MinION Mk1b/GridION/PromethION acquisition; motor-assisted translocation has no sequencing primer. No arbitrary nucleotide read length is assigned; historical leader/adaptor/end geometry and kit version remain unreviewed.'))
    g.add_source_item(SourceItem('historical_kit','Original LSK-SQK108 and NEBNextUltraII documents; finalVOR applicability','Associated mainMethods145 names exact historical kits; exact versions/assembly unreviewed.'))
    return g
