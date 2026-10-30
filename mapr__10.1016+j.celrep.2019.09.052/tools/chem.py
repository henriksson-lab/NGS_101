"""Protocol-owned source-limited model; no finished-publication marker."""
from dataclasses import dataclass
from pathlib import Path
import sys
HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE), str(HERE.parents[1] / 'lib')]
from chemdraw import Scene, Segment, Construct, feature, revcomp, oligo, panel, workflow_from_protocol_graph, workflow_panel
from protocol_model import role, molecule
from protocol_graph import ProtocolGraph, InputRole as R, MaterialKind as K, StepKind as S, Terminal, Disposition, SourceItem
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
    c('Every DNA molecular output continues to an operation or explicit terminal',not any(m.kind in {K.MOLECULAR,K.LIBRARY} and m.id not in g.terminals and not any(port.material==m.id for step in g.steps.values() for port in step.inputs) for m in g.materials.values()))
    c('MNase activation uses already introduced fusion, not a new catalytic aliquot',any(port.material=='target_stock' and port.role is R.CATALYST for port in g.steps['mnase'].inputs))
    c('RNase/release mixture is not physically discarded',all(g.materials[mid].kind is not K.WASTE for mid in g.steps['stop_rna_release'].outputs))
    from checks import Source,have
    from storage import sources_root
    src=Source(sources_root()/HERE.parent.name/'j.celrep.2019.09.052_PMC6870988.xml.txt','Associated original MapR manuscript, standard get_sources')
    if have(c,src,label='internal source values'):
        text=src.path.read_text()
        c('Detailed printed RNaseA token remains source-qualified','5 0μg RNaseA' in text)
        c('Detailed proteinase stock mass differs from short5microgram method','2.5 μL Proteinase K (20 mg/ml)' in text and '5 μg proteinase K' in text)
    c.report()

def unassembled(label):
    return mixture(('source-derived DNA',(role(label),)),('possible carryover',(role('other original DNA/RNA, reaction fragments or protein carryover; quantities and exact partition unassigned'),)))

def operation(g,ident,label,kind,x,state,reagent,*,reagent_kind=K.CONDITION,reagent_role=R.CONDITION):
    stock=g.input(ident+'_stock',reagent,reagent_kind)
    out=g.product(ident+'_product',label,K.MOLECULAR,state=state)
    kw={'reason':'Source operation established; exact molecular geometry/mixture relation not computed.'} if kind in {S.MAGIC,S.GENERIC} else {}
    if kind is S.GENERIC:kw['changes_molecule']=False
    g.react(ident,label,kind,((x,R.TEMPLATE),(stock,reagent_role)),(out,),**kw)
    return out

def ligate(g,ident,label,x,adapter_label):
    adapter_scene=molecule(adapter_label,(role('source adapter; exact historical sequence/modifications and geometry unassembled'),))
    adapter=g.input(ident+'_adapter',adapter_label,K.ADAPTER,state=adapter_scene)
    enzyme=g.input(ident+'_enzyme','source T4DNA ligase',K.CATALYST)
    out=g.product(ident+'_product',label,K.MOLECULAR,state=unassembled('source adapter-ligated product roles; no exact sealed junction claimed'))
    geometry=AdapterLigationGeometry('unavailable','unavailable',(adapter_scene,), 'Historical exact end/junction assembly remains unreviewed; no source-absence claim.')
    g.react(ident,label,S.LIGATION,((x,R.TEMPLATE),(adapter,R.STRUCTURAL),(enzyme,R.CATALYST)),(out,),geometry=geometry)
    return out

def amplify(g,ident,label,x,primer_label):
    f=g.input(ident+'_forward',primer_label+' forward member; exact member unassigned',K.OLIGO,state=molecule('forward primer',(role('source primer identity unassembled'),)))
    r=g.input(ident+'_reverse',primer_label+' reverse member; exact member unassigned',K.OLIGO,state=molecule('reverse primer',(role('source primer identity unassembled'),)))
    enzyme=g.input(ident+'_enzyme','source Q5DNApolymerase',K.CATALYST)
    out=g.product(ident+'_product',label,K.LIBRARY,state=unassembled('source amplified library and residual original templates/cleaved species; no molecular deletion assigned to PCR'))
    g.react(ident,label,S.PCR,((x,R.TEMPLATE),(f,R.STRUCTURAL),(r,R.STRUCTURAL),(enzyme,R.CATALYST)),(out,),geometry=PCRGeometry.unavailable(x.state,reason='Exact original adapters/primer/member/strand geometry unassembled.'))
    return out

PAGE_FILE='mapr.html'
TITLE='MapR'
NOTES='01_mapr.html'
SOURCE='<a href="https://doi.org/10.1016/j.celrep.2019.09.052">Defining article</a>; Associated author manuscript PMC6870988 main58, Methods160–244 and SI published FigureS1–S5 legends'
CAVEAT='Selected original2019 detailed MapR branch, not RHΔCUT&RUN control or later2020 revision. Short Methods162 specify5µgRNaseA and5µgProteinaseK, whereas Step-by-Step214/220 print ambiguous RNaseA token5 0µg and2.5µl20mg/mlProteinaseK(50µg); both values remain visible, detailed branch selected without claiming correction. Exact adapter/read geometry and linked CUT&RUN original remain unreviewed.'
READOUT='Source NextSeq500; exact read role/cycle/index geometry unassembled. RNaseA handles RNA before dsDNA library preparation; no perfect strand reannealing or R-loop purity asserted.'
PREVIEWS=(("Source-qualified operation graph","lm_pcr_product"),)

@dataclass(frozen=True)
class HybridPopulation:
    stage: str
    def __post_init__(self):
        if self.stage not in {'native','cut','rna_processed','bisulfite','second_strand','udg'}:raise ValueError('unsupported selected source population')
    def rows(self):
        sc=Scene()
        sc.strand('template',[role('original hybrid-template genomic DNA; source coordinates and sequence unassigned')],label='original template DNA')
        if self.stage in {'native','cut','rna_processed','bisulfite','second_strand','udg'}:
            sc.strand('displaced',[role('original displaced genomic strand; exposed cytosines selected as C-to-U class' if self.stage in {'bisulfite','second_strand'} else 'UDG-treated displaced-DNA/AP-product class; fragmented state not computed' if self.stage=='udg' else 'original displaced genomic DNA; no U label assigned')],label='displaced DNA, independent strand')
        sc.strand('rna',[role('RNA-fragment population and possible residual hybrid RNA; no physical discard' if self.stage in {'rna_processed','second_strand','udg'} else 'RNA associated noncovalently with template DNA; not a covalent DNA extension')],label='RNA or RNA-derived products')
        if self.stage in {'second_strand','udg'}:
            sc.strand('new_dna',[role('new dUTP-containing DNA second strand' if self.stage=='second_strand' else 'UDG-treated dU-DNA/AP-product class; exact secondary fragmentation unassembled')],label='new DNA, independent polymer')
        if self.stage=='udg':sc.strand('free_bases',[role('excised uracil and other reaction products remain until actual cleanup')],label='reaction products, not sequenced base calls')
        sc.strand('other',[role('other genomic DNA/protein and possible background or digestion fragments; recovery unquantified')],label='other physical material')
        return sc.rows()

def protocol_graph():
    g=ProtocolGraph(TITLE,authority=SOURCE,representation='published_instance')
    current=g.input('cells','unfixed source cells with native R-loop population',K.MOLECULAR,state=HybridPopulation('native'),primary=True)
    current=cleanup(g,'cell_wash','Source wash buffer; centrifuge600g3min and repeat',current,current.state)
    current=operation(g,'immobilize','ConcanavalinA beads, source binding buffer preparation and1h roomtemperature binding',S.BINDING,current,current.state,'source washed ConA-coated beads',reagent_kind=K.SELECTION,reagent_role=R.SELECTION)
    current=cleanup(g,'magnetic_transfer','Magnetic removal of liquid after source cell immobilization',current,current.state)
    current=operation(g,'permeabilize','Source0.02%digitonin Dig-wash permeabilization',S.GENERIC,current,current.state,'source Dig-wash buffer')
    current=operation(g,'target','Source GST-RHdelta-MNase1µM overnight4C; MNase-only separate control',S.MAGIC,current,current.state,'purified source RNaseH-inactive/MNase fusion',reagent_kind=K.CATALYST,reagent_role=R.CATALYST)
    current=cleanup(g,'target_washes','Three source Dig-wash magnetic washes',current,current.state)
    calcium=g.input('calcium','source calcium activation of already introduced fusion',K.CONDITION)
    fragmented=g.product('mnase_product','source fusion-targeted cleavage population before physical release',K.MOLECULAR,state=HybridPopulation('cut'))
    g.react('mnase','Activate already bound MNase with source calcium; do not add another enzyme aliquot',S.DIGESTION,((current,R.TEMPLATE),(g.materials['target_stock'],R.CATALYST),(calcium,R.CONDITION)),(fragmented,))
    current=fragmented
    spike=g.input('yeast_spike','source heterologous yeast DNA in STOP buffer; no exact sequence assigned',K.MOLECULAR,state=unassembled('source spike-in genomic DNA'))
    stop=g.input('stop_buffer','source EDTA/EGTA STOP mixture with RNaseA; detailed printed5 0µg versus short5µg; numerical typography unresolved',K.CONDITION)
    released=g.product('release_mixture','stopped,37C10min RNaseA/release mixture, before physical centrifugation',K.MOLECULAR,state=HybridPopulation('rna_processed'))
    g.react('stop_rna_release','Source STOP/RNaseA and37C10min release; no physical deletion of RNA products',S.MAGIC,((current,R.TEMPLATE),(spike,R.STRUCTURAL),(stop,R.CONDITION)),(released,),reason='Coupled source stop/RNaseA/release envelope; exact RNA digestion and reannealing not predicted.')
    current=cleanup(g,'supernatant','Source16000g5min4C and magnetic separation; transfer supernatant',released,released.state)
    current=operation(g,'proteolyze','SDS/proteinaseK70C10min; detailed2.5µl20mg/ml(50µg) versus short5µg discrepancy',S.DEPROTEINIZATION,current,current.state,'source ProteinaseK withSDS; peptide products retained',reagent_kind=K.CATALYST,reagent_role=R.CATALYST)
    current=cleanup(g,'organic','Source phenol/chloroform/isoamyl extraction and phase-lock aqueous transfer',current,unassembled('recovered MapR DNA; RNA/protein-product carryover unquantified'))
    current=cleanup(g,'ethanol','Source ammoniumacetate/linearacrylamide/ethanol−80C1h; pellet/75%ethanolwash/redissolve',current,current.state)
    current=operation(g,'end_repair','Source End-It DNA end repair',S.MAGIC,current,unassembled('source end-repaired DNA; exact strand ends unassembled'),'source End-It repair reagents')
    current=operation(g,'a_tail','Source Klenowexo-minus dA-tailing',S.MAGIC,current,unassembled('source dA-tailed DNA; exact junction geometry unassembled'),'source Klenowexo-minus and nucleotides')
    current=ligate(g,'adapter_ligation','Source custom adapter/T4DNA ligation',current,'source custom Illumina-compatible adapter')
    current=cleanup(g,'size_select','Source SPRI selection for fragments greater150bp',current,current.state)
    current=amplify(g,'lm_pcr','Source Q5 ligation-mediated PCR; cycle count unassigned',current,'source custom barcoded Illumina primers')
    g.mark_terminal(Terminal(current.id,Disposition.UNRESOLVED,'Source NextSeq500 endpoint; complete molecular/read geometry and named source discrepancies unresolved.'))
    # Display terminal aliases are physical material IDs, not new molecular transformations.
    g.add_source_item(SourceItem('pending_originals','CUT&RUN2018 original and historical custom adapter/run construction','Main Methods162 cites Skene2018 DOI10.1038/nprot.2018.015; exact component/version and adapters unreviewed.'))
    return g
