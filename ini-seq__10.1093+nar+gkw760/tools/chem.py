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
    c('Original nuclei assay and historical kit obligation stays uncovered','pending_originals' in g.source_coverage,False)
    c('HeLa extract contributes physical source material',next(p.role for p in g.steps['initiate'].inputs if p.material=='hela_extract'),R.STRUCTURAL)
    c('Private preview has no finished marker', 'data-protocol-complete=' not in render_page())
    c('Every DNA molecular output continues to an operation or explicit terminal',not any(m.kind in {K.MOLECULAR,K.LIBRARY} and m.id not in g.terminals and not any(port.material==m.id for step in g.steps.values() for port in step.inputs) for m in g.materials.values()))
    c('Source input and IP remain separate downstream libraries',{'library','input_library'}<=set(g.terminals))
    c('RNase products remain molecular until source organic extraction',all(g.materials[mid].kind is not K.WASTE for mid in g.steps['rnase'].outputs))
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

PAGE_FILE='ini-seq.html'
TITLE='ini-seq'
NOTES='01_ini-seq.html'
SOURCE='<a href="https://doi.org/10.1093/nar/gkw760">Defining article</a>; Published PMC5137433 Fig1/main54 and Methods62–80; published File010SI FigS1–S8 and TableS1(qPCR) scope'
CAVEAT='Selected DIG-labelled cell-free origin initiation, not EdU/BrdU live-cell timing. Standard assay conditions/nuclei preparation defer to cited sources; source ChIP-seq sample-prep kit version and molecular library/read geometry remain unreviewed. qPCR validation primers inTableS1 are not NGS adapters. No complete origin specificity or uniform DIG incorporation assumed.'
READOUT='Source HiSeq2000; exact acquired read roles/cycles/index members unassigned. Enrichment near origins is not a direct physical read of replication initiation nucleotide.'
PREVIEWS=(("Source-qualified operation graph","library"),)

@dataclass(frozen=True)
class DigPopulation:
    stage:str
    def __post_init__(self):
        if self.stage not in {'nuclei','replicated','proteolyzed','rna_processed','sheared','captured'}:raise ValueError('unsupported selected ini-seq population')
    def rows(self):
        sc=Scene()
        sc.strand('parent',[role('original chromosomal DNA templates; sequence/coordinates unassigned')],label='parental genomic DNA')
        if self.stage!='nuclei':sc.strand('nascent',[role('new DNA with incorporated DIG-11-dU class; positions/fractions and exact nucleotide chemistry unassembled')],label='nascent DNA; affinity tag is not an identifier')
        sc.strand('rna',[role('RNA-derived digestion products and possible residual RNA' if self.stage in {'rna_processed','sheared','captured'} else 'source RNA and other material; no nonsource removal asserted')],label='other physical nucleic acids')
        sc.strand('protein',[role('proteinase-generated peptides/carryover' if self.stage in {'proteolyzed','rna_processed','sheared','captured'} else 'source nuclear proteins and replication factors')],label='protein-associated material')
        return sc.rows()

def protocol_graph():
    g=ProtocolGraph(TITLE,authority=SOURCE,representation='published_instance')
    current=g.input('nuclei','source mimosine-arrested lateG1EJ30 nuclei; preparation defers to cited predecessor',K.MOLECULAR,state=DigPopulation('nuclei'),primary=True)
    extract=g.input('hela_extract','Actual soluble proliferating HeLa cell extract; source catalytic mixture/composition and standard cell-free conditions unreviewed, not a buffer-only condition',K.MOLECULAR)
    nucleotide=g.input('dig_dutp','source digoxigenin-11-dUTP plus standard exogenous nucleotides',K.MOLECULAR,state=molecule('source nucleotide pool',(role('DIG-11-dUTP and nucleotides; exact standard-assay concentrations unassigned'),)))
    labelled=g.product('labelled','source15min initiation/elongation population, no complete origin specificity asserted',K.MOLECULAR,state=DigPopulation('replicated'))
    g.react('initiate','Cell-free replication initiation and DIG-dUTP incorporation15min',S.MAGIC,((current,R.TEMPLATE),(extract,R.STRUCTURAL),(nucleotide,R.STRUCTURAL)),(labelled,),reason='Source cell-free reaction known; named predecessor standard composition/replication molecular geometry unreviewed.')
    current=operation(g,'stop','Source stop by ice-coldPBS',S.QUENCHING,labelled,labelled.state,'ice-coldPBS')
    current=cleanup(g,'nuclei_recovery','Source centrifugation recovery of nuclei',current,current.state)
    current=operation(g,'proteolyze','Source sarcosinate/proteinaseK2mg/ml55C18h',S.DEPROTEINIZATION,current,DigPopulation('proteolyzed'),'source proteinaseK inDNApurificationbuffer',reagent_kind=K.CATALYST,reagent_role=R.CATALYST)
    current=cleanup(g,'organic_first','Source phenol/chloroform extraction after proteinaseK',current,current.state)
    current=operation(g,'rnase','Source RNaseA200µg/ml37C1h; no physical discard at digestion',S.DIGESTION,current,DigPopulation('rna_processed'),'source RNaseA',reagent_kind=K.CATALYST,reagent_role=R.CATALYST)
    current=cleanup(g,'organic_second','Source subsequent phenol/chloroform extraction',current,current.state)
    current=operation(g,'shear','Source sonication100–1000bp',S.FRAGMENTATION,current,DigPopulation('sheared'),'source sonication conditions; exact sequence/end classes unassembled')
    current=cleanup(g,'debris_spin','Source centrifugation removes physical debris after sonication',current,current.state)
    ip=g.product('ip_input','DNA allocated to source antiDIGIP',K.MOLECULAR,state=current.state)
    control=g.product('input_control','source inputDNA for sequencing control',K.MOLECULAR,state=current.state)
    g.react('allocate_input','Represent source input and IP aliquots separately; amounts unassigned',S.SPLIT,((current,R.TEMPLATE),),(ip,control))
    current=operation(g,'capture','Source antiDIG-coupledProteinGDynabeads2hRT; uncoupled beads separate negativecontrol',S.BINDING,ip,DigPopulation('captured'),'source antiDIG-coupledProteinGDynabeads',reagent_kind=K.SELECTION,reagent_role=R.SELECTION)
    current=cleanup(g,'ip_wash','Source sixLiCl/Igepal/deoxycholate washes then twoTrisNaCl washes; transfer beads',current,current.state)
    current=operation(g,'elute','Source twice15min65C SDS-buffer elution; exact tag/antibody partition unquantified',S.GENERIC,current,current.state,'source50mMTrispH8/1mMEDTA/1%SDS')
    # Recovery includes the actual source aqueous extraction and precipitation, not assumed molecular disappearance.
    for ident,branch in (('ip',current),('input',control)):
        recovered=cleanup(g,ident+'_organic','Source phenol/chloroform extraction of '+ident+'DNA',branch,unassembled('source '+ident+' DNA; residual tag/protein/RNA recovery unquantified'))
        recovered=cleanup(g,ident+'_ethanol','Source ethanol precipitation of '+ident+'DNA',recovered,recovered.state)
        kit=g.input(ident+'_kit','Aggregate actual Illumina ChIP-seq kit reagent allocations; historical version, packaging/stock count and full operations unreviewed; no invented single premix',K.MOLECULAR)
        library=g.product('library' if ident=='ip' else 'input_library','source '+ident+' HiSeq2000 library roles',K.LIBRARY,state=unassembled('source '+ident+' ChIP-kit-derived library; exact ends/PCR/run/index geometry unassembled'))
        g.react(ident+'_library_boundary','Named historical kit boundary for '+ident+'DNA',S.MAGIC,((recovered,R.TEMPLATE),(kit,R.STRUCTURAL)),(library,),reason='Source names kit/readout; no actual constituent library operations inferred from modern kits or matching endpoint.')
        g.mark_terminal(Terminal(library.id,Disposition.UNRESOLVED,'Source HiSeq2000 endpoint; full molecular/read geometry remains unassembled.'))
    # Literal independently reviewed source inventory. The two branch labels are
    # the actual matched input/IP source scopes, not graph-discovered targets.
    g.declare_reviewed_scope(('nuclei','hela_extract','dig_dutp','initiate','stop_stock','stop','nuclei_recovery_condition','nuclei_recovery_step','proteolyze_stock','proteolyze','organic_first_condition','organic_first_step','rnase_stock','rnase','organic_second_condition','organic_second_step','shear_stock','shear','debris_spin_condition','debris_spin_step','allocate_input','capture_stock','capture','ip_wash_condition','ip_wash_step','elute_stock','elute'),g.authority)
    for k in ('ip','input'):
        g.declare_reviewed_scope(tuple(k+t for t in ('_organic_condition','_organic_step','_ethanol_condition','_ethanol_step','_kit','_library_boundary')),g.authority)
    g.add_source_item(SourceItem('pending_originals','Cell-free assay/nuclei preparation and historical IlluminaChIP-seq kit originals','Main Methods64/68 refs27/28 and Methods80 named kit; source-specific complete versions remain unreviewed.'))
    return g
