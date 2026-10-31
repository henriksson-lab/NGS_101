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
    c('T4PNK is an actual catalyst port',next(p.role for p in g.steps['thiophosphate'].inputs if p.material=='thiophosphate_stock'),R.CATALYST)
    c('Maleimide is an actual molecular reagent',next(p.role for p in g.steps['maleimide'].inputs if p.material=='maleimide_stock'),R.STRUCTURAL)
    c('Original historical assembly obligation stays uncovered','unreviewed_assembly' in g.source_coverage,False)
    c('Private preview has no finished marker', 'data-protocol-complete=' not in render_page())

    second=experiment_graph(2)
    c('RNase digestion retains molecular material until physicalcollection',g.steps['rnase'].outputs,('rnase_product',))
    c('ReleasedR derives from physically collected post-RNase mixture',g.steps['collect_r'].inputs[0].material,'rnase_product')
    c('Captured RNA tag remains upstream in residual fraction','biotinRNA' in str(g.materials['residual'].state.rows()),True)
    c('Experiment1 retains original NEBNext librarybranch','NEBNextChIPSeq' in g.steps['library_prep'].label,True)
    c('Experiment2 retains original IntegenX librarybranch','IntegenXPrepX' in second.steps['library_prep'].label,True)
    c('Original yield-onlyKlenow is not an NGSlibrary reaction',all('Klenow' not in s.id for s in g.steps.values()),True)
    c('Original unamplified control split precedesWGA',g.steps['wga'].inputs[0].material,'wga_input')
    c('ColumnCcollection followsalkalireaction',g.steps['collect_c'].inputs[0].material,'alkali_product')
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


def native(g,x,limitation):
    # Source acquisition may be stated without fabricated sequencing primers/cycle counts.
    g.sequence(x,SequencingPlan(limitation=limitation))

def pending(g,label):
    g.add_source_item(SourceItem('unreviewed_assembly',label,'Named source operations are represented; historical end/kit/primer assembly and complete original-version applicability unreviewed.'))

PAGE_FILE='nscr.html'
TITLE='NSCR'
NOTES='01_nscr.html'
SOURCE='Published Methods151–157/Fig1 and Supplemental Methods49–54; 10.1101/gr.176099.114.'
CAVEAT='Original experiment1/2 library branches remain distinct; WGAII chemistry and historical adapters/primers unassembled. RNA-DNA chimeras and brokenDNA background can both enter affinity capture. No ideal RNA-specific recovery or tag transfer to releasedDNA is asserted.'
READOUT='Original HiSeq2000 PE50 experiment1 versus PE72 experiment2. Klenow random-primed duplex conversion is a yield assay explicitly omitted for NGS. The later detailed protocol is unreviewed and not a substitute for this acquisition.'
PREVIEWS=(('Experiment1 preparation','library'),)
def protocol_graph():
    return experiment_graph(1)

def experiment_graph(experiment):
    if experiment not in (1,2):raise ValueError('original experiment1 or2 required')
    g=ProtocolGraph(TITLE,authority=SOURCE,representation='published_instance')
    x=g.input('cells','individual source WT or MCM2-deficient MEF preparation',K.MOLECULAR,state=pop('genomicDNA, nascent RNA-DNA chains and other cellular material'),primary=True)
    x=op(g,'lysis','Source SDS/EDTA/NaCl/Tris/PK37C18h lysis; molecular debris remains',S.MAGIC,x,pop('DNA/RNA plus proteolytic and other molecular remnants'),'source proteinaseK and lysis conditions',K.CATALYST,R.CATALYST)
    x=cleanup(g,'phenol','Source phenol/chloroform extraction',x,x.state)
    x=cleanup(g,'trizol','Source TRIzol DNA recovery and washes; no RNase destruction of RNA leaders',x,pop('recovered DNA with possible RNA-DNA chains/carryover'))
    x=op(g,'denature','Source heat denaturation; original numerical heat parameters unassigned',S.MAGIC,x,pop('denatured genomic fragments and RNA-DNA chimeras'),'source heat condition')
    x=cleanup(g,'gradient','Source16ml sucrose gradient; pooled400–2000nt fractions then ethanol precipitation',x,x.state)
    x=op(g,'thiophosphate','Source T4PNK-mediated5′ thiophosphorylation; RNA leaders and brokenDNA can be labelled',S.MAGIC,x,pop('5′thiophosphorylated RNA-DNA and brokenDNA classes plus unreacted species'),'Source T4PNK catalyst; thiophosphate donor/end and kit geometry remain unreviewed',K.CATALYST,R.CATALYST)
    maleimide=g.input('maleimide_stock','Actual source biotin-maleimide SP1501 molecular supply; dose and whole covalent adduct unassembled',K.MOLECULAR)
    coupling=g.input('maleimide_conditions','Source MB9001 coupling conditions remain under historical-kit review; no numerical defaults',K.CONDITION)
    tagged=g.product('maleimide_product','Source5primebiotin-linked RNA-DNA/brokenDNA classes and possible uncoupled molecules',K.MOLECULAR,state=pop('source labelled chain population; exact phosphate/linker/junction chemistry unassembled'))
    g.react('maleimide','Source biotin-maleimide coupling, not a buffer-only condition',S.CHEMICAL_CONVERSION,((x,R.TEMPLATE),(maleimide,R.STRUCTURAL),(coupling,R.CONDITION)),(tagged,))
    x=tagged
    beads=g.input('beads','source NEBS1420S streptavidin beads',K.SELECTION)
    bound=g.product('bound','captured tagged chains; possible background/carryover',K.MOLECULAR,state=x.state)
    flow=g.product('flowthrough','source collected F flowthrough fraction',K.MOLECULAR,state=pop('unbound molecules and carryover; composition unquantified'))
    g.react('capture','Source1hRT BufferA streptavidin capture',S.AFFINITY_CAPTURE,((x,R.TEMPLATE),(beads,R.SELECTION)),(bound,flow));g.mark_terminal(Terminal(flow.id,Disposition.STORED,'Source F comparison fraction is collected, not physically discarded.'))
    bound=cleanup(g,'wash','Source BufferA/B then4Murea70C5min,20xSSC twice70C10min,TE/RNB washes',bound,bound.state)
    condition=g.input('mock_condition','source200µl RNB15minRT without RNase',K.CONDITION)
    b=g.product('mock_retained','retained beads after mock wash',K.MOLECULAR,state=bound.state)
    w=g.product('mock_w','source collected W mock-wash supernatant',K.MOLECULAR,state=pop('mock released molecules/background; composition unquantified'))
    g.react('mock','Source15minRT mock wash then physical supernatant collection',S.CLEANUP,((bound,R.TEMPLATE),(condition,R.CONDITION)),(b,w));g.mark_terminal(Terminal(w.id,Disposition.STORED,'Source W comparison fraction.'))
    b=op(g,'rnase','Source RNaseI0.02U/µl in200µlRNB15minRT; no physical removal yet',S.MAGIC,b,pop('releasedDNA classes with possible RNA remnants; biotin-bearing cleaved RNA remains upstream; other bound DNA and debris persist'),'source RNaseI',K.CATALYST,R.CATALYST)
    sep=g.input('release_condition','source physical collection of RNase supernatant',K.CONDITION)
    released=g.product('released','source R supernatant',K.MOLECULAR,state=pop('RNase-releasedDNA and possible RNA/DNA/reagent carryover; no transfer of upstream biotin tag toDNA'))
    residual=g.product('residual','post-R retained bead fraction',K.MOLECULAR,state=pop('bead-associated biotinRNA fragments, directly taggedbrokenDNA and possible unreleased chains'))
    g.react('collect_r','Collect source R supernatant after RNase reaction',S.CLEANUP,((b,R.TEMPLATE),(sep,R.CONDITION)),(released,residual))
    residual=cleanup(g,'residual_tewash','Source residual-bead1mlTEwash before alkali',residual,residual.state)
    residual=op(g,'alkali','Source0.4NNaOH5minRT/0.4NHCl; exact alkaline end/remnant fates unresolved',S.MAGIC,residual,pop('alkali-treated column fraction and bead/carryover; no unchangedRNA or specific de-biotinylation assumption'),'source alkali/neutralization conditions')
    ccond=g.input('column_collection','Source post-neutralization supernatant collection',K.CONDITION)
    column=g.product('column_c','source collected C supernatant',K.MOLECULAR,state=pop('alkali-recovered column supernatant; exact DNA/RNA/remnant recovery and bead carryover unquantified'))
    column_beads=g.product('column_beads','retained column beads after C collection',K.MOLECULAR,state=pop('retained bead-associated material after alkali; recovery unquantified'))
    g.react('collect_c','Source physical C-supernatant recovery after NaOH/HCl',S.CLEANUP,((residual,R.TEMPLATE),(ccond,R.CONDITION)),(column,column_beads))
    g.mark_terminal(Terminal(column.id,Disposition.STORED,'Source collected C comparison fraction.'))
    g.mark_terminal(Terminal(column_beads.id,Disposition.UNRESOLVED,'Physical residual beads; downstream disposition unassigned.'))
    splitcond=g.input('aliquot_condition','Source save unamplified fractionR aliquot before WGA',K.CONDITION)
    control=g.product('unamplified_control','saved unamplified R control',K.MOLECULAR,state=released.state)
    wga_input=g.product('wga_input','remaining R selected for WGA',K.MOLECULAR,state=released.state)
    g.react('split_r','Source unamplified R control aliquot separation before amplification',S.GENERIC,((released,R.TEMPLATE),(splitcond,R.CONDITION)),(control,wga_input),changes_molecule=False,reason='Physical source aliquot division, exact fraction unassigned.')
    g.mark_terminal(Terminal(control.id,Disposition.STORED,'Source saved unamplified R control; no synthetic library route assigned.'))
    x=op(g,'wga','Source WGAII '+('SEQX experiment1' if experiment==1 else 'SEQXE experiment2')+'; NGS explicitly omits yield-assay Klenow conversion',S.MAGIC,wga_input,pop('source whole-genome-amplifiedDNA preparation plus remnants; internals unassembled'),'Aggregate actual original WGAII reagent allocations; packaging/stock count unassigned, no inferred Phi29 or PCR primer recipe',K.MOLECULAR,R.STRUCTURAL)
    x=op(g,'library_prep','Source '+('NEBNextChIPSeq experiment1' if experiment==1 else 'IntegenXPrepXDNAChIP experiment2')+' library preparation',S.MAGIC,x,pop('source Illumina library preparation; historical geometry and workup internals unreviewed'),'Aggregate historical library reagent allocations; manufacturing/packaging unreviewed',K.MOLECULAR,R.STRUCTURAL)
    final=g.product('library','source sequencing library handoff',K.LIBRARY,state=x.state)
    cond=g.input('handoff_condition','source library handoff; no invented additional purification',K.CONDITION)
    g.react('handoff','Original source library endpoint',S.GENERIC,((x,R.TEMPLATE),(cond,R.CONDITION)),(final,),changes_molecule=False,reason='Typed acquisition handoff only; not a new chemical operation.')
    native(g,final,'Original HiSeq2000 paired-end'+str(50 if experiment==1 else 72)+' bases; historical priming sites/whole library unreviewed, no primer geometry fabricated.')
    # Actual reviewed six-step/control/WGA scope; exact historical library
    # manufacturing and its supplier applicability remain unaccounted.
    g.declare_reviewed_scope(('cells','lysis_stock','lysis','phenol_condition','phenol_step','trizol_condition','trizol_step','denature_stock','denature','gradient_condition','gradient_step','thiophosphate_stock','thiophosphate','maleimide_stock','maleimide_conditions','maleimide','beads','capture','wash_condition','wash_step','mock_condition','mock','rnase_stock','rnase','release_condition','collect_r','residual_tewash_condition','residual_tewash_step','alkali_stock','alkali','column_collection','collect_c','aliquot_condition','split_r','wga_stock','wga','handoff_condition','handoff','library'),g.authority)
    pending(g,'Original WGAII/NEBNextChIPSeq/IntegenXPrepX documents and detailed cited NSCR method')
    return g
