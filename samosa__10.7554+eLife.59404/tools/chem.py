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
    out = g.product(ident, 'recovered source preparation', K.LIBRARY if template.kind is K.LIBRARY else K.MOLECULAR, state=state)
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

    x=AdenineMarks('AACA',frozenset({0}),frozenset({0,3}))
    c('Native and experimental methylation sharechemicalmark',x.physical,frozenset({0,3}))
    c('Methylation preserves DNAalphabet',x.sequence,'AACA')
    c('Exonuclease reaction retains debris before cleanup',g.steps['exonuclease'].outputs,('exonuclease_product',))
    c('First physical cleanup consumes post-digestion mixture',g.steps['cleanup_one_step'].inputs[0].material,'exonuclease_product')
    c('Actual two MNase templates are allocated before separate processing',tuple(g.steps['mnase_allocations'].outputs),('mnase_first','mnase_second'))
    c('Second processed material derives from second aliquot',g.steps['other_aliquot_processing'].inputs[0].material,'mnase_second')
    c('Naked positive derives from purified nonmethylated control',g.steps['positive_allocation'].inputs[0].material,'nonmeth_purified')
    c('Positive is repurified after its own marking',g.steps['positive_purification'].inputs[0].material,'positive_marked')
    c('SAM initial and replenishing contributions are explicit',tuple(p.material for p in g.steps['mark'].inputs[-3:]),('sam_initial','sam_replenish','cutsmart'))
    c('Final QA allocations do not become production library',tuple(g.terminals[k].disposition for k in ('qubit_1ul','bioanalyzer_1ul')),(Disposition.QC,Disposition.QC))
    c('Secondphysicalcleanup retained',g.steps['library_step'].inputs[0].material,'cleanup_one')
    c('Solubilization does not invent dialysis hardware','membrane/exchange operation invented' in g.steps['solubilize'].label,True)
    c.report()
def pop(label):
    return mixture(('source DNA',(role(label),)),('other physical material',(role('other original DNA/RNA/protein or reaction remnants; exact partition and recovery unquantified'),)))
def op(g,key,label,kind,x,out_state,reagent_label,reagent_kind=K.CONDITION,reagent_role=R.CONDITION,additional=()):
    reagent=g.input(key+'_stock',reagent_label,reagent_kind)
    product=g.product(key+'_product',label,K.MOLECULAR,state=out_state)
    kw={'reason':'Source operation known; exact molecular substrate/product relation and mixture fractions unassembled.'} if kind in {S.MAGIC,S.GENERIC} else {}
    if kind is S.GENERIC:kw['changes_molecule']=False
    g.react(key,label,kind,((x,R.TEMPLATE),(reagent,reagent_role))+tuple(additional),(product,),**kw)
    return product


def native(g,x,limitation):
    # Source acquisition may be stated without fabricated sequencing primers/cycle counts.
    g.sequence(x,SequencingPlan(limitation=limitation))

def pending(g,label):
    g.add_source_item(SourceItem('unreviewed_assembly',label,'Named source operations are represented; historical end/kit/primer assembly and complete original-version applicability unreviewed.'))

PAGE_FILE='samosa.html'
TITLE='SAMOSA'
NOTES='01_samosa.html'
SOURCE='Published eLife.59404 VORv2 (exact manifest-versioned original), Methods260–270; in-vitro array branch244–252 kept separate.'
CAVEAT='Selected original solubleK562 oligonucleosome route. EcoGII physical adenine marking and PacBioIPD/HMM computational calls remain distinct. No complete methylation or perfect nucleosome specificity. Historical SMRTbell assembly/motor/primer geometry unreviewed. INFERRED allocation convention: the matched nonmethylation control is drawn from a representative prepared soluble-chromatin stock; exact experimental stock-sharing is unassigned.'
READOUT='Source SequelII,three30h movies with2h first pre-extension and0.7h later pre-extension. These times are not nucleotide cycle counts. Native m6dA retained; single-moleculeIPD/NN/HMM nucleosome classification is computational.'
PREVIEWS=(('K562 soluble-oligonucleosome route','production_library'),)
@dataclass(frozen=True)
class AdenineMarks:
    sequence:str
    native:frozenset[int]=frozenset()
    installed:frozenset[int]=frozenset()
    def __post_init__(self):
        if not self.sequence or set(self.sequence)-set('ACGT'):raise ValueError('explicitDNA required')
        if any(not isinstance(p,frozenset) or any(type(i)is not int or not 0<=i<len(self.sequence) or self.sequence[i]!='A' for i in p) for p in (self.native,self.installed)):raise ValueError('m6dA must occupy actual adenines')
    @property
    def physical(self):return self.native|self.installed

def control_purification(g, prefix, source):
    rnase=g.input(prefix+'_rnase','Actual source RNaseA control-preparation allocation;10ul/RT10min recipe',K.CATALYST)
    pk=g.input(prefix+'_pk','Actual source control ProteinaseK20ul20mg/ml',K.CATALYST)
    sds=g.input(prefix+'_sds','Actual source control20ul10%SDS',K.MOLECULAR)
    supplies=g.input(prefix+'_workup','Actual phenol/chloroform,NaOAc,GlycoBlue,ethanol and EB workup contributions; packaging/stockcount unassembled',K.MOLECULAR)
    purified=g.product(prefix+'_purified','Source purified control DNA in50ulEB; recovery and native marks unassigned',K.MOLECULAR,state=pop('source control DNA; molecular distributions/recovery unassembled'))
    waste=g.product(prefix+'_waste','Actual control organic/precipitation separated fractions; partition unquantified',K.WASTE)
    g.react(prefix+'_purification','Source control RNase→PK/SDS→organic/ethanol recovery; ordered boundary, molecular substeps unassembled',S.MAGIC,((source,R.TEMPLATE),(rnase,R.CATALYST),(pk,R.CATALYST),(sds,R.STRUCTURAL),(supplies,R.STRUCTURAL)),(purified,waste),reason='Actual source266 shared control purification recipe; no complete chemical/workup audit claimed.')
    g.discard(waste,'Source physical workup fractions, not disappearance at RNase or proteinase alone.')
    return purified


def control_preparations(g, chromatin):
    buffer=g.input('nonmeth_buffer','Actual source nonmethylation-control CutSmart/reaction-buffer contribution, including initialSAM; no EcoGII or replenishingSAM',K.MOLECULAR)
    nonmeth=g.product('nonmeth_reaction','Source nonmethylated-control reaction; native DNA marks unassigned, no EcoGII installation asserted',K.MOLECULAR)
    g.react('nonmeth_incubation','Source control266: same reaction environment minusEcoGII and SAMreplenishment',S.MAGIC,((chromatin,R.TEMPLATE),(buffer,R.STRUCTURAL)),(nonmeth,),reason='Source-positive control incubation; exact buffering/stock schedule unassembled.')
    purified=control_purification(g,'nonmeth',nonmeth)
    positive=g.product('positive_25ul','Actual25ul approximately500ng of post-purified NONmethylated control for naked-DNA positive preparation',K.MOLECULAR)
    remaining=g.product('nonmeth_stock','Nonmethylated-control material remaining after25ul selection; QA/575ng library allocation and full producer unassembled',K.MOLECULAR)
    g.react('positive_allocation','Source naked-DNA positive derives from purified nonmethylated control',S.SPLIT,((purified,R.TEMPLATE),),(positive,remaining))
    g.mark_terminal(Terminal(remaining.id,Disposition.UNRESOLVED,'Source nonmethylated-control SMRTbell continuation/quantitative allocations remain unassembled.'))
    enzyme=g.input('positive_ecogii','Actual source EcoGII contribution for naked-DNA positive, as previously stated; no independent naked-DNA input invented',K.CATALYST)
    sam=g.input('positive_sam','Actual positive-control SAM/buffer supplies per referred recipe; exact allocation packaging unassembled',K.MOLECULAR)
    marked=g.product('positive_marked','Source positive EcoGII-marked naked DNA and unreacted/reagent mixture',K.MOLECULAR)
    g.react('positive_marking','Source post-purification naked-DNA EcoGII positive reaction',S.MAGIC,((positive,R.TEMPLATE),(enzyme,R.CATALYST),(sam,R.STRUCTURAL)),(marked,),reason='Source266 says methylated as previously stated; full occupied DNA geometry unassembled.')
    repurified=control_purification(g,'positive',marked)
    g.mark_terminal(Terminal(repurified.id,Disposition.UNRESOLVED,'Source positive approximately320ng SMRTbell continuation/QA/stock allocations unassembled; independent control preparation retained.'))

def protocol_graph():
    g=ProtocolGraph(TITLE,authority=SOURCE,representation='published_instance')
    x=g.input('cells','source K562 cell preparation',K.MOLECULAR,state=pop('chromatin, RNA and cellular components; methylation positions unassigned'),primary=True)
    x=cleanup(g,'cells_wash','Source300g5min/PBS wash',x,x.state)
    x=op(g,'nuclei','Source HEPES/KCl/MgCl2/Triton/glycerol nuclei isolation,5minice',S.MAGIC,x,x.state,'source nuclei-isolationbuffer')
    x=cleanup(g,'nuclei_wash','Source600g5min4C nuclei recovery/BufferM wash',x,x.state)
    first=g.product('mnase_first','Actual first1ml nuclei aliquot for source MNase reaction',K.MOLECULAR,state=x.state)
    second=g.product('mnase_second','Actual second1ml nuclei aliquot for independently prepared source MNase reaction',K.MOLECULAR,state=x.state)
    g.react('mnase_allocations','Source262 nuclei distributed into two1ml aliquots BEFORE digestion',S.SPLIT,((x,R.TEMPLATE),),(first,second))
    other_reagents=g.input('other_mnase_supplies','Actual second-aliquot MNase1U/50millionnuclei,Ca,EGTA,Tep20/lysolecithin workflow supplies; stock packaging unassembled',K.MOLECULAR)
    other_soluble=g.product('other_soluble','Actual independently processed second-aliquot soluble chromatin; exact molecular recovery/nextallocation unassembled',K.MOLECULAR,state=pop('independent second-aliquot soluble oligonucleosomes'))
    other_pellet=g.product('other_insoluble','Actual second-aliquot insoluble fraction; molecular partition unassigned',K.MOLECULAR)
    g.react('other_aliquot_processing','Source second aliquot MNase/quench/recovery/overnight solubilization/centrifugation; internal chemistry unassembled',S.MAGIC,((second,R.TEMPLATE),(other_reagents,R.STRUCTURAL)),(other_soluble,other_pellet),reason='Source262 actual distinct pre-reaction aliquot; no cloning a reacted product or source-unsupported recombination.')
    g.mark_terminal(Terminal(other_soluble.id,Disposition.UNRESOLVED,'Source does not locate recombination; second soluble preparation retained without invented pooling.'))
    g.mark_terminal(Terminal(other_pellet.id,Disposition.UNRESOLVED,'Actual separated insoluble material; later disposition unassigned.'))
    x=first
    x=op(g,'mnase','Source1UMNase/50millionnuclei,1min37C with1mMCa',S.MAGIC,x,pop('MNase-cut chromatin fragments, original otherDNA and cleavage remnants; exact end geometry unreviewed'),'source S.aureusMNase andCa; variantMNase conditions separate')
    x=op(g,'quench','Source2mMEGTA immediately/onice; no fragment disposal',S.QUENCHING,x,x.state,'source EGTA')
    x=cleanup(g,'nuclei_recovery','Source600g5min4C recovery after MNase',x,x.state)
    x=op(g,'solubilize','Source overnight4C Tep20/300µg/mllysolecithin; heading saysdialysis but no membrane/exchange operation invented',S.MAGIC,x,x.state,'source Tep20/lysolecithin incubation')
    condition=g.input('soluble_condition','source12000g5min4C',K.CONDITION)
    soluble=g.product('soluble','source collected soluble chromatin supernatant',K.MOLECULAR,state=x.state)
    pellet=g.product('insoluble','source separated insoluble pellet',K.MOLECULAR,state=pop('insoluble material; molecular partition unquantified'))
    g.react('soluble_separation','Source soluble chromatin supernatant collection',S.CLEANUP,((x,R.TEMPLATE),(condition,R.CONDITION)),(soluble,pellet));g.mark_terminal(Terminal(pellet.id,Disposition.UNRESOLVED,'Physical separated pellet; downstream disposition not assigned.'))
    methylation_input=g.product('methylation_input','Actual1.5ug soluble chromatin for selected200ul methylation reaction',K.MOLECULAR,state=soluble.state)
    control_input=g.product('nonmeth_input','INFERRED matched soluble-chromatin control allocation; actual stock-sharing/mass unassembled',K.MOLECULAR,state=soluble.state,inferred=True)
    stock=g.product('soluble_stock','Untransferred soluble chromatin stock; amount/nextuse unassigned',K.MOLECULAR,state=soluble.state)
    g.react('mark_allocations','INFERRED — allocation convention: source266 methylation/control contributions; exact stock-sharing unassigned',S.SPLIT,((soluble,R.TEMPLATE),),(methylation_input,control_input,stock))
    g.mark_terminal(Terminal(stock.id,Disposition.UNRESOLVED,'No source-unsupported allocation mass or fate assigned.'))
    control_preparations(g,control_input)
    sam_initial=g.input('sam_initial','Actual initial1mM SAM molecular contribution',K.MOLECULAR)
    sam_replenish=g.input('sam_replenish','Actual SAM replenishment to6.25mM at15min; separate physical contribution',K.MOLECULAR)
    cutsmart=g.input('cutsmart','Actual1xCutSmart reaction-buffer contribution',K.CONDITION)
    x=op(g,'mark','Source1.5µgchromatin/200µl,1xCutSmart/1mMSAM/2.5µlEcoGIIstock2.5e4U/ml,37C30min;SAM6.25mMat15min',S.MAGIC,methylation_input,pop('EcoGII-markedm6dA chromatin and unmarked molecules; original native marks possible; no guaranteed saturation'),'Actual2.5ulEcoGII stock2.5e4U/ml',K.CATALYST,R.CATALYST,additional=((sam_initial,R.STRUCTURAL),(sam_replenish,R.STRUCTURAL),(cutsmart,R.CONDITION)))
    x=op(g,'rnase','Source10µlRNaseART10min; RNAcleavage remnants remain',S.MAGIC,x,pop('markedDNA/protein, RNAcleavage fragments and reagents'),'source RNaseA',K.CATALYST,R.CATALYST)
    x=op(g,'proteolyze','Source20µlPK20mg/ml+20µl10%SDS65C at least2h/overnight; peptide remnants remain',S.DEPROTEINIZATION,x,pop('markedDNA, RNA fragments/peptides and reagents'),'source proteinaseK',K.CATALYST,R.CATALYST)
    x=cleanup(g,'organic','Source phenol:chloroform extraction/aqueous recovery',x,x.state)
    x=cleanup(g,'ethanol','Source NaOAc/GlycoBlue/3volEtOHovernight−20C;spin/70%EtOHwash/EB',x,x.state)
    x=op(g,'shear','Selected original G-tube3380g1min approximately10kb; variantMNase experiments omit',S.FRAGMENTATION,x,x.state,'source G-tube; no in-vitroarray0.8x substitution')
    x=cleanup(g,'concentrate','Source0.45xAMPurePB concentration',x,x.state)
    x=op(g,'repair','Source historical PacBio DNA damage repair and end repair',S.MAGIC,x,x.state,'source manufacturer repair modules; exact polymer/end changes unreviewed')
    ds=mixture(('sourceSMRTbelladaptor',(role('unique source SMRTbellhairpin adapters,100µMstock/20µMannealed; exact polymers and identifier meaning unassigned'),)))
    donor=g.input('adaptor','Source85C5min/RT30s/4Chold annealedSMRTbell adaptors',K.ADAPTER,state=ds)
    enzyme=g.input('ligase','Source manufacturer SMRTbell ligation reagents',K.CATALYST)
    ligated=g.product('ligated','source SMRTbell-ligatedDNA preparation',K.LIBRARY,state=pop('source covalent SMRTbell-linked markedDNA molecules plus unligatedDNA/adaptor/reagents; exact assembly unreviewed'))
    g.react('ligate','Source SMRTbell ligation',S.LIGATION,((x,R.TEMPLATE),(donor,R.STRUCTURAL),(enzyme,R.CATALYST)),(ligated,),geometry=AdapterLigationGeometry('unavailable','unavailable',(ds,),'Historical SMRTbell end/hairpin assembly unreviewed; no source-absence assertion.'))
    x=op(g,'exonuclease','Source exonuclease cleanup reaction; digestedDNA remains until bead workup',S.MAGIC,ligated,pop('SMRTbell library plus exonucleolyticDNA remnants and other carryover; exact susceptibility/recovery unresolved'),'source historical exonuclease module')
    x=cleanup(g,'cleanup_one','First source0.45xAMPurePB cleanup',x,x.state)
    x=cleanup(g,'library','Second source0.45xAMPurePB cleanup;Qubit/Bioanalyzer then acquisition',x,x.state)
    qubit=g.product('qubit_1ul','Actual1ul final-library Qubit analytical allocation',K.QC)
    bioanalyzer=g.product('bioanalyzer_1ul','Actual1ul final-library Bioanalyzer analytical allocation',K.QC)
    production=g.product('production_library','Remaining source SMRTbell library for acquisition; volume/loaded fraction unassembled',K.LIBRARY,state=x.state)
    g.react('final_qa_allocations','Source270 final-library1ulQubit and1ulBioanalyzer allocations',S.SPLIT,((x,R.TEMPLATE),),(qubit,bioanalyzer,production))
    g.mark_terminal(Terminal(qubit.id,Disposition.QC,'Source Qubit analytical material; assay internals unassembled.'))
    g.mark_terminal(Terminal(bioanalyzer.id,Disposition.QC,'Source Bioanalyzer analytical material; assay internals unassembled.'))
    x=production
    native(g,x,'Source SequelII native SMRT acquisition,three30h movies,2h then0.7h pre-extension; historical sequencing-primer/motor/closed-library geometry unreviewed, no cycle count invented.')
    pending(g,'Original historical PacBio SMRTbell/repair/exonuclease/SequelII documents and precise adaptor polymers')
    return g
