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

    x=AdenineMarks('AACA',frozenset({0}),frozenset({0,3}))
    c('Native and experimental methylation sharechemicalmark',x.physical,frozenset({0,3}))
    c('Methylation preserves DNAalphabet',x.sequence,'AACA')
    c('Exonuclease reaction retains debris before cleanup',g.steps['exonuclease'].outputs,('exonuclease_product',))
    c('First physical cleanup consumes post-digestion mixture',g.steps['cleanup_one_step'].inputs[0].material,'exonuclease_product')
    c('Secondphysicalcleanup retained',g.steps['library_step'].inputs[0].material,'cleanup_one')
    c('Solubilization does not invent dialysis hardware','membrane/exchange operation invented' in g.steps['solubilize'].label,True)
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

PAGE_FILE='samosa.html'
TITLE='SAMOSA'
NOTES='01_samosa.html'
SOURCE='Published eLife.59404 VORv2 (14December2020), Methods260–270; in-vitro array branch244–252 kept separate.'
CAVEAT='Selected original solubleK562 oligonucleosome route. EcoGII physical adenine marking and PacBioIPD/HMM computational calls remain distinct. No complete methylation or perfect nucleosome specificity. Historical SMRTbell assembly/motor/primer geometry unreviewed.'
READOUT='Source SequelII,three30h movies with2h first pre-extension and0.7h later pre-extension. These times are not nucleotide cycle counts. Native m6dA retained; single-moleculeIPD/NN/HMM nucleosome classification is computational.'
PREVIEWS=(('K562 soluble-oligonucleosome route','library'),)
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

def protocol_graph():
    g=ProtocolGraph(TITLE,authority=SOURCE,representation='published_instance')
    x=g.input('cells','source K562 cell preparation',K.MOLECULAR,state=pop('chromatin, RNA and cellular components; methylation positions unassigned'),primary=True)
    x=cleanup(g,'cells_wash','Source300g5min/PBS wash',x,x.state)
    x=op(g,'nuclei','Source HEPES/KCl/MgCl2/Triton/glycerol nuclei isolation,5minice',S.MAGIC,x,x.state,'source nuclei-isolationbuffer')
    x=cleanup(g,'nuclei_wash','Source600g5min4C nuclei recovery/BufferM wash',x,x.state)
    x=op(g,'mnase','Source1UMNase/50millionnuclei,1min37C with1mMCa',S.MAGIC,x,pop('MNase-cut chromatin fragments, original otherDNA and cleavage remnants; exact end geometry unreviewed'),'source S.aureusMNase andCa; variantMNase conditions separate')
    x=op(g,'quench','Source2mMEGTA immediately/onice; no fragment disposal',S.QUENCHING,x,x.state,'source EGTA')
    x=cleanup(g,'nuclei_recovery','Source600g5min4C recovery after MNase',x,x.state)
    x=op(g,'solubilize','Source overnight4C Tep20/300µg/mllysolecithin; heading saysdialysis but no membrane/exchange operation invented',S.MAGIC,x,x.state,'source Tep20/lysolecithin incubation')
    condition=g.input('soluble_condition','source12000g5min4C',K.CONDITION)
    soluble=g.product('soluble','source collected soluble chromatin supernatant',K.MOLECULAR,state=x.state)
    pellet=g.product('insoluble','source separated insoluble pellet',K.MOLECULAR,state=pop('insoluble material; molecular partition unquantified'))
    g.react('soluble_separation','Source soluble chromatin supernatant collection',S.CLEANUP,((x,R.TEMPLATE),(condition,R.CONDITION)),(soluble,pellet));g.mark_terminal(Terminal(pellet.id,Disposition.UNRESOLVED,'Physical separated pellet; downstream disposition not assigned.'))
    x=op(g,'mark','Source1.5µgchromatin/200µl,1xCutSmart/1mMSAM/2.5µlEcoGIIstock2.5e4U/ml,37C30min;SAM6.25mMat15min',S.MAGIC,soluble,pop('EcoGII-markedm6dA chromatin and unmarked molecules; original native marks possible; no guaranteed saturation'),'source EcoGII/SAM methylation reaction')
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
    native(g,x,'Source SequelII native SMRT acquisition,three30h movies,2h then0.7h pre-extension; historical sequencing-primer/motor/closed-library geometry unreviewed, no cycle count invented.')
    pending(g,'Original historical PacBio SMRTbell/repair/exonuclease/SequelII documents and precise adaptor polymers')
    return g
