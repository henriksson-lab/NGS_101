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
    c.report()

PAGE_FILE='ap-seq.html'
TITLE='AP-seq — ARP tagging before shearing and recombined PreCR fractions'
NOTES='01_ap-seq.html'
SOURCE='<a href="https://doi.org/10.1186/s13059-018-1582-2">Poetsch et al. (2018)</a>; published PMC6284305 Methods171–183; publisherFigureS1 and synthetic-controlTableS1.'
CAVEAT='DirectAP and OGG1-AP are separate source preparations. ARP labels both intact AP aldehydes and AP beta-elimination intermediates; no single exact lesion/end structure or quantitative conversion is assumed. Tagging precedes fragmentation and the matched input is taken after tagging. PreCR-treated soluble DNA is purified and recombined with bead-associated material; neither fraction is silently discarded. Historical PreCR/KAPA molecular steps and full library/run/index assembly are unresolved. Synthetic qPCR oligos are not genomic library primers.'
READOUT='Source KAPA KK8504 preparation and HiSeq2000 paired-end125bp rapid/high-output runs. Regional enrichment is relative to matched input, not a single-base lesion-encoded sequencing letter. OGG1-AP includes nativeAP signal and is not modeled as perfect subtraction or complete8-oxoG conversion.'
PREVIEWS=(('DirectAP enriched preparation','direct_library'),('DirectAP matched tagged input','direct_input_library'),('OGG1-AP enriched preparation','ogg_library'),('OGG1-AP matched tagged input','ogg_input_library'),('Prepare affinity reagent','washed_beads'))
@dataclass(frozen=True)
class DamageHistory:
    stage:str
    ogg1_exposed:bool=False
    def __post_init__(self):
        if self.stage not in {'input','ogg1_products','tagged','fragmented','captured','precr_unknown','library'}:
            raise ValueError('unknown AP-seq preparation stage')
    def rows(self):
        sc=Scene()
        if self.stage=='library':
            sc.strand('library',[role('source unassembled KAPA end'),role('genomic sequence from selected/input preparation'),role('opposite original kit end')],label='source library roles')
            sc.strand('carryover',[role('possible background/template-derived material; exact preparation purity unassigned')],label='other possible physical classes')
        elif self.stage=='precr_unknown':
            sc.strand('repaired',[role('source PreCR-treated genomic DNA; adducts, ends and tag-retention states unassembled')],label='repair-treated DNA population')
            sc.strand('other',[role('possible repair/ARP-associated products; exact chemistry and fractionation unassigned')],label='physical reaction products, not assumed destroyed')
        else:
            tagged=self.stage in {'tagged','fragmented','captured'}
            sc.strand('ap',[role('intact AP-site DNA'+(';covalent ARP-biotin conjugate' if tagged else ';aldehyde-bearing class'))],label='intact-AP class')
            sc.strand('be',[role('AP beta-elimination intermediate DNA'+(';covalent ARP-biotin conjugate' if tagged else ';reactive carbonyl class'))],label='AP beta-elimination class')
            sc.strand('other',[role('other genomic DNA; complete lesion/selection distribution unassigned')],label='other DNA remains a physical population')
            if self.stage=='ogg1_products':sc.strand('products',[role('independent OGG1 reaction products; retained until source purification')],label='glycosylase reaction products')
        if self.ogg1_exposed:sc.note('library' if self.stage=='library' else 'repaired' if self.stage=='precr_unknown' else 'ap','OGG1 pretreatment recorded; residual8-oxoG/AP/other end-product fractions unassigned. NativeAP is not subtracted away.')
        return sc.rows()

def protocol_graph():
    g=ProtocolGraph(TITLE,authority='Published PMC6284305 Methods171–183,source directAP andOGG1-AP branches; publisherFigureS1 chemical classes')
    raw=g.input('raw_beads','source MyOne Dynabeads before washing',K.MOLECULAR)
    buffer=g.input('bead_buffer','Three1MNaCl/TE washes;resuspend2MNaCl/TE',K.CONDITION)
    beads=g.product('washed_beads','physically prewashed source beads',K.SELECTION)
    bead_waste=g.product('bead_waste','removed bead prewash fractions',K.WASTE)
    g.react('wash_beads','Prepare streptavidin reagent before genomicDNA addition',S.CLEANUP,((raw,R.TEMPLATE),(buffer,R.CONDITION)),(beads,bead_waste));g.discard(bead_waste,'Source prewashing of affinity reagent; no genomicDNA present yet.')
    arp=g.input('arp','5mMARP/PBS37C2h;source direct preparation5.7µgDNA; exact branch scaling not generalized',K.CONDITION)
    kit=g.input('kapa','source KAPA KK8504 library preparation; historical oligos and exact recipe unassembled',K.CONDITION)
    reviewed=['raw_beads','bead_buffer','wash_beads','arp','kapa']
    for prefix,exposed in [('direct',False),('ogg',True)]:
        dna=g.input(prefix+'_dna','independently prepared source genomicDNA for '+('OGG1-AP' if exposed else 'directAP'),K.MOLECULAR,state=DamageHistory('input'),primary=True)
        current=dna
        if exposed:
            ogg=g.input('ogg_enzyme','source recombinantOGG1,0.1µg/µgDNA,NEBbuffer2/BSA37C1h',K.CATALYST)
            treated=g.product('ogg_treated','OGG1-treated genomicDNA and reaction products',K.MOLECULAR,state=DamageHistory('ogg1_products',True))
            g.react('ogg_digest','Source OGG1 pretreatment before ARP labeling',S.DIGESTION,((dna,R.TEMPLATE),(ogg,R.CATALYST)),(treated,))
            current=cleanup(g,'ogg_clean','Source AMPure cleanup after OGG1 treatment',treated,DamageHistory('input',True))
            reviewed.extend(('ogg_enzyme','ogg_digest','ogg_clean_condition','ogg_clean_step'))
        tagged=g.product(prefix+'_tagged','source ARP-biotin conjugated genomicDNA classes',K.MOLECULAR,state=DamageHistory('tagged',exposed))
        g.react(prefix+'_tag','ARP tags AP aldehydes and AP beta-elimination carbonyls BEFORE shearing',S.CHEMICAL_CONVERSION,((current,R.TEMPLATE),(arp,R.CONDITION)),(tagged,))
        clean=cleanup(g,prefix+'_tag_clean','1.8xAMPure/two70%ethanol washes; do not dry beads',tagged,tagged.state)
        covaris=g.input(prefix+'_covaris','source Covaris in130µL,mean300bp',K.CONDITION)
        fragmented=g.product(prefix+'_fragmented','source tagged DNA fragments',K.MOLECULAR,state=DamageHistory('fragmented',exposed))
        g.react(prefix+'_shear','Fragment only after tagging and cleanup',S.FRAGMENTATION,((clean,R.TEMPLATE),(covaris,R.CONDITION)),(fragmented,))
        allocation=g.input(prefix+'_allocation','Reserve30µL/130µL matched input; remainingDNA for capture',K.CONDITION)
        input_dna=g.product(prefix+'_input','matched tagged/fragmented inputDNA',K.MOLECULAR,state=fragmented.state)
        enrich=g.product(prefix+'_enrich','remaining taggedDNA for enrichment',K.MOLECULAR,state=fragmented.state)
        g.react(prefix+'_split','Allocate source matched tagged input and capture preparation',S.SPLIT,((fragmented,R.TEMPLATE),(allocation,R.CONDITION)),(input_dna,enrich))
        capture_cond=g.input(prefix+'_capture_condition','source10h rotation room temperature',K.CONDITION)
        captured=g.product(prefix+'_captured','source bead-associated enriched DNA population',K.MOLECULAR,state=DamageHistory('captured',exposed))
        flow=g.product(prefix+'_flow','physically removed capture supernatant',K.WASTE)
        g.react(prefix+'_capture','Capture ARP-biotin DNA with prewashed beads',S.AFFINITY_CAPTURE,((enrich,R.TEMPLATE),(beads,R.SELECTION),(capture_cond,R.CONDITION)),(captured,flow));g.discard(flow,'Source capture separation; no complete purity/specificity asserted.')
        selected=cleanup(g,prefix+'_capture_wash','Three1MNaCl/TE washes;resuspend50µLTE',captured,captured.state)
        precr=g.input(prefix+'_precr','source PreCR M0309L; exact original repair/adduct/end chemistry unreviewed',K.CONDITION)
        soluble=g.product(prefix+'_repair_soluble','soluble PreCR-treated pulldown DNA fraction',K.MOLECULAR,state=DamageHistory('precr_unknown',exposed))
        bound=g.product(prefix+'_repair_bound','bead-associated PreCR-treated material retained',K.MOLECULAR,state=DamageHistory('precr_unknown',exposed))
        g.react(prefix+'_repair','Source PreCR treatment; retain soluble AND bead-associated material',S.MAGIC,((selected,R.TEMPLATE),(precr,R.CONDITION)),(soluble,bound),reason='Source physical handling retains both supernatant and beads. Exact PreCR/ARP tag/end mechanisms and quantitative partition unassembled; neither DNA fraction is discarded.')
        purified=cleanup(g,prefix+'_soluble_clean','AMPure purification of pulldown supernatant',soluble,soluble.state)
        recombined=g.product(prefix+'_recombined','purified soluble pulldownDNA recombined with source beads',K.MOLECULAR,state=DamageHistory('precr_unknown',exposed))
        g.react(prefix+'_recombine','Recombine recovered supernatantDNA with bead-associated material',S.POOL,((purified,R.TEMPLATE),(bound,R.TEMPLATE)),(recombined,))
        repaired_input=g.product(prefix+'_input_repaired','source PreCR-treated matched input',K.MOLECULAR,state=DamageHistory('precr_unknown',exposed))
        g.react(prefix+'_repair_input','PreCR also treats the matched tagged input',S.MAGIC,((input_dna,R.TEMPLATE),(precr,R.CONDITION)),(repaired_input,),reason='Source names PreCR on input and pulldown; exact original repair chemistry unassembled. No lesion-specific base replacement inferred.')
        clean_input=cleanup(g,prefix+'_input_clean','AMPure purification of PreCR-treated input',repaired_input,repaired_input.state)
        for suffix,x in [('',recombined),('_input',clean_input)]:
            library=g.product(prefix+suffix+'_library','source '+prefix+suffix+' KAPA library roles',K.LIBRARY,state=DamageHistory('library',exposed))
            g.react(prefix+suffix+'_kit','Historical KAPA library preparation boundary',S.MAGIC,((x,R.TEMPLATE),(kit,R.CONDITION)),(library,),reason='Source kit andHiSeqPE125 named; full molecular steps, adapters/indexes/PCR and run-primer geometry not reconstructed from another kit/protocol.')
            g.mark_terminal(Terminal(library.id,Disposition.UNRESOLVED,'Source HiSeq2000 PE125; full KAPA product and physical run/index geometry unassembled.'))
        reviewed.extend(prefix+suffix for suffix in ('_dna','_tag','_tag_clean_condition','_tag_clean_step','_covaris','_shear','_allocation','_split','_capture_condition','_capture','_capture_wash_condition','_capture_wash_step','_precr','_repair','_soluble_clean_condition','_soluble_clean_step','_recombine','_repair_input','_input_clean_condition','_input_clean_step','_kit','_input_kit'))
    g.declare_reviewed_scope(reviewed,g.authority)
    g.add_source_item(SourceItem('historical_reagent_and_library_assembly','Historical PreCR/KAPA original reaction details and full library assembly','Published Methods181–183 names PreCR M0309L andKAPA KK8504; complete historical source/version and molecule/run geometry unreviewed.'))
    return g

def selftest():
    from checks import Check
    c=Check();g=protocol_graph()
    for prefix in ('direct','ogg'):
        producer=g.producer_of(prefix+'_fragmented')
        c(prefix+' input allocated only after tagging/shearing',g.steps[prefix+'_split'].inputs[0].material,prefix+'_fragmented')
        c(prefix+' repaired soluble DNA is recombined with retained bead fraction',tuple(p.material for p in g.steps[prefix+'_recombine'].inputs),(prefix+'_soluble_clean',prefix+'_repair_bound'))
    c('Incomplete chemistry/library geometry blocks finished publication',bool(g.publication_problems()))
    c('Private preview has no finished marker','data-protocol-complete=' not in render_page())
    c.report()
