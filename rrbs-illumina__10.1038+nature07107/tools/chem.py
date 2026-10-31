"""Source-owned bounded model; no inherited historical library recipe."""
from dataclasses import dataclass
from pathlib import Path
import sys
HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE), str(HERE.parents[1] / 'lib')]
from chemdraw import Construct, Scene, Segment, panel, oligo, revcomp
from protocol_graph import ProtocolGraph, InputRole as R, MaterialKind as K, StepKind as S, SequencingPlan, SourceItem, Terminal, Disposition
from protocol_model import role, molecule
from adapter_ligation import AdapterLigationGeometry, reported_adapter_end
from pcr import PCRGeometry
import seqprimers as sp

def mixture(*species):
    sc = Scene()
    for n, (label, parts) in enumerate(species):
        sc.strand('species_' + str(n), list(parts), label=label)
    return sc

def cleanup(g, ident, source, x, state, *, library=False):
    condition = g.input(ident + '_condition', source, K.CONDITION)
    out = g.product(ident + '_dna', 'recovered DNA', K.LIBRARY if library else K.MOLECULAR, state=state)
    waste = g.product(ident + '_waste', 'physically removed fractions', K.WASTE)
    g.react(ident, source, S.CLEANUP, ((x, R.TEMPLATE), (condition, R.CONDITION)), (out, waste))
    g.discard(waste, 'Actual source-described physical separation; no complete-purity assumption.')
    return out

def selected_stock(g,ident,stock,label,*,kind=None):
    selected=g.product(ident+'_selected',label,kind or stock.kind,state=stock.state,chemical_state=stock.chemical_state)
    remaining=g.product(ident+'_remaining','Untransferred stock, if any; source amount/fate unassigned',stock.kind,state=stock.state,chemical_state=stock.chemical_state)
    g.react(ident,label,S.SPLIT,((stock,R.TEMPLATE),),(selected,remaining))
    g.mark_terminal(Terminal(remaining.id,Disposition.UNRESOLVED,'Source selection does not establish total stock or fate/amount of any remainder.'))
    return selected

def qc_stock(g,ident,stock,label):
    measured=g.product(ident+'_qc',label,K.QC,state=stock.state)
    remaining=g.product(ident+'_remaining','Remaining actual source stock after QC aliquot',stock.kind,state=stock.state,chemical_state=stock.chemical_state)
    g.react(ident,label,S.SPLIT,((stock,R.TEMPLATE),),(measured,remaining))
    g.mark_terminal(Terminal(measured.id,Disposition.QC,'Actual source assay aliquot; quantity unassigned, not sequencing material.'))
    return remaining

def render_page():
    from protocol_views import assemble_page,WorkflowSection
    from page import head,info
    from chemdraw import panel,workflow_from_protocol_graph,workflow_panel
    import html
    g=protocol_graph()
    # Intentionally private: the original required source item remains unmet.
    assert not g.publishable
    return assemble_page(g,(head(TITLE),'<main class="wrap"><h1>'+html.escape(TITLE)+'</h1>',
        info('PRIVATE SOURCE-SCOPE HOLD — this page is not a published complete protocol.'),info(SOURCE),'<p>'+SUMMARY+'</p>',info(CAVEAT),
        '<h2>Source-supported preparation</h2><p>'+'MspI-selected mouse genomic fractions are A-extended, ligated to methylated adapters, bisulfite treated, PCR copied and gel recovered. The 3′A substrate class is shown without inventing a complementary adapter end. Chemical U precedes copied/read T; the diagram is a source idealization, not a conversion-efficiency model.'+'</p>',
        '<h2>Readout and limits</h2><p>'+'Source Illumina 1G 36-base 5′end reads. Actual historical adapter/run primers, index architecture, whole occupied construct and kit recipe remain unassembled. No 2005 RRBS or later TruSeq recipe is substituted.'+'</p>',info('Required source hold: the full published VOR comparison and original cited Illumina DNA preparation/EpiTect workup originals remain relevant and unreviewed. The current authority is the associated accepted manuscript plus the stated publisher supplement, not presumed complete final-version equivalence.'),
        '<details open><summary>Full private workflow and unresolved molecular checks</summary>',workflow_panel(workflow_from_protocol_graph(g,next(iter(g.sequencing)))),'</details>',
        '<p><a href="'+NOTES+'">Research notes and remaining sources</a></p></main>'))

def selftest():
    from checks import Check,Source,have
    from storage import sources_root
    c=Check();g=protocol_graph();page=render_page();ports=lambda i:tuple(p.material for p in g.steps[i].inputs)
    c('original required-source item remains unmet',g.publication_problems(),["authoritative source item 'unreviewed_original' is not accounted for"])
    c('no complete marker emitted','data-protocol-complete=' in page,False)
    c('end processing has molecular nucleotide supply','repair_nucleotides' in ports('end_processing'),True)
    c('ligation catalyst supply explicit','ligation_catalysts' in ports('ligate'),True)
    c('bisulfite is actual molecular reagent','bisulfite_reactive' in ports('conversion'),True)
    c('PCR structural nucleotide supply','pcr_nucleotides' in ports('pcr'),True)
    c('actual source final stain contributes','sybr_green' in ports('stain_final_gel'),True)
    c('source stained material precedes final band recovery','stained_library' in ports('final'),True)
    c('library loading remainder retained',g.terminals['loading_remaining'].disposition,Disposition.UNRESOLVED)
    c('unreviewed source not excluded','unreviewed_original' in g.source_coverage,False)
    source=Source(sources_root()/HERE.parent.name/'nature07107_PMC2896277.xml.txt','AssociatedAMMethodsSummary116')
    if have(c,source,label='RRBS selected source'):
        text=' '.join(source.path.read_text().split())
        c('source A extension before adapter ligation','3′-terminal-A extended' in text,True)
        c('source final SYBR stain','SYBR-green-stained gel slices' in text,True)
    c.report()

PAGE_FILE = 'rrbs-illumina.html'
TITLE = 'RRBS — 2008 Illumina reduced representation bisulfite sequencing'
NOTES = '01_rrbs-illumina.html'
SOURCE = '<a href="https://doi.org/10.1038/nature07107">Meissner et al. (2008)</a>; associated accepted manuscript PMC2896277 Methods Summary; publisher Supplementary Figure S2.'
SUMMARY = 'MspI size selection precedes methylated Illumina adapter ligation, bisulfite conversion and PCR. Physical U is distinguished from the T subsequently copied and reported.'
CAVEAT = 'INFERRED — contextual A-end view only; genomic sequence and whole occupied adapter junction remain unassembled. Partial model from the associated accepted manuscript and publisher supplement; complete VOR comparison, historical adapter polymers/end geometry, original kit workup, PCR and run primers remain unreviewed. Genomic base classes are schematic, not assigned genomic coordinates or quantitative conversion predictions. No other RRBS library recipe is imported.'
@dataclass(frozen=True)
class BisulfiteState:
    stage: str
    adapters: bool = False
    debris: bool = False
    def __post_init__(self):
        if self.stage not in {'input', 'converted', 'copied'}:
            raise ValueError('unknown RRBS molecular stage')
    def rows(self):
        sc = Scene()
        classes = [('unmodified C', 'C'), ('5mC', 'C')]
        for n, (label, base) in enumerate(classes):
            value = 'U' if self.stage == 'converted' and n == 0 else 'T' if self.stage == 'copied' and n == 0 else base
            parts = [role('genomic flank'), Segment('physical base ' + value, value), role('genomic flank')]
            if self.adapters:
                parts = [role('copied source adapter-derived end; chemistry unassembled' if self.stage == 'copied' else 'source methylated adapter; exact polymer unassembled')] + parts + [role('opposite source adapter-derived end')]
            sc.strand(str(n), parts, label=label + ' starting class')
            sc.note(str(n), 'Physical product: ' + ('5mC' if n == 1 and self.stage != 'copied' else value) + '; source idealized path, not conversion yield')
        if self.debris:
            sc.strand('carryover', [role('reaction components; not covalently joined to DNA')], label='same tube until cleanup')
        return sc.rows()
FINAL_LIBRARY = Construct((role('unassembled source Illumina PCR end'), role('RRBS genomic copy; U-derived T, retained-C-derived C'), role('opposite source PCR end')), name='RRBS library roles')
FINAL_CAPTION = 'Source-positive library components; exact original polymer assembly and run primer landing sites remain unassigned.'
SEQ_PRIMERS = ()
SEQUENCING_UNAVAILABLE = 'Illumina 1G Genome Analyser, 36 bp end reads; complete original adapter/run-primer/index architecture not assembled in this pass.'
SEQUENCING_ENDING = 'Alignment against a converted reduced-representation reference permits strand-aware C/T interpretation. This graph does not establish perfect conversion or discrimination of every modified cytosine.'
READ_SPECS = (sp.read_spec('Read 1', 36, 'bisulfite-derived genomic copy', 'include'),)
def protocol_graph():
    g = ProtocolGraph(TITLE, representation='published_instance', authority='Meissner 2008 associated manuscript Methods Summary and publisher Supplementary Figure S2')
    dna = g.input('dna', '1–10 µg purified mouse genomic DNA', K.MOLECULAR, state=BisulfiteState('input'), primary=True)
    enzyme = g.input('mspi', '10–100 U MspI', K.CATALYST)
    cut = g.product('cut', 'MspI digested genomic DNA; all fragments still present', K.MOLECULAR, state=BisulfiteState('input', debris=True))
    g.react('digest', 'MspI digestion', S.DIGESTION, ((dna, R.TEMPLATE), (enzyme, R.CATALYST)), (cut,))
    selected = cleanup(g, 'select', 'Size selection: 40–120 bp and 120–220 bp fractions', cut, BisulfiteState('input'))
    selected=selected_stock(g,'fraction_selection',selected,'Selected actual source fraction contribution;40–120/120–220 interval assignment and inter-fraction pooling unassigned')
    repair = g.input('repair', 'source fill-in and 3-prime terminal-A extension reagents; members unreviewed', K.CATALYST)
    repair_nt=g.input('repair_nucleotides','Aggregate nucleotide-containing supplies for source fill-in/3A extension; exact identity/dose/stock packaging unassembled',K.MOLECULAR)
    ended = g.product('ended', 'filled and 3-prime-A-extended DNA', K.MOLECULAR, state=molecule('A-extended genomic strand', (role('filled genomic DNA; opposite end also A-extended'), Segment('terminal A', 'A'))))
    g.react('end_processing', 'Fill in and add terminal 3-prime A', S.END_REPAIR, ((selected, R.TEMPLATE), (repair, R.CATALYST),(repair_nt,R.STRUCTURAL)), (ended,))
    ready = cleanup(g, 'end_cleanup', 'Phenol extraction and ethanol precipitation', ended, ended.state)
    adapter = g.input('adapter', 'pre-annealed Illumina adapter with 5mC instead of C; original polymer unresolved', K.OLIGO, state=molecule('source methylated adapter', (role('exact original methylated adapter; end geometry unreviewed'),)))
    kit = g.input('ligation_kit', 'source Illumina DNA preparation kit and protocol', K.CONDITION)
    ligase=g.input('ligation_catalysts','Aggregate source kit ligating catalyst supplies; identity/dose/packaging unassigned, no fresh named enzyme inferred',K.CATALYST)
    ligated = g.product('ligated', 'methylated adapter-bearing genomic DNA', K.MOLECULAR, state=BisulfiteState('input', adapters=True, debris=True))
    geometry_scene=Scene()
    geometry_scene.strand('substrate',[role('genomic DNA body; occupied sequence unassembled'),Segment('source terminal A','A',inferred=True)],label='INFERRED contextual view of source3A substrate class',mod3='OH')
    geometry_scene.ligation_end('substrate','source terminal A',kind='cohesive',polarity='3-prime',owner='substrate')
    geometry_scene.note('substrate','Source-known terminal A class; opposite strand/end and whole historical adapter/junction remain unassembled. No adapter T is inferred.')
    geometry = AdapterLigationGeometry('unavailable', 'unavailable', (geometry_scene,), 'Source-positive pre-annealed methylated adapter ligation; exact original adapter end and kit revision remain unreviewed, not established absent.')
    g.react('ligate', 'Ligate methylated Illumina adapters', S.LIGATION, ((ready, R.TEMPLATE), (adapter, R.STRUCTURAL), (kit, R.CONDITION),(ligase,R.CATALYST)), (ligated,), geometry=geometry)
    clean = cleanup(g, 'ligation_cleanup', 'QIAquick cleanup', ligated, BisulfiteState('input', adapters=True))
    bis = g.input('bisulfite', 'Qiagen EpiTect Bisulphite kit', K.CONDITION)
    reactive=g.input('bisulfite_reactive','Actual named bisulfite-containing reactive molecular supply; kit formulation/dose/stock count unassigned',K.MOLECULAR)
    converted = g.product('converted', 'chemical U-containing bisulfite products; original 5mC retained in source diagram', K.MOLECULAR, state=BisulfiteState('converted', adapters=True, debris=True))
    g.react('conversion', 'Bisulfite treatment; C to U precedes copying', S.CHEMICAL_CONVERSION, ((clean, R.TEMPLATE), (bis, R.CONDITION),(reactive,R.STRUCTURAL)), (converted,))
    kit_workup = g.input('kit_workup', 'EpiTect source kit workup; original procedural details deferred', K.CONDITION)
    pcr_ready = g.product('pcr_ready', 'source bisulfite-treated DNA supplied to PCR', K.MOLECULAR, state=BisulfiteState('converted', adapters=True, debris=True))
    g.react('workup_boundary', 'Original EpiTect workup boundary', S.MAGIC, ((converted, R.TEMPLATE), (kit_workup, R.CONDITION)), (pcr_ready,), reason='Named kit treatment is source-positive, but original kit workup not read. Carryover is not silently discarded and complete purity is not asserted.')
    forward = g.input('forward', 'original PCR forward primer; member unresolved', K.OLIGO, state=molecule('PCR forward primer', (role('historical PCR forward primer'),)))
    reverse = g.input('reverse', 'original PCR reverse primer; member unresolved', K.OLIGO, state=molecule('PCR reverse primer', (role('historical PCR reverse primer'),)))
    pcr = g.input('pcr_condition', 'preparative-scale PCR in source; exact program/reagents deferred', K.CONDITION)
    pcr_enzyme=g.input('pcr_catalysts','Aggregate source preparative-PCR catalyst supplies; identity/amount/packaging unassigned',K.CATALYST)
    pcr_nt=g.input('pcr_nucleotides','Aggregate source preparative-PCR nucleotide supplies; formulation/amount/packaging unassigned',K.MOLECULAR)
    amplified = g.product('amplified', 'copied RRBS library', K.MOLECULAR, state=BisulfiteState('copied', adapters=True, debris=True))
    g.react('pcr', 'Preparative PCR copies U-derived sequence as T', S.PCR, ((pcr_ready, R.TEMPLATE), (forward, R.STRUCTURAL), (reverse, R.STRUCTURAL), (pcr, R.CONDITION),(pcr_enzyme,R.CATALYST),(pcr_nt,R.STRUCTURAL)), (amplified,), geometry=PCRGeometry.unavailable(molecule('converted DNA', (role('source adapter-bearing bisulfite product'),)), reason='Original PCR primer members and full adapter/strand assembly unassigned.'))
    purified = cleanup(g, 'pcr_cleanup', 'QIAquick purification of PCR products', amplified, BisulfiteState('copied', adapters=True))
    dye=g.input('sybr_green','Actual source SYBR-green stain for final selected agarose gel; amount/carryover unassigned',K.MOLECULAR)
    stained=g.product('stained_library','Source PCR products in SYBR-green-stained gel; dye carryover unquantified',K.MOLECULAR,state=BisulfiteState('copied',adapters=True,debris=True))
    g.react('stain_final_gel','Source SYBR-green staining before band excision/recovery',S.BINDING,((purified,R.TEMPLATE),(dye,R.STRUCTURAL)),(stained,),changes_molecule=False)
    library = cleanup(g, 'final', '4% NuSieve gel: 130–210 or 210–310 bp; QIAquick gel recovery', stained, FINAL_LIBRARY, library=True)
    library=selected_stock(g,'loading',library,'Selected source36bp sequencing contribution; loading amount unassigned',kind=K.LIBRARY)
    g.sequence(library, SequencingPlan(library=FINAL_LIBRARY, limitation=SEQUENCING_UNAVAILABLE, read_lengths={'Read 1': 36}, acquired_roles=('Read 1',)))
    g.declare_reviewed_scope(('dna','mspi','cut','digest','select_condition','select_dna','select_waste','select','repair','ended','end_processing','end_cleanup_condition','end_cleanup_dna','end_cleanup_waste','end_cleanup','adapter','ligation_kit','ligated','ligate','ligation_cleanup_condition','ligation_cleanup_dna','ligation_cleanup_waste','ligation_cleanup','bisulfite','converted','conversion','kit_workup','pcr_ready','workup_boundary','forward','reverse','pcr_condition','amplified','pcr','pcr_cleanup_condition','pcr_cleanup_dna','pcr_cleanup_waste','pcr_cleanup','final_condition','final_dna','final_waste','final'), g.authority)
    g.declare_reviewed_scope(ADDITIONAL_SOURCE_INVENTORY,g.authority)
    g.add_source_item(SourceItem('unreviewed_original', 'Published article and original library kit/workup not fully reviewed', 'Publisher VOR DOI10.1038/nature07107 and original kit instructions cited in Methods Summary'))
    return g

ADDITIONAL_SOURCE_INVENTORY=('fraction_selection','repair_nucleotides','ligation_catalysts','bisulfite_reactive','pcr_catalysts','pcr_nucleotides','sybr_green','stain_final_gel','loading','loading_selected')
