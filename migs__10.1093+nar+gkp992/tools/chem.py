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
        '<h2>Source-supported preparation</h2><p>'+'A His-MBD2/anti-His/protein-G bead stock is physically washed before a selected capture aliquot contacts fragmented DNA. Possible independent blocking-vector DNA is retained as an unresolved carryover species, not a genomic insert or an index. Fragment enrichment is not a site-specific methylation base call.'+'</p>',
        '<h2>Readout and limits</h2><p>'+'Published source GAII 36-base reads and BioAnalyzer insert-distribution observation near 120 bp. Exact historical kit adapters, read-end acquisition/index mask and full PCR/library recipe remain unassigned; the bisulfite-validation primers are not MiGS library primers.'+'</p>',info('Required source hold: the original manufacturer protocol for the cited Illumina ChIP-Seq Sample Prep Kit remains relevant and unreviewed. Its study-used revision and direct original-document URL are unassigned; this is not evidence that instructions were unpublished or inaccessible. The published article itself is available and inspected.'),
        '<details open><summary>Full private workflow and unresolved molecular checks</summary>',workflow_panel(workflow_from_protocol_graph(g,next(iter(g.sequencing)))),'</details>',
        '<p><a href="'+NOTES+'">Research notes and remaining sources</a></p></main>'))

def selftest():
    from checks import Check
    g = protocol_graph(); page = render_page()
    c = Check()
    c('Incomplete source boundary blocks finished publication', bool(g.publication_problems()))
    c('Partial preview does not advertise completion', 'data-protocol-complete=' not in page)
    c('Genomic capture receives physically washed conjugates',
      [port.material for port in g.steps['capture'].inputs if port.role is R.SELECTION],
      ['capture_bead_aliquot_selected'])
    c('Conjugate wash has a physical discard terminal',
      g.terminals['conjugate_wash_waste'].disposition, Disposition.DISCARD)
    for material in ('captured', 'released', 'elution_cleanup_dna', 'library'):
        state = g.materials[material].state
        blocker_rows = [strand for strand in state.strands.values()
                        if any(seg.name == BLOCKER.name for seg in strand.segs)]
        c(material + ' retains possible blocker as a separate molecular species',
          len(blocker_rows) == 1 and len(blocker_rows[0].segs) == 1)
    ports=lambda i:tuple(p.material for p in g.steps[i].inputs)
    c('original required-source item remains unmet',g.publication_problems(),["authoritative source item 'unreviewed_original' is not accounted for"])
    c('capture selected actual DNA stock','capture_dna_aliquot_selected' in ports('capture'),True)
    c('kit receives selected actual10ng contribution','kit_dna_aliquot_selected' in ports('library_prep'),True)
    c('kit oligo inventory explicitly molecular','kit_oligos' in ports('library_prep'),True)
    c('source fragment QC typed',g.materials['fragment_gel_qc_qc'].kind,K.QC)
    c('source library BioAnalyzer QC typed',g.materials['bioanalyzer_qc_qc'].kind,K.QC)
    c('library loading remainder retained',g.terminals['loading_remaining'].disposition,Disposition.UNRESOLVED)
    c('unreviewed source not excluded','unreviewed_original' in g.source_coverage,False)
    from checks import Source,have
    from storage import sources_root
    source=Source(sources_root()/HERE.parent.name/'nar_gkp992_PMC2811030.xml.txt','PublishedMethods44/48')
    if have(c,source,label='MiGS published source'):
        text=' '.join(source.path.read_text().split())
        c('literal published Tris unit retained','20 ml Tris (pH 8)' in text,True)
        c('source actual10ng kit input','Ten nanograms of MBD-isolated DNA' in text,True)
    c.report()

PAGE_FILE = 'migs.html'
TITLE = 'MiGS — MBD2 affinity isolation of methylated genomic fragments'
NOTES = '01_migs.html'
SOURCE = '<a href="https://doi.org/10.1093/nar/gkp992">Serre et al. (online 2009; issue 2010)</a>; published PMC2811030 Materials and Methods and gkp992_1.pdf.'
SUMMARY = 'Sonicated genomic DNA is captured by His-MBD2 domain on anti-His/protein-G beads, washed and released by SDS/Proteinase K, then prepared with the original Illumina ChIP-Seq kit for GAII 36 bp reads.'
CAVEAT = 'Partial source-owned capture model. Exact historical Illumina kit polymers, PCR/run primers and complete read architecture are unassembled. Supplemental bisulfite-validation primers are not MiGS library primers. Possible blocking-vector DNA carryover is retained separately through release and DNA recovery; neither proteolysis nor extraction/precipitation establishes its complete removal. Its recovery and eventual kit/read contribution are unassigned. The published elution buffer prints “20 ml Tris (pH 8)”; this wording is flagged rather than silently replaced with a concentration. Captured-fragment enrichment is not a base-resolution methylation call.'
BLOCKER = role('possible pCR4-TOPO blocking DNA; recovery and eventual kit fate unassigned')
FINAL_LIBRARY = Construct((role('unassembled original ChIP-Seq kit end'), role('MBD-enriched genomic sequence; ordinary DNA letters'), role('opposite source kit end')), name='MiGS GAII library roles')
FINAL_CAPTION = 'Representative source-positive genomic library roles; this panel does not establish removal of blocking DNA or exclude blocker-derived material. No other Illumina library recipe or validation primer is imported.'
SEQ_PRIMERS = ()
SEQUENCING_UNAVAILABLE = 'GAII 36 bp end reads; exact historical kit adapter/run primer and index architecture unassigned.'
SEQUENCING_ENDING = 'Enrichment is interpreted over genomic fragments/windows. Identical mapped starts are computationally collapsed after sequencing, without a declared UMI. A read C is not itself a site-specific methylation call.'
READ_SPECS = (sp.read_spec('Read 1', 36, 'genomic sequence from MBD-enriched fragments', 'include'),)
def protocol_graph():
    g = ProtocolGraph(TITLE, representation='published_instance', authority='Serre et al. published PMC2811030 Materials and Methods38–52; supplement branch classification')
    dna = g.input('dna', 'purified genomic DNA after source Proteinase K/phenol/chloroform/ethanol extraction', K.MOLECULAR, state=molecule('genomic DNA', (role('genomic DNA with heterogeneous methyl-CpG pattern'),)), primary=True)
    sonic = g.input('sonication', 'source sonication in BW buffer to approximately 150–600 bp', K.CONDITION)
    fragmented = g.product('fragmented', 'sonicated genomic fragments and reaction components', K.MOLECULAR, state=mixture(('genomic DNA', (role('sonicated genomic DNA'),)), ('buffer components', (role('independent sonication-buffer components'),))))
    g.react('shear', 'Sonicate genomic DNA', S.FRAGMENTATION, ((dna, R.TEMPLATE), (sonic, R.CONDITION)), (fragmented,))
    fragmented=qc_stock(g,'fragment_gel_qc',fragmented,'Source1percentagarose fragment-size QC; actual assay aliquot unassigned')
    clean = cleanup(g, 'fragment_cleanup', 'QIAquick cleanup, source excludes fragments smaller than 100 bp; 5 µg per capture', fragmented, dna.state)
    clean=selected_stock(g,'capture_dna_aliquot',clean,'Source5µg purified genomicDNA selected per capture; transfer volume unassigned')
    beads = g.input('protein_g_beads', '41.7 µl/ml protein G magnetic bead slurry', K.SELECTION)
    antibody = g.input('antibody', '10 µg/ml anti-His RGS antibody', K.SELECTION)
    mbd = g.input('mbd', '160 nM His-tagged recombinant MBD2 methyl-CpG-binding domain', K.SELECTION)
    blocker = g.input('blocker', '1 µg/ml self-ligated pCR4-TOPO blocking DNA', K.MOLECULAR, state=molecule('blocking vector DNA', (role('self-ligated pCR4-TOPO; not a genomic insert or index'),)))
    conjugation = g.input('conjugation_condition', 'BW buffer, overnight rotation at 4C', K.CONDITION)
    conjugates = g.product('conjugates', 'source antibody/proteinG/MBD bead preparation with blocker', K.MOLECULAR)
    g.react('prepare_capture', 'Prepare the source MBD-domain bead affinity reagent', S.BINDING, ((blocker, R.TEMPLATE), (beads, R.SELECTION), (antibody, R.SELECTION), (mbd, R.SELECTION), (conjugation, R.CONDITION)), (conjugates,), changes_molecule=False)
    conjugate_wash = g.input('conjugate_wash_condition', 'Wash conjugates with one volume cold 1x BW buffer; resuspend in one volume fresh 1x BW buffer', K.CONDITION)
    washed_conjugates = g.product('washed_conjugates', 'Physical washed/resuspended MBD-conjugate stock; residual blocker unquantified', K.MOLECULAR)
    conjugate_waste = g.product('conjugate_wash_waste', 'physically removed conjugate wash; reagent/blocker partition unquantified', K.WASTE)
    g.react('conjugate_wash', 'Wash and resuspend bead conjugates before genomic capture', S.CLEANUP,
            ((conjugates, R.TEMPLATE), (conjugate_wash, R.CONDITION)), (washed_conjugates, conjugate_waste))
    g.discard(conjugate_waste, 'Source physical conjugate wash before DNA addition; no complete blocker removal or quantitative partition asserted.')
    washed_conjugates=selected_stock(g,'capture_bead_aliquot',washed_conjugates,'Source5ml physically washed/resuspended conjugate mixture per capture; source stock total unassigned',kind=K.SELECTION)
    capture_condition = g.input('capture_condition', '5 ml source conjugates per 5 µg fragmented DNA; overnight rotation at 4C, five BW washes', K.CONDITION)
    captured = g.product('captured', 'bead-retained enriched DNA population', K.MOLECULAR, state=mixture(('retained genomic DNA', (role('MBD-affinity-enriched genomic DNA; methylation pattern unspecified'),)), ('possible blocking DNA carryover', (BLOCKER,)), ('source reagent species', (role('independent reagent-associated material; no complete-purity claim'),))))
    flowthrough = g.product('flowthrough', 'physically removed capture/wash fractions', K.WASTE)
    g.react('capture', 'Capture with MBD affinity and physically wash', S.AFFINITY_CAPTURE, ((clean, R.TEMPLATE), (washed_conjugates, R.SELECTION), (capture_condition, R.CONDITION)), (captured, flowthrough))
    g.discard(flowthrough, 'Source physical capture and five washes; unbound does not mean every fragment is unmethylated.')
    pk = g.input('elution', 'SDS/Proteinase K buffer, 60C 2h; printed “20 ml Tris (pH8)” remains unresolved', K.CATALYST)
    released = g.product('released', 'eluted DNA and digestion/reagent products', K.MOLECULAR, state=mixture(('DNA', (role('MBD-enriched genomic DNA'),)), ('possible blocking DNA carryover', (BLOCKER,)), ('reagent debris', (role('Proteinase K digestion products and SDS; present until workup'),))))
    g.react('elute', 'SDS/Proteinase K release without treating digestion as purification', S.DEPROTEINIZATION, ((captured, R.TEMPLATE), (pk, R.CATALYST)), (released,))
    isolated = cleanup(g, 'elution_cleanup', 'Phenol/chloroform extraction; isopropanol precipitation; resuspend in 45 µl water', released, mixture(('recovered genomic DNA', (role('MBD-enriched genomic DNA'),)), ('possible co-recovered blocking DNA', (BLOCKER,))))
    isolated=selected_stock(g,'kit_dna_aliquot',isolated,'Source10ng from45µl recovered DNA stock; transfer volume unassigned')
    kit = g.input('kit', '10 ng DNA; original Illumina ChIP-Seq Sample Prep Kit, manufacturer protocol', K.CONDITION)
    kit_oligos=g.input('kit_oligos','Aggregate actual source kit library-contributing oligo classes; members/counts/whole polymers and packaging unassigned',K.OLIGO,state=molecule('historical kit oligo inventory',(role('source kit polymer components; whole orderings unassembled'),)))
    kit_nt=g.input('kit_nucleotides','Aggregate kit nucleotide-containing molecular supplies; source formulation/doses and stock packaging unassembled',K.MOLECULAR)
    kit_enzyme=g.input('kit_catalysts','Aggregate kit catalyst supplies; actual enzyme identities/doses and packaging unassembled',K.CATALYST)
    library = g.product('library', 'source prepared GAII library; insert distribution near 120 bp', K.LIBRARY, state=mixture(('representative genomic library roles', tuple(FINAL_LIBRARY)), ('blocking DNA or kit-derived material; fate unresolved', (BLOCKER,))))
    g.react('library_prep', 'Original Illumina ChIP-Seq library preparation boundary', S.MAGIC, ((isolated, R.TEMPLATE), (kit, R.CONDITION),(kit_oligos,R.STRUCTURAL),(kit_nt,R.STRUCTURAL),(kit_enzyme,R.CATALYST)), (library,), reason='Named kit and sequencing endpoint source-positive. Exact historical kit processing, adapter members, primer polymers and amplification/run settings unreviewed; no family recipe or validation-primer substitution.')
    library=qc_stock(g,'bioanalyzer_qc',library,'SourceBioAnalyzer library-size QC;~120bp insert distribution is observation, not a fixed synthetic insert')
    library=selected_stock(g,'loading',library,'Selected sourceGAII36bp sequencing contribution; loading amount unassigned',kind=K.LIBRARY)
    g.sequence(library, SequencingPlan(library=FINAL_LIBRARY, limitation=SEQUENCING_UNAVAILABLE, read_lengths={'Read 1':36}, acquired_roles=('Read 1',)))
    g.declare_reviewed_scope(('dna','sonication','fragmented','shear','fragment_cleanup_condition','fragment_cleanup_dna','fragment_cleanup_waste','fragment_cleanup','protein_g_beads','antibody','mbd','blocker','conjugation_condition','conjugates','prepare_capture','conjugate_wash_condition','washed_conjugates','conjugate_wash_waste','conjugate_wash','capture_condition','captured','flowthrough','capture','elution','released','elute','elution_cleanup_condition','elution_cleanup_dna','elution_cleanup_waste','elution_cleanup','kit','library','library_prep'), g.authority)
    g.declare_reviewed_scope(ADDITIONAL_SOURCE_INVENTORY,g.authority)
    g.add_source_item(SourceItem('unreviewed_original', 'Original Illumina ChIP-Seq manufacturer protocol not reviewed', 'Published PMC2811030 Library preparation and sequencing explicitly refers to original manufacturer kit protocol'))
    return g

ADDITIONAL_SOURCE_INVENTORY=('fragment_gel_qc','capture_dna_aliquot','capture_bead_aliquot','kit_dna_aliquot','kit_oligos','kit_nucleotides','kit_catalysts','bioanalyzer_qc','loading','loading_selected')
