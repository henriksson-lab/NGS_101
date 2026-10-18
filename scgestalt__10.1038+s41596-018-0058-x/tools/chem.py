"""scGESTALT native source-owned molecular populations and ordered preparation."""
from types import SimpleNamespace
from chemdraw import Construct,Scene,Segment,feature,revcomp
from protocol_model import role
from protocol_graph import InputRole,MaterialKind,ProtocolGraph,SequencingPlan,StepKind,SourceClaim,SourceDiscrepancy,Terminal,Disposition
from pcr import PCRGeometry
from illumina import P5,P7
import seqprimers as sp
TITLE='scGESTALT — matched lineage and transcriptome libraries'
NOTES='01_scgestalt.html'
SOURCE='<a href="https://doi.org/10.1038/s41596-018-0058-x">Raj andGagnon,NatureProtocols2018</a>, procedure82–105,Fig3 and originalSupplementaryTable2.'
SUMMARY='inDrops cell-indexed cDNA supplies IVT RNA. Physicalaliquots retain all transcript classes:9µl is fragmented for expression;5µl ofunfragmented remainder undergoes randomhexamerRT and three source PCR rounds for lineage readout.'
CAVEAT='The defining protocol names/refers to inDrops source73 for native capture, RT, second-strand and expression finishing components. Known operations are retained but those exact source components/first-cycle geometry remain bounded; no generic adapters or expectedfinal-template substitution. GP12N9 adds sequence diversity, not aUMI. Both aliquots retain transcriptmixtures and originalDNA/copy/offproduct recovery is unresolved except actual cleanup fractions. Published PCRround2 cycle rows2–9 followed16/17 and R1PCR supplementary-table crossreference inconsistencies remain public below.'
CELL1=feature('scgestalt.cell1','cell_barcode','combinatorial',group='cell',part='1')
CELL2=feature('scgestalt.cell2','cell_barcode','combinatorial',group='cell',part='2')
UMI=feature('scgestalt.umi','umi','random')
SAMPLE=feature('scgestalt.sample','sample_index','whitelist',whitelist='publishedTable2ix1–24')
SEQ_PRIMERS=()
FINAL_LIBRARY=None
SHOW_READ_CONTENT=False
_CLEANUP_WASTE='actual nonretained AMPure fraction; complete classpurity not presumed'
def identifiers():return [role('inDropscellpart1; nativepositionunassigned',length_bp=8,molecular_feature=CELL1),role('inDropscellpart2; nativepositionunassigned',length_bp=8,molecular_feature=CELL2),role('nativeUMIlineage',length_bp=6,molecular_feature=UMI)]
def native_state(stage):
 sc=Scene()
 if stage=='cell':
  sc.strand('genome',[role('integrated edited/un-edited lineagecassette in genomicDNA')],label='originalgenomicDNA');sc.strand('RNA',[role('independent lineage/transcriptomeRNA population; physicalU')],label='originalRNA');return sc
 rna=stage in ('RNA','IVTRNA','fragmentedRNA')
 parts=[role(('original' if stage=='RNA' else 'newlytranscribedIVT')+'RNA withphysicalU; mixtureoflineageandcellulartranscripts' if rna else 'actual '+stage+' DNA; full nativegeometry/ends unassigned')]
 if stage!='RNA':parts+=identifiers()
 if stage in ('PCR2','indexed'):parts.append(role('GP12N9 sequence-diversity footprint; notUMI',length_bp=9))
 if stage in ('indexed','expressionPCR'):parts.append(role('installedR1-PCRix sampleindex',length_bp=8,molecular_feature=SAMPLE))
 if stage=='duplex':parts.insert(0,role('nativeT7promoter DNAduplex; bases/extentunassigned'))
 if stage in ('duplex','PCR1','PCR2','indexed','expressionPCR'):
  sc=Scene.duplex(parts,label='sourcepositive'+stage+' DNAduplex; fullnativeendgeometrybounded',unpaired=tuple(p.name+"'" for p in parts if p.is_role_token()))
 else:sc.strand('nativepopulation',parts,label=stage+(' RNA' if rna else ' singleDNA/copyclass; fullnativegeometry bounded'))
 residual='independent originalRNA/genomicDNA/cellularproteins; recoveryunassigned' if stage=='RNA' else 'independent originalRNA/genomicDNA/earlierDNAcopies/oligos; residualIVTtemplateDNA afterIVT; physicalrecoveryunresolved'
 sc.strand('residualorigins',[role(residual)],label='copying or selection eligibility is not physicaldisposal');return sc

READ_SPECS=(sp.read_spec('Read 1','61expression /250lineage','transcript or editedlineagearray','include'),sp.read_spec('Read 2',14,'cellpart2 eightbases followedUMI6','extract'),sp.read_spec('Index 1 (i7)',8,'inDropscellpart1','extract'),sp.read_spec('Index 2 (i5)',8,'librarysampleindex','demultiplex'))
SEQUENCING_ENDING='Defining protocol explicit61/8/8/14 expression and250/8/8/14 lineage; standardIlluminaprimers named, native bindinggeometry deferred.'
DEFINITION=SimpleNamespace(title=TITLE,citation='Raj andGagnon2018',doi='10.1038/s41596-018-0058-x',endpoint=SEQUENCING_ENDING)

def _ampure(graph, ident, label, material, ratio, *, state=None, kind=MaterialKind.MOLECULAR, kept_label=None, source=None):
    """An AMPure XP bead cleanup with the retained eluate and the removed supernatant."""
    beads = graph.input(ident + '_beads', f'AMPure XP beads, {ratio}', MaterialKind.SELECTION, source=source)
    kept = graph.product(ident + '_eluate', kept_label or label + ' (eluate)', kind, state=material.state if state is None else state, source=source)
    waste = graph.product(ident + '_supernatant', _CLEANUP_WASTE, MaterialKind.WASTE, source=source)
    graph.react(ident, label, StepKind.CLEANUP, ((material, InputRole.TEMPLATE), (beads, InputRole.SELECTION)), (kept, waste), changes_molecule=False, source=source)
    graph.discard(waste, 'supernatant removed from the magnet-held AMPure XP beads')
    return kept
def _oligo(name, *parts):
    return Construct([Segment(label, bases, tag, feature=feat) for label, bases, tag, feat in parts], name=name)
GP6 = _oligo('inDrops_GP6', ('DsRed-transgene forward', 'GAGGACTACACCATCGTGGAG', None, None))
PE1SA = _oligo('inDrops_PE1Sa', ('PE1Sa (inDrops adaptor overlap)', 'CTCTTTCCCTACACGACGCTGGGTGTCGGGTGCAG', None, None))
GP12 = _oligo('inDrops_GP12', ("GP12 Nextera-type 5' arm", 'TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG', None, None), ('9 random bases (diversity)', 'NNNNNNNNN', None, None), ('barcode-array forward', 'TCGAGCTCAAGCTTCGGAC', None, None))
PE1SB = _oligo('inDrops_PE1Sb', ('PE1Sb (inDrops adaptor overlap)', 'CTCTTTCCCTACACGACGCT', None, None))
R1_PCR_IX1 = _oligo('R1-PCR-ix1', ('P5', P5, 'p5', None), ('library index ix1', 'CTCTCTAT', 'cbc', SAMPLE), ("R1-PCR 3' end", 'TCGTCGGCAGCGTC', None, None))
R2_PCR = _oligo('R2-PCR', ('P7', P7, 'p7', None), ('PE1 overlap', 'GGGTGTCGGGTGCAG', None, None))
_TABLE='Defining published Reagents primer table'
OLIGO_SOURCES = ((GP6, _TABLE), (PE1SA, _TABLE), (GP12, _TABLE), (PE1SB, _TABLE), (R1_PCR_IX1, 'Supplementary Table 2 (first of 24 R1-PCR-ix members)'), (R2_PCR, _TABLE))
def scgestalt_graph(definition) -> ProtocolGraph:
    """scGESTALT graph with matched transcriptome and lineage-library outcomes."""
    source = f'{definition.citation}; doi:{definition.doi}'
    graph = ProtocolGraph(definition.title, authority=source)
    cell = graph.input('transgenic_cell', 'cell with integrated expressed lineage cassette', MaterialKind.BIOLOGICAL, primary=True, state=native_state('cell'))
    guides = graph.input('guide_pools', 'sgRNAs 1–4 injected at one cell; U6 sgRNAs 5–9 transgenic', MaterialKind.OLIGO)
    cas9 = graph.input('cas9', 'injected Cas9 protein, then heat-shock-induced hsp:Cas9', MaterialKind.CATALYST)
    edited_cell = graph.product('edited_cell', 'cell with cumulatively edited lineage cassette', MaterialKind.MOLECULAR, state=native_state('cell'))
    graph.react('edit_lineage', 'Induce cumulative lineage-cassette editing', StepKind.GENOME_EDITING, ((cell, InputRole.TEMPLATE), (guides, InputRole.STRUCTURAL), (cas9, InputRole.CATALYST)), (edited_cell,))
    host_rnap = graph.input('host_rna_polymerase', 'cellular RNA polymerase (heat-shock-induced cassette)', MaterialKind.CATALYST)
    cell_rna = graph.product('cellular_rna', 'polyadenylated edited-cassette RNA plus cellular mRNA', MaterialKind.MOLECULAR, state=native_state('RNA'))
    graph.react('transcribe_cassette', 'Transcribe and polyadenylate the edited cassette', StepKind.TRANSCRIPTION, ((edited_cell, InputRole.TEMPLATE), (host_rnap, InputRole.CATALYST)), (cell_rna,))
    indrops_primer = graph.input('indrops_capture_primer', 'inDrops V3 hydrogel-bead primer (T7–PE1–cell ID–UMI–oligo(dT))', MaterialKind.OLIGO, source=source + '; Fig. 3b legend; bases in inDrops ref. 73, unassigned in inspected scope')
    rtase = graph.input('reverse_transcriptase', 'inDrops in-droplet reverse transcriptase (ref. 73)', MaterialKind.CATALYST)
    first_strand = graph.product('barcoded_first_strand', 'cell-indexed lineage and transcriptome cDNA', MaterialKind.MOLECULAR, state=native_state('firstss'))
    graph.react('reverse_transcription', 'Capture lineage RNA and transcriptome by inDrops', StepKind.REVERSE_TRANSCRIPTION, ((cell_rna, InputRole.TEMPLATE), (indrops_primer, InputRole.STRUCTURAL), (rtase, InputRole.CATALYST)), (first_strand,))
    second_primer = graph.input('second_strand_primer', 'inDrops second-strand synthesis reagents (ref. 73)', MaterialKind.OLIGO, inferred=True)
    second_polymerase = graph.input('second_strand_polymerase', 'second-strand polymerase', MaterialKind.CATALYST)
    duplex = graph.product('barcoded_cdna_duplex', 'barcoded second-strand cDNA', MaterialKind.MOLECULAR, state=native_state('duplex'))
    graph.react('second_strand', 'Synthesize the second cDNA strand', StepKind.EXTENSION, ((first_strand, InputRole.TEMPLATE), (second_primer, InputRole.STRUCTURAL), (second_polymerase, InputRole.CATALYST)), (duplex,))
    t7 = graph.input('t7_rna_polymerase', 'T7 RNA polymerase', MaterialKind.CATALYST)
    amplified_rna = graph.product('amplified_rna', 'IVT-amplified barcoded RNA', MaterialKind.MOLECULAR, state=native_state('IVTRNA'))
    graph.react('ivt', 'Linearly amplify barcoded cDNA by T7 IVT', StepKind.TRANSCRIPTION, ((duplex, InputRole.TEMPLATE), (t7, InputRole.CATALYST)), (amplified_rna,))
    purified_rna = _ampure(graph, 'post_ivt_cleanup', 'Purify IVT RNA (1.3× AMPure XP; 25 µl RNA elution buffer)', amplified_rna, '1.3×')
    transcriptome_rna = graph.product('transcriptome_rna', '9 µl RNA for transcriptome library', MaterialKind.MOLECULAR, state=native_state('IVTRNA'))
    lineage_rna = graph.product('lineage_rna', 'remaining unfragmented lineage-cassette RNA', MaterialKind.MOLECULAR, state=native_state('IVTRNA'))
    graph.react('split_ivt', 'Split purified IVT RNA: 9 µl to transcriptome, rest kept unfragmented', StepKind.SPLIT, ((purified_rna, InputRole.TEMPLATE),), (transcriptome_rna, lineage_rna))
    fragmentation = graph.input('rna_fragmentation_reagent', 'chemical RNA fragmentation', MaterialKind.CONDITION)
    fragmented = graph.product('fragmented_transcriptome_rna', 'fragmented amplified RNA', MaterialKind.MOLECULAR, state=native_state('fragmentedRNA'))
    graph.react('fragment_transcriptome', 'Chemically fragment the transcriptome branch', StepKind.FRAGMENTATION, ((transcriptome_rna, InputRole.TEMPLATE), (fragmentation, InputRole.CONDITION)), (fragmented,))
    transcriptome_rt_primer = graph.input('transcriptome_rt_primer', 'random-hexamer RT primer with PE2 adaptor (ref. 73)', MaterialKind.OLIGO, inferred=True, source=source + '; Fig. 3b legend; bases unassigned in inspected scope')
    expression_rtase = graph.input('expression_rtase', 'inDrops expressionlibrary RT catalyst; native source component unassigned', MaterialKind.CATALYST)
    expression_cdna = graph.product('expression_cdna', 'PE2-tagged transcriptome ssDNA', MaterialKind.MOLECULAR, inferred=True, state=native_state('expressionSS'))
    graph.react('transcriptome_rt', 'Random-hexamer–PE2 RT of the fragmented RNA', StepKind.REVERSE_TRANSCRIPTION, ((fragmented, InputRole.TEMPLATE), (transcriptome_rt_primer, InputRole.STRUCTURAL), (expression_rtase, InputRole.CATALYST)), (expression_cdna,))
    kapa = graph.input('kapa_hifi', '2× KAPA HiFi HotStart PCR mix', MaterialKind.CATALYST)
    r1 = graph.input('r1_pcrix', 'R1-PCRix (24 indexed members; ix1 shown)', MaterialKind.OLIGO, state=R1_PCR_IX1, source=source + '; Supplementary Table 2')
    r2 = graph.input('r2_pcr', 'R2-PCR', MaterialKind.OLIGO, state=R2_PCR, source=source + '; Reagents primer table')
    expression_library = graph.product('expression_library', 'matched scGESTALT transcriptome library', MaterialKind.LIBRARY, inferred=True, state=native_state('expressionPCR'))
    graph.react('expression_pcr', 'PCR the transcriptome branch with the same R1-PCRix/R2-PCR pair', StepKind.PCR, ((expression_cdna, InputRole.TEMPLATE), (r1, InputRole.STRUCTURAL), (r2, InputRole.STRUCTURAL), (kapa, InputRole.CATALYST)), (expression_library,), geometry=PCRGeometry.unavailable(expression_cdna.state, reason='The PE2 RT-primer and bead PE1 native component geometry is deferred to inDrops ref. 73, not this protocol, so primer landing is not assembled.'))
    graph.sequence(expression_library, SequencingPlan(limitation='NextSeq 75-cycle high output, standard Illumina primers.', read_lengths={'Read 1': 61, 'Index 1': 8, 'Index 2': 8, 'Read 2': 14}))
    lineage_sample_rule=graph.input('lineage_sample_rule','source5µl unfragmentedRNA into lineageRT',MaterialKind.CONDITION)
    lineage_input=graph.product('lineage_input','5µl actualmixedunfragmentedIVTRNA forlineageRT',MaterialKind.MOLECULAR,state=lineage_rna.state)
    lineage_remaining=graph.product('lineage_remaining','unallocatedunfragmentedRNA after5µlsampling',MaterialKind.MOLECULAR,state=lineage_rna.state)
    graph.react('lineage_sample','Take5µl actualunfragmentedRNA mixture forRT',StepKind.SAMPLING,((lineage_rna,InputRole.TEMPLATE),(lineage_sample_rule,InputRole.CONDITION)),(lineage_input,lineage_remaining),changes_molecule=False)
    graph.mark_terminal(Terminal(lineage_remaining.id,Disposition.UNRESOLVED,'Stored/unallocatedRNA not assumedphysicallydestroyed; compositionnotmodalitypurified.'))
    lineage_rt_primer = graph.input('lineage_random_hexamer', 'random hexamers (3.75 µM)', MaterialKind.OLIGO)
    primescript = graph.input('primescript', 'PrimeScript reverse transcriptase + RNaseOut', MaterialKind.CATALYST)
    lineage_cdna = graph.product('lineage_cdna', 'full target-array lineage cDNA', MaterialKind.MOLECULAR, state=native_state('lineageSS'))
    graph.react('lineage_rt', 'Random-hexamer RT of unfragmented RNA (70 °C 3 min; 30 °C 10 min, 42 °C 1 h, 70 °C 15 min)', StepKind.REVERSE_TRANSCRIPTION, ((lineage_input, InputRole.TEMPLATE), (lineage_rt_primer, InputRole.STRUCTURAL), (primescript, InputRole.CATALYST)), (lineage_cdna,))
    lineage_cdna = _ampure(graph, 'lineage_rt_cleanup', 'Purify RT product (1.2× AMPure XP)', lineage_cdna, '1.2×')
    q5 = graph.input('q5', 'Q5 DNA polymerase with Q5 enhancer', MaterialKind.CATALYST)
    nested1_forward = graph.input('lineage_nested1_forward', 'inDrops_GP6 (DsRed)', MaterialKind.OLIGO, state=GP6)
    nested1_reverse = graph.input('lineage_nested1_reverse', 'inDrops_PE1Sa', MaterialKind.OLIGO, state=PE1SA)
    nested = graph.product('lineage_nested_amplicon', 'long DsRed–array–cell-ID amplicon', MaterialKind.MOLECULAR, state=native_state('PCR1'))
    graph.react('lineage_nested_pcr1', 'PCR round 1:7.5µl recoveredRT input; GP6 + PE1Sa,Q5,15cycles (61 °C anneal)', StepKind.PCR, ((lineage_cdna, InputRole.TEMPLATE), (nested1_forward, InputRole.STRUCTURAL), (nested1_reverse, InputRole.STRUCTURAL), (q5, InputRole.CATALYST)), (nested,), geometry=PCRGeometry.unavailable(lineage_cdna.state, reason='The DsRed/array cassette and inDrops PE1 bead bases are unassigned in inspected scope, so primer landing is not assembled.'))
    nested = _ampure(graph, 'pcr1_cleanup', 'Purify round-1 PCR (0.6× AMPure XP)', nested, '0.6×')
    phusion = graph.input('phusion', 'Phusion High-Fidelity DNA polymerase (HF buffer)', MaterialKind.CATALYST)
    nested2_forward = graph.input('lineage_nested2_forward', 'inDrops_GP12 (N9 diversity)', MaterialKind.OLIGO, state=GP12)
    nested2_reverse = graph.input('lineage_nested2_reverse', 'inDrops_PE1Sb', MaterialKind.OLIGO, state=PE1SB)
    short = graph.product('lineage_short_amplicon', 'short array–cell-ID amplicon with GP12 arm', MaterialKind.MOLECULAR, state=native_state('PCR2'))
    graph.react('lineage_nested_pcr2', 'PCR round 2:8µl recoveredround1input; GP12 + PE1Sb,Phusion (60 °C anneal; cycle count ambiguous)', StepKind.PCR, ((nested, InputRole.TEMPLATE), (nested2_forward, InputRole.STRUCTURAL), (nested2_reverse, InputRole.STRUCTURAL), (phusion, InputRole.CATALYST)), (short,), geometry=PCRGeometry.unavailable(nested.state, reason='The barcode-array and PE1 bead bases are unassigned in inspected scope, so primer landing is not assembled.'))
    short = _ampure(graph, 'pcr2_cleanup', 'Purify round-2 PCR (0.6× AMPure XP)', short, '0.6×')
    lineage_library = graph.product('lineage_library', 'scGESTALT lineage library', MaterialKind.LIBRARY, inferred=True, state=native_state('indexed'))
    graph.react('lineage_final_pcr', 'Final PCR on 1:10 diluted eluate: R1-PCRix + R2-PCR, KAPA HiFi (2 cycles 55 °C, 6–8 cycles 65 °C)', StepKind.PCR, ((short, InputRole.TEMPLATE), (r1, InputRole.STRUCTURAL), (r2, InputRole.STRUCTURAL), (kapa, InputRole.CATALYST)), (lineage_library,), geometry=PCRGeometry.unavailable(short.state, reason='R1-PCRix overlaps the GP12 arm and R2-PCR the PE1 overlap, but the bead PE1/cell-ID native component geometry is deferred to inDrops ref. 73.'))
    lineage_library = _ampure(graph, 'final_cleanup', 'Purify final library (0.6× AMPure XP)', lineage_library, '0.6×', kind=MaterialKind.LIBRARY, kept_label='sequencing-ready scGESTALT lineage library')
    graph.sequence(lineage_library, SequencingPlan(limitation=definition.endpoint, read_lengths={'Read 1': 250, 'Index 1': 8, 'Index 2': 8, 'Read 2': 14}))
    graph.annotate(SourceDiscrepancy('r1_pcrix', 'Which supplementary table lists the R1-PCRix primers', (SourceClaim('Reagents primer table, R1-PCRix row (PMC6279253)', 'Index sequences: see Supplementary Table 2.'), SourceClaim("Reagent setup, 'R1-PCRix/R2-PCR primer mix'", 'See Supplementary Table 1 for a list of 24 R1-PCRix primer sequences.'), SourceClaim('Supplementary Information (MOESM1)', 'Table 1 lists sgRNA target sites; Table 2 lists R1-PCR-ix1 to ix24.')), 'The model takes R1-PCRix members from Supplementary Table 2 and shows ix1 as one example.', 'The supplement itself puts the 24 R1-PCR-ix sequences in Table 2.', 'Only the table label is in question; no primer sequence differs between the passages.'))
    graph.annotate(SourceDiscrepancy('lineage_nested_pcr2', 'PCR round-2 cycle numbering', (SourceClaim('Procedure, PCR round 2 thermocycler table', "Cycle rows read '1', '2–9', then '16' (final extension) and '17' (hold)."), SourceClaim('Procedure, PCR round 1 thermocycler table', "The analogous table reads '1', '2–16', '17', '18' (15 amplification cycles).")), 'The model gives no round-2 cycle count; it states only the 60 °C anneal.', "Rows 2–9 imply 8 cycles, but '16' for the next row implies 2–15 (14 cycles).", 'Which count the authors intended is not stated; Step 96 says to lower input or cycles if yield is high.'))
    graph.declare_reviewed_scope(REVIEWED, source)
    return graph.require_publishable()
def oligo_rows():
    """Key-oligo rows rendered from the graph's source-transcribed oligo states."""
    return tuple((f'<b>{c.name}</b> ({where}): 5′-{c.top()}-3′' for c, where in OLIGO_SOURCES))
REVIEWED=('transgenic_cell', 'guide_pools', 'cas9', 'edited_cell', 'edit_lineage', 'host_rna_polymerase', 'cellular_rna', 'transcribe_cassette', 'indrops_capture_primer', 'reverse_transcriptase', 'barcoded_first_strand', 'reverse_transcription', 'second_strand_primer', 'second_strand_polymerase', 'barcoded_cdna_duplex', 'second_strand', 't7_rna_polymerase', 'amplified_rna', 'ivt', 'transcriptome_rna', 'lineage_rna', 'split_ivt', 'rna_fragmentation_reagent', 'fragmented_transcriptome_rna', 'fragment_transcriptome', 'transcriptome_rt_primer', 'expression_rtase', 'expression_cdna', 'transcriptome_rt', 'kapa_hifi', 'r1_pcrix', 'r2_pcr', 'expression_library', 'expression_pcr', 'lineage_random_hexamer', 'primescript', 'lineage_cdna', 'lineage_rt', 'q5', 'lineage_nested1_forward', 'lineage_nested1_reverse', 'lineage_nested_amplicon', 'lineage_nested_pcr1', 'phusion', 'lineage_nested2_forward', 'lineage_nested2_reverse', 'lineage_short_amplicon', 'lineage_nested_pcr2', 'lineage_library', 'lineage_final_pcr', 'post_ivt_cleanup_beads', 'post_ivt_cleanup', 'lineage_rt_cleanup_beads', 'lineage_rt_cleanup', 'pcr1_cleanup_beads', 'pcr1_cleanup', 'pcr2_cleanup_beads', 'pcr2_cleanup', 'final_cleanup_beads', 'final_cleanup')
REVIEWED += ('final_cleanup_eluate','lineage_sample_rule','lineage_sample')
def oligos():return oligo_rows()
def protocol_graph():return scgestalt_graph(DEFINITION)
