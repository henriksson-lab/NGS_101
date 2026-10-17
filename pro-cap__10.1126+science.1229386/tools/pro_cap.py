"""Original2013 PRO-cap source-owned graph; later chemistry is comparator only."""
from adapter_ligation import AdapterLigationGeometry
from chemdraw import Construct, revcomp
import illumina as il
from pcr import PCRGeometry
from protocol_graph import Disposition, InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind, Terminal
from protocol_model import molecule, role
import seqprimers as sp

TITLE = 'PRO-cap — original 2013 capped precision run-on sequencing'
NOTES = '01_pro_cap.html'
SOURCE = '<a href="https://doi.org/10.1126/science.1229386">Kwak et al., Science (2013)</a>; original author-hosted SOM contains both PRO-seq and PRO-cap.'
SUMMARY = 'Run-on with four biotin-NTPs marks nascent RNA. Three affinity enrichments and adapter-first AP/TAP processing preserve the capped 5-prime coordinate.'
CAVEAT = ('The original 2013 branch omits fragmentation and PNK; Antarctic phosphatase precedes TAP. The later 2016 protocol adds an exonuclease branch and is not substituted here. '
          'Exact historical standard Illumina oligos, selected index configuration and loaded run primers are unassigned in this model; no default eight-base index is introduced. '
          'INFERRED — possible uncapped, biotin-labelled 5-prime-OH RNA is a logical outcome class, not a physically separated fraction. Its later fate is unresolved; streptavidin selects biotin rather than cap status.')
FINAL_LIBRARY = Construct((role('original Illumina library end; architecture unassigned'),
                           role('nascent-RNA-derived cDNA; unknown length'),
                           role('opposite library end and indexing configuration unassigned')),
                          name='original 2013 PRO-cap role-level endpoint')
FINAL_CAPTION = 'Original 2013 library roles only. Exact standard adapter/primer identities and index width are unassigned; later comparator oligos are not asserted as historical reagents.'
SEQ_PRIMERS = ()
SEQUENCING_ENDING = 'Original SOM uses the 5-prime ends of non-reverse-complemented PRO-cap reads to locate initiation sites. GenomeAnalyzer 2/HiSeq 2000 are named; acquired cycles and loaded primers are unassigned.'
SEQUENCING_UNAVAILABLE = 'A complete original 2013 loaded sequencing-primer set and corresponding library geometry have not been declared in this model.'
READ_SPECS = (sp.read_spec('Read 1','unassigned','nascent-RNA sense 5-prime boundary; initiation coordinate','include'),)

# Later2016 Table1 comparator declarations, deliberately not primary graph inputs.
RA3_2016 = il.TRUSEQ_SMALL_RNA_RA3.replace('T','U')
RA5_2016 = il.TRUSEQ_SMALL_RNA_RA5.replace('T','U')
RTP_2016 = 'G' + revcomp(il.TRUSEQ_SMALL_RNA_RA3)
RP1_2016 = il.P5 + il.TRUSEQ_SMALL_RNA_RA5[:21]
RPI_2016 = il.P7 + 'N'*6 + 'GTGACTGGAGTTCCTTGGCACCCGAGAATTCCA'

def _ligate_adapter(graph, end, template, state, adapter_label, ligase_label):
    adapter = graph.input(end + '_adapter', adapter_label, MaterialKind.ADAPTER)
    ligase = graph.input(end + '_ligase', ligase_label, MaterialKind.CATALYST)
    product = graph.product(end + '_product', end + '-adapter-ligated nascent RNA',
                            MaterialKind.MOLECULAR, state=state)
    graph.react(end + '_ligation', 'Ligate the ' + end + ' RNA adapter',
                StepKind.LIGATION,
                ((template, InputRole.TEMPLATE), (adapter, InputRole.STRUCTURAL),
                 (ligase, InputRole.CATALYST)), (product,),
                geometry=AdapterLigationGeometry.unavailable(
                    'Historical standard Illumina adapter selection and ligating-end geometry are unassigned in this model; later source comparators do not establish the original selected reagent.'))
    return product


def protocol_graph():
    source = SOURCE
    final = FINAL_LIBRARY
    graph = ProtocolGraph(TITLE, authority=source)
    nuclei = graph.input('nuclei', 'isolated nuclei containing engaged RNA polymerases', MaterialKind.BIOLOGICAL, state=molecule('isolated nuclei', (role('engaged RNA polymerase'), role('nascent RNA'))), primary=True)
    label_name = 'four biotin-NTPs'
    label = graph.input('labelled_ntp', label_name + ' run-on substrates', MaterialKind.MOLECULAR)
    polymerase = graph.input('engaged_polymerase', 'engaged endogenous RNA polymerase', MaterialKind.CATALYST)
    runon = graph.product('runon_rna', 'affinity-labelled nuclear run-on RNA', MaterialKind.MOLECULAR, state=molecule('labelled run-on RNA', (role('capped or uncapped nascent 5-prime end'), role('labelled run-on 3-prime end'))))
    graph.react('run_on', 'Extend engaged polymerases with ' + label_name, StepKind.RUN_ON, ((nuclei, InputRole.TEMPLATE), (label, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)), (runon,))
    preparation = graph.input('rna_preparation','source Trizol RNA extraction and P30 desalting; no fragmentation',MaterialKind.CONDITION)
    extracted = graph.product('extracted_runon_rna','unfragmented extracted run-on RNA',MaterialKind.MOLECULAR,state=runon.state)
    preparation_wash = graph.product('rna_preparation_wash','nonretained extraction/desalting material',MaterialKind.WASTE)
    graph.react('extract_rna','Extract and desalt RNA; omit PRO-seq fragmentation',StepKind.CLEANUP,
        ((runon,InputRole.TEMPLATE),(preparation,InputRole.CONDITION)),(extracted,preparation_wash))
    graph.discard(preparation_wash,'Source RNA extraction/desalting retains RNA while removing other fractions; exact molecular composition unassigned.')
    capture_label = 'streptavidin affinity reagent'
    capture = graph.input('runon_capture', capture_label, MaterialKind.SELECTION)
    enriched = graph.product('enriched_runon', 'affinity-enriched run-on RNA', MaterialKind.MOLECULAR, state=runon.state)
    unlabelled = graph.product('unlabelled_rna', 'RNA not extended in the run-on', MaterialKind.WASTE)
    graph.react('enrich_runon', 'Affinity-enrich labelled run-on RNA', StepKind.AFFINITY_CAPTURE, ((extracted, InputRole.TEMPLATE), (capture, InputRole.SELECTION)), (enriched, unlabelled))
    graph.discard(unlabelled, 'did not acquire the run-on affinity label')
    ligated3 = _ligate_adapter(graph, 'three_prime', enriched,
                              molecule('3-prime-adapter-ligated capped-or-uncapped RNA',
                                       (role('capped or uncapped nascent 5-prime end'), role('nascent RNA body'),
                                        role('3-prime RNA adapter'))),
                              'standard Illumina 3-prime RNA adapter (SOM; sequence/modifications unassigned)',
                              'T4 RNA ligase I (NEB), 6 h 20C in the SOM first-ligation passage; PRO-cap changes only the adapter')
    enriched3 = graph.product('enriched_three_prime_ligated_rna', 'biotin-enriched 3-prime-adapter-ligated RNA', MaterialKind.MOLECULAR, state=ligated3.state)
    rejected3 = graph.product('rejected_three_prime_ligated_rna', 'unlabelled RNA removed after 3-prime ligation', MaterialKind.WASTE)
    graph.react('enrich_three_prime_ligated_rna', 'Re-enrich biotin-labelled RNA after 3-prime ligation', StepKind.AFFINITY_CAPTURE, ((ligated3, InputRole.TEMPLATE), (capture, InputRole.SELECTION)), (enriched3, rejected3))
    graph.discard(rejected3, 'lacks the run-on biotin label')
    phosphatase = graph.input('phosphatase', 'Antarctic phosphatase (original 2013)', MaterialKind.CATALYST)
    dephosphorylated = graph.product('dephosphorylated_rna', '3-prime-ligated RNA after uncapped-end blocking', MaterialKind.MOLECULAR, state=molecule('INFERRED — RNA after AP', (role('protected capped RNA and possible uncapped 5-prime-OH RNA', inferred=True), role('3-prime adapter'))), inferred=True)
    graph.react('dephosphorylate', 'INFERRED — AP removes uncapped 5-prime phosphates while capped ends remain protected', StepKind.DEPHOSPHORYLATION, ((enriched3, InputRole.TEMPLATE), (phosphatase, InputRole.CATALYST)), (dephosphorylated,))
    tap = graph.input('tap', 'tobacco acid pyrophosphatase', MaterialKind.CATALYST)
    cap_selected = graph.product('cap_selected_rna', 'decapped nascent RNA with ligatable 5-prime monophosphate', MaterialKind.MOLECULAR, state=molecule('PRO-cap cap-selected RNA', (role('ligatable decapped 5-prime monophosphate'), role('nascent RNA body'), role('3-prime RNA adapter'))))
    uncapped = graph.product('uncapped_biotin_rna', 'INFERRED — possible uncapped biotin RNA with 5-prime OH', MaterialKind.MOLECULAR,
        state=molecule('uncapped biotin RNA outcome', (role('5-prime OH', inferred=True), role('biotin-labelled RNA body', inferred=True), role('3-prime adapter'))), inferred=True)
    graph.react('decap', 'INFERRED — TAP creates ligatable ends from capped RNA; uncapped 5-prime-OH class remains', StepKind.CHEMICAL_CONVERSION,
        ((dephosphorylated, InputRole.TEMPLATE), (tap, InputRole.CATALYST)), (cap_selected, uncapped))
    graph.mark_terminal(Terminal(uncapped.id, Disposition.UNRESOLVED,
        'Possible chemical outcome, not a physical fractionation. Biotin affinity alone does not remove all uncapped 5-prime-OH RNA; downstream recovery/amplification fate unassigned.'))
    ligated5 = _ligate_adapter(graph, 'five_prime', cap_selected,
                              molecule('PRO-cap RNA with both adapters',
                                       (role('5-prime RNA adapter'), role('nascent RNA body'),
                                        role('3-prime RNA adapter'))),
                              'standard Illumina 5-prime RNA adapters (SOM; sequence unassigned)',
                              '5-prime ligation enzyme; not named in the original SOM')
    enriched5 = graph.product('enriched_five_prime_ligated_rna', 'biotin-enriched RNA with both adapters', MaterialKind.MOLECULAR, state=ligated5.state)
    rejected5 = graph.product('rejected_five_prime_ligated_rna', 'unlabelled RNA removed after 5-prime ligation', MaterialKind.WASTE)
    graph.react('enrich_five_prime_ligated_rna', 'Re-enrich biotin-labelled RNA after 5-prime ligation', StepKind.AFFINITY_CAPTURE, ((ligated5, InputRole.TEMPLATE), (capture, InputRole.SELECTION)), (enriched5, rejected5))
    graph.discard(rejected5, 'lacks the run-on biotin label')
    rt_primer = graph.input('rt_primer', 'original standard Illumina RT primer (25 pmol); exact sequence unassigned; excess also primes PCR', MaterialKind.OLIGO)
    rt = graph.input('reverse_transcriptase', 'source reverse-transcription enzyme', MaterialKind.CATALYST)
    cdna = graph.product('cdna', 'reverse-transcribed PRO-cap RNA', MaterialKind.MOLECULAR,
        state=Construct((role('original primer-derived region unassigned'),role('reverse-transcribed nascent RNA'),role('copied opposite adapter region unassigned')), name='PRO-cap first-strand cDNA'))
    graph.react('reverse_transcribe', 'Reverse-transcribe with original standard Illumina primer', StepKind.REVERSE_TRANSCRIPTION,
        ((enriched5,InputRole.TEMPLATE),(rt_primer,InputRole.STRUCTURAL),(rt,InputRole.CATALYST)),(cdna,))
    pcr_reverse = graph.input('pcr_reverse', 'second PCR primer for PRO-cap; SOM names GX1 or RPI-index for PRO-seq, PRO-cap identity unassigned', MaterialKind.OLIGO)
    polymerase = graph.input('pcr_polymerase', 'Phusion polymerase (original SOM)', MaterialKind.CATALYST)
    library = graph.product('library', final.name, MaterialKind.LIBRARY, state=final)
    graph.react('library_pcr', 'Amplify 12-18 cycles (trial PCR on an RT portion sets cycles); excess RT primer is one PCR primer', StepKind.PCR,
        ((cdna,InputRole.TEMPLATE),(rt_primer,InputRole.STRUCTURAL),(pcr_reverse,InputRole.STRUCTURAL),(polymerase,InputRole.CATALYST)),(library,),
        geometry=PCRGeometry.unavailable(cdna.state,reason='The exact original standard Illumina RT/PCR oligos and selected indexing variant are unassigned; later comparator sequences are not substituted.'))
    page = graph.input('page_selection','source PAGE library purification',MaterialKind.SELECTION)
    purified = graph.product('purified_library',final.name,MaterialKind.LIBRARY,state=final)
    excluded = graph.product('page_excluded','material outside selected PAGE library band',MaterialKind.WASTE)
    graph.react('purify_library','PAGE-purify amplified library',StepKind.CLEANUP,
        ((library,InputRole.TEMPLATE),(page,InputRole.SELECTION)),(purified,excluded))
    graph.discard(excluded,'Original SOM reports PAGE purification; exact excluded molecular species unassigned.')
    graph.sequence(purified,SequencingPlan(library=final,limitation=SEQUENCING_UNAVAILABLE))
    graph.declare_reviewed_scope((
        'nuclei','labelled_ntp','engaged_polymerase','run_on','rna_preparation','extract_rna','runon_capture','enrich_runon',
        'three_prime_adapter','three_prime_ligase','three_prime_ligation','enrich_three_prime_ligated_rna',
        'phosphatase','dephosphorylate','tap','decap','five_prime_adapter','five_prime_ligase','five_prime_ligation',
        'enrich_five_prime_ligated_rna','rt_primer','reverse_transcriptase','reverse_transcribe',
        'pcr_reverse','pcr_polymerase','library_pcr','library','page_selection','purified_library','purify_library'),source)
    return graph.require_publishable()
