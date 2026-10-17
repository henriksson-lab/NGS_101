"""Source-owned K562 total-RNA TimeLapse route with dated-kit scope explicit."""
from chemdraw import Construct, Scene, Segment, feature
from pcr import PCRGeometry
from protocol_graph import Disposition, InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind, Terminal
from protocol_model import molecule, role
import seqprimers as sp

TITLE = 'TimeLapse-seq — K562 total-RNA recoding route'
NOTES = '01_timelapse_seq.html'
SOURCE = '<a href="https://doi.org/10.1038/nmeth.4582">Schofield et al., Nature Methods (2018)</a>; original Supplementary Protocol and Note; dated SMARTer Pico manual (112216) used only as a conditional kit implementation.'
SUMMARY = 'Pulse-labelled RNA contains 4sU. Periodate/TFEA chemistry creates a cytidine analogue C*, whose altered pairing produces apparent U-to-C transitions after RT and sequencing; C* is not ordinary RNA cytidine.'
CAVEAT = ('Published pulse discrepancy: original Supplementary Protocol step 2 prints 100 mM for the K562 four-hour pulse, while the original main article K562-mRNA section and Supplementary Note report 100 µM; the Note gives a successful 10 µM–1 mM range. This page follows the reported 100 µM experiment and leaves pulse concentration unassigned in the graph; the intended Protocol entry remains uncertain. '
          'This page represents the non-enriched K562 total-RNA route. Targeted Actb/Gapdh/MYC assays, MEF experiments and the MTS-biotin-enriched TT-TimeLapse branch are distinct source variants. '
          'INFERRED — kit stages below use the dated 2016 SMARTer Pico manual as a conditional implementation of the named kit. The study-specific revision, input-dependent fragmentation choice, exact RT/TSO/primer oligos and selected indexes remain unassigned. '
          'RNA populations and rRNA-cleavage outcomes are separate logical classes, not claimed physical fractionations or quantitative conversion yields.')
I5 = feature('kit_i5','sample_index','unknown')
I7 = feature('kit_i7','sample_index','unknown')
FINAL_LIBRARY = Construct((role('kit adapter end; exact bases unassigned'),
    role('conditional kit i5 index; selected member/width unassigned',molecular_feature=I5,inferred=True),
    role('RNA-derived cDNA; sequence and length unassigned'),
    role('conditional kit i7 index; selected member/width unassigned',molecular_feature=I7,inferred=True),
    role('opposite kit adapter end; exact bases unassigned')),
    name='INFERRED — dated SMARTer Pico implementation of TimeLapse library')
FINAL_CAPTION = 'INFERRED — conditional dated-kit adapter/index roles, without a default index width. Both unlabelled transcripts and transcripts with pulse-incorporated marks can contribute; recoding is a molecular mark rather than an affinity selection.'
SEQ_PRIMERS = ()
SEQUENCING_ENDING = 'The original K562 four-hour total-RNA experiment used paired-end 150-base sequencing on HiSeq 4000; separate MEF and targeted assays used other runs.'
SEQUENCING_UNAVAILABLE = 'Exact kit oligos, selected indexes and physically loaded run primers are unassigned. No complete sequencing-primer geometry is declared.'
READ_SPECS = (
    sp.read_spec('Read 1','150 in represented K562 experiment','transcript-derived sequence with possible apparent T-to-C conversions','include'),
    sp.read_spec('Read 2','150 in represented K562 experiment','paired transcript-derived sequence','include'),
    sp.read_spec('Index 1 (i7)','unassigned','sample index in conditional dated-kit implementation','demultiplex'),
    sp.read_spec('Index 2 (i5)','unassigned','sample index in conditional dated-kit implementation','demultiplex'),
)

def rna_pool(mark):
    """RNA classes distinguished by pulse incorporation, without a whole-transcript age claim."""
    sc=Scene()
    sc.strand('old',[role('RNA without a pulse-installed mark')],label='RNA without pulse-installed mark')
    sc.strand('new',[role('RNA flank; may predate pulse'),Segment(mark,mark,placeholder=True,length_bp=1),role('RNA flank after marked site')],label=('RNA with representative converted pulse-marked site' if mark=='C*' else 'RNA with representative pulse-labelled site'))
    sc.note('new','A labelled transcript may contain a pre-existing 5-prime tract; no complete labelling or conversion yield is claimed.')
    if mark=='C*':
        sc.strand('possible_unconverted',[role('possible pulse-labelled RNA retaining unconverted sites',inferred=True)],label='INFERRED — possible unconverted pulse-labelled RNA class')
    return sc


def protocol_graph():
    graph=ProtocolGraph(TITLE,authority='Original TimeLapse main K562 four-hour experiment; Supplementary Protocol/Note; conditional kit manual112216')
    cells=graph.input('cells','living K562 cells during pulse',MaterialKind.BIOLOGICAL,state=molecule('cells',(role('cellular RNA populations'),)),primary=True)
    label=graph.input('four_su','4-thiouridine pulse; graph concentration unassigned',MaterialKind.MOLECULAR)
    polymerase=graph.input('cellular_rna_polymerase','cellular RNA polymerases',MaterialKind.CATALYST)
    labelled=graph.product('labelled_rna','RNA pool with and without pulse-incorporated labels',MaterialKind.MOLECULAR,state=rna_pool('4sU'))
    graph.react('pulse_label','Pulse-label newly synthesized RNA',StepKind.TRANSCRIPTION,
        ((cells,InputRole.TEMPLATE),(label,InputRole.STRUCTURAL),(polymerase,InputRole.CATALYST)),(labelled,))
    extraction=graph.input('rna_extraction','source modified RNeasy RNA extraction',MaterialKind.CONDITION)
    isolated=graph.product('isolated_rna','isolated total RNA',MaterialKind.MOLECULAR,state=labelled.state)
    wash=graph.product('extraction_wash','nonretained extraction fractions',MaterialKind.WASTE)
    graph.react('extract_rna','Extract total RNA without affinity enrichment',StepKind.CLEANUP,
        ((labelled,InputRole.TEMPLATE),(extraction,InputRole.CONDITION)),(isolated,wash))
    graph.discard(wash,'Source extraction discards flow-through/washes; RNA molecular recovery is not quantified here.')
    chemistry=graph.input('recoding_chemistry','NaIO4 oxidation with TFEA nucleophilic substitution',MaterialKind.CONDITION)
    recoded=graph.product('recoded_rna','RNA pool after recoding; C* analogue is not canonical C',MaterialKind.MOLECULAR,state=rna_pool('C*'))
    graph.react('recode','Oxidize and recode 4sU sites to cytidine analogues',StepKind.CHEMICAL_CONVERSION,
        ((isolated,InputRole.TEMPLATE),(chemistry,InputRole.CONDITION)),(recoded,))
    beads=graph.input('rnaclean','source RNAClean bead cleanup reagent',MaterialKind.SELECTION)
    cleaned=graph.product('cleaned_recoded_rna','recoded RNA after first bead cleanup',MaterialKind.MOLECULAR,state=recoded.state)
    first_wash=graph.product('first_cleanup_wash','first cleanup supernatant/washes',MaterialKind.WASTE)
    graph.react('cleanup_recoding','RNAClean cleanup after recoding',StepKind.CLEANUP,
        ((recoded,InputRole.TEMPLATE),(beads,InputRole.SELECTION)),(cleaned,first_wash))
    graph.discard(first_wash,'Original Supplementary Protocol discards supernatant/washes before RNA elution.')
    dtt=graph.input('post_recoding_reduction','source post-recoding DTT reducing buffer',MaterialKind.CONDITION)
    reduced=graph.product('reduced_rna','RNA after post-recoding reducing treatment',MaterialKind.MOLECULAR,state=recoded.state)
    graph.react('reduce','Apply source post-recoding reducing treatment',StepKind.CHEMICAL_CONVERSION,
        ((cleaned,InputRole.TEMPLATE),(dtt,InputRole.CONDITION)),(reduced,))
    ready=graph.product('kit_input','recoded RNA prepared for named kit; main article uses 10 ng total RNA per library',MaterialKind.MOLECULAR,state=reduced.state)
    second_wash=graph.product('second_cleanup_wash','second cleanup supernatant/washes',MaterialKind.WASTE)
    graph.react('cleanup_reduction','RNAClean cleanup after reducing treatment',StepKind.CLEANUP,
        ((reduced,InputRole.TEMPLATE),(beads,InputRole.SELECTION)),(ready,second_wash))
    graph.discard(second_wash,'Original Supplementary Protocol repeats bead cleanup before kit entry.')
    conditioned=graph.product('conditioned_kit_input','kit-ready RNA; exact input conditioning unassigned',MaterialKind.MOLECULAR,state=ready.state,inferred=True)
    graph.react('kit_input_conditioning','INFERRED — study-specific kit input conditioning unassigned',StepKind.MAGIC,
        ((ready,InputRole.TEMPLATE),),(conditioned,),reason='The dated manual offers fragmentation-dependent and no-fragmentation routes according to RNA integrity. The actual study revision/input choice is not assigned here.')
    primer=graph.input('random_primer','conditional dated-kit N6/random primer mix; exact oligo unassigned',MaterialKind.OLIGO,inferred=True)
    tso=graph.input('lnatso','conditional dated-kit LNA template-switch oligo mix; exact oligo unassigned',MaterialKind.OLIGO,inferred=True)
    rt=graph.input('rt','SMARTScribe RT in dated-kit implementation',MaterialKind.CATALYST,inferred=True)
    cdna=graph.product('cdna','INFERRED — random-primed/template-switched cDNA',MaterialKind.MOLECULAR,
        state=Construct((role('random-primer-derived end',inferred=True),role('RNA-derived cDNA; C* readout depends on RT',inferred=True),role('template-switch-derived end',inferred=True)),name='conditional kit cDNA'),inferred=True)
    graph.react('reverse_transcribe','INFERRED — dated-kit random priming and LNA template switching',StepKind.REVERSE_TRANSCRIPTION,
        ((conditioned,InputRole.TEMPLATE),(primer,InputRole.STRUCTURAL),(tso,InputRole.STRUCTURAL),(rt,InputRole.CATALYST)),(cdna,))
    f1=graph.input('pcr1_forward','conditional kit indexed forward primer; exact sequence unassigned',MaterialKind.OLIGO,inferred=True)
    r1=graph.input('pcr1_reverse','conditional kit indexed reverse primer; exact sequence unassigned',MaterialKind.OLIGO,inferred=True)
    pol=graph.input('pcr_polymerase','dated-kit PCR enzyme mix',MaterialKind.CATALYST,inferred=True)
    indexed=graph.product('indexed_cdna','INFERRED — PCR1 indexed cDNA pool',MaterialKind.MOLECULAR,state=FINAL_LIBRARY,inferred=True)
    graph.react('pcr1','INFERRED — dated-kit PCR1 adds indexed adapters',StepKind.PCR,
        ((cdna,InputRole.TEMPLATE),(f1,InputRole.STRUCTURAL),(r1,InputRole.STRUCTURAL),(pol,InputRole.CATALYST)),(indexed,),
        geometry=PCRGeometry.unavailable(cdna.state,reason='Exact dated-kit primer/RT/TSO sequences and study-specific revision are not declared.'))
    kit_beads=graph.input('kit_cleanup_beads','dated-kit AMPure cleanup reagent',MaterialKind.SELECTION,inferred=True)
    indexed_clean=graph.product('indexed_clean','INFERRED — purified PCR1 cDNA',MaterialKind.MOLECULAR,state=FINAL_LIBRARY,inferred=True)
    kit_wash1=graph.product('kit_wash1','PCR1 cleanup wash material',MaterialKind.WASTE)
    graph.react('cleanup_pcr1','INFERRED — dated-kit cleanup after PCR1',StepKind.CLEANUP,
        ((indexed,InputRole.TEMPLATE),(kit_beads,InputRole.SELECTION)),(indexed_clean,kit_wash1))
    graph.discard(kit_wash1,'Dated manual removes PCR1 cleanup wash fractions; molecular recovery unquantified.')
    indexed_clean2=graph.product('indexed_clean2','INFERRED — PCR1 cDNA after second consecutive cleanup',MaterialKind.MOLECULAR,state=FINAL_LIBRARY,inferred=True)
    kit_wash1b=graph.product('kit_wash1b','second PCR1 cleanup wash material',MaterialKind.WASTE)
    graph.react('cleanup_pcr1_again','INFERRED — second dated-kit cleanup after PCR1',StepKind.CLEANUP,
        ((indexed_clean,InputRole.TEMPLATE),(kit_beads,InputRole.SELECTION)),(indexed_clean2,kit_wash1b))
    graph.discard(kit_wash1b,'Dated manual SectionV.C explicitly repeats AMPure cleanup before ZapR.')
    zap=graph.input('zapr','conditional dated-kit ZapR cleavage enzyme',MaterialKind.CATALYST,inferred=True)
    probes=graph.input('rprobes','conditional dated-kit mammalian R-Probes',MaterialKind.OLIGO,inferred=True)
    intact=graph.product('nonrrna_cdna','INFERRED — intact non-ribosomal cDNA outcome',MaterialKind.MOLECULAR,state=FINAL_LIBRARY,inferred=True)
    cut=graph.product('rrna_fragments','INFERRED — cleaved ribosomal-cDNA fragments',MaterialKind.MOLECULAR,
        state=molecule('cleaved ribosomal-cDNA outcome',(role('rRNA-derived cDNA cleavage fragments',inferred=True),)),inferred=True)
    graph.react('deplete_rrna','INFERRED — dated-kit ribosomal cDNA cleavage',StepKind.FRAGMENTATION,
        ((indexed_clean2,InputRole.TEMPLATE),(zap,InputRole.CATALYST),(probes,InputRole.STRUCTURAL)),(intact,cut))
    graph.mark_terminal(Terminal(cut.id,Disposition.UNRESOLVED,'Logical cleavage-product class, not a claimed physical separation. Exact fragment recovery/fate through later amplification and cleanup is unassigned.'))
    f2=graph.input('pcr2_forward','dated-kit universal PCR2 forward primer; exact sequence unassigned',MaterialKind.OLIGO,inferred=True)
    r2=graph.input('pcr2_reverse','dated-kit universal PCR2 reverse primer; exact sequence unassigned',MaterialKind.OLIGO,inferred=True)
    library=graph.product('library',FINAL_LIBRARY.name,MaterialKind.LIBRARY,state=FINAL_LIBRARY,inferred=True)
    graph.react('pcr2','INFERRED — dated-kit PCR2 amplifies remaining intact libraries',StepKind.PCR,
        ((intact,InputRole.TEMPLATE),(f2,InputRole.STRUCTURAL),(r2,InputRole.STRUCTURAL),(pol,InputRole.CATALYST)),(library,),
        geometry=PCRGeometry.unavailable(FINAL_LIBRARY,reason='Exact universal PCR2 oligos and binding geometry are unassigned.'))
    final_clean=graph.product('final_clean',FINAL_LIBRARY.name,MaterialKind.LIBRARY,state=FINAL_LIBRARY,inferred=True)
    kit_wash2=graph.product('kit_wash2','final library cleanup wash material',MaterialKind.WASTE)
    graph.react('cleanup_pcr2','INFERRED — dated-kit final cleanup after PCR2',StepKind.CLEANUP,
        ((library,InputRole.TEMPLATE),(kit_beads,InputRole.SELECTION)),(final_clean,kit_wash2))
    graph.discard(kit_wash2,'Dated manual final purification removes wash fractions; exact discarded molecular species unassigned.')
    graph.sequence(final_clean,SequencingPlan(library=FINAL_LIBRARY,limitation=SEQUENCING_UNAVAILABLE))
    graph.declare_reviewed_scope(('cells','four_su','cellular_rna_polymerase','pulse_label','rna_extraction','extract_rna',
        'recoding_chemistry','recode','rnaclean','cleanup_recoding','post_recoding_reduction','reduce','cleanup_reduction',
        'kit_input_conditioning','random_primer','lnatso','rt','reverse_transcribe','pcr1_forward','pcr1_reverse','pcr_polymerase',
        'pcr1','zapr','rprobes','deplete_rrna','pcr2_forward','pcr2_reverse','pcr2','library','kit_cleanup_beads','cleanup_pcr1','cleanup_pcr1_again','cleanup_pcr2','final_clean'),graph.authority)
    return graph.require_publishable()
