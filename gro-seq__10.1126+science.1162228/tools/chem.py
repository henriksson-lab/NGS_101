"""Original Core2008 GRO-seq chemistry, owned by its defining source."""
from types import SimpleNamespace
from adapter_ligation import AdapterLigationGeometry
from chemdraw import Construct
from pcr import PCRGeometry
from protocol_graph import (AnnotationSlot, Highlight, InputRole, MaterialKind,
                            ProtocolGraph, SequencingPlan, StepKind)
from protocol_model import molecule, role
import seqprimers as sp

TITLE = "GRO-seq — global run-on sequencing"
NOTES = "01_gro-seq.html"
SOURCE = ('<a href="https://doi.org/10.1126/science.1162228">Core, Waterfall and Lis, Science (2008)</a>; '
          'original Supporting Online Material, Methods pp.17–19 and Figure S1.')
SUMMARY = ("Already-engaged polymerases incorporate BrUTP during nuclear run-on. After RNA recovery, "
           "base hydrolysis and end repair, the original method ligates the 5′ adapter before the 3′ "
           "adapter, with three anti-BrU capture/elution rounds in total. Reverse transcription, "
           "RNA digestion, PCR and PAGE yield the sequencing library.")
CAVEAT = ("The original 2008 supplement has been recovered and specifies the historical reaction order. "
          "It names Illumina Small RNA Isolation Kit adapters and primers; exact ordered bases and kit "
          "revision remain open source questions in this model. The later poly(A)-primed/circularized "
          "GRO-seq implementation is a separate route. BrU is an affinity label, not an eight-base UMI.")
FINAL_LIBRARY = Construct([role("historical library end derived from the 5′ RNA adapter"),
                           role("nascent-RNA-derived insert"),
                           role("historical library end derived from the 3′ RNA adapter")], name="GRO-seq 2008 library")
SEQ_PRIMERS = ()
READ_LENGTHS = {"Read 1": 33}
ACQUIRED_ROLES = tuple(READ_LENGTHS)
READ_SPECS = (sp.read_spec("Read 1", READ_LENGTHS["Read 1"],
                           "5′ end of fragmented BrU-containing RNA; historical primer site unassigned", "include"),)
FINAL_CAPTION = ("Role-level amplified historical library. No modern canonical arms, sample index or fixed insert "
                 "length are selected. The original supplement reports the amplified library averaging about 100 bases "
                 "(~190-base product minus ~90 bases of adapter), which bounds but does not identify the unassigned "
                 "historical adapter bases.")
SEQUENCING_INTRO = ""
SEQUENCING_ENDING = ("The original study obtained 33-base reads on the Illumina 1G Genome Analyzer, "
                     "sequencing from the 5′ end of fragmented BrU-containing RNA. The read 5′ coordinate "
                     "was used as the polymerase-position measurement. Iterative computational 3′ trimming "
                     "to a minimum of 29 bases and removal of rRNA alignments do not change the physical library. "
                     "Exact historical sequencing-primer binding sites remain unassigned.")
SEQUENCING_UNAVAILABLE = ""


def oligos():
    return ()


def _state(label, *parts):
    return molecule(label, tuple(role(part) for part in parts))


def _extract(graph, prefix, template, label):
    condition = graph.input(prefix + "_condition", label, MaterialKind.CONDITION)
    recovered = graph.product(prefix + "_rna", "recovered RNA", MaterialKind.MOLECULAR, state=template.state)
    removed = graph.product(prefix + "_removed", "RNA extraction removed fraction", MaterialKind.WASTE)
    graph.react(prefix, label, StepKind.EXTRACTION,
                ((template, InputRole.TEMPLATE), (condition, InputRole.CONDITION)), (recovered, removed))
    graph.discard(removed, "is separated from the source RNA fraction during extraction/precipitation")
    return recovered


def _capture_elute(graph, prefix, template):
    """Repeated original2008 anti-BrU binding, washes, elution and extraction."""
    beads = graph.input(prefix + "_beads", "anti-deoxyBrU beads", MaterialKind.SELECTION)
    retained = graph.product(prefix + "_bound", "bead-bound BrU-containing RNA", MaterialKind.MOLECULAR, state=template.state)
    wash = graph.product(prefix + "_wash", "unbound binding/wash fraction", MaterialKind.WASTE)
    graph.react(prefix + "_capture", "Capture BrU-containing RNA (" + prefix + ")", StepKind.AFFINITY_CAPTURE,
                ((template, InputRole.TEMPLATE), (beads, InputRole.SELECTION)), (retained, wash))
    graph.discard(wash, "is removed during source antibody binding/washes")
    elution = graph.input(prefix + "_elution", "source Buffer E elution, four collected fractions", MaterialKind.CONDITION)
    eluate = graph.product(prefix + "_eluate", "combined source RNA eluates", MaterialKind.MOLECULAR, state=template.state)
    spent = graph.product(prefix + "_spent", "beads separated from RNA eluates", MaterialKind.WASTE)
    graph.react(prefix + "_release", "Elute RNA from antibody beads", StepKind.EXTRACTION,
                ((retained, InputRole.TEMPLATE), (elution, InputRole.CONDITION)), (eluate, spent))
    graph.discard(spent, "is separated from the collected RNA eluates")
    return _extract(graph, prefix + "_extract", eluate, "Extract and precipitate antibody-enriched RNA")


def protocol_graph():
    source = "Core et al.2008, original Supporting Online Material Methods629–714/overview381–410; doi:10.1126/science.1162228"
    graph = ProtocolGraph("GRO-seq 2008", authority=source)
    cells = graph.input("cells", "IMR90 cells", MaterialKind.BIOLOGICAL,
                        state=_state("source cells", "nuclei with engaged polymerases and nascent RNA"), primary=True)
    isolation = graph.input("isolation", "cold swelling, Igepal lysis and nuclei washes", MaterialKind.CONDITION)
    nuclei = graph.product("nuclei", "isolated nuclei", MaterialKind.MOLECULAR, state=cells.state)
    supernatant = graph.product("isolation_removed", "lysis and nuclei-wash supernatants", MaterialKind.WASTE)
    graph.react("isolate_nuclei", "Isolate and wash nuclei", StepKind.EXTRACTION,
                ((cells, InputRole.TEMPLATE), (isolation, InputRole.CONDITION)), (nuclei, supernatant))
    graph.discard(supernatant, "is removed during source nuclei preparation")
    ntps = graph.input("run_on_ntps", "500 µM ATP, GTP and BrUTP with limiting 2 µM CTP in the reaction buffer", MaterialKind.MOLECULAR)
    tracer = graph.input("ctp_tracer", "source α-32P-CTP tracer", MaterialKind.MOLECULAR)
    polymerase = graph.input("polymerase", "already-engaged RNA polymerases", MaterialKind.CATALYST)
    run_condition = graph.input("run_condition", "source 1% Sarkosyl, 300 mM KCl run-on buffer; 5 min at 30 °C", MaterialKind.CONDITION)
    labelled = graph.product("labelled_rna", "run-on-labelled nuclear RNA", MaterialKind.MOLECULAR,
                             state=_state("run-on RNA", "pre-existing nascent RNA", "extension containing BrU"))
    graph.react("run_on", "Extend already-engaged transcripts", StepKind.RUN_ON,
                ((nuclei, InputRole.TEMPLATE), (ntps, InputRole.STRUCTURAL), (tracer, InputRole.STRUCTURAL),
                 (polymerase, InputRole.CATALYST), (run_condition, InputRole.CONDITION)), (labelled,))
    dnase = graph.input("dnase_first", "RNase-free DNase I", MaterialKind.CATALYST)
    dna_digested = graph.product("dna_digested", "RNA mixture after DNA digestion", MaterialKind.MOLECULAR, state=labelled.state)
    graph.react("digest_dna_first", "Digest DNA after run-on", StepKind.DIGESTION,
                ((labelled, InputRole.TEMPLATE), (dnase, InputRole.CATALYST)), (dna_digested,))
    protein = graph.input("protein_removal", "SDS and proteinase K", MaterialKind.CONDITION)
    deproteinized = graph.product("deproteinized", "RNA mixture after protein digestion", MaterialKind.MOLECULAR, state=labelled.state)
    graph.react("remove_proteins", "Digest proteins with SDS/proteinase K", StepKind.DEPROTEINIZATION,
                ((dna_digested, InputRole.TEMPLATE), (protein, InputRole.CONDITION)), (deproteinized,))
    recovered = _extract(graph, "recover_rna", deproteinized, "Acid phenol/chloroform extraction and ethanol precipitation")
    hydrolysis = graph.input("hydrolysis", "NaOH hydrolysis and Tris neutralization", MaterialKind.CONDITION)
    fragmented = graph.product("fragmented", "hydrolysed RNA pool", MaterialKind.MOLECULAR,
                               state=_state("hydrolysed RNA", "RNA fragments including hydrolysis-generated 5′-OH ends", "hydrolysis-generated 3′-phosphate ends"))
    graph.react("hydrolyse", "Base-hydrolyse RNA", StepKind.FRAGMENTATION,
                ((recovered, InputRole.TEMPLATE), (hydrolysis, InputRole.CONDITION)), (fragmented,))
    p30 = graph.input("p30", "two P-30 spin-column passes", MaterialKind.SELECTION)
    desalted = graph.product("desalted", "column-recovered RNA", MaterialKind.MOLECULAR, state=fragmented.state)
    column_removed = graph.product("column_removed", "spin-column removed fraction", MaterialKind.WASTE)
    graph.react("p30_cleanup", "Recover RNA through P-30 columns", StepKind.CLEANUP,
                ((fragmented, InputRole.TEMPLATE), (p30, InputRole.SELECTION)), (desalted, column_removed))
    graph.discard(column_removed, "is separated from the collected RNA during source spin-column cleanup")
    dnase2 = graph.input("dnase_second", "second RNase-free DNase I treatment", MaterialKind.CATALYST)
    edta = graph.input("edta", "EDTA stop after second DNA digestion", MaterialKind.CONDITION)
    redigested = graph.product("redigested", "RNA after second DNA digestion", MaterialKind.MOLECULAR, state=desalted.state)
    graph.react("digest_dna_second", "Further remove DNA; stop with EDTA", StepKind.DIGESTION,
                ((desalted, InputRole.TEMPLATE), (dnase2, InputRole.CATALYST), (edta, InputRole.CONDITION)), (redigested,))
    first = _capture_elute(graph, "first", redigested)
    tap = graph.input("tap", "tobacco acid pyrophosphatase", MaterialKind.CATALYST)
    decapped = graph.product("decapped", "decapped RNA pool", MaterialKind.MOLECULAR, state=first.state)
    graph.react("decap", "Remove remaining 5′ caps with TAP", StepKind.DIGESTION,
                ((first, InputRole.TEMPLATE), (tap, InputRole.CATALYST)), (decapped,))
    pnk_low = graph.input("pnk_low", "PNK under source low-pH conditions", MaterialKind.CATALYST)
    oh = graph.product("three_prime_oh", "RNA with repaired 3′-OH ends", MaterialKind.MOLECULAR,
                       state=_state("3′-repaired RNA", "RNA with 3′-OH ends"))
    graph.react("repair3", "Remove hydrolysis-derived 3′ phosphate", StepKind.DEPHOSPHORYLATION,
                ((decapped, InputRole.TEMPLATE), (pnk_low, InputRole.CATALYST)), (oh,))
    pnk_high = graph.input("pnk_high", "PNK with raised pH and ATP", MaterialKind.CATALYST)
    p = graph.product("five_prime_p", "RNA with ligation-ready 5′-phosphate/3′-OH", MaterialKind.MOLECULAR,
                      state=_state("end-repaired RNA", "5′ phosphate", "RNA fragment", "3′ hydroxyl"))
    graph.react("repair5", "Add RNA 5′ phosphate", StepKind.PHOSPHORYLATION,
                ((oh, InputRole.TEMPLATE), (pnk_high, InputRole.CATALYST)), (p,))
    ready = _extract(graph, "repair_extract", p, "Extract and precipitate end-repaired RNA")
    adapter5 = graph.input("adapter5", "5′ adapter, historical Illumina Small RNA Isolation Kit", MaterialKind.ADAPTER)
    ligase5 = graph.input("ligase5", "T4 RNA ligase-1 with source PEG mix", MaterialKind.CATALYST)
    ligated5 = graph.product("ligated5", "5′-adapter-ligated RNA", MaterialKind.MOLECULAR,
                            state=_state("5′-adapted RNA", "5′ RNA adapter", "BrU-containing RNA", "3′ hydroxyl"))
    graph.react("ligate5", "Ligate the 5′ RNA adapter FIRST", StepKind.LIGATION,
                ((ready, InputRole.TEMPLATE), (adapter5, InputRole.STRUCTURAL), (ligase5, InputRole.CATALYST)),
                (ligated5,), geometry=AdapterLigationGeometry.unavailable("Source RNA end repair is specified, but exact historical adapter identity/end geometry is not assigned by this model."))
    second = _capture_elute(graph, "second", ligated5)
    adapter3 = graph.input("adapter3", "3′ adapter, historical Illumina Small RNA Isolation Kit", MaterialKind.ADAPTER)
    ligase3 = graph.input("ligase3", "T4 RNA ligase-1 with source PEG mix", MaterialKind.CATALYST)
    ligated3 = graph.product("ligated3", "RNA bearing both adapters", MaterialKind.MOLECULAR,
                            state=_state("doubly adapted RNA", "5′ RNA adapter", "BrU-containing RNA", "3′ RNA adapter"))
    graph.react("ligate3", "Ligate the 3′ RNA adapter SECOND", StepKind.LIGATION,
                ((second, InputRole.TEMPLATE), (adapter3, InputRole.STRUCTURAL), (ligase3, InputRole.CATALYST)),
                (ligated3,), geometry=AdapterLigationGeometry.unavailable("Historical 3′ adapter sequence/modifications remain open source questions."))
    third = _capture_elute(graph, "third", ligated3)
    rt_primer = graph.input("rt_primer", "RT primer, historical Illumina Small RNA Isolation Kit", MaterialKind.OLIGO)
    rtase = graph.input("rtase", "SuperScript III reverse transcriptase", MaterialKind.CATALYST)
    cdna = graph.product("cdna", "RNA/cDNA after reverse transcription", MaterialKind.MOLECULAR,
                         state=_state("reverse transcription products", "adapter-primed nascent-RNA-derived cDNA", "RNA template"))
    graph.react("reverse_transcribe", "Reverse-transcribe in two separate 10-µL reactions", StepKind.REVERSE_TRANSCRIPTION,
                ((third, InputRole.TEMPLATE), (rt_primer, InputRole.STRUCTURAL), (rtase, InputRole.CATALYST)), (cdna,))
    rnases = graph.input("rnases", "RNase cocktail and RNase H", MaterialKind.CATALYST)
    template = graph.product("pcr_template", "cDNA after RNA digestion", MaterialKind.MOLECULAR,
                             state=_state("cDNA mixture", "adapter-derived handle", "nascent-RNA-derived cDNA", "adapter-derived handle at opposite end"))
    graph.react("digest_rna", "Digest RNA before PCR", StepKind.DIGESTION,
                ((cdna, InputRole.TEMPLATE), (rnases, InputRole.CATALYST)), (template,))
    forward = graph.input("pcr_forward", "PCR forward primer specified by historical Illumina kit", MaterialKind.OLIGO)
    reverse = graph.input("pcr_reverse", "PCR reverse primer specified by historical Illumina kit", MaterialKind.OLIGO)
    phusion = graph.input("phusion", "Phusion high-fidelity DNA polymerase", MaterialKind.CATALYST)
    amplified = graph.product("amplified", "amplified cDNA library pool", MaterialKind.MOLECULAR, state=FINAL_LIBRARY)
    graph.react("pcr", "Amplify the source library for 15 cycles", StepKind.PCR,
                ((template, InputRole.TEMPLATE), (forward, InputRole.STRUCTURAL), (reverse, InputRole.STRUCTURAL),
                 (phusion, InputRole.CATALYST)), (amplified,), geometry=PCRGeometry.unavailable(FINAL_LIBRARY,
                    reason="Historical primer sequences and initial cDNA priming geometry remain unassigned; later 2018 primers are not substituted."))
    gel = graph.input("gel", "native 8% TBE PAGE; excise cDNA greater than 90 nt", MaterialKind.SELECTION)
    library = graph.product("library", "GRO-seq 2008 sequencing library", MaterialKind.LIBRARY, state=FINAL_LIBRARY)
    gel_removed = graph.product("gel_removed", "material outside selected gel fraction", MaterialKind.WASTE)
    graph.react("size_select", "PAGE-select, elute, extract and precipitate library", StepKind.SELECTION,
                ((amplified, InputRole.TEMPLATE), (gel, InputRole.SELECTION)), (library, gel_removed))
    graph.discard(gel_removed, "is outside the source selected PAGE fraction")
    graph.sequence(library, SequencingPlan(library=FINAL_LIBRARY, acquired_roles=ACQUIRED_ROLES,
                                          read_lengths=READ_LENGTHS, limitation=SEQUENCING_ENDING))
    graph.annotate(Highlight("ligate5", "Original 2008 adapter order",
                            "The original SOM overview and Figure S1 specify 5′ adapter first, then 3′ adapter, with immunopurification after each. No later poly(A)/circularization recipe is imported."))
    graph.annotate(Highlight("size_select", "Reported library length",
                            "Original SOM overview p.10: after amplification and PAGE purification the library averages "
                            "about 100 bases (~190-base product minus ~90 bases of adapter), with a pre-amplification "
                            "complexity of 1x10^9 molecules estimated by re-amplification of a known aliquot. Those "
                            "quality-control measurements are not separate library-construction reactions here."))
    graph.annotate(Highlight("library", "Physical library versus computational read filtering",
                            "Original reads were 33 bases. Iterative 3′ trimming to 29 and rRNA alignment removal occur after sequencing; they do not discard molecules during this preparation.", slot=AnnotationSlot.SEQUENCING))
    targets = ["cells", "isolation", "isolate_nuclei", "run_on_ntps", "ctp_tracer", "polymerase", "run_condition", "run_on", "dnase_first", "digest_dna_first", "protein_removal", "remove_proteins", "hydrolysis", "hydrolyse", "p30", "p30_cleanup", "dnase_second", "edta", "digest_dna_second", "tap", "decap", "pnk_low", "repair3", "pnk_high", "repair5", "adapter5", "ligase5", "ligate5", "adapter3", "ligase3", "ligate3", "rt_primer", "rtase", "reverse_transcribe", "rnases", "digest_rna", "pcr_forward", "pcr_reverse", "phusion", "pcr", "gel", "size_select", "library"]
    for prefix in ("recover_rna", "repair_extract"):
        targets.extend((prefix, prefix + "_condition"))
    for prefix in ("first", "second", "third"):
        targets.extend(prefix + suffix for suffix in ("_beads", "_capture", "_elution", "_release", "_extract", "_extract_condition"))
    graph.declare_reviewed_scope(targets, source)
    return graph.require_publishable()


SPEC = SimpleNamespace(**{name: globals()[name] for name in (
    "TITLE", "NOTES", "SOURCE", "SUMMARY", "CAVEAT", "FINAL_LIBRARY", "SEQ_PRIMERS", "READ_LENGTHS",
    "ACQUIRED_ROLES", "READ_SPECS", "FINAL_CAPTION", "SEQUENCING_INTRO", "SEQUENCING_ENDING", "SEQUENCING_UNAVAILABLE", "oligos", "protocol_graph")})
