"""Source-owned ChRO-seq route; leChRO and adapter variants are explicit boundaries."""
from types import SimpleNamespace
from adapter_ligation import AdapterLigationGeometry
from chemdraw import Construct
from pcr import PCRGeometry
from protocol_graph import (AnnotationSlot, Highlight, InputRole, MaterialKind,
                            ProtocolGraph, SequencingPlan, StepKind)
from protocol_model import molecule, role

TITLE = "ChRO-seq — chromatin run-on sequencing"
NOTES = "01_chro-seq.html"
SOURCE = ('<a href="https://doi.org/10.1038/s41588-018-0244-3">Chu et al., Nature Genetics (2018)</a>; '
          'cited library predecessor: <a href="https://doi.org/10.1038/nprot.2016.086">Mahat et al. (2016), Table 1</a>.')
SUMMARY = ("Isolated chromatin retains engaged polymerases. ATP/GTP plus biotin-CTP/biotin-UTP "
           "support run-on, followed by RNA recovery, base hydrolysis, three streptavidin captures, "
           "two RNA-adapter ligations, reverse transcription and PCR. This graph shows the ChRO branch.")
CAVEAT = ("The cited PRO-seq protocol prints VRA3/VRA5 and RP1/RPI-n; those source oligos are listed "
          "as predecessor reagents. They are not assigned unchanged to every ChRO library. The article "
          "and Supplementary Table 1 identify sample-specific six-base UMI variants. Exact study "
          "oligo orders and run masks remain open model questions. The biotin run-on label is a "
          "chemical affinity label, not a barcode or an eight-base UMI.")
FINAL_LIBRARY = Construct([role("Illumina-compatible left end"),
                           role("nascent-RNA-derived cDNA"),
                           role("Illumina-compatible right end")], name="ChRO-seq library")
SEQ_PRIMERS = ()
READ_LENGTHS = None
FINAL_CAPTION = "Role-level amplified ChRO library; exact study adapter arms, index architecture and primer sites are not selected."
SEQUENCING_INTRO = ""
SEQUENCING_ENDING = ("The source names NextSeq 500. For ChRO-seq, the nascent RNA 3′ end is represented "
                     "by the 5′ end of the sequenced read; leChRO-seq instead represents the RNA 5′ end "
                     "by the read 3′ end. The first six sequenced bases are a molecule identifier only "
                     "for libraries marked in Supplementary Table 1; this does not set an instrument-wide run mask.")
SEQUENCING_UNAVAILABLE = ""
OLIGOS = (
    "Cited PRO-seq predecessor, 3′ RNA adapter VRA3: 5′-pGAUCGUCGGACUGUAGAACUCUGAAC-/Inverted dT/-3′",
    "Cited PRO-seq predecessor, 5′ RNA adapter VRA5: 5′-CCUUGGCACCCGAGAAUUCCA-3′",
    "Cited PRO-seq predecessor, RT/PCR primer RP1: 5′-AATGATACGGCGACCACCGAGATCTACACGTTCAGAGTTCTACAGTCCGA-3′",
    "Cited PRO-seq predecessor, indexed PCR primer RPI-n: 5′-CAAGCAGAAGACGGCATACGAGATNNNNNNGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA-3′",
)

def oligos():
    return OLIGOS


def _state(label, *parts):
    return molecule(label, tuple(role(part) for part in parts))


def _capture_release(graph, prefix, template):
    """One locally declared source streptavidin capture and Trizol recovery."""
    beads = graph.input(prefix + "_beads", "streptavidin beads, source S1421S", MaterialKind.SELECTION)
    bound = graph.product(prefix + "_bound", "biotin-RNA retained on streptavidin beads",
                          MaterialKind.MOLECULAR, state=template.state)
    wash = graph.product(prefix + "_wash", "unbound capture/wash fraction", MaterialKind.WASTE)
    graph.react(prefix + "_capture", "Capture biotin-containing RNA (" + prefix + ")",
                StepKind.AFFINITY_CAPTURE,
                ((template, InputRole.TEMPLATE), (beads, InputRole.SELECTION)), (bound, wash))
    graph.discard(wash, "is removed during source streptavidin binding/washes")
    trizol = graph.input(prefix + "_trizol", "Trizol RNA recovery", MaterialKind.CONDITION)
    released = graph.product(prefix + "_rna", "recovered biotin-RNA",
                             MaterialKind.MOLECULAR, state=template.state)
    removed = graph.product(prefix + "_removed", "beads and extraction removed fraction", MaterialKind.WASTE)
    graph.react(prefix + "_release", "Recover RNA from beads using Trizol",
                StepKind.EXTRACTION, ((bound, InputRole.TEMPLATE), (trizol, InputRole.CONDITION)),
                (released, removed))
    graph.discard(removed, "is separated from recovered RNA during source bead release/extraction")
    return released


def protocol_graph():
    source = "Chu et al.2018, original Methods chromatin isolation and ChRO library preparation; doi:10.1038/s41588-018-0244-3"
    graph = ProtocolGraph("ChRO-seq", authority=source)
    sample = graph.input("sample", "cultured cells or prepared tissue", MaterialKind.BIOLOGICAL,
                         state=_state("source specimen", "chromatin containing engaged polymerases and nascent RNA"), primary=True)
    nun = graph.input("nun", "NUN isolation and wash conditions with RNase inhibitor", MaterialKind.CONDITION)
    chromatin = graph.product("chromatin", "washed chromatin pellet", MaterialKind.MOLECULAR,
                              state=sample.state)
    supernatant = graph.product("isolation_removed", "isolation and wash supernatants", MaterialKind.WASTE)
    graph.react("isolate", "Isolate and wash chromatin in NUN buffer", StepKind.EXTRACTION,
                ((sample, InputRole.TEMPLATE), (nun, InputRole.CONDITION)), (chromatin, supernatant))
    graph.discard(supernatant, "is removed from the chromatin pellet during isolation/washes")
    sonication = graph.input("sonication", "sonication in chromatin storage buffer", MaterialKind.CONDITION)
    suspended = graph.product("suspended_chromatin", "chromatin suspension", MaterialKind.MOLECULAR, state=chromatin.state)
    graph.react("suspend", "Sonicate to suspend the chromatin pellet", StepKind.GENERIC,
                ((chromatin, InputRole.TEMPLATE), (sonication, InputRole.CONDITION)), (suspended,), changes_molecule=False)
    polymerase = graph.input("engaged_polymerase", "already-engaged RNA polymerase", MaterialKind.CATALYST)
    ntps = graph.input("run_on_ntps", "ATP, GTP, biotin-11-CTP and biotin-11-UTP", MaterialKind.MOLECULAR)
    run_condition = graph.input("run_condition", "2× chromatin run-on buffer (300 mM KCl, 1% Sarkosyl, SUPERase In), 37 °C 5 min", MaterialKind.CONDITION)
    labelled = graph.product("labelled_rna", "run-on-labelled chromatin-associated RNA", MaterialKind.MOLECULAR,
                             state=_state("run-on RNA", "pre-existing nascent RNA", "run-on extension containing biotin nucleotides"))
    graph.react("run_on", "Extend already-engaged nascent transcripts", StepKind.RUN_ON,
                ((suspended, InputRole.TEMPLATE), (polymerase, InputRole.CATALYST),
                 (ntps, InputRole.STRUCTURAL), (run_condition, InputRole.CONDITION)), (labelled,))
    trizol = graph.input("stop_trizol", "Trizol LS stop and GlycoBlue-assisted RNA recovery", MaterialKind.CONDITION)
    recovered = graph.product("recovered_rna", "recovered run-on RNA", MaterialKind.MOLECULAR, state=labelled.state)
    removed = graph.product("run_on_removed", "run-on extraction removed fraction", MaterialKind.WASTE)
    graph.react("recover_rna", "Stop run-on and recover RNA", StepKind.EXTRACTION,
                ((labelled, InputRole.TEMPLATE), (trizol, InputRole.CONDITION)), (recovered, removed))
    graph.discard(removed, "is removed during source RNA extraction/precipitation")
    base = graph.input("hydrolysis_base", "65 °C 40 s heat denaturation, then 0.2 N NaOH on ice 8 min", MaterialKind.CONDITION)
    fragments = graph.product("rna_fragments", "hydrolysed RNA pool", MaterialKind.MOLECULAR,
                              state=_state("RNA fragments", "RNA fragments, including biotin-containing run-on ends"))
    graph.react("hydrolyse", "Fragment RNA by base hydrolysis (ChRO branch)", StepKind.FRAGMENTATION,
                ((recovered, InputRole.TEMPLATE), (base, InputRole.CONDITION)), (fragments,))
    first = _capture_release(graph, "first", fragments)
    adapter3 = graph.input("adapter3", "study 3′ RNA adapter (variant under review)", MaterialKind.ADAPTER)
    ligase3 = graph.input("ligase3", "T4 RNA ligase I, source M0204L", MaterialKind.CATALYST)
    ligated3 = graph.product("ligated3", "3′-adapter-ligated RNA", MaterialKind.MOLECULAR,
                            state=_state("3′-adapted RNA", "nascent RNA with internal biotin near run-on 3′ end", "3′ RNA adapter"))
    graph.react("ligate3", "Ligate the 3′ RNA adapter", StepKind.LIGATION,
                ((first, InputRole.TEMPLATE), (adapter3, InputRole.STRUCTURAL), (ligase3, InputRole.CATALYST)),
                (ligated3,), geometry=AdapterLigationGeometry.unavailable("Study-specific adapter variant and exact source end geometry remain under review; predecessor Table1 is not assigned unchanged."))
    second = _capture_release(graph, "second", ligated3)
    rpph = graph.input("rpph", "RppH, source M0356S", MaterialKind.CATALYST)
    decapped = graph.product("decapped", "5′-decapped RNA pool", MaterialKind.MOLECULAR, state=second.state)
    graph.react("decap", "Remove remaining 5′ caps with RppH", StepKind.DIGESTION,
                ((second, InputRole.TEMPLATE), (rpph, InputRole.CATALYST)), (decapped,))
    pnk = graph.input("pnk", "polynucleotide kinase, source M0201L", MaterialKind.CATALYST)
    phosphorylated = graph.product("phosphorylated", "5′-phosphorylated RNA", MaterialKind.MOLECULAR, state=decapped.state)
    graph.react("phosphorylate", "Phosphorylate RNA 5′ ends", StepKind.PHOSPHORYLATION,
                ((decapped, InputRole.TEMPLATE), (pnk, InputRole.CATALYST)), (phosphorylated,))
    cleanup = graph.input("repair_trizol", "Trizol purification after 5′ end preparation", MaterialKind.CONDITION)
    ready = graph.product("ready5", "recovered 5′-prepared RNA", MaterialKind.MOLECULAR, state=phosphorylated.state)
    waste = graph.product("repair_removed", "5′ preparation purification removed fraction", MaterialKind.WASTE)
    graph.react("purify5", "Purify RNA after RppH and PNK", StepKind.EXTRACTION,
                ((phosphorylated, InputRole.TEMPLATE), (cleanup, InputRole.CONDITION)), (ready, waste))
    graph.discard(waste, "is removed during source Trizol purification")
    adapter5 = graph.input("adapter5", "study 5′ RNA adapter (variant under review)", MaterialKind.ADAPTER)
    ligase5 = graph.input("ligase5", "RNA ligase for inherited 5′-adapter step", MaterialKind.CATALYST)
    ligated5 = graph.product("ligated5", "RNA bearing both adapters", MaterialKind.MOLECULAR,
                            state=_state("doubly adapted RNA", "5′ RNA adapter", "biotin-containing nascent RNA", "3′ RNA adapter"))
    graph.react("ligate5", "Ligate the 5′ RNA adapter", StepKind.LIGATION,
                ((ready, InputRole.TEMPLATE), (adapter5, InputRole.STRUCTURAL), (ligase5, InputRole.CATALYST)),
                (ligated5,), geometry=AdapterLigationGeometry.unavailable("Study-specific 5′-adapter variant remains under review."))
    third = _capture_release(graph, "third", ligated5)
    rt_primer = graph.input("rt_primer", "study adapter-complementary RT primer", MaterialKind.OLIGO)
    rtase = graph.input("rtase", "reverse transcription reagent, source 18080-044", MaterialKind.CATALYST)
    cdna = graph.product("cdna", "nascent-RNA-derived cDNA", MaterialKind.MOLECULAR,
                         state=_state("cDNA", "RNA-adapter-derived handles", "nascent-RNA-derived cDNA"))
    graph.react("reverse_transcribe", "Reverse-transcribe adapter-ligated RNA", StepKind.REVERSE_TRANSCRIPTION,
                ((third, InputRole.TEMPLATE), (rt_primer, InputRole.STRUCTURAL), (rtase, InputRole.CATALYST)), (cdna,))
    forward = graph.input("pcr_forward", "study PCR forward primer", MaterialKind.OLIGO)
    reverse = graph.input("pcr_reverse", "study PCR reverse primer", MaterialKind.OLIGO)
    pcr = graph.input("pcr_enzyme", "PCR reagent, source M0491L", MaterialKind.CATALYST)
    final = graph.product("library", "ChRO-seq sequencing library", MaterialKind.LIBRARY, state=FINAL_LIBRARY)
    graph.react("library_pcr", "Amplify to an Illumina-compatible library", StepKind.PCR,
                ((cdna, InputRole.TEMPLATE), (forward, InputRole.STRUCTURAL),
                 (reverse, InputRole.STRUCTURAL), (pcr, InputRole.CATALYST)), (final,),
                geometry=PCRGeometry.unavailable(FINAL_LIBRARY, reason="Exact study primer variants and initial cDNA geometry remain unresolved."))
    graph.sequence(final, SequencingPlan(library=FINAL_LIBRARY, limitation=SEQUENCING_ENDING))
    graph.annotate(Highlight("hydrolyse", "ChRO and leChRO are different source branches",
                            "leChRO is length-extension ChRO: RNase cocktail replaces inhibitor during initial isolation and post-run-on hydrolysis is omitted. It is not a linear-amplification branch."))
    graph.annotate(Highlight("library", "Sample-specific six-base molecule identifiers",
                            "The article says first six sequenced bases only for marked libraries. Original Supplementary Table 1 marks UMU94-13_leChr, GBM-88-04_Chr and PDX-88-04_P57 among others. Blank aggregate/subsample cells are not a universal oligo assignment.", slot=AnnotationSlot.SEQUENCING))
    targets = ["sample", "nun", "isolate", "sonication", "suspend", "engaged_polymerase", "run_on_ntps", "run_condition", "run_on", "stop_trizol", "recover_rna", "hydrolysis_base", "hydrolyse", "adapter3", "ligase3", "ligate3", "rpph", "decap", "pnk", "phosphorylate", "repair_trizol", "purify5", "adapter5", "ligase5", "ligate5", "rt_primer", "rtase", "reverse_transcribe", "pcr_forward", "pcr_reverse", "pcr_enzyme", "library_pcr", "library"]
    for prefix in ("first", "second", "third"):
        targets.extend(prefix + suffix for suffix in ("_beads", "_capture", "_trizol", "_release"))
    graph.declare_reviewed_scope(targets, source)
    return graph.require_publishable()


SPEC = SimpleNamespace(**{name: globals()[name] for name in (
    "TITLE", "NOTES", "SOURCE", "SUMMARY", "CAVEAT", "FINAL_LIBRARY", "SEQ_PRIMERS",
    "READ_LENGTHS", "FINAL_CAPTION", "SEQUENCING_INTRO", "SEQUENCING_ENDING", "SEQUENCING_UNAVAILABLE",
    "oligos", "protocol_graph")})
