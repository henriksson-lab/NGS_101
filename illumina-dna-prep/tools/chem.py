"""Local source and reaction model for illumina-dna-prep.

Authority: Illumina DNA Prep Reference Guide 1000000025416 v12 (TAG/PTC/BLT PCR
programs pp. 23-27, double-sided IPB cleanup pp. 28-30) and the Illumina Prep and
Nextera adapter-sequence table for the public adapter and index-PCR primer arms.
"""
import illumina as il
import nextera as nx
import seqprimers as sp
from chemdraw import Construct, Scene, Segment, Workflow, feature
from pcr import PCRGeometry
from protocol_graph import (AnnotationSlot, Highlight, InputRole, MaterialKind,
                            ProtocolGraph, SequencingPlan, StepKind)
from illumina_vendor_workflows import (Definition, _duplex, _single, _tagged_intermediate,
                                       _workflow, nextera_final_rows)
from types import SimpleNamespace

# Guide p.26: IDT for Illumina DNA/RNA UD and Nextera DNA UD index codes are 10 bp;
# Nextera DNA CD indexes, also listed as compatible on p.5, use 8 bp.  The page draws
# the UD (10 bp) variant and keeps the well-specific bases symbolic.
UD_INDEX_NT = 10

TAG_PROGRAM = "TAG program 55°C 15 min (BLT + TB1, 50 µl)"
PTC_PROGRAM = "TSB added; PTC program 37°C 15 min"
PCR_CYCLES = {"1–9 ng": 12, "10–24 ng": 8, "25–49 ng": 6, "50–99 ng": 5,
              "100–500 ng": 5, "blood/saliva": 5}
PCR_PROGRAM = ("BLT PCR: 68°C 3 min, 98°C 3 min, X × (98°C 45 s, 62°C 30 s, 68°C 2 min), "
               "68°C 1 min; X = " + ", ".join(f"{k}: {v}" for k, v in PCR_CYCLES.items()))


def dna_prep_library(insert_name: str, title: str):
    """Public Nextera endpoint with UD-length (10 bp) dual sample-index windows."""
    lib = Construct([
        Segment("P5", il.P5, "p5"),
        Segment("i5", "J" * UD_INDEX_NT, "cbc", placeholder=True,
                feature=feature("sample_index_i5", "sample_index", "unknown")),
        Segment("S5", nx.S5, "s5"),
        Segment("left mosaic end", nx.ME, "me"),
        Segment(insert_name, "N" * 42, placeholder=True),
        Segment("right mosaic end reverse complement", nx.ME_RC, "me"),
        Segment("S7 reverse complement", nx.S7_RC, "s7"),
        Segment("i7 reverse complement", "I" * UD_INDEX_NT, "cbc", placeholder=True,
                feature=feature("sample_index_i7", "sample_index", "unknown")),
        Segment("P7 reverse complement", il.P7_RC, "p7")], name=title)
    primers = (sp.NEXTERA["R1"], sp.NEXTERA["I1"], sp.NEXTERA["I2"], sp.NEXTERA["R2"])
    problems = sp.verify(lib, primers, required_roles=tuple(p.role for p in primers))
    if problems:
        raise ValueError(f"{title}: " + "; ".join(problems))
    return lib, primers


def _dna_prep_workflow(final: Construct) -> Workflow:
    tagged = _tagged_intermediate("genomic DNA insert")
    bead = _single("bead-bound", [Segment("bead-linked transposome", "[BLT]",
                                                   placeholder=True),
                                           Segment("genomic DNA", "N" * 36,
                                                   placeholder=True)],
                   label="BLT-bound genomic DNA")
    return _workflow("genomic DNA", _duplex("genomic DNA"), (
        ("Bind genomic DNA to bead-linked transposomes", "BLT-bound DNA", bead,
         "At inputs of 100 ng or more, the finite bead-linked complexes support yield normalization."),
        ("On-bead tagmentation", "bead-bound tagged fragments", tagged,
         TAG_PROGRAM + "; BLT fragments a set number of DNA molecules and transfers partial adapters."),
        ("Post-tagmentation wash", "washed tagged fragments", tagged,
         PTC_PROGRAM + "; supernatant discarded, then TWB washes on the magnet."),
        ("Input-dependent limited-cycle index PCR", "dual-index sequencing library",
         Scene.duplex(list(final), label="Illumina DNA Prep library"),
         "PCR adds the i7/i5 index adapters and cluster-generation sequences; lower input requires additional cycles."),
    ))

def _dna_prep_graph(d, library, primers):
    source = d.source
    graph = ProtocolGraph(d.title, authority=source)
    dna = graph.input("dna", "genomic DNA (1–500 ng in 30 µl)", MaterialKind.BIOLOGICAL,
                      state=_duplex("genomic DNA"), primary=True)
    beads = graph.input("blt_beads", "bead-linked transposome complexes (BLT)",
                        MaterialKind.SELECTION, inferred=True,
                        source=source + " Reagent-level BLT composition is proprietary.")
    tb1 = graph.input("tb1", "Tagmentation Buffer 1 (TB1)", MaterialKind.CONDITION)
    bound = graph.product("bound", "genomic DNA bound to bead-linked transposomes",
                          MaterialKind.MOLECULAR,
                          state=_single("bead-bound", [
                              Segment("bead-linked transposome", "[BLT]", placeholder=True),
                              Segment("genomic DNA", "N" * 36, placeholder=True)],
                              label="BLT-bound genomic DNA"))
    # Guide p.1: the bead complex fragments a set number of DNA molecules.  That
    # saturation limits yield; it is not a physical removal.  Excess DNA leaves in the
    # post-tagmentation supernatant/wash below.
    graph.react("bind_blt", "Bind genomic DNA to bead-linked transposomes in TB1",
                StepKind.BINDING,
                ((dna, InputRole.TEMPLATE), (beads, InputRole.SELECTION),
                 (tb1, InputRole.CONDITION)),
                (bound,))
    adapter = graph.input("transposome_adapter", "BLT partial Nextera adapter",
                          MaterialKind.ADAPTER, inferred=True,
                          source=source + " Public endpoint architecture is used; BLT reagent composition is proprietary.")
    tn5 = graph.input("tn5", "bead-linked transposase", MaterialKind.CATALYST)
    tag_condition = graph.input("tag_program", TAG_PROGRAM, MaterialKind.CONDITION)
    tagged = graph.product("tagged", "bead-bound tagged genomic fragments",
                           MaterialKind.MOLECULAR,
                           state=_tagged_intermediate("genomic DNA insert"))
    tn5_waste = graph.product("tn5_waste", "non-library or same-end tagmentation products",
                              MaterialKind.WASTE)
    graph.react("tagment", "Tagment genomic DNA on bead (55°C 15 min)",
                StepKind.TAGMENTATION,
                ((bound, InputRole.TEMPLATE), (adapter, InputRole.STRUCTURAL),
                 (tn5, InputRole.CATALYST), (tag_condition, InputRole.CONDITION)),
                (tagged, tn5_waste))
    graph.discard(tn5_waste, "does not form a heterologous amplifiable library fragment")
    tsb = graph.input("tsb", "Tagment Stop Buffer (TSB), " + PTC_PROGRAM,
                      MaterialKind.CONDITION)
    stopped = graph.product("stopped", "stopped bead-bound tagged fragments",
                            MaterialKind.MOLECULAR, state=tagged.state)
    graph.react("stop_tagmentation", "Stop tagmentation with TSB (37°C 15 min)",
                StepKind.QUENCHING, ((tagged, InputRole.TEMPLATE), (tsb, InputRole.CONDITION)),
                (stopped,), changes_molecule=False)
    wash = graph.input("wash_condition",
                       "magnet: discard supernatant; Tagment Wash Buffer (TWB) 100 µl × 3, last TWB kept until PCR",
                       MaterialKind.CONDITION)
    washed = graph.product("washed", "washed bead-bound tagged fragments",
                           MaterialKind.MOLECULAR, state=tagged.state)
    wash_waste = graph.product("wash_waste",
                               "supernatant and TWB: unbound excess DNA and reaction components",
                               MaterialKind.WASTE)
    graph.react("wash", "Wash bead-bound tagged fragments", StepKind.CLEANUP,
                ((stopped, InputRole.TEMPLATE), (wash, InputRole.CONDITION)),
                (washed, wash_waste), changes_molecule=False)
    graph.discard(wash_waste, "is removed with the post-tagmentation supernatant and TWB washes")
    # Guide p.27: 10 µl pre-paired i7 and i5 index adapters per sample from one well.
    i7_adapter = graph.input("i7_adapter", "i7 index adapter (pre-paired well, UD set)",
                             MaterialKind.OLIGO)
    i5_adapter = graph.input("i5_adapter", "i5 index adapter (same pre-paired well)",
                             MaterialKind.OLIGO)
    epm = graph.input("epm", "Enhanced PCR Mix (EPM)", MaterialKind.CATALYST)
    pcr_program = graph.input("pcr_program", PCR_PROGRAM, MaterialKind.CONDITION)
    amplified = graph.product("amplified", "amplified dual-index library (unpurified)",
                              MaterialKind.MOLECULAR, state=library)
    graph.react("index_pcr", "Amplify tagmented DNA with index adapters (BLT PCR)",
                StepKind.PCR,
                ((washed, InputRole.TEMPLATE), (i7_adapter, InputRole.STRUCTURAL),
                 (i5_adapter, InputRole.STRUCTURAL),
                 (epm, InputRole.CATALYST), (pcr_program, InputRole.CONDITION)),
                (amplified,), geometry=PCRGeometry.unavailable(
                    library,
                    reason="The guide refers index-adapter bases to Illumina Adapter Sequences and leaves the well choice to the user; the public final Nextera sites are exact, the pre-PCR bead-bound population and chosen well are not."))
    ipb = graph.input("ipb", "Illumina Purification Beads (IPB), double-sided: 45 µl IPB + 40 µl water, then 15 µl IPB to the supernatant",
                      MaterialKind.SELECTION)
    final = graph.product("library", "Illumina DNA Prep sequencing library",
                          MaterialKind.LIBRARY, state=library)
    large = graph.product("large_fraction", "first-bead fraction (larger fragments)",
                          MaterialKind.WASTE)
    small = graph.product("small_fraction", "second supernatant (short products and reaction components)",
                          MaterialKind.WASTE)
    graph.react("ipb_cleanup", "Double-sided IPB cleanup, 80% EtOH washes, elute in RSB",
                StepKind.CLEANUP, ((amplified, InputRole.TEMPLATE), (ipb, InputRole.SELECTION)),
                (final, large, small), changes_molecule=False)
    graph.discard(large, "is discarded with the first bead plate")
    graph.discard(small, "is discarded with the second-bind supernatant")
    graph.sequence(final, SequencingPlan(library=library, primers=tuple(primers)))
    graph.annotate(Highlight("i7_adapter", "Index kit variant",
        "Guide p.5 lists IDT for Illumina DNA/RNA UD, IDT for Illumina Nextera DNA UD and "
        "Nextera DNA CD index plates; p.26 says UD index codes are 10 bp while Nextera DNA CD "
        "codes are 8 bp. The schematic draws the 10 bp UD variant with symbolic well-specific "
        "bases; a CD-indexed library has 8 bp index windows."))
    graph.annotate(Highlight("index_pcr", "Input-dependent PCR cycles",
        "Guide p.26-27 BLT PCR table: " + ", ".join(f"{k} → {v} cycles" for k, v in PCR_CYCLES.items())
        + ". Libraries from <100 ng are not yield-normalized and must be quantified individually."))
    reviewed = ["dna", "blt_beads", "tb1", "bind_blt", "transposome_adapter", "tn5",
                "tag_program", "tagment", "tsb", "stop_tagmentation", "wash_condition",
                "wash", "i7_adapter", "i5_adapter", "epm", "pcr_program", "index_pcr", "ipb",
                "ipb_cleanup", "library"]
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()

DEFINITION = Definition(
        "Illumina DNA Prep", "01_illumina_dna_prep.html",
        'Illumina <a href="https://support.illumina.com/content/dam/illumina-support/documents/documentation/chemistry_documentation/illumina_prep/1000000025416_12_idp.pdf">Illumina DNA Prep Reference Guide 1000000025416 v12</a>.',
        "Bead-linked transposomes fragment and tag genomic DNA, TSB stops the reaction and washes remove unbound material, input-dependent limited-cycle PCR adds the index adapters, and a double-sided bead cleanup purifies unique-dual-index libraries.",
        "Illumina does not disclose the reagent-level bead/transposome composition. The molecular endpoint uses the public Nextera adapter and primer architecture with 10 bp UD index windows.",
        "A representative Illumina DNA Prep fragment after bead-linked tagmentation, wash, dual-index PCR and IPB cleanup.",
        "Standard Nextera primers read both genomic-insert ends and the unique dual sample indexes (10 bp index reads for UD index adapters). The guide does not fix the paired read lengths.",
        "nextera", "genomic DNA insert", _dna_prep_workflow)


def protocol_graph():
    return _dna_prep_graph(DEFINITION, FINAL_LIBRARY, SEQ_PRIMERS)

FINAL_LIBRARY, SEQ_PRIMERS = dna_prep_library(DEFINITION.insert_name, DEFINITION.title + " library")
d = DEFINITION
SPEC = SimpleNamespace(
    TITLE=d.title, NOTES=d.notes, SOURCE=d.source, SUMMARY=d.summary,
    CAVEAT=d.caveat, FINAL_LIBRARY=FINAL_LIBRARY, SEQ_PRIMERS=SEQ_PRIMERS,
    FINAL_CAPTION=d.final_caption, SEQUENCING_INTRO=d.sequencing_intro,
    workflow=lambda: _dna_prep_workflow(FINAL_LIBRARY),
    final_rows=lambda: nextera_final_rows(FINAL_LIBRARY),
    protocol_graph=protocol_graph,
)
globals().update(vars(SPEC))
