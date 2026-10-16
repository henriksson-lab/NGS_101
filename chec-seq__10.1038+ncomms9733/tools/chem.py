"""Protocol-owned source, oligos, workflow and reaction graph."""
from adapter_ligation import AdapterLigationGeometry
from chemdraw import Construct, Scene, Segment, feature
import illumina as il
import seqprimers as sp
from pcr import PCRGeometry
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            SourceClaim, SourceDiscrepancy, StepKind)
from protocol_model import da_tailed, molecule, role

CHEC_UNIVERSAL_ADAPTER = il.P5[:-4] + il.TRUSEQ_READ1

CHEC_INDEXED_ADAPTER = il.CANONICAL_ARM_READ2 + "N" * 6 + il.P7_RC

CHEC_P5_PCR = il.P5[:21]

CHEC_P7_PCR = il.P7[:22]

SAMPLE_INDEX = feature("orsi_illustrative_i7", "sample_index", "unknown")


def indexed_adapter_segments():
    return [Segment("Orsi Read 2 adapter arm", CHEC_INDEXED_ADAPTER[:33], "r2"),
            Segment("six-position sample index; oligo orientation", "N" * 6, "i7",
                    placeholder=True, feature=SAMPLE_INDEX),
            Segment("Orsi P7 complement", CHEC_INDEXED_ADAPTER[39:], "p7")]


def orsi_adapter_scene():
    # The universal terminal T is a one-base 3′ overhang, not part of the paired core.
    scene = Scene()
    scene.strand("universal", [Segment("Orsi P5/Read1 Y arm", CHEC_UNIVERSAL_ADAPTER[:-13]),
                                Segment("universal paired core", CHEC_UNIVERSAL_ADAPTER[-13:-1]),
                                Segment("adapter terminal dT", CHEC_UNIVERSAL_ADAPTER[-1])],
                 label="Orsi Universal Adapter", mod3="terminal phosphorothioate")
    bottom = [Segment("indexed paired core", CHEC_INDEXED_ADAPTER[:12]),
              Segment("indexed Read2 Y arm", CHEC_INDEXED_ADAPTER[12:33]),
              *indexed_adapter_segments()[1:]]
    scene.anneal("indexed", bottom, to="universal",
                 pair=("indexed paired core", "universal paired core"),
                 unpaired=("indexed Read2 Y arm", "six-position sample index; oligo orientation", "Orsi P7 complement"),
                 label="Orsi Indexed Adapter", mod5="/5Phos/", mod3="terminal phosphorothioate")
    scene.overhang("universal", "adapter terminal dT", polarity="3-prime", owner="adapter",
                   text="unpaired source terminal dT")
    return scene


ORSI_DA_LIGATION = AdapterLigationGeometry.known("cohesive", orsi_adapter_scene(),
    da_tailed("Orsi illustrative dA-tailed insert", "ChEC insert"),
    statement="Separate Orsi Basic Protocol3 supplies both protected Y-adapter oligos and the dA-tailed insert workflow. The universal adapter's terminal T is unpaired and complements the fragment dA; these exact oligos are an illustrated variant, not a demonstrated 2015 sample-index assignment.")


def chec_library():
    library = Construct([Segment("complete Orsi Universal Adapter; P5/Read1 overlap", CHEC_UNIVERSAL_ADAPTER, "r1"),
                         role("DNA released by two nearby protein-tethered MNase cuts"),
                         Segment("fragment-derived dA junction", "A"), *indexed_adapter_segments()],
                        name="ChEC readout illustrated with the separately documented Orsi adapter variant")
    primers = (sp.TRUSEQ["R1"], sp.TRUSEQ["I1"], sp.TRUSEQ["R2"])
    errors = sp.verify(library, primers, required_roles=tuple(p.role for p in primers))
    if errors:
        raise ValueError("Orsi-derived ChEC readout: " + "; ".join(errors))
    return library, primers


def oligos():
    return (
        f"TruSeq Universal Adapter — 5′-{CHEC_UNIVERSAL_ADAPTER[:-1]}*{CHEC_UNIVERSAL_ADAPTER[-1]}-3′",
        f"TruSeq Indexed Adapter — 5′-/5Phos/{CHEC_INDEXED_ADAPTER[:-1]}*{CHEC_INDEXED_ADAPTER[-1]}-3′",
        f"P5 library-enrichment primer — 5′-{CHEC_P5_PCR}-3′",
        f"P7 library-enrichment primer — 5′-{CHEC_P7_PCR}-3′",
    )


FINAL_LIBRARY, SEQ_PRIMERS = chec_library()
TITLE = 'ChEC-seq — calcium-triggered chromatin endogenous cleavage sequencing'
NOTES = '01_chec_seq.html'
SOURCE = '<a href="https://doi.org/10.1038/ncomms9733">Zentner et al., Nature Communications (2015)</a>; the exact displayed library oligos come from the Henikoff laboratory’s separate <a href="https://doi.org/10.1002/0471142727.mb2131s110">Orsi et al. stepwise protocol (2015)</a>.'
SUMMARY = 'A protein of interest is genetically fused to MNase. Digitonin admits calcium into living yeast, activating cleavage beside occupied binding sites; short two-cut fragments are selected and converted into a paired-end Illumina library.'
CAVEAT = 'ChEC-seq reports calcium-triggered cleavage by a mobile protein–MNase fusion, not an antibody-captured footprint. A free nuclear MNase control and early kinetic samples are essential because slower cleavage contains accessibility, enzyme-bias and transient-collision signal. The exact oligo strings shown here come from a separate same-laboratory Orsi procedure; the defining ChEC article does not itself print a sample-index assignment.'
FINAL_CAPTION = 'The selected insert is bounded by two MNase-created ends near a protein-binding event. End repair and dA-tailing erase the original overhang geometry but retain both genomic boundary coordinates; the illustrated Orsi variant has one six-base i7 index and no i5 index. This does not identify the exact adapters/index member in the defining 2015 ChEC experiment.'
READ_LENGTHS = {'Read 1': 25, 'Index 1 (i7)': 6, 'Read 2': 25}
SEQUENCING_INTRO = 'The defining data set used 25-base paired-end reads. Standard TruSeq Read 1 and Read 2 sites are computed here on the separately illustrated Orsi-derived construct. Its six-position index is shown with six Index1 cycles; this is not a recovered 2015 instrument mask or sample-index assignment.'
SEQUENCING_ENDING = ''
SEQUENCING_UNAVAILABLE = ''

def protocol_graph() -> ProtocolGraph:
    final = FINAL_LIBRARY
    primers = SEQ_PRIMERS
    source = ("Zentner et al., Nature Communications (2015); doi:10.1038/ncomms9733; "
              "separate illustrated library variant: Orsi et al., Current Protocols (2015)")
    graph = ProtocolGraph("ChEC-seq", authority=source)
    cells = graph.input(
        "fusion_cells", "living yeast expressing the chromatin-protein–MNase fusion",
        MaterialKind.BIOLOGICAL,
        state=molecule("living budding yeast expressing a protein–MNase fusion", (
            role("chromatin"), role("DNA-bound protein fused to 3xFLAG-MNase"))),
        primary=True)
    digitonin = graph.input("digitonin", "0.1% digitonin permeabilization condition",
                            MaterialKind.CONDITION)
    chromatin = graph.product(
        "permeabilized_chromatin", "accessible intracellular chromatin",
        MaterialKind.MOLECULAR,
        state=molecule("permeabilized cells", (role("accessible intracellular chromatin"),)))
    tethered_mnase = graph.product("tethered_mnase", "inactive protein-tethered MNase",
                                   MaterialKind.CATALYST)
    graph.react("permeabilize", "Permeabilize cells with digitonin",
                StepKind.PERMEABILIZATION,
                ((cells, InputRole.TEMPLATE), (digitonin, InputRole.CONDITION)),
                (chromatin, tethered_mnase))
    calcium = graph.input("calcium", "2 mM calcium activation condition",
                          MaterialKind.CONDITION)
    cut = graph.product(
        "cut_chromatin", "chromatin cut near occupied protein sites",
        MaterialKind.MOLECULAR,
        state=molecule("chromatin cut near occupied protein sites", (
            role("distal chromatin"), role("MNase cleavage end"),
            role("protein-protected binding region"), role("MNase cleavage end"),
            role("distal chromatin"))))
    graph.react("timed_digest", "Activate tethered MNase and remove timed aliquots",
                StepKind.DIGESTION,
                ((chromatin, InputRole.TEMPLATE), (tethered_mnase, InputRole.CATALYST),
                 (calcium, InputRole.CONDITION)), (cut,))
    quench = graph.input("quench_reagents", "EDTA/EGTA, SDS and high-salt quench",
                         MaterialKind.CONDITION)
    stopped = graph.product(
        "stopped_digest", "stopped ChEC digest", MaterialKind.MOLECULAR,
        state=molecule("stopped ChEC digest", (
            role("large background genomic DNA"), role("short two-cut protein-site fragment"))))
    graph.react("quench", "Quench each timed digestion aliquot", StepKind.QUENCHING,
                ((cut, InputRole.TEMPLATE), (quench, InputRole.CONDITION)), (stopped,))
    extraction = graph.input("extraction_condition", "nucleic-acid extraction and RNase condition",
                             MaterialKind.CONDITION)
    extracted = graph.product("extracted_dna", "extracted ChEC DNA",
                              MaterialKind.MOLECULAR, state=stopped.state)
    graph.react("extract", "Extract nucleic acid and digest RNA", StepKind.EXTRACTION,
                ((stopped, InputRole.TEMPLATE), (extraction, InputRole.CONDITION)),
                (extracted,))
    beads = graph.input("ampure_selection", "2.5:1 AMPure size-selection rule",
                        MaterialKind.SELECTION)
    selected = graph.product(
        "selected_fragments", "short ChEC DNA fragments", MaterialKind.MOLECULAR,
        state=molecule("size-selected ChEC DNA", (
            role("short fragment bounded by two MNase cleavage ends"),)))
    large = graph.product("large_genomic_dna", "large bead-bound genomic DNA",
                          MaterialKind.WASTE)
    graph.react("size_select", "Retain the small-DNA supernatant", StepKind.SELECTION,
                ((extracted, InputRole.TEMPLATE), (beads, InputRole.SELECTION)),
                (selected, large))
    graph.discard(large, "binds the AMPure beads while the short ChEC fragments remain in supernatant")
    reviewed = [
        "fusion_cells", "digitonin", "permeabilize", "calcium", "timed_digest",
        "quench_reagents", "quench", "extraction_condition", "extract",
        "ampure_selection", "size_select",
    ]
    # Source-owned Orsi variant. The 2015 ChEC short-fragment selection above and
    # these positive downstream cleanups retain different physical fractions.
    polish = graph.input("chec_repair_polymerases", "Orsi end-polish T4 DNA polymerase and Klenow fragment", MaterialKind.CATALYST)
    kinase = graph.input("chec_kinase", "Orsi T4 PNK in the same end-polish mix", MaterialKind.CATALYST)
    repaired = graph.product("chec_repaired", "End-polished/phosphorylated DNA; illustrated Orsi variant", MaterialKind.MOLECULAR)
    graph.react("chec_repair", "Combined end polishing and phosphorylation (Orsi variant)", StepKind.END_REPAIR,
                ((selected, InputRole.TEMPLATE), (polish, InputRole.CATALYST), (kinase, InputRole.CATALYST)), (repaired,))
    organic = graph.input("chec_organic_cleanup", "Orsi phenol/chloroform and S-300 cleanup", MaterialKind.CONDITION)
    clean = graph.product("chec_repair_clean", "Cleaned end-polished DNA", MaterialKind.MOLECULAR)
    repair_waste = graph.product("chec_repair_waste", "Material removed by Orsi end-polish cleanup", MaterialKind.WASTE)
    graph.react("chec_clean_repair", "Clean up after end polishing", StepKind.CLEANUP,
                ((repaired, InputRole.TEMPLATE), (organic, InputRole.CONDITION)), (clean, repair_waste))
    graph.discard(repair_waste, "Reported phenol/chloroform and S-300 cleanup; no exact molecular identities assigned.")
    tailing = graph.input("chec_tailing_enzyme", "Orsi Klenow exo−, source materials M0212M", MaterialKind.CATALYST)
    datp = graph.input("chec_datp", "Orsi dATP tailing condition", MaterialKind.CONDITION)
    tailed = graph.product("chec_tailed", "Orsi dA-tailed ChEC fragment", MaterialKind.MOLECULAR,
                           state=da_tailed("dA-tailed DNA", "ChEC insert"))
    graph.react("chec_tail", "Add fragment 3′ dA tails (Orsi variant)", StepKind.TAILING,
                ((clean, InputRole.TEMPLATE), (tailing, InputRole.CATALYST), (datp, InputRole.CONDITION)), (tailed,))
    tail_clean = graph.product("chec_tail_clean", "Cleaned dA-tailed DNA", MaterialKind.MOLECULAR, state=tailed.state)
    tail_waste = graph.product("chec_tail_waste", "Material removed by Orsi A-tail cleanup", MaterialKind.WASTE)
    graph.react("chec_clean_tail", "Clean up after A-tailing", StepKind.CLEANUP,
                ((tailed, InputRole.TEMPLATE), (organic, InputRole.CONDITION)), (tail_clean, tail_waste))
    graph.discard(tail_waste, "Reported phenol/chloroform and S-300 cleanup.")
    adapters = graph.input("chec_adapters", "Separately documented protected Orsi Y adapters", MaterialKind.ADAPTER,
                           state=orsi_adapter_scene().rows())
    ligase = graph.input("chec_adapter_ligase", "Orsi Rapid DNA ligase, Enzymatics L6030-HC-L", MaterialKind.CATALYST)
    adapted = graph.product("chec_adapted", "Y-adapter-ligated DNA; physical Y arms not fully constructed", MaterialKind.MOLECULAR)
    graph.react("chec_ligation", "Ligate the illustrated Orsi adapters", StepKind.LIGATION,
                ((tail_clean, InputRole.TEMPLATE), (adapters, InputRole.STRUCTURAL), (ligase, InputRole.CATALYST)),
                (adapted,), geometry=ORSI_DA_LIGATION)
    positive_beads = graph.input("chec_positive_beads", "Orsi positive AMPure cleanup; retain bead-bound DNA", MaterialKind.SELECTION)
    ligation_clean = graph.product("chec_ligation_clean", "Bead-purified adapter-ligated DNA", MaterialKind.MOLECULAR)
    ligation_waste = graph.product("chec_ligation_waste", "Unbound material from post-ligation cleanup", MaterialKind.WASTE)
    graph.react("chec_clean_ligation", "Retain bead-bound ligated DNA", StepKind.CLEANUP,
                ((adapted, InputRole.TEMPLATE), (positive_beads, InputRole.SELECTION)), (ligation_clean, ligation_waste))
    graph.discard(ligation_waste, "Source positive cleanup discards supernatant; distinct from initial ChEC negative selection.")
    p5 = graph.input("chec_p5_primer", "Orsi P5 enrichment primer", MaterialKind.OLIGO,
                     state=Construct([Segment("P5 primer", CHEC_P5_PCR)]))
    p7 = graph.input("chec_p7_primer", "Orsi P7 enrichment primer", MaterialKind.OLIGO,
                     state=Construct([Segment("P7 primer", CHEC_P7_PCR)]))
    polymerase = graph.input("chec_pcr_polymerase", "KAPA HiFi HotStart polymerase; defining ChEC and Orsi sources", MaterialKind.CATALYST)
    amplified = graph.product("chec_amplified", "PCR-amplified illustrated Orsi library", MaterialKind.MOLECULAR, state=final)
    graph.react("chec_pcr", "PCR-enrich the illustrated Orsi library", StepKind.PCR,
                ((ligation_clean, InputRole.TEMPLATE), (p5, InputRole.STRUCTURAL),
                 (p7, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)), (amplified,),
                geometry=PCRGeometry.unavailable(molecule("initial Y-adapter-ligated DNA", (role("Y-adapter-ligated DNA"),)),
                    reason="Both Orsi enrichment primers are source-known, but the complete initial Y-arm strand state and first PCR copies are not constructed here. A fully complementary amplified endpoint is not substituted as the physical pre-PCR template."))
    library = graph.product("chec_library", final.name, MaterialKind.LIBRARY, state=final)
    pcr_waste = graph.product("chec_pcr_waste", "Unbound material from post-PCR cleanup", MaterialKind.WASTE)
    graph.react("chec_clean_pcr", "Retain bead-bound amplified library", StepKind.CLEANUP,
                ((amplified, InputRole.TEMPLATE), (positive_beads, InputRole.SELECTION)), (library, pcr_waste))
    graph.discard(pcr_waste, "Source post-PCR positive AMPure cleanup; no consensus or mapping filter is represented as physical removal.")
    graph.sequence(library, SequencingPlan(final, tuple(primers), tuple(p.role for p in primers), READ_LENGTHS))
    reviewed.extend(("chec_repair_polymerases", "chec_kinase", "chec_repair", "chec_organic_cleanup", "chec_clean_repair",
        "chec_tailing_enzyme", "chec_datp", "chec_tail", "chec_clean_tail", "chec_adapters", "chec_adapter_ligase",
        "chec_ligation", "chec_positive_beads", "chec_clean_ligation", "chec_p5_primer", "chec_p7_primer",
        "chec_pcr_polymerase", "chec_pcr", "chec_clean_pcr", "chec_library"))
    graph.declare_reviewed_scope(reviewed, source)
    graph.annotate(SourceDiscrepancy(
        "fusion_cells", "Published correction: original Reb1 result table",
        (SourceClaim("Original Supplementary Data 1, Reb1 sheet",
                     "The initially published workbook sorted columns E, F, G, H and K independently of the other columns."),
         SourceClaim("2015 corrigendum, doi:10.1038/ncomms10264, sole correction paragraph",
                     "The authors replaced that workbook and state that the incorrect version was not used in the paper's analyses.")),
        "This page describes molecular preparation and uses no peak values from the original result table.",
        "The defining article links the original corrigendum, which identifies the affected sheet and columns.",
        "A historical downloaded result workbook needs its own version check; the correction does not identify an adapter or reaction change."))
    graph.annotate(SourceDiscrepancy(
        "timed_digest", "Published correction: DNA-shape interpretation",
        (SourceClaim("Defining article Fig. 7 and DNA-shape Methods",
                     "High- and low-scoring sites were motif-centered, whereas random controls were not."),
         SourceClaim("2017 corrigendum, doi:10.1038/ncomms15723, correction paragraphs",
                     "With motif-centered random controls, the reported analysis cannot establish shape-dependent, sequence-independent binding at low-scoring sites.")),
        "The schematic retains calcium-triggered MNase cleavage and makes no sequence-independent DNA-shape binding claim.",
        "The authors' correction withdraws that interpretation and states that it does not invalidate ChEC-seq as a method to identify binding sites.",
        "Cleavage kinetics alone do not establish the molecular basis of a particular binding event; the linked correspondence and reply are separate analysis sources."))
    return graph.require_publishable()
