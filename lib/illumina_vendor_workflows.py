"""Authoritative commercial Illumina library-preparation workflows.

Product pages stay separate because their sample selection and adapter-installation
chemistries differ.  This module centralizes the recurring molecular states and final
TruSeq/Nextera constructs so those pages cannot silently drift apart.
"""
from __future__ import annotations

from dataclasses import dataclass
from types import SimpleNamespace

import nextera as nx
import illumina as il
from batch_ngs import (nextera_library, seg, truseq_full_length_adapter, truseq_library,
                       truseq_pcr_free_library, truseq_pcr_free_rows)
from chemdraw import Construct, MolecularState, Scene, Segment, Workflow
from endprep import (ThreePrimeOverhangAdapter, cohesive_ligation_geometry,
                     dA_tailed_scene, ligation_product_scene, repair_and_dA_tail)
from pcr import PCRGeometry
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind)


TRUSEQ_Y_END = ThreePrimeOverhangAdapter(
    "TruSeq Y adapter", il.STEM + "T", il.STEM_COMPLEMENT)

# The PCR-free workflow ligates the whole indexed adaptor, so its reagent geometry has to
# show both forked arms rather than only the stem and its 3′ dT.
TRUSEQ_FULL_LENGTH = truseq_full_length_adapter()
PCR_FREE_LIGATION_STATEMENT = (
    "The reference guide's ligation control joins “3′-T overhang adapters to 3′-A "
    "overhang inserts”, and Illumina Adapter Sequences prints both full-length arms, so "
    "the ligated adaptor already carries P5, i5, i7 and P7.")


def _duplex(name: str, *, bases: str = "N" * 36) -> Scene:
    return Scene.duplex([Segment(name, bases, placeholder=True)], label=name)


def _single(name: str, parts: list[Segment], *, label: str) -> Scene:
    sc = Scene()
    sc.strand(name, parts, label=label)
    sc.labels(name)
    return sc


def _workflow(initial_name: str, initial: Scene,
              steps: tuple[tuple[str, str, Scene, str], ...]) -> Workflow:
    w = Workflow(MolecularState(initial_name, tuple(initial.rows())))
    for action, name, state, note in steps:
        w.react(action, state.rows(), name=name, note=note)
    return w


def _truseq_final(insert_name: str, title: str, *, dual_index: bool = True):
    return truseq_library([seg(insert_name, "N" * 42, placeholder=True)], title,
                          dual_index=dual_index)


def _truseq_pcr_free_final(insert_name: str, title: str):
    return truseq_pcr_free_library([seg(insert_name, "N" * 42, placeholder=True)], title,
                                   adapter=TRUSEQ_FULL_LENGTH)


def _nextera_final(insert_name: str, title: str):
    return nextera_library([seg(insert_name, "N" * 42, placeholder=True)], title)


def truseq_final_rows(library: Construct) -> list:
    insert = next(s.name for s in library if "insert" in s.name or "cDNA" in s.name)
    return ligation_product_scene(
        library, insert, label=library.name, adapter_group="TruSeq adapter",
        highlight_components=False).rows()


def nextera_final_rows(library: Construct) -> list:
    sc = Scene.duplex(list(library), label=library.name)
    insert = next(s.name for s in library if "insert" in s.name or "DNA" in s.name)
    sc.junction("top", "left mosaic end", insert, "transposition boundary")
    sc.junction("top", insert, "right mosaic end reverse complement",
                "transposition boundary")
    return sc.rows()


def _pcr_free_workflow(final: Construct) -> Workflow:
    sheared = _duplex("Covaris-sheared genomic insert")
    repaired = _duplex("blunt size-selected genomic insert")
    end_prep = repair_and_dA_tail(Construct([
        Segment("size-selected genomic insert", "N" * 36, placeholder=True)]))
    tailed = dA_tailed_scene(end_prep, "dA-tailed genomic insert")
    w = _workflow("high-quality genomic DNA", _duplex("high-molecular-weight genomic DNA"), (
        ("Covaris shear", "fragmented DNA", sheared,
         "The guide specifies 350-bp or 550-bp insert workflows."),
        ("End-repair and bead size-selection", "blunt size-selected DNA", repaired,
         "End repair removes the 3-prime and 5-prime overhangs made by shearing."),
        ("Add one 3-prime dA to each end", "dA-tailed DNA", tailed,
         "The complementary adapter overhang is a single 3-prime T."),
    ))
    w.ligate_adapter(
        "Ligate full-length indexed TruSeq adapters",
        TRUSEQ_FULL_LENGTH.geometry(end_prep, statement=PCR_FREE_LIGATION_STATEMENT),
        truseq_pcr_free_rows(final, "genomic DNA insert",
                             label="PCR-free sequencing library"),
        name="PCR-free sequencing library",
        note="Sealing the adaptor finishes the library: no enrichment PCR follows, and qPCR "
             "measures the cluster-forming molecules.")
    return w


def _stranded_mrna_workflow(final: Construct) -> Workflow:
    selected = _single("RNA", [Segment("fragmented poly(A) RNA", "R" * 36,
                                       placeholder=True)], label="fragmented poly(A) RNA")
    first = _single("first strand", [Segment("antisense first-strand cDNA", "N" * 36,
                                              placeholder=True)], label="first-strand cDNA")
    marked = Scene.duplex([
        Segment("antisense first-strand cDNA", "N" * 36, placeholder=True)
    ], label="dUTP-marked cDNA duplex")
    marked.mark("bottom", "antisense first-strand cDNA'", "second strand contains dUTP")
    end_prep = repair_and_dA_tail(Construct([
        Segment("stranded cDNA insert", "N" * 36, placeholder=True)]))
    tailed = dA_tailed_scene(end_prep, "dA-tailed stranded cDNA")
    w = _workflow("total RNA", _single("RNA", [Segment("total RNA", "R" * 40,
                                                        placeholder=True)], label="total RNA"), (
        ("Select poly(A) RNA twice and fragment at elevated temperature",
         "fragmented poly(A) RNA", selected,
         "Poly-T magnetic beads enrich polyadenylated RNA; divalent cations fragment it."),
        ("Random-prime first-strand reverse transcription",
         "antisense first-strand cDNA", first,
         "Actinomycin D suppresses DNA-dependent synthesis during first-strand synthesis."),
        ("RNase H / DNA polymerase I second-strand synthesis with dUTP",
         "dUTP-marked duplex", marked,
         "Replacing dTTP with dUTP marks the second strand for exclusion during enrichment."),
        ("End-repair and add one 3-prime dA", "dA-tailed cDNA", tailed,
         "A-tailed cDNA accepts the complementary T-overhang adapters."),
    ))
    w.ligate_adapter(
        "Ligate TruSeq adapters and selectively PCR-enrich",
        cohesive_ligation_geometry(end_prep, TRUSEQ_Y_END,
            statement="The public TruSeq adapter has a 3′-T end complementary to dA-tailed cDNA."),
        ligation_product_scene(final, "antisense cDNA insert", label="stranded library",
                               adapter_group="TruSeq adapter").rows(),
        name="stranded sequencing library",
        note="The enrichment polymerase does not copy past dUTP, so Read 1 maps antisense and Read 2 sense.")
    return w


def _tagged_intermediate(label: str) -> Scene:
    parts = [Segment("S5", nx.S5, "s5"), Segment("left mosaic end", nx.ME, "me"),
             Segment(label, "N" * 36, placeholder=True),
             Segment("right mosaic end reverse complement", nx.ME_RC, "me"),
             Segment("S7 reverse complement", nx.S7_RC, "s7")]
    sc = Scene.duplex(parts, label="tagmented DNA")
    sc.junction("top", "left mosaic end", label, "transposition boundary")
    sc.junction("top", label, "right mosaic end reverse complement",
                "transposition boundary")
    return sc


def _nextera_xt_workflow(final: Construct) -> Workflow:
    tagged = _tagged_intermediate("small-genome / amplicon insert")
    return _workflow("1 ng double-stranded DNA", _duplex("small-genome or amplicon DNA"), (
        ("Solution-transposome tagmentation", "tagged fragments", tagged,
         "The engineered transposome fragments DNA and transfers partial adapters in one reaction."),
        ("Neutralize transposase", "neutralized tagged fragments", tagged,
         "Neutralization stops further tagmentation before amplification."),
        ("Limited-cycle PCR with two index primers", "dual-index sequencing library",
         Scene.duplex(list(final), label="Nextera XT library"),
         "PCR both amplifies the insert and completes P5, P7, i5 and i7 sequences."),
    ))


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
         "BLT fragments a set number of DNA molecules and transfers partial adapters."),
        ("Post-tagmentation wash", "washed tagged fragments", tagged,
         "The solid phase allows removal of unbound DNA and reaction components."),
        ("Input-dependent limited-cycle index PCR", "dual-index sequencing library",
         Scene.duplex(list(final), label="Illumina DNA Prep library"),
         "PCR completes the platform arms; lower input requires additional cycles."),
    ))


@dataclass(frozen=True)
class Definition:
    title: str
    notes: str
    source: str
    summary: str
    caveat: str
    final_caption: str
    sequencing_intro: str
    final_kind: str
    insert_name: str
    workflow_builder: object


DEFINITIONS = {
    "truseq-dna-pcr-free": Definition(
        "TruSeq DNA PCR-Free", "01_truseq_dna_pcr_free.html",
        'Illumina <a href="https://emea.support.illumina.com/content/dam/illumina-support/documents/documentation/chemistry_documentation/samplepreps_truseq/truseq-dna-pcr-free-workflow/truseq-dna-pcr-free-workflow-reference-1000000039279-00.pdf">TruSeq DNA PCR-Free Reference Guide 1000000039279 v00</a>.',
        "Covaris-sheared genomic DNA is end-repaired, bead-size-selected, dA-tailed and ligated directly to indexed TruSeq adapters without library enrichment PCR.",
        "The kit offers single, combinatorial-dual and unique-dual index components. The schematic shows the dual-index form; index identities are placeholders, not invented sequences.",
        "A cluster-capable dual-index TruSeq molecule made entirely by adapter ligation; there is no library PCR.",
        "Standard TruSeq primers read both genomic-insert ends and the two sample indexes.",
        "truseq_pcr_free", "genomic DNA insert", _pcr_free_workflow),
    "truseq-stranded-mrna": Definition(
        "TruSeq Stranded mRNA", "01_truseq_stranded_mrna.html",
        'Illumina <a href="https://emea.support.illumina.com/content/dam/illumina-support/documents/documentation/chemistry_documentation/samplepreps_truseq/truseq-stranded-mrna-workflow/truseq-stranded-mrna-workflow-reference-1000000040498-00.pdf">TruSeq Stranded mRNA Reference Guide 1000000040498 v00</a>.',
        "Poly(A) RNA is bead-selected and fragmented, random-prime copied, marked with dUTP in the second strand, adapter-ligated and selectively PCR-enriched.",
        "The guide supports single- and dual-index variants. The schematic uses the dual-index HS architecture and does not assign a particular index identity.",
        "Directional TruSeq cDNA library: Read 1 derives from antisense cDNA and Read 2 from the sense transcript strand.",
        "Standard TruSeq primers sequence the antisense-derived insert in Read 1, the sense-derived mate in Read 2 and both sample indexes.",
        "truseq", "antisense cDNA insert", _stranded_mrna_workflow),
    "nextera-xt-dna-library-prep": Definition(
        "Nextera XT DNA Library Prep", "01_nextera_xt.html",
        'Illumina <a href="https://support.illumina.com/content/dam/illumina-support/documents/documentation/chemistry_documentation/samplepreps_nextera/nextera-xt/15031942_07_nextera-xt-product-documentation.pdf">Nextera XT DNA Library Prep Product Documentation 15031942 v07</a>.',
        "A solution transposome fragments and tags 1 ng DNA; neutralization and limited-cycle PCR complete a dual-index Nextera library.",
        "The commercial transposome composition is proprietary. Public Nextera mosaic-end, primer-entry and sequencing-primer sequences are imported from the central canonical model.",
        "A representative heterologous S5/S7 product; same-end tagmentation products cannot bridge P5 to P7 and are not the sequenced library.",
        "Standard Nextera primers read the insert from both mosaic-end boundaries and read i5/i7 sample indexes.",
        "nextera", "small-genome / amplicon insert", _nextera_xt_workflow),
    "illumina-dna-prep": Definition(
        "Illumina DNA Prep", "01_illumina_dna_prep.html",
        'Illumina <a href="https://support.illumina.com/content/dam/illumina-support/documents/documentation/chemistry_documentation/illumina_prep/1000000025416_12_idp.pdf">Illumina DNA Prep Reference Guide 1000000025416 v12</a>.',
        "Bead-linked transposomes fragment and tag genomic DNA, a wash removes unbound material, and input-dependent limited-cycle PCR completes unique-dual-index libraries.",
        "Illumina does not disclose the reagent-level bead/transposome composition. The molecular endpoint uses the public Nextera adapter and primer architecture.",
        "A representative Illumina DNA Prep fragment after bead-linked tagmentation, wash and dual-index PCR.",
        "Standard Nextera primers read both genomic-insert ends and the unique dual sample indexes.",
        "nextera", "genomic DNA insert", _dna_prep_workflow),
}


def _sequence(graph, material, library, primers, endpoint):
    graph.sequence(material, SequencingPlan(
        library=library, primers=tuple(primers), limitation=""))


def _finish_pcr(graph, prefix, template, library, primers, endpoint, reviewed, *, label):
    forward = graph.input(prefix + "_forward", label + " index PCR forward primer",
                          MaterialKind.OLIGO)
    reverse = graph.input(prefix + "_reverse", label + " index PCR reverse primer",
                          MaterialKind.OLIGO)
    polymerase = graph.input(prefix + "_polymerase", label + " PCR polymerase",
                             MaterialKind.CATALYST)
    product = graph.product(prefix + "_library", label + " sequencing library",
                            MaterialKind.LIBRARY, state=library)
    graph.react(prefix + "_pcr", "Complete and index the " + label + " library",
                StepKind.PCR,
                ((template, InputRole.TEMPLATE), (forward, InputRole.STRUCTURAL),
                 (reverse, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)),
                (product,), geometry=PCRGeometry.unavailable(
                    library,
                    reason="The vendor workflow supports several index-component revisions; the canonical final library and sequencing sites are exact, while one universal PCR-primer pair is not."))
    _sequence(graph, product, library, primers, endpoint)
    reviewed.extend((prefix + "_forward", prefix + "_reverse", prefix + "_polymerase",
                     prefix + "_pcr", prefix + "_library"))
    return product


def _pcr_free_graph(d, library, primers):
    source = d.source
    graph = ProtocolGraph(d.title, authority=source)
    dna = graph.input("dna", "high-quality high-molecular-weight genomic DNA",
                      MaterialKind.BIOLOGICAL,
                      state=_duplex("high-molecular-weight genomic DNA"), primary=True)
    shear = graph.input("shear_condition", "Covaris shearing condition",
                        MaterialKind.CONDITION)
    fragments = graph.product("fragments", "Covaris-sheared genomic DNA",
                              MaterialKind.MOLECULAR,
                              state=_duplex("Covaris-sheared genomic insert"))
    graph.react("shear", "Covaris-shear genomic DNA", StepKind.FRAGMENTATION,
                ((dna, InputRole.TEMPLATE), (shear, InputRole.CONDITION)), (fragments,))
    repair = graph.input("repair_enzyme", "DNA end-repair enzyme mix",
                         MaterialKind.CATALYST)
    repaired = graph.product("repaired", "blunt repaired genomic fragments",
                             MaterialKind.MOLECULAR,
                             state=_duplex("blunt size-selected genomic insert"))
    graph.react("end_repair", "End-repair sheared DNA", StepKind.END_REPAIR,
                ((fragments, InputRole.TEMPLATE), (repair, InputRole.CATALYST)), (repaired,))
    beads = graph.input("size_selection_beads", "size-selection magnetic beads",
                        MaterialKind.SELECTION)
    selected = graph.product("selected", "size-selected blunt genomic inserts",
                             MaterialKind.MOLECULAR, state=repaired.state)
    rejected = graph.product("rejected_sizes", "off-size DNA fragments", MaterialKind.WASTE)
    graph.react("size_select", "Select the requested insert-size range",
                StepKind.SELECTION,
                ((repaired, InputRole.TEMPLATE), (beads, InputRole.SELECTION)),
                (selected, rejected))
    graph.discard(rejected, "falls outside the selected library insert range")
    tailing = graph.input("tailing_enzyme", "single-dA tailing enzyme",
                          MaterialKind.CATALYST)
    end_prep = repair_and_dA_tail(Construct([
        Segment("size-selected genomic insert", "N" * 36, placeholder=True)]))
    tailed = graph.product("tailed", "single-dA-tailed genomic inserts",
                           MaterialKind.MOLECULAR,
                           state=dA_tailed_scene(end_prep, "dA-tailed genomic insert"))
    graph.react("tail", "Add one 3-prime dA to each end", StepKind.TAILING,
                ((selected, InputRole.TEMPLATE), (tailing, InputRole.CATALYST)), (tailed,))
    adapter = graph.input("adapter", "full-length indexed TruSeq 3-prime-T adapters",
                          MaterialKind.ADAPTER)
    ligase = graph.input("ligase", "TruSeq adapter DNA ligase", MaterialKind.CATALYST)
    final_product = graph.product("library", "PCR-free indexed TruSeq library",
                                  MaterialKind.LIBRARY, state=library)
    geometry = TRUSEQ_FULL_LENGTH.geometry(
        end_prep, statement=PCR_FREE_LIGATION_STATEMENT)
    graph.react("ligate", "Ligate full-length indexed TruSeq adapters without enrichment PCR",
                StepKind.LIGATION,
                ((tailed, InputRole.TEMPLATE), (adapter, InputRole.STRUCTURAL),
                 (ligase, InputRole.CATALYST)), (final_product,), geometry=geometry)
    _sequence(graph, final_product, library, primers, d.sequencing_intro)
    reviewed = ["dna", "shear_condition", "shear", "repair_enzyme", "end_repair",
                "size_selection_beads", "size_select", "tailing_enzyme", "tail",
                "adapter", "ligase", "ligate", "library"]
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()


def _stranded_mrna_graph(d, library, primers):
    source = d.source
    graph = ProtocolGraph(d.title, authority=source)
    rna = graph.input("rna", "total RNA", MaterialKind.BIOLOGICAL,
                      state=_single("RNA", [Segment("total RNA", "R" * 40,
                                                   placeholder=True)], label="total RNA"),
                      primary=True)
    beads = graph.input("polyt_beads", "poly-T RNA-selection beads",
                        MaterialKind.SELECTION)
    selected = graph.product("selected", "poly(A)-selected RNA", MaterialKind.MOLECULAR,
                             state=_single("RNA", [Segment("poly(A) RNA", "R" * 36,
                                                         placeholder=True)], label="poly(A) RNA"))
    rejected = graph.product("non_polya", "non-polyadenylated RNA", MaterialKind.WASTE)
    graph.react("select_polya", "Select poly(A) RNA twice", StepKind.AFFINITY_CAPTURE,
                ((rna, InputRole.TEMPLATE), (beads, InputRole.SELECTION)),
                (selected, rejected))
    graph.discard(rejected, "is depleted by the poly-T selection")
    fragmentation = graph.input("fragmentation_condition", "elevated-temperature divalent-cation fragmentation condition",
                                MaterialKind.CONDITION)
    fragmented = graph.product("fragmented", "fragmented poly(A) RNA",
                               MaterialKind.MOLECULAR, state=selected.state)
    graph.react("fragment", "Fragment selected poly(A) RNA", StepKind.FRAGMENTATION,
                ((selected, InputRole.TEMPLATE), (fragmentation, InputRole.CONDITION)),
                (fragmented,))
    random_primer = graph.input("rt_primer", "random first-strand RT primer",
                                MaterialKind.OLIGO)
    rt = graph.input("reverse_transcriptase", "reverse transcriptase with actinomycin D",
                     MaterialKind.CATALYST)
    first = graph.product("first_strand", "antisense first-strand cDNA",
                          MaterialKind.MOLECULAR,
                          state=_single("first strand", [Segment("antisense first-strand cDNA",
                              "N" * 36, placeholder=True)], label="first-strand cDNA"))
    graph.react("reverse_transcribe", "Random-prime first-strand reverse transcription",
                StepKind.REVERSE_TRANSCRIPTION,
                ((fragmented, InputRole.TEMPLATE), (random_primer, InputRole.STRUCTURAL),
                 (rt, InputRole.CATALYST)), (first,))
    second_primer = graph.input("second_strand_primer", "second-strand synthesis primer",
                                MaterialKind.OLIGO)
    second_polymerase = graph.input("second_strand_polymerase", "RNase H / DNA polymerase I mix",
                                    MaterialKind.CATALYST)
    marked = graph.product("marked_duplex", "dUTP-marked cDNA duplex",
                           MaterialKind.MOLECULAR, state=first.state)
    graph.react("second_strand", "Synthesize the dUTP-marked second strand",
                StepKind.EXTENSION,
                ((first, InputRole.TEMPLATE), (second_primer, InputRole.STRUCTURAL),
                 (second_polymerase, InputRole.CATALYST)), (marked,))
    repair = graph.input("repair_enzyme", "cDNA end-repair enzyme mix",
                         MaterialKind.CATALYST)
    repaired = graph.product("repaired", "blunt stranded-cDNA duplexes",
                             MaterialKind.MOLECULAR, state=marked.state)
    graph.react("end_repair", "End-repair the cDNA duplex", StepKind.END_REPAIR,
                ((marked, InputRole.TEMPLATE), (repair, InputRole.CATALYST)), (repaired,))
    tailing = graph.input("tailing_enzyme", "single-dA tailing enzyme",
                          MaterialKind.CATALYST)
    end_prep = repair_and_dA_tail(Construct([
        Segment("stranded cDNA insert", "N" * 36, placeholder=True)]))
    tailed = graph.product("tailed", "dA-tailed stranded cDNA", MaterialKind.MOLECULAR,
                           state=dA_tailed_scene(end_prep, "dA-tailed stranded cDNA"))
    graph.react("tail", "Add one 3-prime dA to the cDNA ends", StepKind.TAILING,
                ((repaired, InputRole.TEMPLATE), (tailing, InputRole.CATALYST)), (tailed,))
    adapter = graph.input("adapter", "TruSeq 3-prime-T Y adapters", MaterialKind.ADAPTER)
    ligase = graph.input("ligase", "TruSeq adapter DNA ligase", MaterialKind.CATALYST)
    adapted = graph.product("adapted", "adapter-ligated directional cDNA",
                            MaterialKind.MOLECULAR, state=tailed.state)
    geometry = cohesive_ligation_geometry(
        end_prep, TRUSEQ_Y_END,
        statement="The public TruSeq adapter has a 3′-T end complementary to dA-tailed cDNA.")
    graph.react("ligate", "Ligate TruSeq adapters", StepKind.LIGATION,
                ((tailed, InputRole.TEMPLATE), (adapter, InputRole.STRUCTURAL),
                 (ligase, InputRole.CATALYST)), (adapted,), geometry=geometry)
    reviewed = ["rna", "polyt_beads", "select_polya", "fragmentation_condition", "fragment",
                "rt_primer", "reverse_transcriptase", "reverse_transcribe",
                "second_strand_primer", "second_strand_polymerase", "second_strand",
                "repair_enzyme", "end_repair", "tailing_enzyme", "tail", "adapter",
                "ligase", "ligate"]
    _finish_pcr(graph, "library", adapted, library, primers, d.sequencing_intro,
                reviewed, label="TruSeq stranded mRNA")
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()


def _nextera_xt_graph(d, library, primers):
    source = d.source
    graph = ProtocolGraph(d.title, authority=source)
    dna = graph.input("dna", "1-ng double-stranded small-genome or amplicon DNA",
                      MaterialKind.BIOLOGICAL,
                      state=_duplex("small-genome or amplicon DNA"), primary=True)
    transposome = graph.input("transposome", "solution transposome carrying partial Nextera adapters",
                              MaterialKind.ADAPTER, inferred=True,
                              source=source + " Commercial transposome composition is proprietary.")
    tn5 = graph.input("tn5", "engineered transposase", MaterialKind.CATALYST)
    tagged = graph.product("tagged", "heterologously tagged DNA fragments",
                           MaterialKind.MOLECULAR,
                           state=_tagged_intermediate("small-genome / amplicon insert"))
    unusable = graph.product("unusable_tn5_products", "same-end or non-library tagmentation products",
                             MaterialKind.WASTE)
    graph.react("tagment", "Fragment and tag DNA with the solution transposome",
                StepKind.TAGMENTATION,
                ((dna, InputRole.TEMPLATE), (transposome, InputRole.STRUCTURAL),
                 (tn5, InputRole.CATALYST)), (tagged, unusable))
    graph.discard(unusable, "cannot bridge the heterologous P5/P7 library arms")
    neutralizer = graph.input("neutralizer", "transposase-neutralization condition",
                              MaterialKind.CONDITION)
    neutral = graph.product("neutral", "neutralized tagged fragments",
                            MaterialKind.MOLECULAR, state=tagged.state)
    graph.react("neutralize", "Neutralize transposase before PCR", StepKind.QUENCHING,
                ((tagged, InputRole.TEMPLATE), (neutralizer, InputRole.CONDITION)),
                (neutral,))
    reviewed = ["dna", "transposome", "tn5", "tagment", "neutralizer", "neutralize"]
    _finish_pcr(graph, "library", neutral, library, primers, d.sequencing_intro,
                reviewed, label="Nextera XT")
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()


def _dna_prep_graph(d, library, primers):
    source = d.source
    graph = ProtocolGraph(d.title, authority=source)
    dna = graph.input("dna", "genomic DNA", MaterialKind.BIOLOGICAL,
                      state=_duplex("genomic DNA"), primary=True)
    beads = graph.input("blt_beads", "bead-linked transposome complexes",
                        MaterialKind.SELECTION, inferred=True,
                        source=source + " Reagent-level BLT composition is proprietary.")
    bound = graph.product("bound", "genomic DNA bound to bead-linked transposomes",
                          MaterialKind.MOLECULAR,
                          state=_single("bead-bound", [
                              Segment("bead-linked transposome", "[BLT]", placeholder=True),
                              Segment("genomic DNA", "N" * 36, placeholder=True)],
                              label="BLT-bound genomic DNA"))
    unbound = graph.product("unbound", "unbound excess genomic DNA", MaterialKind.WASTE)
    graph.react("bind_blt", "Bind genomic DNA to bead-linked transposomes",
                StepKind.BINDING,
                ((dna, InputRole.TEMPLATE), (beads, InputRole.SELECTION)),
                (bound, unbound))
    graph.discard(unbound, "is removed by the finite bead-linked transposome capacity")
    adapter = graph.input("transposome_adapter", "BLT partial Nextera adapter",
                          MaterialKind.ADAPTER, inferred=True,
                          source=source + " Public endpoint architecture is used; BLT reagent composition is proprietary.")
    tn5 = graph.input("tn5", "bead-linked transposase", MaterialKind.CATALYST)
    tagged = graph.product("tagged", "bead-bound tagged genomic fragments",
                           MaterialKind.MOLECULAR,
                           state=_tagged_intermediate("genomic DNA insert"))
    tn5_waste = graph.product("tn5_waste", "non-library or same-end tagmentation products",
                              MaterialKind.WASTE)
    graph.react("tagment", "Tagment genomic DNA on bead", StepKind.TAGMENTATION,
                ((bound, InputRole.TEMPLATE), (adapter, InputRole.STRUCTURAL),
                 (tn5, InputRole.CATALYST)), (tagged, tn5_waste))
    graph.discard(tn5_waste, "does not form a heterologous amplifiable library fragment")
    wash = graph.input("wash_condition", "post-tagmentation wash condition",
                       MaterialKind.CONDITION)
    washed = graph.product("washed", "washed bead-bound tagged fragments",
                           MaterialKind.MOLECULAR, state=tagged.state)
    wash_waste = graph.product("wash_waste", "removed unbound DNA and reaction components",
                               MaterialKind.WASTE)
    graph.react("wash", "Wash bead-bound tagged fragments", StepKind.CLEANUP,
                ((tagged, InputRole.TEMPLATE), (wash, InputRole.CONDITION)),
                (washed, wash_waste))
    graph.discard(wash_waste, "is removed during the post-tagmentation wash")
    reviewed = ["dna", "blt_beads", "bind_blt", "transposome_adapter", "tn5",
                "tagment", "wash_condition", "wash"]
    _finish_pcr(graph, "library", washed, library, primers, d.sequencing_intro,
                reviewed, label="Illumina DNA Prep")
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()


def get(key: str):
    d = DEFINITIONS[key]
    if d.final_kind == "truseq_pcr_free":
        library, primers = _truseq_pcr_free_final(d.insert_name, d.title + " library")
        rows = lambda: truseq_pcr_free_rows(library, d.insert_name,
                                            highlight_components=False)
    elif d.final_kind == "truseq":
        library, primers = _truseq_final(d.insert_name, d.title + " library")
        rows = lambda: truseq_final_rows(library)
    else:
        library, primers = _nextera_final(d.insert_name, d.title + " library")
        rows = lambda: nextera_final_rows(library)
    spec = SimpleNamespace(
        TITLE=d.title, NOTES=d.notes, SOURCE=d.source, SUMMARY=d.summary,
        CAVEAT=d.caveat, FINAL_LIBRARY=library, SEQ_PRIMERS=primers,
        FINAL_CAPTION=d.final_caption, SEQUENCING_INTRO=d.sequencing_intro,
        workflow=lambda: d.workflow_builder(library), final_rows=rows,
    )
    factories = {
        "truseq-dna-pcr-free": lambda: _pcr_free_graph(d, library, primers),
        "truseq-stranded-mrna": lambda: _stranded_mrna_graph(d, library, primers),
        "nextera-xt-dna-library-prep": lambda: _nextera_xt_graph(d, library, primers),
        "illumina-dna-prep": lambda: _dna_prep_graph(d, library, primers),
    }
    spec.protocol_graph = factories[key]
    return spec
