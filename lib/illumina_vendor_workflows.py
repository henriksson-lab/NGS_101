"""Authoritative commercial Illumina library-preparation workflows.

Product pages stay separate because their sample selection and adapter-installation
chemistries differ.  This module centralizes the recurring molecular states and final
TruSeq/Nextera constructs so those pages cannot silently drift apart.
"""
from __future__ import annotations

from dataclasses import dataclass
from types import SimpleNamespace

import nextera as nx
from batch_ngs import nextera_library, seg, truseq_library
from chemdraw import Construct, MolecularState, Scene, Segment, Workflow
from endprep import dA_tailed_scene, ligation_product_scene, repair_and_dA_tail


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
    tailed = dA_tailed_scene(repair_and_dA_tail(
        Construct([Segment("size-selected genomic insert", "N" * 36,
                           placeholder=True)])), "dA-tailed genomic insert")
    return _workflow("high-quality genomic DNA", _duplex("high-molecular-weight genomic DNA"), (
        ("Covaris shear", "fragmented DNA", sheared,
         "The guide specifies 350-bp or 550-bp insert workflows."),
        ("End-repair and bead size-selection", "blunt size-selected DNA", repaired,
         "End repair removes the 3-prime and 5-prime overhangs made by shearing."),
        ("Add one 3-prime dA to each end", "dA-tailed DNA", tailed,
         "The complementary adapter overhang is a single 3-prime T."),
        ("Ligate indexed TruSeq adapters", "PCR-free sequencing library",
         ligation_product_scene(final, "genomic DNA insert",
                                label="PCR-free sequencing library",
                                adapter_group="TruSeq adapter"),
         "No enrichment PCR follows ligation; qPCR measures cluster-forming molecules."),
    ))


def _stranded_mrna_workflow(final: Construct) -> Workflow:
    selected = _single("RNA", [Segment("fragmented poly(A) RNA", "R" * 36,
                                       placeholder=True)], label="fragmented poly(A) RNA")
    first = _single("first strand", [Segment("antisense first-strand cDNA", "N" * 36,
                                              placeholder=True)], label="first-strand cDNA")
    marked = Scene.duplex([
        Segment("antisense first-strand cDNA", "N" * 36, placeholder=True)
    ], label="dUTP-marked cDNA duplex")
    marked.mark("bottom", "antisense first-strand cDNA'", "second strand contains dUTP")
    tailed = dA_tailed_scene(repair_and_dA_tail(
        Construct([Segment("stranded cDNA insert", "N" * 36, placeholder=True)])),
        "dA-tailed stranded cDNA")
    return _workflow("total RNA", _single("RNA", [Segment("total RNA", "R" * 40,
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
        ("Ligate TruSeq adapters and selectively PCR-enrich",
         "stranded sequencing library",
         ligation_product_scene(final, "antisense cDNA insert", label="stranded library",
                                adapter_group="TruSeq adapter"),
         "The enrichment polymerase does not copy past dUTP, so Read 1 maps antisense and Read 2 sense."),
    ))


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
        "truseq", "genomic DNA insert", _pcr_free_workflow),
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


def get(key: str):
    d = DEFINITIONS[key]
    if d.final_kind == "truseq":
        library, primers = _truseq_final(d.insert_name, d.title + " library")
        rows = lambda: truseq_final_rows(library)
    else:
        library, primers = _nextera_final(d.insert_name, d.title + " library")
        rows = lambda: nextera_final_rows(library)
    return SimpleNamespace(
        TITLE=d.title, NOTES=d.notes, SOURCE=d.source, SUMMARY=d.summary,
        CAVEAT=d.caveat, FINAL_LIBRARY=library, SEQ_PRIMERS=primers,
        FINAL_CAPTION=d.final_caption, SEQUENCING_INTRO=d.sequencing_intro,
        workflow=lambda: d.workflow_builder(library), final_rows=rows,
    )
