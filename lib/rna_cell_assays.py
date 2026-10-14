"""Shared molecular constructors for run-on, RNA probing and cell-state assays.

The public methods remain separate protocol pages.  These constructors centralize the
relationships that must not drift between them: a run-on label belongs at the nascent
RNA 3' end, an RT stop is upstream of a chemical adduct, genotype/lineage branches keep
their cellular identifiers, and a proximity report exists only after the two probe
oligos have been covalently joined.
"""
from __future__ import annotations

from dataclasses import dataclass
from types import SimpleNamespace

import seqprimers as sp
from adapter_ligation import AdapterLigationGeometry, reported_adapter_end
from batch_ngs import nextera_library, seg
from chemdraw import Construct, MolecularState, Scene, Segment, Workflow, feature, revcomp
from circular import CircularizationGeometry
from pcr import PCRGeometry
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            SourceCoverage, SourceItem, StepKind)


def token(name: str, *, tag: str | None = None, length: int = 18,
          molecular_feature=None, inferred: bool = False) -> Segment:
    return Segment(name, "X" * length, tag, placeholder=True,
                   feature=molecular_feature, inferred=inferred)


def unknown_region(name: str, *, tag: str | None = None, molecular_feature=None) -> Segment:
    """A source-established role whose exact sequence or length is not established."""
    return Segment(name, f"[{name}]", tag, placeholder=True,
                   feature=molecular_feature)


def single(label: str, parts: tuple[Segment, ...], *, note: str = "",
           junctions: tuple[tuple[str, str, str], ...] = ()) -> Scene:
    """One selectable strand with checked, named covalent junctions."""
    sc = Scene()
    sc.strand("molecule", list(parts), label=label)
    for left, right, caption in junctions:
        sc.junction("molecule", left, right, caption)
    if note:
        sc.note("molecule", note)
    sc.labels("molecule")
    return sc


def branches(label: str, rows: tuple[tuple[str, tuple[Segment, ...]], ...]) -> Scene:
    sc = Scene()
    for row, parts in rows:
        sc.strand(row, list(parts), label=label)
        sc.labels(row)
    return sc


def chain(initial_name: str, initial: Scene,
          transitions: tuple[tuple[str, Scene, str], ...]) -> Workflow:
    w = Workflow(MolecularState(initial_name, tuple(initial.rows())))
    for action, state, note in transitions:
        w.react(action, state.rows(), name=action, note=note)
    return w


def run_on(*, chromatin: bool) -> Workflow:
    """Nascent-RNA run-on; the label is structurally fixed at the extended 3' end."""
    substrate = "isolated chromatin" if chromatin else "isolated nuclei"
    labelled = "biotin-NTP run-on end" if chromatin else "BrU run-on end"
    selection = ("streptavidin selects biotin-labelled nascent RNA" if chromatin else
                 "anti-BrdU immunopurification is repeated three times")
    return chain(substrate, single(substrate, (token("engaged RNA polymerase"),
                                               token("nascent RNA"))), (
        ("Allow engaged polymerases to run on; prevent new initiation",
         single("extended nascent RNA", (token("pre-existing nascent RNA"),
                                          token(labelled, tag="umi", length=8))),
         "The nucleotide analogue is incorporated only at the RNA 3′ end extended in vitro."),
        ("Fragment RNA and select analogue-containing molecules",
         single("selected nascent RNA", (token("nascent RNA fragment"),
                                          token(labelled, tag="umi", length=8))), selection),
        ("Install the RNA-library adapter and reverse-transcribe",
         single("nascent-RNA cDNA", (token("library handle"),
                                      token("nascent-RNA-derived cDNA"))),
         "The genomic end of the read reports the position and orientation of engaged polymerase."),
        ("PCR to a sequencing-compatible library",
         single("run-on library", (token("left sequencing arm"),
                                    token("nascent-RNA-derived cDNA"),
                                    token("right sequencing arm"))),
         "The defining source does not expose one universal modern primer set for this library."),
    ))


def run_on_graph(*, chromatin: bool, definition, final: Construct) -> ProtocolGraph:
    """Typed run-on graph shared by GRO-seq and ChRO-seq."""
    source = f"{definition.citation}; doi:{definition.doi}"
    graph = ProtocolGraph(definition.title, authority=source)
    substrate_label = "isolated chromatin" if chromatin else "isolated nuclei"
    label = "biotin-NTP" if chromatin else "BrUTP"
    substrate = graph.input(
        "run_on_substrate", substrate_label, MaterialKind.BIOLOGICAL,
        state=single(substrate_label, (token("engaged RNA polymerase"),
                                       token("nascent RNA"))), primary=True)
    polymerase = graph.input("engaged_polymerase", "already-engaged RNA polymerase",
                             MaterialKind.CATALYST)
    analogue = graph.input("labelled_ntp", label + " run-on nucleotide",
                           MaterialKind.MOLECULAR)
    labelled = graph.product(
        "labelled_nascent_rna", "run-on-labelled nascent RNA",
        MaterialKind.MOLECULAR,
        state=single("extended nascent RNA", (
            token("pre-existing nascent RNA"),
            token(label + " run-on end", tag="umi", length=8))))
    graph.react("run_on", "Extend already-engaged nascent transcripts",
                StepKind.TRANSCRIPTION,
                ((substrate, InputRole.TEMPLATE),
                 (polymerase, InputRole.CATALYST),
                 (analogue, InputRole.STRUCTURAL)), (labelled,))

    fragmentation = graph.input("rna_fragmentation_reagent",
                                "RNA fragmentation reagent", MaterialKind.CONDITION)
    fragments = graph.product("labelled_rna_fragments", "fragmented labelled nascent RNA",
                              MaterialKind.MOLECULAR, state=labelled.state)
    graph.react("fragment_rna", "Fragment run-on-labelled RNA", StepKind.FRAGMENTATION,
                ((labelled, InputRole.TEMPLATE),
                 (fragmentation, InputRole.CONDITION)), (fragments,))
    selector = graph.input(
        "label_selector",
        "streptavidin" if chromatin else "anti-BrdU immunopurification reagent",
        MaterialKind.SELECTION)
    selected = graph.product(
        "selected_nascent_rna", "selected analogue-containing nascent RNA",
        MaterialKind.MOLECULAR,
        state=single("selected nascent RNA", (
            token("nascent RNA fragment"),
            token(label + " run-on end", tag="umi", length=8))))
    rejected = graph.product("unlabelled_rna", "unselected RNA", MaterialKind.WASTE)
    graph.react("select_labelled_rna", "Select analogue-containing RNA",
                StepKind.AFFINITY_CAPTURE,
                ((fragments, InputRole.TEMPLATE),
                 (selector, InputRole.SELECTION)), (selected, rejected))
    graph.discard(rejected, "lacks the run-on label required for affinity selection")

    adapter = graph.input("rna_adapter", "run-on RNA-library adapter",
                          MaterialKind.ADAPTER, inferred=True)
    ligase = graph.input("rna_ligase", "RNA ligase", MaterialKind.CATALYST)
    ligated = graph.product("adapter_ligated_rna", "adapter-ligated nascent RNA",
                            MaterialKind.MOLECULAR, state=selected.state,
                            inferred=True)
    graph.react("adapter_ligation", "Install the RNA-library adapter",
                StepKind.LIGATION,
                ((selected, InputRole.TEMPLATE), (adapter, InputRole.STRUCTURAL),
                 (ligase, InputRole.CATALYST)), (ligated,),
                geometry=AdapterLigationGeometry.unavailable(
                    "The defining historical implementation does not expose one universal modern RNA-adapter end."))
    rt_primer = graph.input("rt_primer", "adapter-complementary RT primer",
                            MaterialKind.OLIGO, inferred=True)
    rtase = graph.input("reverse_transcriptase", "reverse transcriptase",
                        MaterialKind.CATALYST)
    cdna = graph.product(
        "nascent_cdna", "nascent-RNA-derived cDNA", MaterialKind.MOLECULAR,
        state=single("nascent-RNA cDNA", (
            token("library handle"), token("nascent-RNA-derived cDNA"))),
        inferred=True)
    graph.react("reverse_transcription", "Reverse-transcribe selected nascent RNA",
                StepKind.REVERSE_TRANSCRIPTION,
                ((ligated, InputRole.TEMPLATE), (rt_primer, InputRole.STRUCTURAL),
                 (rtase, InputRole.CATALYST)), (cdna,))
    forward = graph.input("pcr_forward_primer", "library PCR forward primer",
                          MaterialKind.OLIGO, inferred=True)
    reverse = graph.input("pcr_reverse_primer", "library PCR reverse primer",
                          MaterialKind.OLIGO, inferred=True)
    pcr_enzyme = graph.input("pcr_polymerase", "library PCR polymerase",
                             MaterialKind.CATALYST)
    library = graph.product("sequencing_library", definition.title + " library",
                            MaterialKind.LIBRARY, state=final, inferred=True)
    graph.react("library_pcr", "PCR to a sequencing-compatible library", StepKind.PCR,
                ((cdna, InputRole.TEMPLATE), (forward, InputRole.STRUCTURAL),
                 (reverse, InputRole.STRUCTURAL),
                 (pcr_enzyme, InputRole.CATALYST)), (library,),
                geometry=PCRGeometry.unavailable(
                    final, reason="The defining source does not expose one universal modern primer pair."))
    graph.sequence(library, SequencingPlan(limitation=definition.endpoint))
    graph.declare_reviewed_scope((
        "run_on_substrate", "engaged_polymerase", "labelled_ntp", "run_on",
        "rna_fragmentation_reagent", "fragment_rna", "label_selector",
        "select_labelled_rna", "rna_adapter", "rna_ligase", "adapter_ligation",
        "rt_primer", "reverse_transcriptase", "reverse_transcription",
        "pcr_forward_primer", "pcr_reverse_primer", "pcr_polymerase",
        "library_pcr", "sequencing_library",
    ), source)
    return graph.require_publishable()


def icshape() -> Workflow:
    modified = token("NAI-N3 2-prime-O adduct", tag="umi", length=6)
    return chain("RNA in living cells", single("folded cellular RNA", (token("RNA"),)), (
        ("Acylate flexible RNA 2′-hydroxyls with NAI-N3",
         single("chemically probed RNA", (token("RNA upstream"), modified,
                                           token("RNA downstream"))),
         "NAI-N3 supplies both the structure-dependent adduct and an azide enrichment handle."),
        ("Click DIBO-biotin to the azide, fragment, and repair RNA 3′ ends",
         single("biotinylated RNA fragment", (token("RNA fragment"), modified,
                                               token("DIBO-biotin", tag="cbc", length=6))),
         "Copper-free click chemistry makes only NAI-N3-modified fragments selectable."),
        ("Ligate the preadenylated 3′ RNA adapter and enrich on streptavidin",
         single("adapter-ligated selected RNA", (token("RNA fragment"),
                                                  token("3-prime RNA adapter")),
                junctions=(("RNA fragment", "3-prime RNA adapter", "RNA ligation"),)),
         "T4 RNL2 truncated K227Q joins the repaired RNA end to the blocked adapter."),
        ("Reverse-transcribe until the adduct stops synthesis",
         single("barcoded truncated cDNA", (token("13-nt experimental barcode + UMI", tag="cbc", length=13),
                                             token("cDNA ending one base before the adduct"))),
         "The RT primer carries 9 random bases and a 4-base experimental barcode."),
        ("Circularize cDNA with CircLigase II and PCR",
         single("icSHAPE library", (token("left Solexa arm"),
                                     token("13-nt experimental barcode + UMI", tag="cbc", length=13),
                                     token("RT-stop cDNA"), token("right Solexa arm"))),
         "A read's first transcript base marks the RT stop; analysis shifts one nucleotide to the adduct."),
    ))


def icshape_graph(definition, final: Construct) -> ProtocolGraph:
    """Typed icSHAPE graph through click enrichment, RT stop and circularization."""
    source = f"{definition.citation}; doi:{definition.doi}"
    graph = ProtocolGraph(definition.title, authority=source)
    rna = graph.input("cellular_rna", "folded RNA in living cells",
                      MaterialKind.BIOLOGICAL,
                      state=single("folded cellular RNA", (token("RNA"),)),
                      primary=True)
    nai = graph.input("nai_n3", "NAI-N3 acylation reagent", MaterialKind.CONDITION)
    modified_state = single("chemically probed RNA", (
        token("RNA upstream"), token("NAI-N3 2-prime-O adduct", tag="umi", length=6),
        token("RNA downstream")))
    modified = graph.product("modified_rna", "NAI-N3-modified RNA",
                             MaterialKind.MOLECULAR, state=modified_state)
    graph.react("acylate", "Acylate flexible RNA 2-prime-hydroxyls",
                StepKind.CHEMICAL_CONVERSION,
                ((rna, InputRole.TEMPLATE), (nai, InputRole.CONDITION)), (modified,))

    dibo = graph.input("dibo_biotin", "DIBO-biotin click reagent",
                       MaterialKind.MOLECULAR)
    click_condition = graph.input("click_condition", "copper-free click condition",
                                  MaterialKind.CONDITION)
    biotinylated = graph.product(
        "biotinylated_rna", "DIBO-biotin-tagged modified RNA",
        MaterialKind.MOLECULAR,
        state=single("biotinylated RNA", (
            token("RNA"), token("NAI-N3 adduct", tag="umi", length=6),
            token("DIBO-biotin", tag="cbc", length=6))))
    graph.react("click_biotin", "Click DIBO-biotin to the NAI-N3 azide",
                StepKind.CHEMICAL_CONVERSION,
                ((modified, InputRole.TEMPLATE), (dibo, InputRole.STRUCTURAL),
                 (click_condition, InputRole.CONDITION)), (biotinylated,))
    fragmentation = graph.input("rna_fragmentation_reagent", "RNA fragmentation reagent",
                                MaterialKind.CONDITION)
    fragmented = graph.product("fragmented_rna", "fragmented biotinylated RNA",
                               MaterialKind.MOLECULAR, state=biotinylated.state)
    graph.react("fragment_rna", "Fragment the biotinylated RNA",
                StepKind.FRAGMENTATION,
                ((biotinylated, InputRole.TEMPLATE),
                 (fragmentation, InputRole.CONDITION)), (fragmented,))
    repair = graph.input("rna_end_repair_enzyme", "RNA 3-prime-end repair enzyme",
                         MaterialKind.CATALYST)
    repaired = graph.product("repaired_rna", "3-prime-end-repaired RNA fragments",
                             MaterialKind.MOLECULAR, state=fragmented.state)
    graph.react("repair_rna_end", "Repair RNA 3-prime ends", StepKind.END_REPAIR,
                ((fragmented, InputRole.TEMPLATE), (repair, InputRole.CATALYST)),
                (repaired,))

    adapter = graph.input("preadenylated_adapter",
                          "preadenylated blocked 3-prime RNA adapter",
                          MaterialKind.ADAPTER)
    ligase = graph.input("t4_rnl2_truncated", "T4 RNL2 truncated K227Q",
                         MaterialKind.CATALYST)
    ligated = graph.product(
        "adapter_ligated_rna", "3-prime-adapter-ligated RNA",
        MaterialKind.MOLECULAR,
        state=single("adapter-ligated RNA", (
            token("RNA fragment"), token("3-prime RNA adapter")),
            junctions=(("RNA fragment", "3-prime RNA adapter", "RNA ligation"),)))
    adapter_geometry = reported_adapter_end(
        "icSHAPE 3-prime RNA adapter", kind="single-strand", polarity="5-prime",
        terminal_sequence="AGATCGGAAGAGCGGTTCAG",
        statement="The protocol prints the 5′ preadenylated, 3′-blocked RNA adapter sequence used by truncated T4 RNL2.")
    graph.react("adapter_ligation", "Ligate the preadenylated 3-prime RNA adapter",
                StepKind.LIGATION,
                ((repaired, InputRole.TEMPLATE), (adapter, InputRole.STRUCTURAL),
                 (ligase, InputRole.CATALYST)), (ligated,), geometry=adapter_geometry)
    streptavidin = graph.input("streptavidin", "streptavidin capture reagent",
                               MaterialKind.SELECTION)
    selected = graph.product("selected_rna", "streptavidin-selected modified RNA",
                             MaterialKind.MOLECULAR, state=ligated.state)
    unselected = graph.product("unselected_rna", "unmodified or uncaptured RNA",
                               MaterialKind.WASTE)
    graph.react("streptavidin_capture", "Enrich biotinylated RNA",
                StepKind.AFFINITY_CAPTURE,
                ((ligated, InputRole.TEMPLATE),
                 (streptavidin, InputRole.SELECTION)), (selected, unselected))
    graph.discard(unselected, "not retained by streptavidin enrichment")

    rt_primer = graph.input("rt_primer", "barcoded RT-primer-1",
                            MaterialKind.OLIGO)
    rtase = graph.input("reverse_transcriptase", "reverse transcriptase",
                        MaterialKind.CATALYST)
    truncated = graph.product(
        "truncated_cdna", "barcoded cDNA ending before the adduct",
        MaterialKind.MOLECULAR,
        state=single("barcoded truncated cDNA", (
            token("13-nt experimental barcode + UMI", tag="cbc", length=13),
            token("cDNA ending one base before the adduct"))))
    graph.react("reverse_transcription", "Reverse-transcribe until the adduct stops synthesis",
                StepKind.REVERSE_TRANSCRIPTION,
                ((selected, InputRole.TEMPLATE), (rt_primer, InputRole.STRUCTURAL),
                 (rtase, InputRole.CATALYST)), (truncated,))
    circligase = graph.input("circligase_ii", "CircLigase II", MaterialKind.CATALYST)
    circular = graph.product("circular_cdna", "circularized RT-stop cDNA",
                             MaterialKind.MOLECULAR, state=truncated.state)
    graph.react("circularize", "Circularize cDNA", StepKind.CIRCULARIZATION,
                ((truncated, InputRole.TEMPLATE),
                 (circligase, InputRole.CATALYST)), (circular,),
                geometry=CircularizationGeometry.known(
                    "single-strand intramolecular ligation",
                    statement="The icSHAPE workflow circularizes RT-stop cDNA with CircLigase II."))
    forward = graph.input("pcr_forward_primer", "Solexa PCR forward primer",
                          MaterialKind.OLIGO, inferred=True)
    reverse = graph.input("pcr_reverse_primer", "Solexa PCR reverse primer",
                          MaterialKind.OLIGO, inferred=True)
    polymerase = graph.input("pcr_polymerase", "PCR polymerase", MaterialKind.CATALYST)
    library = graph.product("sequencing_library", "icSHAPE sequencing library",
                            MaterialKind.LIBRARY, state=final, inferred=True)
    graph.react("library_pcr", "PCR-amplify the circular cDNA", StepKind.PCR,
                ((circular, InputRole.TEMPLATE), (forward, InputRole.STRUCTURAL),
                 (reverse, InputRole.STRUCTURAL),
                 (polymerase, InputRole.CATALYST)), (library,),
                geometry=PCRGeometry.unavailable(
                    final, reason="The exact selected Solexa PCR primer pair is not fixed by the accessible primary text."))
    graph.sequence(library, SequencingPlan(limitation=definition.endpoint))
    graph.declare_reviewed_scope((
        "cellular_rna", "nai_n3", "acylate", "dibo_biotin", "click_condition",
        "click_biotin", "rna_fragmentation_reagent", "fragment_rna",
        "rna_end_repair_enzyme", "repair_rna_end", "preadenylated_adapter",
        "t4_rnl2_truncated", "adapter_ligation", "streptavidin",
        "streptavidin_capture", "rt_primer", "reverse_transcriptase",
        "reverse_transcription", "circligase_ii", "circularize",
        "pcr_forward_primer", "pcr_reverse_primer", "pcr_polymerase",
        "library_pcr", "sequencing_library",
    ), source)
    return graph.require_publishable()


def genotyping(kind: str) -> Workflow:
    if kind == "got":
        initial = single("10x first-strand cDNA", (
            token("target transcript cDNA"),
            token("10x cell barcode", tag="cbc", length=16,
                  molecular_feature=feature("cell_barcode", "cell_barcode", "whitelist",
                                            whitelist="10x Chromium chemistry-specific whitelist")),
            unknown_region("chemistry-specific UMI", tag="umi",
                  molecular_feature=feature("umi", "umi", "random")),))
        return chain("barcoded 10x cDNA", initial, (
            ("Optionally spike a target primer into RT or cDNA amplification", initial,
             "The paper specifies 5′ RT/additive spike-ins and a 3′ cDNA-amplification spike-in when needed."),
            ("Split amplified cDNA into expression and genotyping branches",
             branches("same cell barcode and UMI", (
                 ("expression branch", (token("whole-transcriptome cDNA"),)),
                 ("genotype branch", (
                     token("target cDNA"),
                     token("cell barcode", tag="cbc", length=16),
                     unknown_region("chemistry-specific UMI", tag="umi",
                                    molecular_feature=feature("umi", "umi", "random")),
                 )),
             )),
             "Only a small cDNA aliquot enters targeted genotyping; the remainder follows standard 10x preparation."),
            ("Target-PCR from upstream of the variant toward the 10x end",
             single("GoT amplicon", (token("variant-containing target"),
                                      token("cell barcode", tag="cbc", length=16),
                                      unknown_region("chemistry-specific UMI", tag="umi",
                                          molecular_feature=feature("umi", "umi", "random")))),
             "The amplicon deliberately retains the cell barcode and UMI."),
            ("Add sequencing handles and sequence with the expression library",
             single("GoT library", (token("sequencing arm"), token("variant-containing target"),
                                     token("cell barcode", tag="cbc", length=16),
                                     unknown_region("chemistry-specific UMI", tag="umi",
                                         molecular_feature=feature("umi", "umi", "random")),
                                     token("indexed sequencing arm"))),
             "Read 1 reports cellular identifiers and Read 2 reaches the selected variant."),
        ))
    if kind == "target":
        start = branches("one FACS-indexed cell", (
            ("genomic molecules", (token("genomic DNA"),)),
            ("RNA molecules", (token("poly(A) RNA"),))))
        preamp = branches("one-tube targeted RT-PCR", (
            ("gDNA amplicons", (token("targeted genomic loci"),)),
            ("cDNA amplicons", (token("targeted coding loci"),)),
            ("transcriptome", (token("Smart-seq+ cDNA"),))))
        return chain("single cell", start, (
            ("Lyse mildly, reverse-transcribe, and multiplex-pre-amplify", preamp,
             "Genomic and transcript-specific primers amplify mutation loci alongside whole-transcriptome cDNA."),
            ("Split the common pre-amplification product", branches("matched cell", (
                ("whole-transcriptome branch", (token("Smart-seq+ cDNA"),)),
                ("targeted genotype branch", (token("gDNA and cDNA target amplicons"),)))),
             "Both branches originate from the same physical well and FACS record."),
            ("Tagment the transcriptome; perform two nested PCRs on genotype targets",
             branches("two sequencing libraries", (
                 ("Nextera transcriptome", (token("Nextera arm"), token("transcript cDNA"), token("Nextera arm"))),
                 ("CS1/CS2 genotype", (token("CS1 site"), token("mutation locus"), token("CS2 site"))))),
             "PCR1 adds CS1/CS2 tails; PCR2 adds the cell-specific Illumina-compatible index."),
        ))
    raise ValueError(kind)


def lineage() -> Workflow:
    return chain("integrated scGESTALT transgene", single("lineage cassette", (
        token("DsRed transcript"), token("CRISPR target array"), token("poly(A)"))), (
        ("Induce Cas9 and guide pools at successive developmental times",
         single("edited lineage cassette", (token("DsRed transcript"),
                                              token("cumulative target-array indels"), token("poly(A)"))),
         "Repair-derived insertions and deletions form a permanent, cumulative lineage record."),
        ("Capture cassette RNA and transcriptome by inDrops",
         single("cell-indexed cDNA", (token("transcript or lineage cDNA"),
                                      token("cell-barcode part 1", tag="cbc", length=8),
                                      token("cell-barcode part 2", tag="cbc", length=8),
                                      token("UMI", tag="umi", length=6))),
         "The inDrops cell identifier is the join key between expression and lineage libraries."),
        ("IVT-amplify, then split fragmented and unfragmented RNA",
         branches("same indexed amplified RNA", (
             ("transcriptome", (token("fragmented amplified RNA"),)),
             ("lineage", (token("unfragmented lineage RNA"),)))),
         "Fragmented RNA makes the expression library; unfragmented RNA preserves the full target array."),
        ("Nested-PCR the target array and add matched Illumina indexes",
         single("scGESTALT lineage library", (token("lineage barcode array"),
                                               token("cell-barcode part 1", tag="cbc", length=8),
                                               token("cell-barcode part 2", tag="cbc", length=8),
                                               token("UMI", tag="umi", length=6), token("library index", tag="cbc", length=8))),
         "The lineage and transcriptome libraries use the same multiplexing index so cells can be matched."),
    ))


def proximity() -> Workflow:
    a = unknown_region("antibody A barcode", tag="cbc",
              molecular_feature=feature("protein_barcode_a", "feature_barcode", "fixed"))
    b = unknown_region("antibody B barcode", tag="cbc",
              molecular_feature=feature("protein_barcode_b", "feature_barcode", "fixed"))
    umi = unknown_region("PLA UMI", tag="umi",
                         molecular_feature=feature("pla_umi", "umi", "random"))
    return chain("cell surface", branches("DNA-conjugated antibodies", (
        ("probe A", (a, token("A connector half"))),
        ("probe B", (token("B connector half"), b, umi, token("poly(A)"))))), (
        ("Bind paired antibody–oligo probes to extracellular proteins",
         branches("two bound targets", (("protein A", (a, token("A connector half"))),
                                         ("protein B", (token("B connector half"), b, umi,
                                                        token("poly(A)"))))),
         "Probe identity is encoded in one fixed barcode on each antibody conjugate."),
        ("Hybridize the universal connector and ligate only proximal probes",
         single("119-base PLA reporter", (a, token("A connector half"),
                                           token("B connector half"), b, umi, token("poly(A)")),
                junctions=(("A connector half", "B connector half",
                            "connector-splinted proximity ligation"),)),
         "The universal connector splints the A/B nick but is not inserted. Only a "
         "ligated A–B product contains both primer sites and is amplifiable."),
        ("Co-capture polyadenylated PLA reporters and cellular mRNA",
         branches("same cell barcode", (("PLA branch", (a, b, umi)),
                                         ("RNA branch", (token("transcript cDNA"),)))),
         "Drop-seq, Smart-seq2 or 10x can supply the single-cell partition and RNA branch."),
        ("PCR each branch and add sequencing indexes",
         single("Prox-seq PLA library", (token("sequencing arm"), a, b, umi,
                                          token("indexed sequencing arm"))),
         "The PLA read reports both protein barcodes and a molecule identifier."),
    ))


def scgestalt_graph(definition, final: Construct) -> ProtocolGraph:
    """scGESTALT graph with matched transcriptome and lineage-library outcomes."""
    source = f"{definition.citation}; doi:{definition.doi}"
    graph = ProtocolGraph(definition.title, authority=source)
    cell = graph.input(
        "transgenic_cell", "cell with integrated expressed lineage cassette",
        MaterialKind.BIOLOGICAL,
        state=single("lineage cassette and transcriptome", (
            token("DsRed transcript"), token("CRISPR target array"),
            token("poly(A)"), token("cellular transcriptome"))), primary=True)
    guides = graph.input("guide_pools", "successive CRISPR guide pools",
                         MaterialKind.OLIGO)
    cas9 = graph.input("cas9", "inducible Cas9", MaterialKind.CATALYST)
    edited_cell = graph.product(
        "edited_cell", "cell with cumulatively edited lineage cassette",
        MaterialKind.MOLECULAR,
        state=single("edited lineage cassette and transcriptome", (
            token("DsRed transcript"), token("cumulative target-array indels"),
            token("poly(A)"), token("cellular transcriptome"))))
    graph.react("edit_lineage", "Induce cumulative lineage-cassette editing",
                StepKind.GENOME_EDITING,
                ((cell, InputRole.TEMPLATE), (guides, InputRole.STRUCTURAL),
                 (cas9, InputRole.CATALYST)), (edited_cell,))
    host_rnap = graph.input("host_rna_polymerase", "cellular RNA polymerase",
                            MaterialKind.CATALYST)
    cell_rna = graph.product(
        "cellular_rna", "polyadenylated edited-cassette RNA plus cellular mRNA",
        MaterialKind.MOLECULAR,
        state=single("polyadenylated cellular RNA", (
            token("edited lineage transcript"), token("cellular transcripts"))))
    graph.react("transcribe_cassette", "Transcribe and polyadenylate the edited cassette",
                StepKind.TRANSCRIPTION,
                ((edited_cell, InputRole.TEMPLATE),
                 (host_rnap, InputRole.CATALYST)), (cell_rna,))

    indrops_primer = graph.input("indrops_capture_primer",
                                 "inDrops cell-barcode/UMI oligo(dT) primer",
                                 MaterialKind.OLIGO)
    rtase = graph.input("reverse_transcriptase", "reverse transcriptase",
                        MaterialKind.CATALYST)
    first_strand = graph.product(
        "barcoded_first_strand", "cell-indexed lineage and transcriptome cDNA",
        MaterialKind.MOLECULAR,
        state=single("cell-indexed cDNA", (
            token("transcript or lineage cDNA"),
            token("cell-barcode part 1", tag="cbc", length=8),
            token("cell-barcode part 2", tag="cbc", length=8),
            token("UMI", tag="umi", length=6))))
    graph.react("reverse_transcription", "Capture lineage RNA and transcriptome by inDrops",
                StepKind.REVERSE_TRANSCRIPTION,
                ((cell_rna, InputRole.TEMPLATE),
                 (indrops_primer, InputRole.STRUCTURAL),
                 (rtase, InputRole.CATALYST)), (first_strand,))
    second_primer = graph.input("second_strand_primer", "inDrops second-strand primer",
                                MaterialKind.OLIGO, inferred=True)
    second_polymerase = graph.input("second_strand_polymerase", "second-strand polymerase",
                                    MaterialKind.CATALYST)
    duplex = graph.product("barcoded_cdna_duplex", "barcoded second-strand cDNA",
                           MaterialKind.MOLECULAR, state=first_strand.state)
    graph.react("second_strand", "Synthesize the second cDNA strand", StepKind.EXTENSION,
                ((first_strand, InputRole.TEMPLATE),
                 (second_primer, InputRole.STRUCTURAL),
                 (second_polymerase, InputRole.CATALYST)), (duplex,))
    t7 = graph.input("t7_rna_polymerase", "T7 RNA polymerase", MaterialKind.CATALYST)
    amplified_rna = graph.product("amplified_rna", "IVT-amplified barcoded RNA",
                                  MaterialKind.MOLECULAR, state=duplex.state)
    graph.react("ivt", "Linearly amplify barcoded cDNA by IVT", StepKind.TRANSCRIPTION,
                ((duplex, InputRole.TEMPLATE), (t7, InputRole.CATALYST)),
                (amplified_rna,))
    transcriptome_rna = graph.product("transcriptome_rna", "RNA for transcriptome library",
                                      MaterialKind.MOLECULAR,
                                      state=single("transcriptome branch", (
                                          token("amplified transcriptome RNA"),)))
    lineage_rna = graph.product("lineage_rna", "unfragmented lineage-cassette RNA",
                                MaterialKind.MOLECULAR,
                                state=single("lineage branch", (
                                    token("unfragmented lineage RNA"),)))
    graph.react("split_ivt", "Split amplified RNA into transcriptome and lineage branches",
                StepKind.SPLIT, ((amplified_rna, InputRole.TEMPLATE),),
                (transcriptome_rna, lineage_rna))

    fragmentation = graph.input("rna_fragmentation_reagent", "RNA fragmentation reagent",
                                MaterialKind.CONDITION)
    fragmented = graph.product("fragmented_transcriptome_rna", "fragmented amplified RNA",
                               MaterialKind.MOLECULAR, state=transcriptome_rna.state)
    graph.react("fragment_transcriptome", "Fragment the transcriptome branch",
                StepKind.FRAGMENTATION,
                ((transcriptome_rna, InputRole.TEMPLATE),
                 (fragmentation, InputRole.CONDITION)), (fragmented,))
    rna_adapter = graph.input("transcriptome_rna_adapter", "inDrops RNA-library adapter",
                              MaterialKind.ADAPTER, inferred=True)
    rna_ligase = graph.input("rna_ligase", "RNA ligase", MaterialKind.CATALYST)
    ligated_rna = graph.product("adapter_ligated_transcriptome_rna",
                                "adapter-ligated transcriptome RNA",
                                MaterialKind.MOLECULAR, state=fragmented.state,
                                inferred=True)
    graph.react("transcriptome_adapter_ligation", "Ligate the transcriptome RNA adapter",
                StepKind.LIGATION,
                ((fragmented, InputRole.TEMPLATE),
                 (rna_adapter, InputRole.STRUCTURAL),
                 (rna_ligase, InputRole.CATALYST)), (ligated_rna,),
                geometry=AdapterLigationGeometry.unavailable(
                    "The accessible scGESTALT protocol establishes the inDrops branch but not one exact adapter-end sequence."))
    transcriptome_rt_primer = graph.input("transcriptome_rt_primer",
                                          "transcriptome RT primer",
                                          MaterialKind.OLIGO, inferred=True)
    expression_cdna = graph.product("expression_cdna", "transcriptome-library cDNA",
                                    MaterialKind.MOLECULAR, state=ligated_rna.state,
                                    inferred=True)
    graph.react("transcriptome_rt", "Reverse-transcribe the fragmented transcriptome RNA",
                StepKind.REVERSE_TRANSCRIPTION,
                ((ligated_rna, InputRole.TEMPLATE),
                 (transcriptome_rt_primer, InputRole.STRUCTURAL),
                 (rtase, InputRole.CATALYST)), (expression_cdna,))
    pcr_polymerase = graph.input("pcr_polymerase", "library PCR polymerase",
                                 MaterialKind.CATALYST)
    expression_forward = graph.input("expression_forward_primer",
                                     "inDrops expression PCR forward primer",
                                     MaterialKind.OLIGO, inferred=True)
    expression_reverse = graph.input("expression_reverse_primer",
                                     "inDrops expression PCR reverse primer",
                                     MaterialKind.OLIGO, inferred=True)
    expression_final = custom_library("scGESTALT matched transcriptome library", (
        token("cell barcode", tag="cbc", length=16), token("UMI", tag="umi", length=6),
        token("transcript cDNA", length=44)))
    expression_library = graph.product(
        "expression_library", "matched scGESTALT transcriptome library",
        MaterialKind.LIBRARY, state=expression_final, inferred=True)
    graph.react("expression_pcr", "PCR the transcriptome branch",
                StepKind.PCR,
                ((expression_cdna, InputRole.TEMPLATE),
                 (expression_forward, InputRole.STRUCTURAL),
                 (expression_reverse, InputRole.STRUCTURAL),
                 (pcr_polymerase, InputRole.CATALYST)), (expression_library,),
                geometry=PCRGeometry.unavailable(
                    expression_final,
                    reason="The accessible protocol states the matched inDrops run but not one assembled primer pair."))
    graph.sequence(expression_library, SequencingPlan(
        limitation="The transcriptome branch follows the matched inDrops library implementation."))

    lineage_rt_primer = graph.input("lineage_random_hexamer", "lineage random-hexamer RT primer",
                                    MaterialKind.OLIGO)
    lineage_cdna = graph.product("lineage_cdna", "full target-array lineage cDNA",
                                 MaterialKind.MOLECULAR, state=lineage_rna.state)
    graph.react("lineage_rt", "Reverse-transcribe unfragmented lineage RNA",
                StepKind.REVERSE_TRANSCRIPTION,
                ((lineage_rna, InputRole.TEMPLATE),
                 (lineage_rt_primer, InputRole.STRUCTURAL),
                 (rtase, InputRole.CATALYST)), (lineage_cdna,))
    nested1_forward = graph.input("lineage_nested1_forward", "lineage nested-PCR-1 forward primer",
                                  MaterialKind.OLIGO)
    nested1_reverse = graph.input("lineage_nested1_reverse", "lineage nested-PCR-1 reverse primer",
                                  MaterialKind.OLIGO)
    nested = graph.product("lineage_nested_amplicon", "first nested lineage amplicon",
                           MaterialKind.MOLECULAR, state=lineage_cdna.state)
    graph.react("lineage_nested_pcr1", "First nested PCR of the lineage array",
                StepKind.PCR,
                ((lineage_cdna, InputRole.TEMPLATE),
                 (nested1_forward, InputRole.STRUCTURAL),
                 (nested1_reverse, InputRole.STRUCTURAL),
                 (pcr_polymerase, InputRole.CATALYST)), (nested,),
                geometry=PCRGeometry.unavailable(
                    final, reason="The page preserves published role geometry; exact assembled primer sites remain unavailable."))
    nested2_forward = graph.input("lineage_nested2_forward", "indexed lineage PCR forward primer",
                                  MaterialKind.OLIGO)
    nested2_reverse = graph.input("lineage_nested2_reverse", "indexed lineage PCR reverse primer",
                                  MaterialKind.OLIGO)
    lineage_library = graph.product("lineage_library", "scGESTALT lineage library",
                                    MaterialKind.LIBRARY, state=final, inferred=True)
    graph.react("lineage_nested_pcr2", "Second nested PCR adds matched Illumina indexes",
                StepKind.PCR,
                ((nested, InputRole.TEMPLATE),
                 (nested2_forward, InputRole.STRUCTURAL),
                 (nested2_reverse, InputRole.STRUCTURAL),
                 (pcr_polymerase, InputRole.CATALYST)), (lineage_library,),
                geometry=PCRGeometry.unavailable(
                    final, reason="Exact inDrops primer bases are role-labelled in the accessible protocol."))
    graph.sequence(lineage_library, SequencingPlan(
        limitation=definition.endpoint,
        read_lengths={"Read 1": 250, "Index 1": 8, "Index 2": 8, "Read 2": 14}))
    graph.declare_reviewed_scope((
        "transgenic_cell", "guide_pools", "cas9", "edit_lineage",
        "host_rna_polymerase", "transcribe_cassette", "indrops_capture_primer",
        "reverse_transcriptase", "reverse_transcription", "second_strand_primer",
        "second_strand_polymerase", "second_strand", "t7_rna_polymerase", "ivt",
        "split_ivt", "rna_fragmentation_reagent", "fragment_transcriptome",
        "transcriptome_rna_adapter", "rna_ligase", "transcriptome_adapter_ligation",
        "transcriptome_rt_primer", "transcriptome_rt", "pcr_polymerase",
        "expression_forward_primer", "expression_reverse_primer", "expression_pcr",
        "expression_library", "lineage_random_hexamer", "lineage_rt",
        "lineage_nested1_forward", "lineage_nested1_reverse", "lineage_nested_pcr1",
        "lineage_nested2_forward", "lineage_nested2_reverse", "lineage_nested_pcr2",
        "lineage_library",
    ), source)
    return graph.require_publishable()


def prox_seq_graph(definition, final: Construct) -> ProtocolGraph:
    """Prox-seq graph with covalent reporter and matched RNA-library branches."""
    source = f"{definition.citation}; doi:{definition.doi}"
    graph = ProtocolGraph(definition.title, authority=source)
    cell = graph.input("single_cell", "single cell with surface proteins and mRNA",
                       MaterialKind.BIOLOGICAL, state=proximity().initial, primary=True)
    antibodies = graph.input("antibody_probe_pairs", "antibody-conjugated A/B probe pairs",
                             MaterialKind.SELECTION)
    probe_oligos = graph.input("probe_oligos", "A/B protein-barcode probe oligos",
                               MaterialKind.OLIGO)
    bound = graph.product("probe_bound_cell", "cell with bound A/B antibody–oligo probes",
                          MaterialKind.MOLECULAR,
                          state=branches("two bound targets", (
                              ("protein A", (unknown_region("antibody A barcode", tag="cbc"),
                                             token("A connector half"))),
                              ("protein B", (token("B connector half"),
                                             unknown_region("antibody B barcode", tag="cbc"),
                                             unknown_region("PLA UMI", tag="umi"),
                                             token("poly(A)"))))))
    unbound = graph.product("unbound_probes", "unbound antibody–oligo probes",
                            MaterialKind.WASTE)
    graph.react("bind_probes", "Bind paired antibody–oligo probes to surface proteins",
                StepKind.BINDING,
                ((cell, InputRole.TEMPLATE), (antibodies, InputRole.SELECTION),
                 (probe_oligos, InputRole.STRUCTURAL)), (bound, unbound))
    graph.discard(unbound, "removed before proximity-reporter construction")
    connector = graph.input("universal_connector", "universal 20-base connector",
                            MaterialKind.OLIGO)
    hybridized = graph.product("connector_hybridized_probes",
                               "connector-hybridized proximal probe pair",
                               MaterialKind.MOLECULAR, state=bound.state)
    graph.react("hybridize_connector", "Hybridize the universal connector",
                StepKind.HYBRIDIZATION,
                ((bound, InputRole.TEMPLATE),
                 (connector, InputRole.STRUCTURAL)), (hybridized,))
    ligase = graph.input("dna_ligase", "DNA ligase", MaterialKind.CATALYST)
    reporter_state = single("119-base PLA reporter", (
        unknown_region("antibody A barcode", tag="cbc"),
        token("A connector half"), token("B connector half"),
        unknown_region("antibody B barcode", tag="cbc"),
        unknown_region("PLA UMI", tag="umi"), token("poly(A)")),
        junctions=(("A connector half", "B connector half",
                    "connector-splinted proximity ligation"),))
    ligated = graph.product("ligated_reporter_cell", "cell bearing ligated PLA reporters",
                            MaterialKind.MOLECULAR, state=reporter_state)
    graph.react("proximity_ligation", "Ligate only connector-bridged proximal probes",
                StepKind.LIGATION,
                ((hybridized, InputRole.TEMPLATE),
                 (connector, InputRole.STRUCTURAL),
                 (ligase, InputRole.CATALYST)), (ligated,),
                geometry=AdapterLigationGeometry.unavailable(
                    "The sources establish a 20-base connector-hybridized ligation region but no universal probe sequence across protein targets."))
    capture = graph.input("oligodt_capture_primer",
                          "platform-specific oligo(dT) cell-capture primer",
                          MaterialKind.OLIGO, inferred=True)
    rtase = graph.input("reverse_transcriptase", "reverse transcriptase",
                        MaterialKind.CATALYST)
    cell_cdna = graph.product("cell_cdna", "co-captured PLA-reporter and transcript cDNA",
                              MaterialKind.MOLECULAR,
                              state=branches("same cell identity", (
                                  ("PLA", (token("PLA reporter cDNA"),)),
                                  ("RNA", (token("transcript cDNA"),)))))
    graph.react("reverse_transcription", "Co-capture PLA reporters and cellular mRNA",
                StepKind.REVERSE_TRANSCRIPTION,
                ((ligated, InputRole.TEMPLATE), (capture, InputRole.STRUCTURAL),
                 (rtase, InputRole.CATALYST)), (cell_cdna,))
    pla_cdna = graph.product("pla_cdna", "PLA-reporter cDNA branch",
                             MaterialKind.MOLECULAR, state=reporter_state)
    expression_cdna = graph.product("expression_cdna", "matched cellular cDNA branch",
                                    MaterialKind.MOLECULAR,
                                    state=single("RNA branch", (token("transcript cDNA"),)))
    graph.react("split_cdna", "Split PLA and transcriptome cDNA branches",
                StepKind.SPLIT, ((cell_cdna, InputRole.TEMPLATE),),
                (pla_cdna, expression_cdna))
    polymerase = graph.input("pcr_polymerase", "library PCR polymerase",
                             MaterialKind.CATALYST)
    pla_forward = graph.input("pla_forward_primer", "PLA reporter PCR forward primer",
                              MaterialKind.OLIGO)
    pla_reverse = graph.input("pla_reverse_primer", "indexed PLA reporter PCR reverse primer",
                              MaterialKind.OLIGO)
    pla_library = graph.product("pla_library", "Prox-seq PLA library",
                                MaterialKind.LIBRARY, state=final, inferred=True)
    graph.react("pla_pcr", "PCR and index the PLA reporter branch", StepKind.PCR,
                ((pla_cdna, InputRole.TEMPLATE), (pla_forward, InputRole.STRUCTURAL),
                 (pla_reverse, InputRole.STRUCTURAL),
                 (polymerase, InputRole.CATALYST)), (pla_library,),
                geometry=PCRGeometry.unavailable(
                    final, reason="Read/index geometry depends on the selected single-cell back end."))
    graph.sequence(pla_library, SequencingPlan(limitation=definition.endpoint))
    expression_forward = graph.input("expression_forward_primer",
                                     "platform-specific expression PCR forward primer",
                                     MaterialKind.OLIGO, inferred=True)
    expression_reverse = graph.input("expression_reverse_primer",
                                     "platform-specific indexed expression PCR reverse primer",
                                     MaterialKind.OLIGO, inferred=True)
    expression_final = custom_library("Prox-seq matched RNA library", (
        token("cell identity"), token("transcript cDNA", length=44)))
    expression_library = graph.product("expression_library", "matched RNA library",
                                       MaterialKind.LIBRARY, state=expression_final,
                                       inferred=True)
    graph.react("expression_pcr", "PCR and index the matched RNA branch", StepKind.PCR,
                ((expression_cdna, InputRole.TEMPLATE),
                 (expression_forward, InputRole.STRUCTURAL),
                 (expression_reverse, InputRole.STRUCTURAL),
                 (polymerase, InputRole.CATALYST)), (expression_library,),
                geometry=PCRGeometry.unavailable(
                    expression_final,
                    reason="The method deliberately supports Drop-seq, Smart-seq2 and 10x back ends."))
    graph.sequence(expression_library, SequencingPlan(
        limitation="The matched RNA endpoint follows the chosen Drop-seq, Smart-seq2 or 10x implementation."))
    graph.declare_reviewed_scope((
        "single_cell", "antibody_probe_pairs", "probe_oligos", "bind_probes",
        "universal_connector", "hybridize_connector", "dna_ligase",
        "proximity_ligation", "oligodt_capture_primer", "reverse_transcriptase",
        "reverse_transcription", "split_cdna", "pcr_polymerase", "pla_forward_primer",
        "pla_reverse_primer", "pla_pcr", "pla_library", "expression_forward_primer",
        "expression_reverse_primer", "expression_pcr", "expression_library",
    ), source)
    return graph.require_publishable()


def find_seq() -> Workflow:
    return chain("heterogeneous cell suspension", single("one cell", (token("genome"),
                                                                        token("poly(A) RNA"))), (
        ("Co-encapsulate one cell, molten allyl agarose, and acrydite oligo(dT)",
         single("hydrogel bead precursor", (token("cell genome"), token("cell RNA"),
                                             token("covalent oligo(dT) matrix"))),
         "APS/TEMED links acrydite oligo(dT) to allyl agarose."),
        ("Lyse strongly and reverse-transcribe mRNA on the hydrogel",
         single("genome-capture bead", (token("retained genome"), token("matrix-bound cDNA"))),
         "Proteinase K/LiDS exposes DNA while the covalent capture matrix retains the transcriptome."),
        ("Re-encapsulate beads, run marker-specific TaqMan PCR, and sort droplets",
         single("fluorescent selected bead", (token("retained genome"), token("matrix-bound cDNA"),
                                               token("TaqMan marker signal"))),
         "Marker PCR is a selection measurement; it is not the transcriptome library."),
        ("Recover selected cDNA and build Smart-seq2 or Smart-seq3 libraries",
         single("FIND-seq library", (token("Nextera left arm"), token("selected-cell cDNA"),
                                      token("Nextera right arm"))),
         "The 2024 protocol supports bulk bins or individually sorted beads and UMI-aware Smart-seq3."),
    ))


def got_graph(definition, final: Construct) -> ProtocolGraph:
    """GoT graph retaining both expression and targeted-genotype outputs."""
    source = f"{definition.citation}; doi:{definition.doi}"
    graph = ProtocolGraph(definition.title, authority=source)
    initial = genotyping("got").initial
    cdna = graph.input("barcoded_10x_cdna", "amplified barcoded 10x cDNA",
                       MaterialKind.MOLECULAR, state=initial, primary=True)
    expression_cdna = graph.product("expression_cdna", "whole-transcriptome cDNA aliquot",
                                    MaterialKind.MOLECULAR,
                                    state=single("expression branch", (
                                        token("whole-transcriptome cDNA"),)))
    genotype_cdna = graph.product("genotype_cdna", "targeted-genotyping cDNA aliquot",
                                  MaterialKind.MOLECULAR,
                                  state=single("genotype branch", (
                                      token("target cDNA"),
                                      token("cell barcode", tag="cbc", length=16),
                                      unknown_region("chemistry-specific UMI", tag="umi"))))
    graph.react("split_cdna", "Split cDNA into expression and genotyping branches",
                StepKind.SPLIT, ((cdna, InputRole.TEMPLATE),),
                (expression_cdna, genotype_cdna))

    expression_forward = graph.input("expression_forward_primer",
                                     "10x expression-library forward primer",
                                     MaterialKind.OLIGO, inferred=True)
    expression_reverse = graph.input("expression_reverse_primer",
                                     "10x expression-library indexed reverse primer",
                                     MaterialKind.OLIGO, inferred=True)
    polymerase = graph.input("library_polymerase", "library PCR polymerase",
                             MaterialKind.CATALYST)
    expression_final = custom_library("GoT matched expression library", (
        token("whole-transcriptome cDNA", length=44),))
    expression_library = graph.product(
        "expression_library", "matched 10x expression library",
        MaterialKind.LIBRARY, state=expression_final, inferred=True)
    graph.react("expression_library_pcr", "Complete the standard expression library",
                StepKind.PCR,
                ((expression_cdna, InputRole.TEMPLATE),
                 (expression_forward, InputRole.STRUCTURAL),
                 (expression_reverse, InputRole.STRUCTURAL),
                 (polymerase, InputRole.CATALYST)), (expression_library,),
                geometry=PCRGeometry.unavailable(
                    expression_final,
                    reason="The exact primers depend on the selected 10x 3-prime or 5-prime chemistry."))
    graph.sequence(expression_library, SequencingPlan(
        limitation="The matched expression endpoint follows the chosen 10x chemistry."))

    target_forward = graph.input("target_forward_primer", "variant-upstream target primer",
                                 MaterialKind.OLIGO)
    target_reverse = graph.input("target_reverse_primer", "10x-end target primer",
                                 MaterialKind.OLIGO)
    amplicon = graph.product(
        "target_amplicon", "cell-barcode- and UMI-retaining GoT amplicon",
        MaterialKind.MOLECULAR,
        state=single("GoT amplicon", (
            token("variant-containing target"),
            token("cell barcode", tag="cbc", length=16),
            unknown_region("chemistry-specific UMI", tag="umi"))))
    graph.react("target_pcr", "Amplify from the variant toward the 10x end",
                StepKind.PCR,
                ((genotype_cdna, InputRole.TEMPLATE),
                 (target_forward, InputRole.STRUCTURAL),
                 (target_reverse, InputRole.STRUCTURAL),
                 (polymerase, InputRole.CATALYST)), (amplicon,),
                geometry=PCRGeometry.unavailable(
                    final, reason="Target-specific primer sequence varies with the selected locus."))
    handle_forward = graph.input("handle_forward_primer", "GoT handle PCR forward primer",
                                 MaterialKind.OLIGO, inferred=True)
    handle_reverse = graph.input("handle_reverse_primer", "GoT indexed handle PCR reverse primer",
                                 MaterialKind.OLIGO, inferred=True)
    genotype_library = graph.product("genotype_library", "GoT targeted genotype library",
                                     MaterialKind.LIBRARY, state=final, inferred=True)
    graph.react("handle_pcr", "Add sequencing handles and sample index", StepKind.PCR,
                ((amplicon, InputRole.TEMPLATE),
                 (handle_forward, InputRole.STRUCTURAL),
                 (handle_reverse, InputRole.STRUCTURAL),
                 (polymerase, InputRole.CATALYST)), (genotype_library,),
                geometry=PCRGeometry.unavailable(
                    final, reason="The endpoint varies with the selected 10x chemistry and target."))
    graph.sequence(genotype_library, SequencingPlan(limitation=definition.endpoint))

    graph.declare_reviewed_scope((
        "barcoded_10x_cdna", "split_cdna", "expression_forward_primer",
        "expression_reverse_primer", "library_polymerase", "expression_library_pcr",
        "expression_library", "target_forward_primer", "target_reverse_primer",
        "target_pcr", "handle_forward_primer", "handle_reverse_primer", "handle_pcr",
        "genotype_library",
    ), source)
    optional = SourceItem("optional_rt_spikein", "optional target-primer spike-in", source)
    graph.add_source_item(optional)
    graph.cover_source_item(SourceCoverage(
        optional.id, "excluded",
        reason="Optional target-primer spike-ins are chemistry- and locus-dependent; this graph begins at amplified barcoded 10x cDNA."))
    return graph.require_publishable()


def target_seq_graph(definition, final: Construct,
                     primers: tuple[sp.SeqPrimer, ...]) -> ProtocolGraph:
    """TARGET-seq graph with matched transcriptome and genotype libraries."""
    source = f"{definition.citation}; doi:{definition.doi}"
    graph = ProtocolGraph(definition.title, authority=source)
    cell = graph.input("single_cell", "one FACS-indexed cell", MaterialKind.BIOLOGICAL,
                       state=genotyping("target").initial, primary=True)
    rt_primers = graph.input("rt_primer_pool", "targeted and whole-transcriptome RT primers",
                             MaterialKind.OLIGO)
    rtase = graph.input("reverse_transcriptase", "reverse transcriptase",
                        MaterialKind.CATALYST)
    copied = graph.product("copied_cell", "cell lysate with cDNA and genomic targets",
                           MaterialKind.MOLECULAR,
                           state=branches("copied single-cell material", (
                               ("genomic", (token("genomic DNA"),)),
                               ("cDNA", (token("target and transcriptome cDNA"),)))))
    graph.react("reverse_transcription", "Reverse-transcribe cellular RNA",
                StepKind.REVERSE_TRANSCRIPTION,
                ((cell, InputRole.TEMPLATE), (rt_primers, InputRole.STRUCTURAL),
                 (rtase, InputRole.CATALYST)), (copied,))
    preamp_forward = graph.input("preamp_forward_pool", "multiplex preamplification forward pool",
                                 MaterialKind.OLIGO)
    preamp_reverse = graph.input("preamp_reverse_pool", "multiplex preamplification reverse pool",
                                 MaterialKind.OLIGO)
    polymerase = graph.input("pcr_polymerase", "preamplification/PCR polymerase",
                             MaterialKind.CATALYST)
    common = graph.product("common_preamp", "one-tube targeted and transcriptome preamplification",
                           MaterialKind.MOLECULAR,
                           state=branches("one-tube targeted RT-PCR", (
                               ("gDNA", (token("targeted genomic loci"),)),
                               ("cDNA", (token("targeted coding loci"),)),
                               ("transcriptome", (token("Smart-seq+ cDNA"),)))))
    graph.react("multiplex_preamp", "Multiplex-pre-amplify genomic and cDNA targets",
                StepKind.PCR,
                ((copied, InputRole.TEMPLATE),
                 (preamp_forward, InputRole.STRUCTURAL),
                 (preamp_reverse, InputRole.STRUCTURAL),
                 (polymerase, InputRole.CATALYST)), (common,),
                geometry=PCRGeometry.unavailable(
                    final, reason="The multiplex pool contains locus-specific pairs rather than one immutable amplicon."))
    transcriptome = graph.product("transcriptome_branch", "Smart-seq+ transcriptome branch",
                                  MaterialKind.MOLECULAR,
                                  state=single("transcriptome branch", (token("Smart-seq+ cDNA"),)))
    genotype = graph.product("genotype_branch", "targeted gDNA/cDNA amplicon branch",
                             MaterialKind.MOLECULAR,
                             state=single("genotype branch", (
                                 token("gDNA and cDNA target amplicons"),)))
    graph.react("split_preamp", "Split the common preamplification product",
                StepKind.SPLIT, ((common, InputRole.TEMPLATE),),
                (transcriptome, genotype))

    loaded = graph.input("loaded_nextera_adapters", "loaded Nextera XT adapters",
                         MaterialKind.ADAPTER)
    tn5 = graph.input("tn5", "Nextera XT transposase", MaterialKind.CATALYST)
    tagged = graph.product("tagged_transcriptome", "tagmented transcriptome cDNA",
                           MaterialKind.MOLECULAR, state=transcriptome.state)
    tn5_waste = graph.product("non_library_tn5_products", "non-library Tn5 products",
                              MaterialKind.WASTE)
    graph.react("tagment_transcriptome", "Tagment the transcriptome branch",
                StepKind.TAGMENTATION,
                ((transcriptome, InputRole.TEMPLATE),
                 (loaded, InputRole.STRUCTURAL), (tn5, InputRole.CATALYST)),
                (tagged, tn5_waste))
    graph.discard(tn5_waste, "not the productive asymmetric Nextera product")
    gap_fill = graph.input("gap_fill_polymerase", "Nextera gap-fill polymerase",
                           MaterialKind.CATALYST)
    filled = graph.product("gap_filled_transcriptome", "gap-filled tagmented transcriptome",
                           MaterialKind.MOLECULAR, state=tagged.state)
    graph.react("fill_tagmentation_gaps", "Fill the staggered tagmentation gaps",
                StepKind.EXTENSION,
                ((tagged, InputRole.TEMPLATE),
                 (gap_fill, InputRole.CATALYST)), (filled,))
    nx_forward = graph.input("nextera_forward_primer", "Nextera XT indexed forward primer",
                             MaterialKind.OLIGO)
    nx_reverse = graph.input("nextera_reverse_primer", "Nextera XT indexed reverse primer",
                             MaterialKind.OLIGO)
    transcriptome_final, transcriptome_primers = nextera_library(
        [token("transcript cDNA", length=44)], "TARGET-seq transcriptome library")
    transcriptome_library = graph.product(
        "transcriptome_library", "TARGET-seq transcriptome library",
        MaterialKind.LIBRARY, state=transcriptome_final)
    graph.react("transcriptome_pcr", "Index and amplify the transcriptome library",
                StepKind.PCR,
                ((filled, InputRole.TEMPLATE), (nx_forward, InputRole.STRUCTURAL),
                 (nx_reverse, InputRole.STRUCTURAL),
                 (polymerase, InputRole.CATALYST)), (transcriptome_library,),
                geometry=PCRGeometry.unavailable(
                    transcriptome_final,
                    reason="The Nextera XT index choice varies by cell; the endpoint geometry is canonical."))
    graph.sequence(transcriptome_library, SequencingPlan(
        transcriptome_final, tuple(transcriptome_primers),
        tuple(p.role for p in transcriptome_primers)))

    nested1_forward = graph.input("nested1_forward_pool", "first nested target-PCR forward pool",
                                  MaterialKind.OLIGO)
    nested1_reverse = graph.input("nested1_reverse_pool", "first nested target-PCR reverse pool",
                                  MaterialKind.OLIGO)
    nested = graph.product("nested_target_amplicons", "CS1/CS2-tailed target amplicons",
                           MaterialKind.MOLECULAR, state=genotype.state)
    graph.react("nested_pcr1", "First nested PCR adds CS1/CS2 tails", StepKind.PCR,
                ((genotype, InputRole.TEMPLATE),
                 (nested1_forward, InputRole.STRUCTURAL),
                 (nested1_reverse, InputRole.STRUCTURAL),
                 (polymerase, InputRole.CATALYST)), (nested,),
                geometry=PCRGeometry.unavailable(
                    final, reason="The multiplex target pool contains locus-specific nested primer pairs."))
    cs1 = graph.input("cs1_index_primer", "cell-indexed CS1 PCR primer",
                      MaterialKind.OLIGO)
    cs2 = graph.input("cs2_index_primer", "cell-indexed CS2 PCR primer",
                      MaterialKind.OLIGO)
    genotype_library = graph.product("genotype_library", "TARGET-seq genotype library",
                                     MaterialKind.LIBRARY, state=final)
    graph.react("nested_pcr2", "Second nested PCR adds the cell-specific index",
                StepKind.PCR,
                ((nested, InputRole.TEMPLATE), (cs1, InputRole.STRUCTURAL),
                 (cs2, InputRole.STRUCTURAL),
                 (polymerase, InputRole.CATALYST)), (genotype_library,),
                geometry=PCRGeometry.unavailable(
                    final, reason="The indexed CS1/CS2 amplification family varies by cell."))
    graph.sequence(genotype_library, SequencingPlan(
        final, tuple(primers), tuple(p.role for p in primers),
        read_lengths={"Read 1": 151, "Read 2": 151}))
    graph.declare_reviewed_scope((
        "single_cell", "rt_primer_pool", "reverse_transcriptase", "reverse_transcription",
        "preamp_forward_pool", "preamp_reverse_pool", "pcr_polymerase",
        "multiplex_preamp", "split_preamp", "loaded_nextera_adapters", "tn5",
        "tagment_transcriptome", "gap_fill_polymerase", "fill_tagmentation_gaps",
        "nextera_forward_primer", "nextera_reverse_primer",
        "transcriptome_pcr", "transcriptome_library", "nested1_forward_pool",
        "nested1_reverse_pool", "nested_pcr1", "cs1_index_primer", "cs2_index_primer",
        "nested_pcr2", "genotype_library",
    ), source)
    return graph.require_publishable()


def find_seq_graph(definition, final: Construct,
                   primers: tuple[sp.SeqPrimer, ...]) -> ProtocolGraph:
    """FIND-seq graph from hydrogel partition through selected-cell library."""
    source = f"{definition.citation}; doi:{definition.doi}"
    graph = ProtocolGraph(definition.title, authority=source)
    cells = graph.input("cell_suspension", "heterogeneous cell suspension",
                        MaterialKind.BIOLOGICAL, state=find_seq().initial, primary=True)
    matrix = graph.input("allyl_agarose", "molten allyl agarose matrix",
                         MaterialKind.CONDITION)
    capture_oligo = graph.input("acrydite_oligodt", "acrydite oligo(dT)",
                                MaterialKind.OLIGO)
    bead = graph.product("cell_hydrogel", "one-cell hydrogel bead precursor",
                         MaterialKind.MOLECULAR,
                         state=single("hydrogel bead precursor", (
                             token("cell genome"), token("cell RNA"),
                             token("covalent oligo(dT) matrix"))))
    empty = graph.product("empty_partition_material", "empty or multiplet partitions",
                          MaterialKind.WASTE)
    graph.react("encapsulate", "Encapsulate cells in oligo(dT) hydrogels",
                StepKind.PARTITION,
                ((cells, InputRole.TEMPLATE), (matrix, InputRole.CONDITION),
                 (capture_oligo, InputRole.STRUCTURAL)), (bead, empty))
    graph.discard(empty, "partitions without one usable cell are not selected for sequencing")
    proteinase = graph.input("proteinase_k", "proteinase K", MaterialKind.CATALYST)
    lids = graph.input("lids", "LiDS lysis condition", MaterialKind.CONDITION)
    lysed = graph.product("lysed_bead", "lysed genome-retaining hydrogel bead",
                          MaterialKind.MOLECULAR, state=bead.state)
    graph.react("lyse", "Lyse cells strongly inside the hydrogel",
                StepKind.DEPROTEINIZATION,
                ((bead, InputRole.TEMPLATE), (proteinase, InputRole.CATALYST),
                 (lids, InputRole.CONDITION)), (lysed,))
    rtase = graph.input("reverse_transcriptase", "reverse transcriptase",
                        MaterialKind.CATALYST)
    copied = graph.product(
        "genome_capture_bead", "genome-retaining bead with matrix-bound cDNA",
        MaterialKind.MOLECULAR,
        state=single("genome-capture bead", (
            token("retained genome"), token("matrix-bound cDNA"))))
    graph.react("reverse_transcription", "Reverse-transcribe mRNA on the hydrogel",
                StepKind.REVERSE_TRANSCRIPTION,
                ((lysed, InputRole.TEMPLATE),
                 (capture_oligo, InputRole.STRUCTURAL),
                 (rtase, InputRole.CATALYST)), (copied,))
    taqman = graph.input("taqman_assay", "marker-specific TaqMan assay",
                         MaterialKind.SELECTION)
    selected = graph.product(
        "selected_bead", "fluorescence-selected genome-capture bead",
        MaterialKind.MOLECULAR,
        state=single("fluorescent selected bead", (
            token("retained genome"), token("matrix-bound cDNA"),
            token("TaqMan marker signal"))))
    rejected = graph.product("rejected_beads", "marker-negative beads",
                             MaterialKind.WASTE)
    graph.react("marker_sort", "Run marker-specific TaqMan PCR and sort droplets",
                StepKind.AFFINITY_CAPTURE,
                ((copied, InputRole.TEMPLATE), (taqman, InputRole.SELECTION)),
                (selected, rejected))
    graph.discard(rejected, "marker-negative droplets are excluded by fluorescence sorting")
    recovery = graph.input("cdna_recovery_condition", "selected-bead cDNA recovery",
                           MaterialKind.CONDITION)
    recovered = graph.product("recovered_cdna", "selected-cell cDNA",
                              MaterialKind.MOLECULAR,
                              state=single("selected-cell cDNA", (token("full-length cDNA"),)))
    graph.react("recover_cdna", "Recover selected-cell cDNA", StepKind.EXTRACTION,
                ((selected, InputRole.TEMPLATE),
                 (recovery, InputRole.CONDITION)), (recovered,))
    preamp_forward = graph.input("smart_forward_primer", "Smart-seq amplification forward primer",
                                 MaterialKind.OLIGO, inferred=True)
    preamp_reverse = graph.input("smart_reverse_primer", "Smart-seq amplification reverse primer",
                                 MaterialKind.OLIGO, inferred=True)
    polymerase = graph.input("library_polymerase", "library amplification polymerase",
                             MaterialKind.CATALYST)
    amplified = graph.product("amplified_cdna", "Smart-seq-amplified selected-cell cDNA",
                              MaterialKind.MOLECULAR, state=recovered.state,
                              inferred=True)
    graph.react("smartseq_pcr", "Amplify selected cDNA by Smart-seq2 or Smart-seq3",
                StepKind.PCR,
                ((recovered, InputRole.TEMPLATE),
                 (preamp_forward, InputRole.STRUCTURAL),
                 (preamp_reverse, InputRole.STRUCTURAL),
                 (polymerase, InputRole.CATALYST)), (amplified,),
                geometry=PCRGeometry.unavailable(
                    final, reason="The protocol permits Smart-seq2 or UMI-aware Smart-seq3."))
    loaded = graph.input("loaded_nextera_adapters", "loaded Nextera adapters",
                         MaterialKind.ADAPTER)
    tn5 = graph.input("tn5", "Tn5 transposase", MaterialKind.CATALYST)
    tagged = graph.product("tagged_cdna", "tagmented selected-cell cDNA",
                           MaterialKind.MOLECULAR, state=amplified.state)
    tn5_waste = graph.product("non_library_tn5_products", "non-library Tn5 products",
                              MaterialKind.WASTE)
    graph.react("tagment", "Tagment selected-cell cDNA", StepKind.TAGMENTATION,
                ((amplified, InputRole.TEMPLATE), (loaded, InputRole.STRUCTURAL),
                 (tn5, InputRole.CATALYST)), (tagged, tn5_waste))
    graph.discard(tn5_waste, "not the productive asymmetric Nextera product")
    gap_fill = graph.input("gap_fill_polymerase", "Nextera gap-fill polymerase",
                           MaterialKind.CATALYST)
    filled = graph.product("gap_filled_cdna", "gap-filled tagmented selected-cell cDNA",
                           MaterialKind.MOLECULAR, state=tagged.state)
    graph.react("fill_tagmentation_gaps", "Fill the staggered tagmentation gaps",
                StepKind.EXTENSION,
                ((tagged, InputRole.TEMPLATE),
                 (gap_fill, InputRole.CATALYST)), (filled,))
    index_forward = graph.input("index_forward_primer", "Nextera indexed forward primer",
                                MaterialKind.OLIGO)
    index_reverse = graph.input("index_reverse_primer", "Nextera indexed reverse primer",
                                MaterialKind.OLIGO)
    library = graph.product("sequencing_library", "FIND-seq transcriptome library",
                            MaterialKind.LIBRARY, state=final)
    graph.react("index_pcr", "Index and amplify the FIND-seq library", StepKind.PCR,
                ((filled, InputRole.TEMPLATE),
                 (index_forward, InputRole.STRUCTURAL),
                 (index_reverse, InputRole.STRUCTURAL),
                 (polymerase, InputRole.CATALYST)), (library,),
                geometry=PCRGeometry.unavailable(
                    final, reason="Index choice varies; the final Nextera endpoint is validated."))
    graph.sequence(library, SequencingPlan(
        final, tuple(primers), tuple(p.role for p in primers)))
    graph.declare_reviewed_scope((
        "cell_suspension", "allyl_agarose", "acrydite_oligodt", "encapsulate",
        "proteinase_k", "lids", "lyse", "reverse_transcriptase",
        "reverse_transcription", "taqman_assay", "marker_sort",
        "cdna_recovery_condition", "recover_cdna", "smart_forward_primer",
        "smart_reverse_primer", "library_polymerase", "smartseq_pcr",
        "loaded_nextera_adapters", "tn5", "tagment", "gap_fill_polymerase",
        "fill_tagmentation_gaps", "index_forward_primer",
        "index_reverse_primer", "index_pcr", "sequencing_library",
    ), source)
    return graph.require_publishable()


def custom_library(name: str, payload: tuple[Segment, ...]) -> Construct:
    return Construct([token("Illumina-compatible left arm"), *payload,
                      token("Illumina-compatible right arm")], name=name)


TARGET_R1 = "ACACTGACGACATGGTTCTACA"
TARGET_R2 = "TACGGTAGCAGAGACTTGGTCT"


def target_genotype_library() -> tuple[Construct, tuple[sp.SeqPrimer, ...]]:
    lib = Construct([token("P5 and cell index"), Segment("CS1", TARGET_R1, "r1"),
                     token("targeted gDNA or cDNA locus", length=40),
                     Segment("CS2 reverse complement", revcomp(TARGET_R2), "r2"),
                     token("P7")], name="TARGET-seq genotyping library")
    primers = (sp.custom("Read 1", "TARGET-seq CS1-seq", TARGET_R1,
                         "Rodriguez-Meira et al. 2019 Key Resources"),
               sp.custom("Read 2", "TARGET-seq CS2-seq", TARGET_R2,
                         "Rodriguez-Meira et al. 2019 Key Resources"))
    problems = sp.verify(lib, primers, required_roles=("Read 1", "Read 2"))
    if problems:
        raise ValueError("TARGET-seq primer geometry: " + "; ".join(problems))
    return lib, primers


@dataclass(frozen=True)
class Definition:
    title: str
    doi: str
    citation: str
    notes: str
    summary: str
    endpoint: str
    caveat: str = ""


DEFS = {
    "gro-seq": Definition("GRO-seq — global run-on sequencing", "10.1126/science.1162228",
        "Core, Waterfall &amp; Lis, <i>Science</i> (2008)", "01_gro-seq.html",
        "A nuclear run-on incorporates BrUTP into short extensions made by already-engaged polymerases; repeated anti-BrdU purification selects those nascent RNAs for strand-specific sequencing.",
        "The original paper used 33-base single-end Illumina reads from the 5′ end of fragmented, BrU-labelled RNA; its historical library does not expose a current canonical primer set."),
    "chro-seq": Definition("ChRO-seq — chromatin run-on sequencing", "10.1038/s41588-018-0244-3",
        "Chu et al., <i>Nature Genetics</i> (2018)", "01_chro-seq.html",
        "Run-on transcription in isolated chromatin incorporates biotinylated nucleotides, allowing active polymerases to be mapped even in solid or degraded tissue.",
        "The final strand-specific Illumina read begins beside the run-on-labelled RNA 3′ end; exact sequencing oligos are implementation-specific in the defining paper."),
    "icshape": Definition("icSHAPE — in vivo click SHAPE", "10.1038/nprot.2016.011",
        "Flynn et al., <i>Nature Protocols</i> (2016)", "01_icshape.html",
        "NAI-N3 marks flexible RNA 2′-hydroxyls in cells, DIBO-biotin enables enrichment, and reverse-transcriptase stops encode nucleotide-resolution RNA structure.",
        "A 75-cycle single-end Illumina read begins with the 13-nt inline experimental barcode/UMI and then reports cDNA whose first mapped base locates the RT stop."),
    "got": Definition("GoT — Genotyping of Transcriptomes", "10.1038/s41586-019-1367-0",
        "Nam et al., <i>Nature</i> (2019)", "01_got.html",
        "A targeted amplicon branch is made from 10x cDNA while retaining its cell barcode and UMI, linking a selected transcribed variant to the cell's whole transcriptome.",
        "For 3′ v3, Read 1 uses 28 cycles for cell barcode and UMI, the 8-cycle index selects the library, and the 98- or 130-cycle Read 2 reaches the targeted variant.",
        "The defining paper supports 10x 3′ v2/v3 and 5′ inputs. The schematic shows the common identifier-preserving branch; exact UMI length follows the chosen 10x chemistry."),
    "target-seq": Definition("TARGET-seq — targeted genotyping and transcriptomics", "10.1016/j.molcel.2019.01.009",
        "Rodriguez-Meira et al., <i>Molecular Cell</i> (2019)", "01_target-seq.html",
        "Targeted gDNA and coding-DNA loci are co-amplified with Smart-seq+ whole-transcriptome cDNA from one FACS-indexed cell, then split into genotype and expression libraries.",
        "The genotype library uses published CS1/CS2 custom Read 1/Read 2 primers (151 cycles each) plus a 10-cycle cell index; its matched transcriptome branch is Nextera XT."),
    "scgestalt": Definition("scGESTALT — single-cell lineage recording", "10.1038/s41596-018-0058-x",
        "Raj &amp; Gagnon, <i>Nature Protocols</i> (2018)", "01_scgestalt.html",
        "Cumulative CRISPR indels in an expressed target array record ancestry; inDrops captures that lineage transcript alongside the cell transcriptome and supplies their shared cell identifier.",
        "The lineage run uses Read 1 = 250 cycles for the edited array, Index 1 = 8 cycles for cell-barcode part 1, Index 2 = 8 cycles for library index, and Read 2 = 14 cycles for cell-barcode part 2 plus 6-nt UMI."),
    "prox-seq": Definition("Prox-seq — proximity sequencing", "10.1038/s41592-022-01684-z",
        "Vistain et al., <i>Nature Methods</i> (2022)", "01_prox-seq.html",
        "Pairs of antibody-conjugated DNA probes ligate only when their extracellular protein targets are close; two protein barcodes and a UMI are recovered with the same cell's mRNA.",
        "The complete 119-base PLA reporter is amplified as a distinct library; the exact read/index layout follows the chosen Drop-seq, Smart-seq2 or 10x single-cell implementation.",
        "Prox-seq intentionally supports several scRNA-seq back ends. The schematic therefore records reporter structure and identifier roles without claiming one universal sequencing-primer geometry."),
    "find-seq": Definition("FIND-seq — nucleic-acid cytometry for rare cells", "10.1038/s41596-024-01021-y",
        "Shin et al., <i>Nature Protocols</i> (2024)", "01_find-seq.html",
        "A hydrogel bead retains one cell's genome and oligo(dT)-captured transcriptome while marker-specific TaqMan PCR drives fluorescent droplet sorting before whole-transcriptome library preparation.",
        "Recovered selected-cell cDNA is prepared by Smart-seq2 or UMI-aware Smart-seq3 and sequenced through its validated Nextera library ends."),
}


def get(key: str):
    d = DEFS[key]
    if key == "target-seq":
        final, primers = target_genotype_library()
        reads = {"Read 1": 151, "Read 2": 151}
    elif key == "find-seq":
        final, primers = nextera_library([token("selected-cell transcript cDNA", length=44)],
                                         "FIND-seq transcriptome library")
        reads = None
    elif key == "icshape":
        final = custom_library("icSHAPE library", (
            token("experimental barcode", tag="cbc", length=4,
                  molecular_feature=feature("experimental_barcode", "sample_index", "fixed")),
            token("UMI", tag="umi", length=9,
                  molecular_feature=feature("umi", "umi", "random")),
            token("RT-stop cDNA", length=40),))
        primers = (); reads = None
    elif key == "got":
        final = custom_library("GoT targeted library", (
            token("cell barcode", tag="cbc", length=16,
                  molecular_feature=feature("cell_barcode", "cell_barcode", "whitelist",
                                            whitelist="10x Chromium chemistry-specific whitelist")),
            unknown_region("chemistry-specific UMI", tag="umi",
                  molecular_feature=feature("umi", "umi", "random")),
            token("variant-containing cDNA", length=44),))
        primers = (); reads = None
    elif key == "scgestalt":
        final = custom_library("scGESTALT lineage library", (
            token("lineage barcode array", length=44),
            token("cell-barcode part 1", tag="cbc", length=8,
                  molecular_feature=feature("cell_barcode_1", "cell_barcode", "combinatorial",
                                            group="cell_barcode", part="part 1")),
            token("cell-barcode part 2", tag="cbc", length=8,
                  molecular_feature=feature("cell_barcode_2", "cell_barcode", "combinatorial",
                                            group="cell_barcode", part="part 2")),
            token("UMI", tag="umi", length=6,
                  molecular_feature=feature("umi", "umi", "random")),
            token("library index", tag="cbc", length=8,
                  molecular_feature=feature("sample_index", "sample_index", "fixed")),))
        primers = (); reads = None
    elif key == "prox-seq":
        final = custom_library("Prox-seq PLA library", (
            unknown_region("protein barcode A", tag="cbc",
                  molecular_feature=feature("protein_barcode_a", "feature_barcode", "fixed")),
            unknown_region("protein barcode B", tag="cbc",
                  molecular_feature=feature("protein_barcode_b", "feature_barcode", "fixed")),
            unknown_region("PLA UMI", tag="umi",
                  molecular_feature=feature("pla_umi", "umi", "random")),))
        primers = (); reads = None
    else:
        final = custom_library(d.title + " library", (token("nascent-RNA-derived cDNA", length=44),))
        primers = (); reads = None

    workflow_fn = {"gro-seq": lambda: run_on(chromatin=False),
                   "chro-seq": lambda: run_on(chromatin=True),
                   "icshape": icshape,
                   "got": lambda: genotyping("got"),
                   "target-seq": lambda: genotyping("target"),
                   "scgestalt": lineage,
                   "prox-seq": proximity,
                   "find-seq": find_seq}[key]
    oligos = ()
    if key == "icshape":
        oligos = (
            "<b>3′ adapter</b> — 5′-/5rApp/AGATCGGAAGAGCGGTTCAG/3ddC/-3′",
            "<b>RT-primer-1 architecture</b> — 5′-/5Phos/DDDNN AACC NNNN AGATCGGAAGAGCGTCGTGGA /iSp18/ GGATCC /iSp18/ TACTGAACCGC-3′",
        )
    elif key == "target-seq":
        oligos = ("<b>CS1-seq custom Read 1 primer</b> — 5′-A+CA+CTG+ACGACATGGTTCTACA-3′",
                   "<b>CS2-seq custom Read 2 primer</b> — 5′-T+AC+GGT+AGCAGAGACTTGGTCT-3′",
                   "The published ordering notation uses + before an LNA base; the model uses the corresponding unmodified binding sequences.")
    elif key == "find-seq":
        oligos = (
            "<b>Smart-seq2 bead capture</b> — 5′-acrydite-TTTTTTAAGCAGTGGTATCAACGCAGAGTAC(T)₃₀VN-3′",
            "<b>Smart-seq2 template switch</b> — 5′-AAGCAGTGGTATCAACGCAGAGTGAATrGrGrG-3′",
            "<b>Smart-seq3 bead capture</b> — 5′-acrydite-TTTTTTACGAGCATCAGCAGCATACGA(T)₃₀V-3′",
            "<b>Smart-seq3 template switch</b> — 5′-biotin-AGAGACAGATTGCGCAATGNNNNNNNNrGrGrG-3′; N₈ is the UMI",
        )
    unavailable = "" if primers else "Exact sequencing-primer sequences are not fixed by the defining protocol or were not exposed in its accessible primary text; read entry and cycle roles are stated explicitly."
    spec = SimpleNamespace(
        TITLE=d.title, NOTES=d.notes,
        SOURCE=f'Defining source: <a href="https://doi.org/{d.doi}">{d.citation}</a>.',
        SUMMARY=d.summary, CAVEAT=d.caveat, FINAL_LIBRARY=final,
        SEQ_PRIMERS=primers, FINAL_CAPTION=d.endpoint, SEQUENCING_INTRO=d.endpoint,
        READ_LENGTHS=reads, SEQUENCING_ENDING=d.endpoint,
        SEQUENCING_UNAVAILABLE=unavailable, workflow=workflow_fn,
        oligos=lambda: oligos,
    )
    if key == "gro-seq":
        spec.protocol_graph = lambda: run_on_graph(
            chromatin=False, definition=d, final=final)
    elif key == "chro-seq":
        spec.protocol_graph = lambda: run_on_graph(
            chromatin=True, definition=d, final=final)
    elif key == "icshape":
        spec.protocol_graph = lambda: icshape_graph(d, final)
    elif key == "got":
        spec.protocol_graph = lambda: got_graph(d, final)
    elif key == "target-seq":
        spec.protocol_graph = lambda: target_seq_graph(d, final, primers)
    elif key == "find-seq":
        spec.protocol_graph = lambda: find_seq_graph(d, final, primers)
    elif key == "scgestalt":
        spec.protocol_graph = lambda: scgestalt_graph(d, final)
    elif key == "prox-seq":
        spec.protocol_graph = lambda: prox_seq_graph(d, final)
    return spec
