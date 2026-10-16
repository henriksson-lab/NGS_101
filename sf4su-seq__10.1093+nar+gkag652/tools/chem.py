"""sf4sU-Seq source definition and size-selection route."""
from types import SimpleNamespace
import seqprimers as sp
from adapter_ligation import AdapterLigationGeometry
from chemdraw import Construct, Scene, Segment, feature, revcomp
from pcr import PCRGeometry
from protocol_graph import AnnotationSlot, Disposition, Highlight, SourceClaim, SourceDiscrepancy, Terminal, InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind
from protocol_model import molecule, role

# Supplementary Methods, "Kits and Reagents", identifies E7300S. The official
# NEB E7300 manual v9.0_07/25 prints these molecules on page 16.
# Study-specific index/revision assignment remains an open source question.
KIT_SOURCE = "NEBNext Small RNA Library Prep Kit E7300S; NEB manual v9.0_07/25, p. 16"
ADAPTER_3P = "AGATCGGAAGAGCACACGTCT"      # 5′-rApp…-NH2-3′
ADAPTER_5P = "GUUCAGAGUUCUACAGUCCGACGAUC"  # RNA, 5′→3′
RT_PRIMER = "AGACGTGTGCTCTTCCGATCT"
PCR_FORWARD = "AATGATACGGCGACCACCGAGATCTACACGTTCAGAGTTCTACAGTCCGA"
PCR_REVERSE_5P = "CAAGCAGAAGACGGCATACGAGAT"
PCR_REVERSE_3P = "GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT"


def kit_oligos():
    return (
        "NEBNext E7300 3′ SR adapter — 5′-rApp" + ADAPTER_3P + "-NH2-3′",
        "NEBNext E7300 5′ SR adapter — 5′-" + "".join("r" + base for base in ADAPTER_5P) + "-3′",
        "NEBNext E7300 SR RT primer — 5′-" + RT_PRIMER + "-3′",
        "NEBNext E7300 SR PCR primer — 5′-" + PCR_FORWARD[:-1] + "*" + PCR_FORWARD[-1] + "-3′",
        "NEBNext E7300 index PCR primer family — 5′-" + PCR_REVERSE_5P +
        "[six-base index]" + PCR_REVERSE_3P[:-1] + "*" + PCR_REVERSE_3P[-1] +
        "-3′; study-specific index/revision assignment remains open",
    )

DEFINITION = SimpleNamespace(
    title="sf4sU-Seq — size-fractionated 4sU sequencing", doi="10.1093/nar/gkag652",
    citation="Gerlevik et al., Nucleic Acids Research (2026)", notes="01_sf4su_seq.html",
    summary="Intact thio-labelled RNA is affinity-enriched, size-selected by denaturing gel and dialysis, decapped, and processed with the source-named small-RNA kit.",
    endpoint="For deposited HEK293 library GSM5452295 (SRX11439177), original SRA metadata declares SINGLE on NextSeq 500; all four archived lane runs contain 75-base R1 reads. This is a sample-specific genomic read length, not a complete instrument/index mask. RNA alignment 5′/3′ coordinates are computed from the reads.")
FINAL_LIBRARY = Construct([role("unassembled small-RNA-library end"), role("size-selected RNA-derived cDNA"),
                           role("unassembled opposite library end")], name="sf4sU library roles")
SPEC = SimpleNamespace(TITLE=DEFINITION.title, NOTES=DEFINITION.notes,
    SOURCE='<a href="https://doi.org/10.1093/nar/gkag652">Gerlevik et al., NAR (2026)</a>; original Supplementary Methods and NEB E7300 manual Version9.0_07/25.',
    SUMMARY=DEFINITION.summary, FINAL_LIBRARY=FINAL_LIBRARY, SEQ_PRIMERS=(), READ_LENGTHS={"Read 1": 75}, ACQUIRED_ROLES=("Read 1",),
    FINAL_CAPTION="Source-established library roles; the kit-guide oligos below are version-qualified, not an assigned study index or complete final assembly.",
    SEQUENCING_INTRO="", SEQUENCING_ENDING=DEFINITION.endpoint,
    SEQUENCING_UNAVAILABLE="Exact index member, technical/index-cycle acquisition and complete primer geometry remain open model questions; the represented deposited genomic read length is known.",
    READ_SPECS=(sp.read_spec("Read 1", 75, "Deposited GSM5452295 single genomic read; alignment endpoints analysed with bedtools -5 and -3. Technical/index cycles unassigned.", "include"),),
    oligos=kit_oligos,
    CAVEAT="The 2026 supplement names E7300S. The available official manual is Version9.0_07/25; its positive adapter/primer family is shown without claiming verified study revision or index member. This is an INFERRED study mapping of a source-positive dated kit reference, not proof that every enriched RNA has compatible ends or is recovered. The sf4sU branch does not inherit SNU bPAP/QuantSeq or TT sonication/directional-library choices. Intact label-containing RNA may contain pre-existing sequence; capture does not establish that every base was synthesized during the pulse.")


def rna_pool(mark):
    sc = Scene()
    sc.strand("label_containing", [role("RNA flank"), role(mark), role("RNA flank after incorporated site")], label="RNA with pulse-incorporated 4sU")
    sc.strand("unlabelled", [role("RNA without pulse-installed label")], label="RNA without pulse-installed label")
    sc.note("label_containing", "A label-containing molecule can include pre-existing RNA sequence; no whole-transcript age or complete incorporation is asserted.")
    return sc


def _recover(graph, prefix, template, condition_name, *, label="recovered RNA fraction"):

    condition = graph.input(prefix + "_condition", condition_name, MaterialKind.CONDITION)
    recovered = graph.product(prefix + "_rna", label, MaterialKind.MOLECULAR, state=template.state)
    removed = graph.product(prefix + "_removed", "RNA recovery removed fraction", MaterialKind.WASTE)
    graph.react(prefix, condition_name, StepKind.EXTRACTION,
                ((template, InputRole.TEMPLATE), (condition, InputRole.CONDITION)), (recovered, removed))
    graph.discard(removed, "is separated from recovered RNA during source purification")
    return recovered


def blocked_adapter_scene(*, inferred=True):
    sc = Scene()
    sc.strand("free_adapter", [Segment("3′ adapter", ADAPTER_3P, inferred=inferred)], label="free 3′ adapter", mod5="rApp", mod3="NH2")
    sc.anneal("rt_primer", [Segment("RT primer", RT_PRIMER, inferred=inferred)], to="free_adapter", pair=("RT primer", "3′ adapter"), label="SR RT primer")
    return sc


def _rna_ligation_geometry(adapter: str, substrate: Construct, *,
                           adapter_end: str, substrate_end: str,
                           statement: str, inferred: bool = True) -> AdapterLigationGeometry:
    scene = Scene()
    scene.strand("substrate", list(substrate), label="selected 4sU RNA")
    scene.strand("adapter", [Segment("kit adapter", adapter, inferred=inferred)],
                 label="NEBNext E7300 adapter",
                 mod5="rApp" if adapter_end == "5-prime" else "",
                 mod3="NH2" if adapter_end == "5-prime" else "")
    scene.ligation_end("substrate", substrate_end, kind="single-strand",
                       polarity="3-prime" if adapter_end == "5-prime" else "5-prime",
                       owner="substrate")
    scene.ligation_end("adapter", "kit adapter", kind="single-strand",
                       polarity=adapter_end, owner="adapter")
    return AdapterLigationGeometry.known("single-strand", scene, statement=statement)


def protocol_graph():
    source = f"{DEFINITION.citation}; doi:{DEFINITION.doi}"
    graph = ProtocolGraph(DEFINITION.title, authority=source)
    cells = graph.input(
        "cells", "HEK293 cells during the 4sU pulse (deposited readout GSM5452295)", MaterialKind.BIOLOGICAL,
        state=molecule("living cells", (role("RNA synthesis during pulse"),)),
        primary=True)
    four_su = graph.input("four_su", "4-thiouridine pulse-label nucleotide",
                          MaterialKind.MOLECULAR)
    polymerase = graph.input("cellular_rna_polymerase", "cellular RNA polymerases",
                             MaterialKind.CATALYST)
    labelled = graph.product("labelled_rna", "4sU pulse-labelled RNA",
                             MaterialKind.MOLECULAR, state=rna_pool("representative pulse-incorporated 4sU"))
    graph.react("pulse_label", "Pulse-label newly synthesized RNA with 4sU",
                StepKind.TRANSCRIPTION,
                ((cells, InputRole.TEMPLATE), (four_su, InputRole.STRUCTURAL),
                 (polymerase, InputRole.CATALYST)), (labelled,))
    reviewed = ["cells", "four_su", "cellular_rna_polymerase", "pulse_label"]

    extracted = _recover(graph, "extract_total_rna", labelled, "QIAzol/chloroform extraction and isopropanol precipitation, without SNU sonication")
    reviewed.extend(("extract_total_rna", "extract_total_rna_condition"))
    biotin = graph.input("thiol_biotin", "EZ-link HPDP-Biotin",
                         MaterialKind.MOLECULAR)
    condition = graph.input("biotinylation_condition", "HPDP-biotinylation in dimethylformamide",
                            MaterialKind.CONDITION)
    biotinylated = graph.product("biotinylated_rna", "biotinylated 4sU RNA",
                                MaterialKind.MOLECULAR, state=rna_pool("4sU conjugated to biotin through the reducible HPDP linkage"))
    graph.react("biotinylate", "Thiol-biotinylate incorporated 4sU",
                StepKind.CHEMICAL_CONVERSION,
                ((extracted, InputRole.TEMPLATE), (biotin, InputRole.STRUCTURAL),
                 (condition, InputRole.CONDITION)), (biotinylated,))
    cleaned_biotin = _recover(graph, "biotin_cleanup", biotinylated, "Phase-Lock/chloroform cleanup and isopropanol precipitation")
    reviewed.extend(("biotin_cleanup", "biotin_cleanup_condition"))
    streptavidin = graph.input("streptavidin", "µMACS streptavidin RNA-capture reagent",
                              MaterialKind.SELECTION)
    enriched = graph.product("enriched_rna", "streptavidin-bound label-containing RNA",
                             MaterialKind.MOLECULAR, state=molecule("bound RNA", (role("RNA retained through biotin-streptavidin binding"),)))
    unlabelled = graph.product("unlabelled_rna", "unbound capture/wash fraction",
                               MaterialKind.WASTE)
    graph.react("enrich", "Affinity-purify biotinylated 4sU RNA",
                StepKind.AFFINITY_CAPTURE,
                ((cleaned_biotin, InputRole.TEMPLATE),
                 (streptavidin, InputRole.SELECTION)), (enriched, unlabelled))
    graph.discard(unlabelled, "is removed during source binding and washes; no perfect biochemical separation is claimed")
    reviewed.extend(("thiol_biotin", "biotinylation_condition", "biotinylate",
                     "streptavidin", "enrich"))

    dtt = graph.input("dtt", "DTT reductive elution", MaterialKind.CONDITION)
    released = graph.product("released_rna", "released label-containing RNA", MaterialKind.MOLECULAR,
                             state=molecule("released RNA", (role("RNA containing incorporated 4sU sites"),)))
    beads_removed = graph.product("capture_beads_removed", "streptavidin support separated from eluate", MaterialKind.WASTE)
    graph.react("release", "Reduce the HPDP linkage and elute RNA with DTT", StepKind.CHEMICAL_CONVERSION,
                ((enriched, InputRole.TEMPLATE), (dtt, InputRole.CONDITION)), (released, beads_removed))
    graph.discard(beads_removed, "is separated from the collected source DTT eluates")
    recovered = _recover(graph, "post_capture", released, "miRNeasy column binding/washes before on-column DNase", label="column-bound RNA")
    dnase = graph.input("capture_dnase_enzyme", "on-column DNase I in source miRNeasy recovery", MaterialKind.CATALYST)
    cleaned = graph.product("cleaned_thio_rna", "DNA-digested recovered thio-RNA", MaterialKind.MOLECULAR, state=recovered.state)
    graph.react("capture_dnase", "Digest DNA on the source recovery column", StepKind.DIGESTION,
                ((recovered, InputRole.TEMPLATE), (dnase, InputRole.CATALYST)), (cleaned,))
    reviewed.extend(("dtt", "release", "post_capture", "post_capture_condition", "capture_dnase", "capture_dnase_enzyme"))
    cleaned = _recover(graph, "column_elution", cleaned, "Elute RNA after miRNeasy on-column DNA digestion")
    reviewed.extend(("column_elution", "column_elution_condition"))

    size_rule = graph.input("size_selection", "source 3.5% TBE-urea gel small-RNA selection",
                            MaterialKind.SELECTION)
    selected_rna = Construct((Segment("selected short RNA", "[size-selected 4sU RNA]",
                                      placeholder=True),), name="size-selected 4sU RNA")
    selected = graph.product("selected_rna", "gel-selected short enriched 4sU RNA",
                             MaterialKind.MOLECULAR, state=cleaned.state,
                             chemical_state=selected_rna)
    rejected = graph.product("rejected_sizes", "enriched RNA outside the selected size",
                             MaterialKind.WASTE)
    graph.react("select_size", "Gel-select short thio-labelled RNA",
                StepKind.SELECTION,
                ((cleaned, InputRole.TEMPLATE), (size_rule, InputRole.SELECTION)),
                (selected, rejected))
    graph.discard(rejected, "falls outside the size-resolved endpoint population")
    reviewed.extend(("size_selection", "select_size"))

    gel_recovered = _recover(graph, "gel_recovery", selected, "Dialyse RNA from excised gel and precipitate with isopropanol")
    reviewed.extend(("gel_recovery", "gel_recovery_condition"))
    rpph = graph.input("rpph", "RNA 5-prime pyrophosphohydrolase", MaterialKind.CATALYST)
    decapped = graph.product("decapped_rna", "RppH-treated short-RNA pool; compatible 5′-phosphate class represented",
                            MaterialKind.MOLECULAR, state=selected_rna, chemical_state=selected_rna)
    graph.react("decap", "Remove the RNA 5-prime cap with RppH", StepKind.DEPHOSPHORYLATION,
                ((gel_recovered, InputRole.TEMPLATE), (rpph, InputRole.CATALYST)), (decapped,))
    reviewed.extend(("rpph", "decap"))

    decap_stop = graph.input("decap_stop", "EDTA and heat inactivation after RppH", MaterialKind.CONDITION)
    stopped = graph.product("stopped_decap", "stopped decapping reaction", MaterialKind.MOLECULAR, state=selected_rna)
    graph.react("stop_decap", "Stop RppH with EDTA and heat", StepKind.QUENCHING,
                ((decapped, InputRole.TEMPLATE), (decap_stop, InputRole.CONDITION)), (stopped,))
    ready = _recover(graph, "decap_recovery", stopped, "Isopropanol recovery after decapping")
    reviewed.extend(("decap_stop", "stop_decap", "decap_recovery", "decap_recovery_condition"))
    adapter3 = graph.input("adapter_3p", "NEBNext E7300 preadenylated, 3-prime-blocked SR adapter",
                           MaterialKind.ADAPTER,
                           state=Construct((Segment("3-prime SR adapter", ADAPTER_3P, inferred=True),)),
                           source=KIT_SOURCE, inferred=True)
    ligase3 = graph.input("ligase_3p", "NEBNext 3-prime RNA-adapter ligase",
                          MaterialKind.CATALYST)
    three_prime = Construct((*selected_rna, Segment("3-prime SR adapter", ADAPTER_3P, inferred=True)),
                            name="3-prime-adapter-ligated short RNA")
    three_ligated = graph.product("three_prime_ligated", "3-prime-adapter-ligated RNA",
                                  MaterialKind.MOLECULAR, chemical_state=three_prime)
    geometry3 = _rna_ligation_geometry(
        ADAPTER_3P, selected_rna, adapter_end="5-prime", substrate_end="selected short RNA",
        statement=("NEB E7300 manual p. 16 prints the 5′-rApp and 3′-NH2 adapter; "
                   "its activated 5′ end joins compatible RNA 3′-OH; this does not establish recovery or compatibility of every source RNA."))
    free3 = graph.product("free_adapter3", "excess unligated 3′ adapter", MaterialKind.MOLECULAR,
                          state=Construct([Segment("free 3′ adapter", ADAPTER_3P, inferred=True)]), source=KIT_SOURCE, inferred=True)
    graph.react("ligate_3p", "INFERRED kit reference — ligate compatible 3′-OH RNA to the 3′ SR adapter",
                StepKind.LIGATION,
                ((ready, InputRole.TEMPLATE), (adapter3, InputRole.STRUCTURAL),
                 (ligase3, InputRole.CATALYST)), (three_ligated, free3), geometry=geometry3)
    reviewed.extend(("adapter_3p", "ligase_3p", "ligate_3p"))

    rt_primer = graph.input("library_rt_primer", "NEBNext E7300 SR RT primer",
                            MaterialKind.OLIGO,
                            state=Construct((Segment("SR RT primer", RT_PRIMER, inferred=True),)),
                            source=KIT_SOURCE, inferred=True)
    blocked = graph.product("blocked_adapter3", "RT-primer-blocked free 3′ adapter duplex", MaterialKind.MOLECULAR,
                            state=blocked_adapter_scene(), source=KIT_SOURCE, inferred=True)
    graph.react("block_adapter", "INFERRED kit reference — hybridize RT primer before 5′ ligation", StepKind.HYBRIDIZATION,
                ((free3, InputRole.TEMPLATE), (rt_primer, InputRole.STRUCTURAL)), (blocked,))
    graph.mark_terminal(Terminal(blocked.id, Disposition.UNRESOLVED,
                                 "The guide describes blocking excess free adapter; its study-specific downstream removal or carryover is unassigned."))
    reviewed.extend(("block_adapter", "blocked_adapter3"))

    adapter5 = graph.input("adapter_5p", "NEBNext E7300 5-prime RNA SR adapter",
                           MaterialKind.ADAPTER,
                           state=Construct((Segment("5-prime SR adapter", ADAPTER_5P, inferred=True),)),
                           source=KIT_SOURCE, inferred=True)
    ligase5 = graph.input("ligase_5p", "NEBNext 5-prime RNA-adapter ligase",
                          MaterialKind.CATALYST)
    doubly_adapted = Construct((Segment("5-prime SR adapter", ADAPTER_5P, inferred=True),
                               *selected_rna, Segment("3-prime SR adapter", ADAPTER_3P, inferred=True)),
                              name="doubly-adapter-ligated short RNA")
    ligated = graph.product("library_adapter_product", "doubly-adapter-ligated RNA",
                            MaterialKind.MOLECULAR, chemical_state=doubly_adapted)
    geometry5 = _rna_ligation_geometry(
        ADAPTER_5P, three_prime, adapter_end="3-prime", substrate_end="selected short RNA",
        statement=("NEB E7300 manual p. 16 prints the 5′ RNA adapter; "
                   "its 3′ end joins compatible monophosphorylated RNA 5′ ends; total recovery is not asserted."))
    graph.react("ligate_5p", "INFERRED kit reference — ligate the 5′ adapter to compatible 5′-phosphate RNA",
                StepKind.LIGATION,
                ((three_ligated, InputRole.TEMPLATE), (adapter5, InputRole.STRUCTURAL),
                 (ligase5, InputRole.CATALYST)), (ligated,), geometry=geometry5)
    reviewed.extend(("adapter_5p", "ligase_5p", "ligate_5p"))

    rt = graph.input("library_reverse_transcriptase", "reverse transcriptase",
                     MaterialKind.CATALYST)
    cdna = graph.product("library_cdna", "adapter-primed first-strand cDNA",
                         MaterialKind.MOLECULAR)
    graph.react("library_rt", "INFERRED kit reference — reverse-transcribe doubly-adapted RNA",
                StepKind.REVERSE_TRANSCRIPTION,
                ((ligated, InputRole.TEMPLATE), (rt_primer, InputRole.STRUCTURAL),
                 (rt, InputRole.CATALYST)), (cdna,))
    reviewed.extend(("library_rt_primer", "library_reverse_transcriptase", "library_rt"))

    forward = graph.input("library_pcr_forward", "NEBNext E7300 SR PCR primer",
                          MaterialKind.OLIGO,
                          state=Construct((Segment("SR PCR primer", PCR_FORWARD, inferred=True),)),
                          source=KIT_SOURCE, inferred=True)
    reverse = graph.input("library_pcr_reverse", "NEBNext E7300 indexed PCR primer family",
                          MaterialKind.OLIGO,
                          state=Construct((Segment("P7 and index landing", PCR_REVERSE_5P, inferred=True),
                                           Segment("unassigned six-base index of this manual family", "N" * 6,
                                                   placeholder=True, inferred=True, feature=feature("e7300_family_index", "sample_index", "unknown")),
                                           Segment("SR reverse annealing arm", PCR_REVERSE_3P, inferred=True))),
                          source=KIT_SOURCE, inferred=True)
    pcr_polymerase = graph.input("library_pcr_polymerase", "library PCR polymerase",
                                 MaterialKind.CATALYST)
    library = graph.product("library_library", "indexed sequencing library",
                            MaterialKind.LIBRARY, state=SPEC.FINAL_LIBRARY)
    dimer = graph.product("adapter_dimer", "source-observed adapter-dimer product class", MaterialKind.MOLECULAR, state=molecule("adapter dimer", (role("adapter-derived short DNA product"),)))
    graph.react("library_pcr", "INFERRED kit reference — index and amplify cDNA", StepKind.PCR,
                ((cdna, InputRole.TEMPLATE), (forward, InputRole.STRUCTURAL),
                 (reverse, InputRole.STRUCTURAL), (pcr_polymerase, InputRole.CATALYST)),
                (library, dimer), geometry=PCRGeometry.unavailable(
                    Construct([role("unassembled adapter-primed first-strand cDNA")], name="initial PCR template"), reason="The E7300 primer family is published; study-specific index/revision assignment and pre-PCR cDNA assembly remain open model questions."))
    graph.mark_terminal(Terminal(dimer.id, Disposition.UNRESOLVED, "Supplementary Figure4C shows an adapter-dimer band; its experiment-specific removal or carryover is not assigned by the inspected source."))
    graph.sequence(library, SequencingPlan(library=FINAL_LIBRARY, acquired_roles=SPEC.ACQUIRED_ROLES, read_lengths=SPEC.READ_LENGTHS, limitation=SPEC.SEQUENCING_UNAVAILABLE))
    reviewed.extend(("library_pcr_forward", "library_pcr_reverse", "library_pcr_polymerase",
                     "library_pcr", "library_library", "adapter_dimer"))
    graph.annotate(Highlight("library_pcr", "Source-observed adapter dimers",
                            "Supplementary Figure 4C shows a separate adapter-dimer band. The graph retains that product class with an unresolved study-specific fate, rather than inventing complete removal."))
    graph.annotate(Highlight("ligate_3p", "Conditional dated kit-guide implementation",
                            "The 2025 guide defines compatible 5′-phosphate/3′-OH substrates and anti-dimer handling. Its oligos are published; study lot/revision, exact index and quantitative recovery remain unassigned."))
    graph.annotate(SourceDiscrepancy("enrich", "Published wash-buffer difference",
                    (SourceClaim("2026 Supplementary Methods p.18, text338", "1 mM NaCl in the capture wash buffer"),
                     SourceClaim("Dölken et al.2008 original Methods p.1970/PDFp13, text966–967", "1 M NaCl in the earlier capture wash buffer")),
                    "Leave the numerical wash salt unassigned in this model.",
                    "Both readers inspected the actual original PDFs; each value is visibly printed, not a text-extraction substitution.",
                    "The 2026 method explicitly uses modifications. Whether the difference is intentional or a transcription error is unresolved; other ancestor conditions are not transferred."))
    graph.annotate(SourceDiscrepancy("select_size", "Size-selection windows differ across source records",
                    (SourceClaim("GSE179306 series Overall Design, line11", "50–100 bp"),
                     SourceClaim("GSM5452295 extraction protocol, line70; SRX11439177 construction protocol", "70–100 bp"),
                     SourceClaim("bioRxiv version1, 14 July2021, Methods1370–1385", "70–100 bp"),
                     SourceClaim("bioRxiv version2, 24 October2025, Methods1664–1680", "50–100 bp")),
                    "Keep the exact excision interval unassigned; the 2026 method describes the small-RNA region.",
                    "Original versioned methods and deposited series/sample records contain these distinct values.",
                    "Whether this reflects protocol revision, specific-library choices or a transcription difference remains unresolved; no typo or universal interval is assumed."))
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()


SPEC.protocol_graph = protocol_graph
