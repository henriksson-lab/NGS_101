"""Protocol-owned source, molecular workflow and reaction graph.

Authority: Sziraki et al., Nature Genetics 2023 (doi:10.1038/s41588-023-01572-y),
Supplementary Protocol 1 (publisher MOESM1, "EasySci-RNA protocol") and the
published oligo workbook MOESM3 Supplementary Tables 16-19.  The constant oligo
handles below are verbatim from those sources; they happen to equal canonical
Illumina/Nextera strings, so they are composed from ``illumina``/``nextera``.
"""
from adapter_ligation import AdapterLigationGeometry
from chemdraw import Construct, Scene, Segment, feature, oligo
from illumina import INDEX2_PRIMER_RC, P5, P7, P7_RC, TRUSEQ_READ1
from nextera import ME, ME_RC, S7, S7_RC
from pcr import PCRGeometry, PCRPrimer
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            SourceClaim, SourceDiscrepancy, StepKind)
from protocol_model import chain, molecule, role
from extended_protocols import Definition

DEFINITION = Definition('EasySci-RNA', '10.1038/s41588-023-01572-y', 'Sziraki et al., Nature Genetics (2023)', 'Barcoded oligo(dT) and random-hexamer RT, splinted barcode ligation and indexed PCR create a scalable three-round combinatorial single-nucleus RNA library.', 'RT-well, ligation-well and PCR-well (i7) barcodes jointly identify each nucleus; an 8-nt UMI identifies each cDNA molecule.', '01_easysci.html')

# Supplementary Protocol 1 materials list / MOESM3 Tables 16-19 (verbatim handles).
RT_HANDLE = TRUSEQ_READ1[15:]          # ACGACGCTCTTCCGATCT, 5'-phosphorylated RT primer head
LIGATION_ARM = TRUSEQ_READ1[:15]       # ACACTCTTTCCCTAC, 3' end of every ligation primer
LIGATION_SPLINT = INDEX2_PRIMER_RC     # "Ligation Adaptor Primer", all-PS bonds + /3ddC/
TN5_ME_B = S7 + ME                     # Tn5-ME-B, step 17
TN5_ME_REV = ME_RC                     # Tn5-ME-rev, 5'-phosphorylated, step 17
SHORT_DT_TAIL = "T" * 15               # Table 16, every one of 192 rows (see discrepancy)
RANDOM_TAIL = "N" * 6                  # Table 17


def _umi():
    return Segment("8-nt UMI", "N" * 8, "umi", placeholder=True,
                   feature=feature("easysci_umi", "umi", "random"))


def _rt_bc():
    return Segment("10-nt RT-well barcode", "N" * 10, "cbc", placeholder=True,
                   feature=feature("easysci_rt_well", "cell_barcode", "whitelist",
                                   group="easysci_cell", part="rt",
                                   whitelist="MOESM3 Supplementary Tables 16-17 (192+192, matched wells)"))


def _lig_bc():
    return Segment("10-nt ligation-well barcode", "N" * 10, "cbc", placeholder=True,
                   feature=feature("easysci_ligation_well", "cell_barcode", "whitelist",
                                   group="easysci_cell", part="ligation",
                                   whitelist="MOESM3 Supplementary Table 18 (384)"))


def _pcr_bc(name="10-nt PCR-well i7 barcode (as read)"):
    return Segment(name, "N" * 10, "cbc", placeholder=True,
                   feature=feature("easysci_pcr_well", "cell_barcode", "whitelist",
                                   group="easysci_cell", part="pcr",
                                   whitelist="MOESM3 Supplementary Table 19 (384)"))


def _priming_end():
    return role("short-dT T15 or randomN N6 priming end")


def _insert():
    return role("nuclear cDNA insert")


def _cdna_head():
    """First-strand cDNA 5' end: the 5'-phosphorylated RT primer and its extension."""
    return [Segment("RT primer handle (5′-Phos)", RT_HANDLE, "r1"), _umi(), _rt_bc(),
            _priming_end(), _insert()]


def _p5_side():
    return [Segment("P5", P5, "p5"), _lig_bc(), Segment("ligation-primer arm", LIGATION_ARM, "r1")]


def tagmented_template() -> Construct:
    """Gap-filled P5-ligated / Tn5-ME-B-ended molecule presented to PCR (derived)."""
    return Construct([*_p5_side(), *_cdna_head(), Segment("Tn5 ME (rev)", ME_RC, "me"),
                      Segment("Tn5-ME-B S7 handle (complement)", S7_RC, "s7")],
                     name="EasySci-RNA PCR template")


FINAL_LIBRARY = Construct([*_p5_side(), *_cdna_head(),
                           Segment("Tn5 ME (rev)", ME_RC, "me"),
                           Segment("Tn5-ME-B S7 handle (complement)", S7_RC, "s7"),
                           _pcr_bc(), Segment("P7 (complement)", P7_RC, "p7")],
                          name="EasySci-RNA final library")

P5_PRIMER = PCRPrimer("Universal P5 primer", P5)
P7_PRIMER = PCRPrimer("indexed P7 primer (i7 shown as N10)", P7 + "N" * 10 + S7,
                      anneal_nt=len(S7))

TITLE = DEFINITION.title
NOTES = DEFINITION.notes
SOURCE = (f'Defining source: <a href="https://doi.org/{DEFINITION.doi}">{DEFINITION.citation}</a>, '
          'Supplementary Protocol 1 and Supplementary Tables 16–19.')
SUMMARY = DEFINITION.summary
CAVEAT = ("The ligation junction, Tn5 gap fill and read content below are derived from the "
          "published oligo sequences and step order; the protocol does not draw them. "
          "Barcodes, UMI, priming end and insert are placeholders.")
SEQ_PRIMERS = ()
FINAL_CAPTION = DEFINITION.endpoint
READ_LENGTHS = {"Read 1": 100, "Index 1 (i7)": 10, "Index 2 (i5)": 10, "Read 2": 100}
ACQUIRED_ROLES = ("Read 1", "Index 1 (i7)", "Index 2 (i5)", "Read 2")
SEQUENCING_INTRO = (
    'Supplementary Protocol 1 step 71: NovaSeq, Read 1 100 bp, Read 2 100 bp, Index 1 10 bp, '
    'Index 2 10 bp. The <a href="https://github.com/JunyueCaoLab/EasySci">authors\' analysis '
    'README</a> labels the demultiplexed files R1 = genomic Read 1, R2 = P5 barcode, R3 = '
    'genomic Read 2 and I1 = P7 barcode, i.e. the i5 read carries the ligation-well barcode '
    'and the i7 read the PCR-well barcode. Derived, not stated: Read 1 would start at the '
    'TruSeq Read 1 site formed by the ligation-primer arm plus RT handle and cover UMI, '
    'RT barcode and priming end before cDNA; Read 2 would start in the Tn5-ME-B handle.')
SEQUENCING_ENDING = SEQUENCING_INTRO
SEQUENCING_UNAVAILABLE = ("The source gives run lengths but names no sequencing-primer reagent; "
                          "primer landing on the TruSeq-Read-1-like and Nextera-Read-2-like "
                          "handles is a derivation, so no binding oligo is asserted.")


def ligation_geometry() -> AdapterLigationGeometry:
    sc = Scene()
    sc.strand("ligation primer", _p5_side(), label="barcoded ligation primer (3′-OH)")
    sc.anneal("splint", [Segment("adapter: RT-handle arm", LIGATION_SPLINT[:18], "r1"),
                         Segment("adapter: ligation-arm pair", LIGATION_SPLINT[18:], "r1")],
              to="ligation primer", pair=("adapter: ligation-arm pair", "ligation-primer arm"),
              label="Ligation Adaptor Primer (all phosphorothioate)", mod3="/3ddC/")
    sc.anneal("cDNA", _cdna_head(), to="splint",
              pair=("RT primer handle (5′-Phos)", "adapter: RT-handle arm"),
              label="first-strand cDNA (5′-Phos RT primer)", above=True,
              unpaired=("8-nt UMI", "10-nt RT-well barcode",
                        "short-dT T15 or randomN N6 priming end", "nuclear cDNA insert"))
    sc.ligation_end("ligation primer", "ligation-primer arm", kind="splinted",
                    polarity="3-prime", owner="adapter",
                    text="ligation-primer 3′-OH presented by the annealed adapter")
    sc.ligation_end("cDNA", "RT primer handle (5′-Phos)", kind="splinted",
                    polarity="5-prime", owner="substrate",
                    text="5′-phosphorylated RT-primer end of the cDNA")
    return AdapterLigationGeometry.known(
        "splinted", sc,
        statement=("Supplementary Protocol 1 steps 13–16 anneal each barcoded ligation primer "
                   "with the 3′-ddC, phosphorothioate Ligation Adaptor Primer; the RT primers "
                   "carry /5Phos/. That the adaptor bridges the ligation-primer 3′ end to the "
                   "RT-primer 5′ end is derived from the sequences (Scene checks the pairing)."))


def workflow():
    return chain(molecule("fixed nuclei distributed into RT wells", (role("nuclear RNA"),)), (
        ("Reverse-transcribe with matched barcoded short-dT and randomN primers",
         molecule("first-strand cDNA in nuclei", tuple(_cdna_head())),
         "Both primer types are added to each RT well (step 31)."),
        ("Pool, redistribute and ligate the barcoded ligation primer",
         molecule("two-round indexed cDNA", (*_p5_side(), *_cdna_head()),
                  junctions=(("ligation-primer arm", "RT primer handle (5′-Phos)", "ligation"),)),
         "Splinted ligation appends P5 and the ligation-well barcode."),
        ("Second strand, tagment and index by PCR",
         molecule("EasySci-RNA library", tuple(FINAL_LIBRARY.segments)),
         "Three combinatorial rounds identify a nucleus."),
    ))


def oligos():
    yield oligo("Short-dT RT primer (Table 16; 192, A1 barcode TTCTCGCATG)",
                [Segment("RT handle", RT_HANDLE, "r1"), _umi(), _rt_bc(),
                 Segment("T15", SHORT_DT_TAIL)], mods="/5Phos/")
    yield oligo("RandomN RT primer (Table 17; 192, A1 barcode CGGTCAAGAA)",
                [Segment("RT handle", RT_HANDLE, "r1"), _umi(), _rt_bc(),
                 Segment("N6", RANDOM_TAIL, placeholder=True)], mods="/5Phos/")
    yield oligo("Ligation primer (Table 18; 384, A1 barcode CCGCGGCTCA)", _p5_side())
    yield oligo("Ligation Adaptor Primer (every bond phosphorothioate; 3′ dideoxy-C)",
                [Segment("adapter", LIGATION_SPLINT, "r1")], three="/3ddC/-3'")
    yield oligo("Tn5-ME-B", [Segment("S7", S7, "s7"), Segment("ME", ME, "me")])
    yield oligo("Tn5-ME-rev", [Segment("ME rev", ME_RC, "me")], mods="/5Phos/")
    yield oligo("Universal P5 primer", [Segment("P5", P5, "p5")])
    yield oligo("Indexed P7 primer (Table 19; 384; i7 printed lowercase)",
                [Segment("P7", P7, "p7"), _pcr_bc("10-nt i7"), Segment("S7", S7, "s7")])


def protocol_graph() -> ProtocolGraph:
    d = DEFINITION
    source = f"{d.citation}; doi:{d.doi}; Supplementary Protocol 1"
    graph = ProtocolGraph(d.title, authority=source)
    nuclei = graph.input(
        "rt_well_nuclei", "fixed nuclei (0.1% formaldehyde) in 2×96 RT wells",
        MaterialKind.BIOLOGICAL,
        state=molecule("fixed nuclei in RT wells", (role("nuclear RNA"),)), primary=True)
    dt = graph.input("barcoded_dt", "/5Phos/ indexed short-dT RT primer (Table 16)",
                     MaterialKind.OLIGO)
    random = graph.input("barcoded_random", "/5Phos/ indexed randomN RT primer (Table 17)",
                         MaterialKind.OLIGO)
    rt = graph.input("reverse_transcriptase", "Maxima H Minus reverse transcriptase",
                     MaterialKind.CATALYST)
    first = graph.product(
        "first_index_cdna", "RT-well-indexed first-strand cDNA in nuclei",
        MaterialKind.MOLECULAR,
        state=molecule("first-strand cDNA in nuclei", tuple(_cdna_head())))
    graph.react("reverse_transcription",
                "Anneal primers 55 °C 5 min, then RT 4→55 °C ramp, 55 °C 15 min",
                StepKind.REVERSE_TRANSCRIPTION,
                ((nuclei, InputRole.TEMPLATE), (dt, InputRole.STRUCTURAL),
                 (random, InputRole.STRUCTURAL), (rt, InputRole.CATALYST)), (first,))
    repartition = graph.input("repartition_condition", "pool, wash in NBB, redistribute into 4×96 ligation wells",
                              MaterialKind.CONDITION)
    redistributed = graph.product("redistributed_nuclei", "nuclei redistributed into ligation wells",
                                  MaterialKind.MOLECULAR, state=first.state)
    collisions = graph.product("repartition_losses", "empty or collision-prone redistributed wells",
                               MaterialKind.WASTE)
    graph.react("repartition", "Pool and redistribute nuclei into ligation wells",
                StepKind.PARTITION,
                ((first, InputRole.TEMPLATE), (repartition, InputRole.CONDITION)),
                (redistributed, collisions), changes_molecule=False)
    graph.discard(collisions, "does not provide an interpretable combinatorial cell identity")
    lig_primer = graph.input("ligation_barcode", "barcoded ligation primer (Table 18)",
                             MaterialKind.ADAPTER)
    splint = graph.input("ligation_splint", "Ligation Adaptor Primer (phosphorothioate, /3ddC/)",
                         MaterialKind.OLIGO,
                         state=Construct([Segment("Ligation Adaptor Primer", LIGATION_SPLINT)],
                                         name="Ligation Adaptor Primer"))
    ligase = graph.input("barcode_ligase", "T4 DNA ligase (RT, 30 min; EDTA stop)", MaterialKind.CATALYST)
    indexed = graph.product(
        "two_round_nuclei", "P5-ligated two-round-indexed cDNA in nuclei",
        MaterialKind.MOLECULAR,
        state=molecule("two-round indexed cDNA", (*_p5_side(), *_cdna_head()),
                       junctions=(("ligation-primer arm", "RT primer handle (5′-Phos)", "ligation"),)))
    graph.react("barcode_ligation", "Splinted ligation of the barcoded ligation primer",
                StepKind.LIGATION,
                ((redistributed, InputRole.TEMPLATE), (lig_primer, InputRole.STRUCTURAL),
                 (splint, InputRole.STRUCTURAL), (ligase, InputRole.CATALYST)), (indexed,),
                geometry=ligation_geometry())
    second_mix = graph.input("second_strand_polymerase",
                             "NEBNext Ultra II Non-Directional RNA Second Strand Synthesis Module",
                             MaterialKind.CATALYST)
    double = graph.product("double_stranded_cdna", "double-stranded two-index cDNA",
                           MaterialKind.MOLECULAR, state=indexed.state)
    graph.react("second_strand", "Second-strand synthesis, 16 °C 1 h (10,000 nuclei/well)",
                StepKind.EXTENSION,
                ((indexed, InputRole.TEMPLATE), (second_mix, InputRole.CATALYST)), (double,),
                reason="Module mechanism and its priming are vendor-internal; not drawn.")
    beads = graph.input("ampure_08x", "DNA binding buffer + 0.8× AMPure XP", MaterialKind.SELECTION)
    cleaned = graph.product("purified_cdna", "bead-purified double-stranded cDNA",
                            MaterialKind.MOLECULAR, state=double.state)
    supernatant = graph.product("cleanup_supernatant", "bead supernatant and ethanol washes",
                                MaterialKind.WASTE)
    graph.react("cdna_cleanup", "0.8× AMPure purification in each well", StepKind.CLEANUP,
                ((double, InputRole.TEMPLATE), (beads, InputRole.SELECTION)), (cleaned, supernatant))
    graph.discard(supernatant, "removed from the magnet-held beads")
    me_b = graph.input("tn5_me_b", "Tn5-ME-B oligo", MaterialKind.OLIGO,
                       state=Construct([Segment("Tn5-ME-B", TN5_ME_B)], name="Tn5-ME-B"))
    me_rev = graph.input("tn5_me_rev", "/5Phos/ Tn5-ME-rev oligo", MaterialKind.OLIGO,
                         state=Construct([Segment("Tn5-ME-rev", TN5_ME_REV)], name="Tn5-ME-rev"))
    tn5 = graph.input("tn5", "in-house Tn5 (Hennig et al. 2018)", MaterialKind.CATALYST)
    tagmented = graph.product("tagmented_cdna", "tagmented cDNA, P5-ligated end + Tn5-ME-B end",
                              MaterialKind.MOLECULAR,
                              state=molecule("tagmented cDNA", tuple(tagmented_template().segments)))
    losses = graph.product("tn5_losses",
                           "candidate internal fragments with Tn5-ME-B at both ends (no P5 end)",
                           MaterialKind.WASTE)
    graph.react("tagment", "Tagment with ME-B-only Tn5, 55 °C 5 min", StepKind.TAGMENTATION,
                ((cleaned, InputRole.TEMPLATE), (me_b, InputRole.STRUCTURAL),
                 (me_rev, InputRole.STRUCTURAL), (tn5, InputRole.CATALYST)), (tagmented, losses))
    graph.discard(losses, "candidate only: lacks the P5 end, so it cannot form a P5/P7 "
                          "cluster; no physical separation is described")
    p5 = graph.input("library_forward", "Universal P5 primer (added with SDS/BSA)", MaterialKind.OLIGO,
                     state=Construct([Segment("P5", P5)], name="Universal P5 primer"))
    p7 = graph.input("library_reverse", "indexed P7 primer (Table 19)", MaterialKind.OLIGO)
    pol = graph.input("library_polymerase", "NEBNext High-Fidelity 2X PCR Master Mix",
                      MaterialKind.CATALYST)
    library = graph.product("library_library", "EasySci-RNA PCR-well-indexed library",
                            MaterialKind.LIBRARY, state=FINAL_LIBRARY)
    graph.react("library_pcr", "SDS 55 °C 15 min, Tween-20; 72 °C 5 min, 12–15 cycles",
                StepKind.PCR,
                ((tagmented, InputRole.TEMPLATE), (p5, InputRole.STRUCTURAL),
                 (p7, InputRole.STRUCTURAL), (pol, InputRole.CATALYST)), (library,),
                geometry=PCRGeometry.exact(
                    tagmented_template(), P5_PRIMER, P7_PRIMER,
                    statement="Universal P5 and indexed P7 primers (Supplementary Protocol 1 "
                              "materials); the initial 72 °C 5 min is taken as Tn5 gap fill (derived)."))
    graph.sequence(library, SequencingPlan(limitation=SEQUENCING_UNAVAILABLE,
                                           acquired_roles=ACQUIRED_ROLES,
                                           read_lengths=READ_LENGTHS))
    graph.annotate(SourceDiscrepancy(
        "reverse_transcription", "Short-dT length differs between protocol and oligo table",
        (SourceClaim("Supplementary Protocol 1, Materials used (MOESM1 p.3, text lines 98–99)",
                     "…NNNNNNNN[10bp barcode] followed by 16 T"),
         SourceClaim("MOESM3 Supplementary Table 16, all 192 rows; Fig. 1d legend",
                     "15 T after the 10-nt barcode; legend: 'a stretch of 15 thymine nucleotides'")),
        "The model uses T15 from the per-well ordering table.",
        "Both values are printed in the published supplements; the table and figure legend agree.",
        "Which length was ordered for any given experiment is not stated beyond the table."))
    graph.annotate(SourceDiscrepancy(
        "reverse_transcription", "RandomN primer description",
        (SourceClaim("Fig. 1d legend", "RandomN, oligonucleotides composed of 8 random nucleotides"),
         SourceClaim("Supplementary Protocol 1 materials; MOESM3 Table 17",
                     "N8 UMI, 10-nt barcode, then a 3′ NNNNNN (6 random nucleotides)")),
        "The model uses the Table 17 sequence (N8 UMI and 3′ N6 random priming end).",
        "The main text also calls this primer a random hexamer.",
        "The legend's 8 may refer to the UMI; this is not stated."))
    reviewed = [
        "rt_well_nuclei", "barcoded_dt", "barcoded_random", "reverse_transcriptase",
        "reverse_transcription", "repartition_condition", "repartition",
        "ligation_barcode", "ligation_splint", "barcode_ligase", "barcode_ligation",
        "second_strand_polymerase", "second_strand", "ampure_08x", "cdna_cleanup",
        "tn5_me_b", "tn5_me_rev", "tn5", "tagment",
        "library_forward", "library_reverse", "library_polymerase", "library_pcr",
        "library_library",
    ]
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()
