"""scRCAT-seq: single-cell RNA cap-and-tail end capture after full-length cDNA.

Source scope (round 2026-10-17-breadth-first-17): published Nat. Commun. article
PMC7555861 (JATS Results line 65, Fig. 1a legend line 71, Methods lines 147-149
library construction, 169-171 TSS/TES read processing; rendered Fig. 1a on PDF p. 3).
The cited Smart-seq2 Nature Protocols original (ref 21) is not open and was not read,
so the TSO and preamplification primer stay unresolved.
"""
from chemdraw import Construct, Segment, feature, revcomp
import nextera as nx
import seqprimers as sp
from pcr import PCRGeometry, PCRPrimer
from protocol_graph import (Disposition, InputRole, MaterialKind, ProtocolGraph,
                            SequencingPlan, SourceClaim, SourceDiscrepancy, StepKind,
                            Terminal)
from protocol_model import molecule, role

TITLE = "scRCAT-seq — single-cell RNA cap and tail sequencing"
NOTES = "01_scrcat-seq.html"
DOI = "10.1038/s41467-020-18976-7"
AUTHORITY = f"Liu et al., Nature Communications (2020); doi:{DOI}"
METHODS = AUTHORITY + ", Methods ‘scRCAT-seq library construction for a single cell’ (PMC7555861)"
SOURCE = f'Defining source: <a href="https://doi.org/{DOI}">Liu et al., <i>Nature Communications</i> (2020)</a>.'
SUMMARY = ("Smart-seq2-style full-length cDNA is made with a barcoded, UMI-tagged oligo(dT) RT "
           "primer and template switching, then tagmented with Nextera XT. One PCR with a "
           "tag-specific P5 index primer and a Nextera P7 index primer captures fragments "
           "carrying either cDNA end in the same library; TSS and TES reads are separated "
           "computationally.")

# Methods line 149, 5'->3'.  Bracketed fields stay symbolic.
RT_HANDLE = "AAGCAGTGGTATCAACGCAGAGT"
P5_OUTER = "AATGATACGGCGACCACCGAGATCTACAC"
P5_INNER = "TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGGTGGTATCAACGCAGAGT"   # after [8 bps of index]
P7_OUTER = "CAAGCAGAAGACGGCATACGAGAT"
P7_INNER = "GTCTCGTGGGCTCGG"
CAPTURE_TAG = P5_INNER[-18:]          # GTGGTATCAACGCAGAGT, shared with the RT handle
# Methods line 169: TSS reads begin with this tag; the TSO part is GTGGTATCAACGCAGAGTACAT.
TSS_READ_TAG = "GTGGTATCAACGCAGAGTACATGGG"
TSO_READ_PART = "GTGGTATCAACGCAGAGTACAT"
assert RT_HANDLE.endswith(CAPTURE_TAG) and TSS_READ_TAG.startswith(TSO_READ_PART)
assert TSO_READ_PART.startswith(CAPTURE_TAG)
assert P5_INNER == nx.S5 + nx.ME + CAPTURE_TAG and P7_INNER == nx.S7


def oligos():
    return [
        "<b>RT primer</b> — 5′-" + RT_HANDLE + "N₈[16 bps of cell barcode]T₃₀VN-3′ (Methods)",
        "<b>P5 index primer</b> — 5′-" + P5_OUTER + "[8 bps of index]" + P5_INNER + "-3′ (Methods)",
        "<b>P7 index primer</b> — 5′-" + P7_OUTER + "[8 bps of index]" + P7_INNER + "-3′ (Methods)",
        "<b>TSS read tag</b> — 5′-" + TSS_READ_TAG + "-3′, of which " + TSO_READ_PART +
        " is trimmed as TSO sequence (read processing; the full TSO is not printed)",
    ]


UMI = feature("rt_umi", "umi", "random", note="N8 of the RT primer")
CELL = feature("rt_cell_barcode", "cell_barcode", "unknown",
               note="16-bp cell barcode field of the RT primer; sequences not printed")
I5 = feature("i5", "sample_index", "unknown",
             note="8-bp P5 index; Results says the cell barcodes are the same as sequencing indexes")
I7 = feature("i7", "sample_index", "unknown", note="8-bp P7 index")


def _p5_end():
    return [Segment("P5", P5_OUTER, "p5"),
            Segment("i5", "J" * 8, "i5", placeholder=True, feature=I5),
            Segment("S5", nx.S5, "r1"), Segment("ME (P5 primer)", nx.ME, "me"),
            Segment("capture tag", CAPTURE_TAG, "r2")]


def _p7_end():
    # The mosaic end on this side comes from the Nextera XT kit adaptor, whose original
    # was not inspected; the P7 primer itself is printed only from S7 outwards.
    return [Segment("ME′ (Nextera XT, inferred)", nx.ME_RC, "me", inferred=True),
            Segment("S7′", revcomp(P7_INNER), "r2"),
            Segment("i7", "I" * 8, "i7", placeholder=True, feature=I7),
            Segment("P7′", revcomp(P7_OUTER), "p7")]


TES_LIBRARY = Construct(_p5_end() + [
    Segment("UMI", "NNNNNNNN", "umi", feature=UMI),
    Segment("cell barcode", "B" * 16, "cbc", placeholder=True, feature=CELL),
    Segment("oligo(dT)", "T" * 30),
    Segment("VN", "VN", placeholder=True),
    Segment("transcript 3′ end (antisense)", "X" * 24, placeholder=True),
] + _p7_end(), name="scRCAT-seq TES-end molecule")

TSS_LIBRARY = Construct(_p5_end() + [
    Segment("TSO ACAT", TSS_READ_TAG[len(CAPTURE_TAG):len(TSO_READ_PART)], "tso"),
    Segment("GGG", TSS_READ_TAG[len(TSO_READ_PART):], "tso"),
    Segment("transcript 5′ end (sense)", "X" * 24, placeholder=True),
] + _p7_end(), name="scRCAT-seq TSS-end molecule")

READ_LENGTHS = {"Read 1": 150, "Read 2": 150}
TES_PRIMERS = sp.standard_primers(TES_LIBRARY)
TSS_PRIMERS = sp.standard_primers(TSS_LIBRARY)
FINAL_LIBRARIES = (
    ("TES-end molecules", TES_LIBRARY, TES_PRIMERS,
     "Derived from the printed RT and index primers. Read 1 crosses the capture tag, UMI, "
     "cell barcode and oligo(dT); Read 2 enters the transcript from the Tn5 end toward the poly(A).",
     "Standard Nextera primer sites are computed on the derived molecule; the paper states only "
     "Illumina X Ten PE150."),
    ("TSS-end molecules (same library)", TSS_LIBRARY, TSS_PRIMERS,
     "Read 1 begins with the printed TSS tag (capture tag, ACAT, GGG) and then the transcript 5′ end. "
     "The Methods tag carries no UMI or cell barcode (see discrepancy notice).",
     "Same physical library and primers as the TES-end molecules."),
)
CAVEAT = ("Both end molecules are assembled from the printed primers and TSS read tag; the mosaic "
          "end on the P7 side is inferred from the Nextera XT kit (not inspected). The Smart-seq2 "
          "TSO and preamplification primer (ref 21, not openly available) are not reconstructed.")
SEQUENCING_UNAVAILABLE = ""


def protocol_graph():
    g = ProtocolGraph(TITLE, authority=AUTHORITY)
    cell = g.input("cell", "single cell with poly(A) RNA and ERCC spike-in", MaterialKind.BIOLOGICAL,
                   primary=True, state=molecule("single-cell RNA", (role("poly(A) RNA"), role("ERCC spike-in"))))
    rt_primer = g.input("rt_primer", "RT primer: handle–N8–16-bp cell barcode–T30VN", MaterialKind.OLIGO,
                        source=METHODS)
    ssiii = g.input("superscript_iii", "SuperScript III (replacing Smart-seq2 SuperScript II)",
                    MaterialKind.CATALYST)
    first = g.product("first_strand", "first-strand cDNA with handle/UMI/cell barcode/T30 at its 5′ end",
                      MaterialKind.MOLECULAR,
                      state=molecule("first strand", (role("RT handle"), role("UMI"), role("cell barcode"),
                                                      role("T30VN"), role("antisense cDNA"))))
    g.react("reverse_transcribe", "Oligo(dT)-primed reverse transcription", StepKind.REVERSE_TRANSCRIPTION,
            ((cell, InputRole.TEMPLATE), (rt_primer, InputRole.STRUCTURAL), (ssiii, InputRole.CATALYST)),
            (first,), source=METHODS)
    tso = g.input("tso", "Smart-seq2 template-switch oligo (full sequence not printed; read part "
                  + TSO_READ_PART + ")", MaterialKind.OLIGO, inferred=True,
                  source=AUTHORITY + ", Methods read processing (partial); Smart-seq2 Nat. Protoc. ref 21 not read")
    switched = g.product("switched", "first strand extended across the TSO (…ACAT+GGG at the 5′ cap end)",
                         MaterialKind.MOLECULAR,
                         state=molecule("template-switched cDNA", (role("RT handle/UMI/cell barcode/T30VN"),
                                                                   role("antisense cDNA"), role("CCC"),
                                                                   role("TSO complement"))))
    g.react("template_switch", "Template switching at the transcript 5′ end", StepKind.TEMPLATE_SWITCHING,
            ((first, InputRole.TEMPLATE), (tso, InputRole.STRUCTURAL), (ssiii, InputRole.CATALYST)),
            (switched,), source=AUTHORITY + ", Results line 65; Methods line 149")
    smart_pcr = g.input("smart_seq2_pcr", "Smart-seq2 preamplification reagents (ref 21, not read)",
                        MaterialKind.CONDITION)
    amplified = g.product("amplified", "full-length cDNA carrying the capture tag at both ends",
                          MaterialKind.MOLECULAR,
                          state=molecule("amplified cDNA", (role("TSO end"), role("cDNA"),
                                                            role("T30/cell barcode/UMI/RT-handle end"))))
    g.react("preamplify", "Amplify full-length cDNA per Smart-seq2", StepKind.MAGIC,
            ((switched, InputRole.TEMPLATE), (smart_pcr, InputRole.CONDITION)), (amplified,),
            source=METHODS,
            reason="The article cites the Smart-seq2 protocol for PCR and does not print its primer or cycling; ref 21 was not available.")
    purification = g.input("purification", "cDNA purification (method not stated)", MaterialKind.SELECTION)
    purified = g.product("purified_cdna", "purified cDNA (0.1 ng taken for tagmentation)",
                         MaterialKind.MOLECULAR, state=amplified.state)
    pur_waste = g.product("purification_waste", "material removed by purification", MaterialKind.WASTE)
    g.react("purify_cdna", "Purify amplified cDNA", StepKind.CLEANUP,
            ((amplified, InputRole.TEMPLATE), (purification, InputRole.SELECTION)), (purified, pur_waste),
            source=METHODS)
    g.discard(pur_waste, "removed by the unspecified purification")

    nxt = g.input("nextera_xt_adaptors", "Nextera XT transposome adaptors (kit; original not inspected)",
                  MaterialKind.ADAPTER, inferred=True,
                  source="Illumina Nextera XT DNA sample preparation kit (named in Methods; not inspected)")
    tn5 = g.input("tn5", "Nextera XT Tn5 transposase", MaterialKind.CATALYST)
    tagmented = g.product("tagmented", "tagmented cDNA: end fragments with the capture tag plus internal fragments",
                          MaterialKind.MOLECULAR,
                          state=molecule("tagmented cDNA", (role("capture tag + cDNA end"), role("insert"),
                                                            role("Nextera adaptor end"))))
    g.react("tagment", "Tagment 0.1 ng cDNA with Nextera XT", StepKind.TAGMENTATION,
            ((purified, InputRole.TEMPLATE), (nxt, InputRole.STRUCTURAL), (tn5, InputRole.CATALYST)),
            (tagmented,), source=METHODS)

    p5 = g.input("p5_index_primer", "tag-specific P5 index primer ([8 bps of index])", MaterialKind.OLIGO,
                 source=METHODS)
    p7 = g.input("p7_index_primer", "Nextera P7 index primer ([8 bps of index])", MaterialKind.OLIGO,
                 source=METHODS)
    polymerase = g.input("capture_polymerase", "end-capture PCR polymerase (not named)", MaterialKind.CATALYST)
    captured = g.product("captured_library", "one indexed library holding TSS-end and TES-end molecules",
                         MaterialKind.LIBRARY, state=TES_LIBRARY, inferred=True,
                         source=METHODS + " (assembly derived; P7-side mosaic end inferred)")
    unselected = g.product("unselected_fragments",
                           "candidate fragments without a capture-tag end (internal Tn5 fragments)",
                           MaterialKind.MOLECULAR, inferred=True)
    g.react("end_capture_pcr", "One end-capture index PCR amplifies both cDNA ends", StepKind.PCR,
            ((tagmented, InputRole.TEMPLATE), (p5, InputRole.STRUCTURAL), (p7, InputRole.STRUCTURAL),
             (polymerase, InputRole.CATALYST)), (captured, unselected), source=METHODS,
            geometry=PCRGeometry.exact(
                TES_LIBRARY, PCRPrimer("P5 index primer after [index]", P5_INNER),
                PCRPrimer("P7 index primer after [index]", P7_INNER, anneal_nt=len(P7_INNER)),
                statement="Primer portions 3′ of the symbolic index windows located on the derived TES-end molecule; the TSS-end molecule shares both sites."))
    g.mark_terminal(Terminal("unselected_fragments", Disposition.UNRESOLVED,
                             "Derived candidate: no capture-tag site for the P5 index primer, so not exponentially "
                             "amplified; the 200–300 bp selection is the only separation described."))

    ampure = g.input("ampure", "1.8× AMPure XP beads", MaterialKind.SELECTION)
    cleaned = g.product("cleaned_library", "bead-purified library", MaterialKind.LIBRARY,
                        state=TES_LIBRARY, inferred=True)
    sup = g.product("ampure_supernatant", "bead supernatant and washes", MaterialKind.WASTE)
    g.react("ampure_cleanup", "Purify with 1.8× AMPure XP", StepKind.CLEANUP,
            ((captured, InputRole.TEMPLATE), (ampure, InputRole.SELECTION)), (cleaned, sup), source=METHODS)
    g.discard(sup, "removed with the bead supernatant")
    egel = g.input("egel", "E-Gel 2% SizeSelect, 200–300 bp", MaterialKind.SELECTION)
    selected = g.product("size_selected_library", "200–300 bp scRCAT-seq library", MaterialKind.LIBRARY,
                         state=TES_LIBRARY, inferred=True)
    rejected = g.product("rejected_sizes", "fragments outside 200–300 bp", MaterialKind.WASTE)
    g.react("size_select", "Select 200–300 bp fragments on E-Gel SizeSelect", StepKind.SELECTION,
            ((cleaned, InputRole.TEMPLATE), (egel, InputRole.SELECTION)), (selected, rejected), source=METHODS)
    g.discard(rejected, "left outside the collected size window")
    g.sequence(selected, SequencingPlan(
        library=TES_LIBRARY, primers=TES_PRIMERS,
        acquired_roles=tuple(p.role for p in TES_PRIMERS), read_lengths=READ_LENGTHS,
        limitation="Illumina X Ten PE150 (Methods); the primer set and index-read cycles are not stated, "
                   "so standard Nextera primer landing is computed on the derived molecule."))

    g.annotate(SourceDiscrepancy(
        "cell", "Identifier content of the TSO-end tag",
        (SourceClaim("Results line 65 and Fig. 1a legend line 71",
                     "a tag containing the UMI and cell identifier is added to both cDNA ends; Fig. 1a draws CS2+TAG on the TSO and T30+TAG+CS1 on the RT primer"),
         SourceClaim("Methods line 169 (TSS read processing)",
                     "TSS reads are selected by " + TSS_READ_TAG + " and " + TSO_READ_PART +
                     " is trimmed as TSO sequence — no random or barcode positions"),
         SourceClaim("Methods line 149",
                     "only the RT primer is printed with N8 and a 16-bp cell barcode; one P5 index primer ending "
                     + CAPTURE_TAG + " captures both ends")),
        "The TSS-end molecule carries only the printed TSS tag; UMI and cell barcode are drawn on the TES end only.",
        "All statements are in the published article and its rendered Fig. 1a.",
        "The full TSO sequence (Smart-seq2 ref 21) is not printed, and the article does not say how a TSO-end UMI or cell identifier would be read."))
    g.declare_reviewed_scope((
        "cell", "rt_primer", "superscript_iii", "reverse_transcribe", "tso", "template_switch",
        "smart_seq2_pcr", "preamplify", "purification", "purify_cdna", "nextera_xt_adaptors", "tn5",
        "tagment", "p5_index_primer", "p7_index_primer", "capture_polymerase", "end_capture_pcr",
        "ampure", "ampure_cleanup", "egel", "size_select", "size_selected_library"), AUTHORITY)
    return g.require_publishable()
