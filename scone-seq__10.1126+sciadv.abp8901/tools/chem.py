"""scONE-seq: one-tube single-cell genomic-DNA and total-RNA preamplification.

Source scope (round 2026-10-17-breadth-first-17): published Sci. Adv. article
PMC9812385 (JATS Results line 111; Materials and Methods lines 211-229) and the
published supplement (Supplementary Text "The molecular design of scONE-seq",
Fig. S1C legend, Table S2).  The preprint is not used.
"""
from types import SimpleNamespace

from chemdraw import Construct, Scene, Segment, feature, revcomp
import seqprimers as sp
from pcr import PCRGeometry
from protocol_graph import (Disposition, InputRole, MaterialKind, ProtocolGraph,
                            SequencingPlan, SourceClaim, SourceDiscrepancy, StepKind,
                            Terminal)
from protocol_model import chain, molecule, role

TITLE = "scONE-seq — one-tube single-cell DNA and RNA sequencing"
DOI = "10.1126/sciadv.abp8901"
CITATION = "Yu et al., Science Advances (2023)"
NOTES = "01_scone_seq.html"
AUTHORITY = f"{CITATION}; doi:{DOI}"
SOURCE = f'Defining source: <a href="https://doi.org/{DOI}">{CITATION}</a>.'
TABLE_S2 = AUTHORITY + ", Supplementary Table S2"
METHODS = AUTHORITY + ", Materials and Methods (PMC9812385)"
SUMMARY = ("In one tube, a custom Tn5 adaptor labels genomic DNA with a DNA-branch tag and "
           "UMI, and MATQ-seq-style RT primers label total RNA with an RNA-branch tag and "
           "UMI; after poly(C) tailing and 3′-adaptor second-strand synthesis both branches "
           "share one amplification primer. A second, P5-only tagmentation and index PCR "
           "keep the tagged end for sequencing.")

# Table S2, 5'->3', spaces removed (the table prints the handle, tag, N8 and priming
# sequence as separate groups).
SCONE_HANDLE = "GTCTCGTGGGCTCGG"           # shared amplification handle
SCONE_DNA_TAG = "TCATG"                    # printed DNA-branch tag group
SCONE_RNA_TAG = "AACTT"                    # printed RNA-branch tag group
SCONE_ME = "AGATGTGTATAAGAGACAG"           # Tn5 mosaic end in both Tn5 adaptors
SCONE_ONE_TN5 = SCONE_HANDLE + SCONE_DNA_TAG + "NNNNNNNN" + SCONE_ME
SCONE_P5_HANDLE = "TCGTCGGCAGCGTC"
SCONE_P5_TN5 = SCONE_P5_HANDLE + SCONE_ME
SCONE_MOSAIC = "CTGTCTCTTATACACATCT"
SCONE_RT_POLYT = SCONE_HANDLE + SCONE_RNA_TAG + "NNNNNNNN" + "T" * 20 + "VN"
SCONE_RT_GGG = SCONE_HANDLE + SCONE_RNA_TAG + "NNNNNNNN" + "GGGHN"
SCONE_RT_TTT = SCONE_HANDLE + SCONE_RNA_TAG + "NNNNNNNN" + "TTTVN"
SCONE_3P_ADAPTER = SCONE_RT_GGG             # printed identically, without biotin
SCONE_AMP_PRIMER = "GATGTGTGGAG" + SCONE_HANDLE
SCONE_P7_OUTER = "CAAGCAGAAGACGGCATACGAGAT"
SCONE_P7_INNER = "AGATGTGTATAAGAGACAG" + SCONE_HANDLE   # P7 index primer after [i7]
SCONE_P5_OUTER = "AATGATACGGCGACCACCGAGATCTACAC"
SCONE_ILLU_R1 = "TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG"
SCONE_ILLU_R2 = "CCGAGCCCACGAGACCTGTCTCTTATACACATCT"   # row labelled "ONE-Illu Read 2"
SCONE_ILLU_I7 = "AGATGTGTATAAGAGACAGGTCTCGTGGGCTCGG"   # row labelled "ONE-Illu I7 Index Primer"

# The two printed custom rows are exact reverse complements; the page discrepancy
# below depends on this relation, so it is enforced at import rather than restated.
assert revcomp(SCONE_ILLU_R2) == SCONE_ILLU_I7 == SCONE_P7_INNER


def oligos():
    """Exact oligos printed in Supplementary Table S2 (index windows stay symbolic)."""
    return (
        "<b>One-Tn5 adaptor</b> — 5′-[biotin]" + SCONE_ONE_TN5 + "-3′",
        "<b>P5-Tn5 adaptor</b> — 5′-" + SCONE_P5_TN5 + "-3′",
        "<b>Mosaic</b> — 5′-[phos]" + SCONE_MOSAIC + "-3′",
        "<b>RT-polyT</b> — 5′-[biotin]" + SCONE_RT_POLYT + "-3′ (mixed 1.5:5:5 with RT-GGG and RT-TTT)",
        "<b>RT-GGG</b> — 5′-[biotin]" + SCONE_RT_GGG + "-3′",
        "<b>RT-TTT</b> — 5′-[biotin]" + SCONE_RT_TTT + "-3′",
        "<b>3′ adaptor</b> — 5′-" + SCONE_3P_ADAPTER + "-3′",
        "<b>Amp-primer</b> — 5′-" + SCONE_AMP_PRIMER + "-3′",
        "<b>P7 index-primer family</b> — 5′-" + SCONE_P7_OUTER + "[i7]" + SCONE_P7_INNER + "-3′",
        "<b>P5 index-primer family</b> — 5′-" + SCONE_P5_OUTER + "[i5]" + SCONE_P5_HANDLE + "-3′",
        "<b>ONE-Illu Read 1</b> — 5′-" + SCONE_ILLU_R1 + "-3′",
        "<b>ONE-Illu Read 2</b> (as labelled) — 5′-" + SCONE_ILLU_R2 + "-3′",
        "<b>ONE-Illu I7 Index Primer</b> (as labelled) — 5′-" + SCONE_ILLU_I7 + "-3′",
    )


_TAG_NOTE = ("Table S2 prints this as a five-base group; the main text calls it a "
             "6-nucleotide barcode (see discrepancy notice).")
DNA_TAG = feature("dna_branch", "inline_barcode", "fixed",
                  note="DNA-branch tag on the One-Tn5 adaptor. " + _TAG_NOTE)
RNA_TAG = feature("rna_branch", "inline_barcode", "fixed",
                  note="RNA-branch tag on RT primers and 3′ adaptor. " + _TAG_NOTE)
DNA_UMI = feature("dna_umi", "umi", "random", note="N8 on the One-Tn5 adaptor")
RNA_UMI = feature("rna_umi", "umi", "random", note="N8 on the RT primers")
I5 = feature("i5", "sample_index", "unknown", note="Table S2 [i5] window; length not printed")
I7 = feature("i7", "sample_index", "unknown", note="Table S2 [i7] window; length not printed")


def one_tn5_adaptor():
    """One-Tn5 adaptor annealed to the phosphorylated Mosaic strand (Methods line 211)."""
    sc = Scene()
    sc.strand("One-Tn5", [Segment("handle", SCONE_HANDLE, "r2"),
                          Segment("DNA tag", SCONE_DNA_TAG, "cbc", feature=DNA_TAG),
                          Segment("DNA UMI", "NNNNNNNN", "umi", feature=DNA_UMI),
                          Segment("ME", SCONE_ME, "me")],
              label="One-Tn5 adaptor (5′ biotin)")
    sc.anneal("Mosaic", [Segment("Mosaic", SCONE_MOSAIC, "me")], to="One-Tn5",
              pair=("Mosaic", "ME"), label="Mosaic", mod5="phos")
    return sc


def p5_tn5_adaptor():
    sc = Scene()
    sc.strand("P5-Tn5", [Segment("P5 handle", SCONE_P5_HANDLE, "r1"),
                         Segment("ME", SCONE_ME, "me")], label="P5-Tn5 adaptor")
    sc.anneal("Mosaic", [Segment("Mosaic", SCONE_MOSAIC, "me")], to="P5-Tn5",
              pair=("Mosaic", "ME"), label="Mosaic", mod5="phos")
    return sc


def amp_primer():
    sc = Scene()
    sc.strand("Amp-primer", [Segment("5′ tail", SCONE_AMP_PRIMER[:11]),
                             Segment("handle", SCONE_HANDLE, "r2")], label="Amp-primer")
    return sc


def rt_primer_mix():
    sc = Scene()
    # V and H are printed IUPAC degenerate positions; they are kept as stand-ins.
    for name, prime, degenerate in (("RT-polyT", "T" * 20, "VN"), ("RT-GGG", "GGG", "HN"),
                                    ("RT-TTT", "TTT", "VN")):
        sc.strand(name, [Segment("handle", SCONE_HANDLE, "r2"),
                         Segment("RNA tag", SCONE_RNA_TAG, "cbc", feature=RNA_TAG),
                         Segment("RNA UMI", "NNNNNNNN", "umi", feature=RNA_UMI),
                         Segment("priming bases", prime),
                         Segment("degenerate 3′ end", degenerate, placeholder=True)],
                  label=name + " (5′ biotin)")
    return sc


def three_prime_adaptor():
    """Table S2 3' adaptor: the RT-GGG bases without the biotin."""
    sc = Scene()
    sc.strand("3′ adaptor", [Segment("handle", SCONE_HANDLE, "r2"),
                             Segment("RNA tag", SCONE_RNA_TAG, "cbc", feature=RNA_TAG),
                             Segment("RNA UMI", "NNNNNNNN", "umi", feature=RNA_UMI),
                             Segment("poly(C)-pairing bases", "GGG"),
                             Segment("degenerate 3′ end", "HN", placeholder=True)],
              label="3′ adaptor")
    return sc


# Derived final library, DNA branch, P5 -> P7 top strand.  Every base comes from a
# Table S2 oligo; the assembly (P5-Tn5 end + insert + tagged preamplification end +
# P7 index primer) follows Fig. S1C legend and the supplement text, and is derived.
FINAL_LIBRARY = Construct([
    Segment("P5", SCONE_P5_OUTER, "p5"),
    role("i5", tag="i5", molecular_feature=I5),
    Segment("P5-Tn5 handle", SCONE_P5_HANDLE, "r1"),
    Segment("ME (P5-Tn5)", SCONE_ME, "me"),
    role("genomic insert (second-tagmentation end)"),
    Segment("ME′ (One-Tn5)", revcomp(SCONE_ME), "me"),
    Segment("DNA UMI′", "NNNNNNNN", "umi", feature=DNA_UMI),
    Segment("DNA tag′", revcomp(SCONE_DNA_TAG), "cbc", feature=DNA_TAG),
    Segment("handle′", revcomp(SCONE_HANDLE), "r2"),
    Segment("P7-primer ME′", revcomp(SCONE_P7_INNER[:19]), "me"),
    role("i7", tag="i7", molecular_feature=I7),
    Segment("P7′", revcomp(SCONE_P7_OUTER), "p7"),
], name="scONE-seq DNA-branch library (derived)")

SEQ_PRIMERS = (
    sp.custom("Read 1", "ONE-Illu Read 1", SCONE_ILLU_R1, TABLE_S2),
    sp.custom("Index 1 (i7)", "row labelled ‘ONE-Illu Read 2’", SCONE_ILLU_R2, TABLE_S2,
              "Computed landing: extends from the handle′/ME′ junction into i7; see discrepancy notice."),
    sp.custom("Read 2", "row labelled ‘ONE-Illu I7 Index Primer’", SCONE_ILLU_I7, TABLE_S2,
              "Computed landing: extends into the branch tag, UMI and insert; see discrepancy notice."),
)
SEQUENCING_ENDING = ("NextSeq 500 with the custom primers of Table S2 (Methods line 229); Table S3 "
                     "costs a High Output 150-cycle kit, but the read configuration is not stated. "
                     "The i5 index primer is left to the instrument manufacturer (Table S2 asterisk). "
                     "An MGI conversion route with ONE-MGI primers is also listed and not drawn.")


def protocol_graph():
    g = ProtocolGraph(TITLE, authority=AUTHORITY)
    cell = g.input("cell", "FACS-sorted single cell or nucleus in 2.3 µl lysis buffer",
                   MaterialKind.BIOLOGICAL, primary=True,
                   state=molecule("single cell or nucleus", (role("genomic DNA"), role("total RNA"))))
    pk = g.input("proteinase_k", "proteinase K (55 °C 15 min; 75 °C 15 min inactivation)",
                 MaterialKind.CATALYST)
    lysate = g.product("lysate", "lysate with genomic DNA and total RNA in one tube",
                       MaterialKind.MOLECULAR, state=cell.state)
    g.react("lyse", "Lyse with proteinase K, then heat-inactivate it",
            StepKind.DEPROTEINIZATION, ((cell, InputRole.TEMPLATE), (pk, InputRole.CATALYST)),
            (lysate,), source=METHODS + " line 215")

    transposome = g.input("one_tn5", "One-Tn5 adaptor/Mosaic loaded on Tn5 (Novoprotein), 24 °C 45 min",
                          MaterialKind.ADAPTER, state=one_tn5_adaptor(),
                          source=TABLE_S2 + "; Methods line 211")
    tn5 = g.input("tn5", "Tn5 transposase (Novoprotein)", MaterialKind.CATALYST)
    kapa = g.input("kapa_tagmentation", "KAPA polymerase in the tagmentation mix (PEG-8000, TAPS)",
                   MaterialKind.CATALYST)
    tagged = g.product("tagged_lysate", "tagmented genomic fragments carrying handle/DNA tag/N8/ME at both ends, plus RNA",
                       MaterialKind.MOLECULAR,
                       state=molecule("tagged DNA in lysate", (
                           role("handle"), role("DNA tag"), role("DNA UMI"), role("ME"),
                           role("genomic fragment"), role("ME′/UMI′/tag′/handle′"),
                           role("total RNA (untagged)"))))
    g.react("tagment", "Tagment genomic DNA in the lysate; 55 °C 10 min then 72 °C 10 min",
            StepKind.TAGMENTATION,
            ((lysate, InputRole.TEMPLATE), (transposome, InputRole.STRUCTURAL),
             (tn5, InputRole.CATALYST), (kapa, InputRole.CATALYST)), (tagged,), source=METHODS + " line 215",
            reason="The 72 °C KAPA step is read as gap filling across the transposition junction; the source lists the reagents and temperatures but not the mechanism.")
    pk2 = g.input("proteinase_k_2", "proteinase K (55 °C 15 min; 75 °C 15 min)", MaterialKind.CATALYST)
    released = g.product("released", "deproteinized tagged DNA and RNA", MaterialKind.MOLECULAR,
                         state=tagged.state)
    g.react("release_tn5", "Proteinase K removes Tn5 and inactivates enzymes",
            StepKind.DEPROTEINIZATION, ((tagged, InputRole.TEMPLATE), (pk2, InputRole.CATALYST)),
            (released,), changes_molecule=False, source=METHODS + " line 215")

    rt_primers = g.input("rt_primers", "RT-polyT : RT-GGG : RT-TTT (1.5:5:5) internal RT primer mix",
                         MaterialKind.OLIGO, state=rt_primer_mix(), source=TABLE_S2)
    ssiii = g.input("superscript_iii", "SuperScript III; 12→42 °C gradient, 50 °C 50 min, 55 °C 50 min",
                    MaterialKind.CATALYST)
    rt_mix = g.product("rt_products", "first-strand cDNA (handle/RNA tag/N8 at 5′) on RNA, plus tagged DNA",
                       MaterialKind.MOLECULAR,
                       state=molecule("first-strand cDNA", (
                           role("handle"), role("RNA tag"), role("RNA UMI"),
                           role("internal or poly(A) priming end"), role("cDNA"))))
    g.react("reverse_transcribe", "Internally prime total RNA (polyA and non-polyA)",
            StepKind.REVERSE_TRANSCRIPTION,
            ((released, InputRole.TEMPLATE), (rt_primers, InputRole.STRUCTURAL),
             (ssiii, InputRole.CATALYST)), (rt_mix,), source=METHODS + " line 215; Results line 111")

    nucleases = g.input("exo_rnase", "thermolabile Exo I, RNase If and RNase H", MaterialKind.CATALYST)
    single_cdna = g.product("cleared_cdna", "single-stranded cDNA and tagged DNA after primer/RNA removal",
                            MaterialKind.MOLECULAR, state=rt_mix.state)
    digested1 = g.product("digested_primers_rna", "residual RT primers and RNA digested in the tube",
                          MaterialKind.WASTE)
    g.react("clear_primers_rna", "Digest residual RT primers (Exo I) and RNA (RNase If, RNase H)",
            StepKind.DIGESTION, ((rt_mix, InputRole.TEMPLATE), (nucleases, InputRole.CATALYST)),
            (single_cdna, digested1), source=METHODS + " line 215")
    g.discard(digested1, "enzymatically degraded in the same tube; not physically separated")

    tdt = g.input("tdt_dctp", "terminal transferase + dCTP, 37 °C 5 min; proteinase K stop",
                  MaterialKind.CATALYST)
    tailed = g.product("c_tailed", "poly(C)-tailed cDNA", MaterialKind.MOLECULAR,
                       state=molecule("C-tailed cDNA", (
                           role("handle"), role("RNA tag"), role("RNA UMI"), role("cDNA"),
                           role("poly(C) tail"))))
    g.react("c_tail", "Poly(C)-tail cDNA 3′ ends", StepKind.TAILING,
            ((single_cdna, InputRole.TEMPLATE), (tdt, InputRole.CATALYST)), (tailed,),
            source=METHODS + " line 215")

    adaptor3 = g.input("adaptor_3p", "3′ adaptor (handle/RNA tag/N8/GGGHN)", MaterialKind.OLIGO,
                       state=three_prime_adaptor(), source=TABLE_S2)
    kapa2 = g.input("kapa_second_strand", "KAPA polymerase; 72 °C 5 min, 10 × (48 °C 1 min, 72 °C 1 min)",
                    MaterialKind.CATALYST)
    ds_cdna = g.product("ds_cdna", "double-stranded cDNA with the handle/RNA tag/N8 at both ends",
                        MaterialKind.MOLECULAR,
                        state=molecule("ds cDNA", (
                            role("handle/RNA tag/UMI/priming end"), role("cDNA"),
                            role("poly(G)/UMI′/tag′/handle′"))))
    g.react("second_strand", "Prime the poly(C) tail with the 3′ adaptor (non-denaturing cycles)",
            StepKind.EXTENSION,
            ((tailed, InputRole.TEMPLATE), (adaptor3, InputRole.STRUCTURAL),
             (kapa2, InputRole.CATALYST)), (ds_cdna,), source=METHODS + " line 215")

    exo2 = g.input("exo_i_2", "Exo I (37 °C 30 min; 55 °C 15 min; 65 °C 15 min)", MaterialKind.CATALYST)
    cleared = g.product("co_templates", "tagged genomic DNA and ds cDNA sharing the handle",
                        MaterialKind.MOLECULAR,
                        state=molecule("co-amplifiable templates", (
                            role("handle"), role("branch tag"), role("UMI"),
                            role("genomic or cDNA insert"), role("UMI′/tag′/handle′"))))
    digested2 = g.product("digested_adaptor", "residual 3′ adaptor digested in the tube",
                          MaterialKind.WASTE)
    g.react("clear_adaptor", "Digest residual 3′ adaptor with Exo I", StepKind.DIGESTION,
            ((ds_cdna, InputRole.TEMPLATE), (exo2, InputRole.CATALYST)), (cleared, digested2),
            source=METHODS + " lines 215, 221")
    g.discard(digested2, "enzymatically degraded in the same tube; not physically separated")

    amp = g.input("amp_primer", "Amp-primer (single primer for both branches)", MaterialKind.OLIGO,
                  state=amp_primer(), source=TABLE_S2)
    hotstart = g.input("kapa_hotstart", "KAPA HotStart ReadyMix, DMSO; 18–22 cycles", MaterialKind.CATALYST)
    preamp = g.product("preamplified", "co-amplified DNA and cDNA with Amp-primer ends",
                       MaterialKind.MOLECULAR,
                       state=molecule("preamplified product", (
                           role("Amp-primer 5′ tail"), role("handle"), role("branch tag"),
                           role("UMI"), role("insert"), role("UMI′/tag′/handle′/tail′"))))
    g.react("preamplify", "Amplify DNA and cDNA together with the Amp-primer", StepKind.PCR,
            ((cleared, InputRole.TEMPLATE), (amp, InputRole.STRUCTURAL),
             (hotstart, InputRole.CATALYST)), (preamp,), source=METHODS + " line 215",
            geometry=PCRGeometry.unavailable(
                Construct([role("handle"), role("branch tag"), role("UMI"), role("insert"),
                           role("UMI′/tag′/handle′")], name="co-amplified templates"),
                symmetric_species=True,
                reason="One Amp-primer species primes the shared handle at both ends; templates are role-level mixtures of genomic and cDNA inserts."))

    ampure = g.input("ampure_1", "AMPure XP beads", MaterialKind.SELECTION)
    purified = g.product("purified_preamp", "purified preamplified product (diluted to 0.1 ng/µl)",
                         MaterialKind.MOLECULAR, state=preamp.state)
    sup1 = g.product("ampure_1_supernatant", "bead supernatant and washes", MaterialKind.WASTE)
    g.react("purify_preamp", "Purify the preamplified product on AMPure XP", StepKind.CLEANUP,
            ((preamp, InputRole.TEMPLATE), (ampure, InputRole.SELECTION)), (purified, sup1),
            source=METHODS + " line 229")
    g.discard(sup1, "removed with the bead supernatant")

    p5tn5 = g.input("p5_tn5", "P5-Tn5 adaptor/Mosaic transposome, 55 °C 15 min",
                    MaterialKind.ADAPTER, state=p5_tn5_adaptor(), source=TABLE_S2 + "; Methods line 229")
    tn5_2 = g.input("tn5_library", "Tn5 transposase loaded with P5-Tn5 (TAPS, PEG-8000)", MaterialKind.CATALYST)
    end_frag = g.product("tagged_end_fragments", "fragments with the preamplified tagged end and one P5-Tn5 end",
                         MaterialKind.MOLECULAR,
                         state=molecule("tagged-end fragment", (
                             role("P5-Tn5 end"), role("insert"), role("UMI′/tag′/handle′/tail′"))))
    internal = g.product("internal_p5_fragments",
                         "candidate internal fragments with P5-Tn5 ends only", MaterialKind.MOLECULAR,
                         inferred=True)
    g.react("retagment", "Second tagmentation adds the Read 1 (P5-Tn5) end", StepKind.TAGMENTATION,
            ((purified, InputRole.TEMPLATE), (p5tn5, InputRole.STRUCTURAL),
             (tn5_2, InputRole.CATALYST)), (end_frag, internal),
            source=METHODS + " line 229; Supplementary Text; Fig. S1C legend")
    g.mark_terminal(Terminal("internal_p5_fragments", Disposition.UNRESOLVED,
                             "Derived candidate: lacks the handle site for the P7 index primer; the "
                             "source states only that the 5′ tagged end is kept and describes no separation."))

    p5_idx = g.input("p5_index_primers", "P5 index-primer family ([i5])", MaterialKind.OLIGO, source=TABLE_S2)
    p7_idx = g.input("p7_index_primers", "P7 index-primer family ([i7], 3′ handle)", MaterialKind.OLIGO,
                     source=TABLE_S2)
    hifi = g.input("kapa_hifi", "KAPA HiFi HotStart; 10–11 cycles", MaterialKind.CATALYST)
    indexed = g.product("indexed_library", "indexed scONE-seq library", MaterialKind.LIBRARY,
                        state=FINAL_LIBRARY, inferred=True,
                        source=TABLE_S2 + " (assembly derived)")
    g.react("index_pcr", "Index PCR with P5 and modified P7 index primers", StepKind.PCR,
            ((end_frag, InputRole.TEMPLATE), (p5_idx, InputRole.STRUCTURAL),
             (p7_idx, InputRole.STRUCTURAL), (hifi, InputRole.CATALYST)), (indexed,),
            source=METHODS + " line 229",
            geometry=PCRGeometry.unavailable(
                FINAL_LIBRARY, reason="Index windows are symbolic and the template is role-level."))

    ampure2 = g.input("ampure_2", "pooling and AMPure XP beads", MaterialKind.SELECTION)
    pooled = g.product("pooled_library", "pooled purified library", MaterialKind.LIBRARY,
                       state=FINAL_LIBRARY, inferred=True)
    sup2 = g.product("ampure_2_supernatant", "bead supernatant and washes", MaterialKind.WASTE)
    g.react("pool_purify", "Pool samples and purify on AMPure XP", StepKind.CLEANUP,
            ((indexed, InputRole.TEMPLATE), (ampure2, InputRole.SELECTION)), (pooled, sup2),
            source=METHODS + " line 229")
    g.discard(sup2, "removed with the bead supernatant")

    scdash = g.input("scdash", "scDASH reagents (refs 33, 58; not inspected)", MaterialKind.CONDITION)
    depleted = g.product("depleted_library", "library depleted of abundant ribosomal and mitochondrial sequences",
                         MaterialKind.LIBRARY, state=FINAL_LIBRARY, inferred=True)
    g.react("scdash_deplete", "scDASH depletion of rRNA and mitochondrial sequences; optional double size selection",
            StepKind.MAGIC, ((pooled, InputRole.TEMPLATE), (scdash, InputRole.CONDITION)), (depleted,),
            source=METHODS + " line 229",
            reason="The scDASH predecessor method is cited, not described; its reagents and the fate of removed molecules are not modelled.")
    g.sequence(depleted, SequencingPlan(
        library=FINAL_LIBRARY, primers=SEQ_PRIMERS,
        acquired_roles=tuple(p.role for p in SEQ_PRIMERS), limitation=SEQUENCING_ENDING))

    g.annotate(SourceDiscrepancy(
        "depleted_library", "Custom Read 2 and i7 primer labels appear swapped",
        (SourceClaim("Table S2, P7 index-primer row",
                     "5′-…[i7]" + SCONE_P7_INNER + "-3′: the 3′ handle primes the tagged end, so this sequence lies between i7 and the branch tag"),
         SourceClaim("Supplementary Text, ‘The molecular design of scONE-seq’",
                     "the Read 2 sequence includes the DNA and RNA barcodes"),
         SourceClaim("Table S2, ONE-Illu rows",
                     "‘ONE-Illu Read 2’ = " + SCONE_ILLU_R2 + "; ‘ONE-Illu I7 Index Primer’ = " + SCONE_ILLU_I7 + " (exact reverse complements)")),
        "The sequencing panel places both printed sequences by computation: the row labelled I7 extends into the branch tag, UMI and insert (Read 2); the row labelled Read 2 extends into i7 (Index 1). Labels are kept as printed.",
        "All three statements are in the published supplement; the derived library is assembled from Table S2 oligos.",
        "The source does not say how the custom primers were loaded on the NextSeq 500; whether the labels are a typographical swap is not established."))
    g.annotate(SourceDiscrepancy(
        "cell", "Branch-barcode length differs between main text and Table S2",
        (SourceClaim("Results, first section (PMC9812385 JATS line 111)",
                     "a 6-nucleotide ‘DNA barcode’ and a 6-nucleotide ‘RNA barcode’"),
         SourceClaim("Table S2 and Methods line 215",
                     "the tag group is printed as " + SCONE_DNA_TAG + " (One-Tn5) and " + SCONE_RNA_TAG + " (RT primers, 3′ adaptor), five bases each")),
        "The model uses the printed bases (five) and does not add a sixth.",
        "Both values are in the published article and its supplement.",
        "Whether the authors count an adjacent handle base or a base is missing from the printed oligos is not stated."))
    g.declare_reviewed_scope((
        "cell", "proteinase_k", "lyse", "one_tn5", "tn5", "kapa_tagmentation", "tagment",
        "proteinase_k_2", "release_tn5", "rt_primers", "superscript_iii", "reverse_transcribe",
        "exo_rnase", "clear_primers_rna", "tdt_dctp", "c_tail", "adaptor_3p",
        "kapa_second_strand", "second_strand", "exo_i_2", "clear_adaptor", "amp_primer",
        "kapa_hotstart", "preamplify", "ampure_1", "purify_preamp", "p5_tn5", "tn5_library", "retagment",
        "p5_index_primers", "p7_index_primers", "kapa_hifi", "index_pcr", "ampure_2",
        "pool_purify", "scdash", "scdash_deplete", "depleted_library"), AUTHORITY)
    return g.require_publishable()


def workflow():
    """Role-level summary chain (the page renders the graph)."""
    return chain(molecule("lysed single cell or nucleus", (role("genomic DNA"), role("total RNA"))), (
        ("Tagment DNA with the One-Tn5 adaptor", molecule("tagged DNA", (
            role("handle"), role("DNA tag"), role("DNA UMI"), role("genomic fragment"))), ""),
        ("Internally prime RNA, C-tail, 3′-adaptor second strand", molecule("ds cDNA", (
            role("handle"), role("RNA tag"), role("RNA UMI"), role("cDNA"))), ""),
        ("Co-amplify, retagment with P5-Tn5, index PCR",
         molecule("scONE-seq library", tuple(FINAL_LIBRARY)), ""),
    ))


FINAL_CAPTION = ("Derived DNA-branch library (P5 → P7). The RNA branch is identical except that "
                 "the tag is " + SCONE_RNA_TAG + " and the UMI is followed by the RT-primer or "
                 "3′-adaptor priming end instead of the mosaic end.")
SEQUENCING_INTRO = SEQUENCING_ENDING
READ_LENGTHS = None
CAVEAT = ("The final library is assembled from Table S2 oligos following the Fig. S1C legend; the "
          "assembly is derived. Index windows, the i5 primer and acquired read lengths are not printed.")
SEQUENCING_UNAVAILABLE = ""
