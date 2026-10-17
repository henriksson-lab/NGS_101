"""Protocol-owned reaction graph for ChAIR (Chai et al., Nature Methods 2025).

Sources inspected (2026-10-17 bounded pass):
- published article, PMC12955729 JATS, Methods "ChAIR" (xml.txt line 179) and
  "ChAIR data processing" (line 187);
- the article's Supplementary Information (Springer MOESM1, "A step-by-step ChAIR
  protocol", dated 9.30.2024), steps 1-82 and the Addendum "Bridge linker preparation";
- the 10x kit it names (Chromium Next GEM Single Cell Multiome ATAC + Gene Expression,
  1000283; user guide CG000338, no revision cited by ChAIR): the Rev G appendix
  "Sequences" (p. 80) supplies the ATAC-library shell drawn below.
The bioRxiv preprint (10.1101/2024.05.03.592322) is superseded and not used.
"""
import illumina as il
import nextera as nx
from chemdraw import Construct, Segment, complement, feature, revcomp
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            SourceClaim, SourceDiscrepancy, StepKind)
from pcr import PCRGeometry
from protocol_model import molecule, role

TITLE = "ChAIR — chromatin accessibility, interaction and RNA"
DOI = "10.1038/s41592-025-02658-7"
CITATION = "Chai et al., Nature Methods (2025)"
SUMMARY = ("In-situ AluI (or AluI + HpyCH4V) digestion, Klenow dA-tailing and ligation of a "
           "biotinylated T-overhang bridge linker record chromatin contacts; Tn5 then tags "
           "accessible chromatin and 10x Multiome droplets barcode DNA and RNA together. "
           "Streptavidin capture separates the biotinylated DNA component from the cDNA.")
NOTES = "01_chair.html"
SOURCE = (f'Defining source: <a href="https://doi.org/{DOI}">{CITATION}</a>, Methods and '
          'Supplementary Information (step-by-step protocol).')
CAVEAT = ("10x kit components (ATAC Buffer B/Enzyme B, Multiome Gel Beads, Pre-Amp Primers, "
          "SI-PCR Primer B, Sample Index N, cDNA Primers and the GEX library reagents) are "
          "kit-defined. ChAIR names kit 1000283 and guide CG000338 but no revision; the drawn "
          "ATAC shell is from the CG000338 Rev G appendix and the sample-index and cell-barcode "
          "values stay unassigned. Read lengths are not stated (\"sequenced in accordance with "
          "the manufacturer's guidelines\").")

# Supplementary Information, Addendum "Bridge linker preparation" (MOESM1 lines 375-376)
# and the published Methods (JATS line 179): 5'-/5Phos/CGTGATATT/iBIOdT/CACGACTCT-3' and
# 5'-/5Phos/GAGTCGTGAAATATCACGT-3'. The biotin is carried on the T at position 10.
BRIDGE_TOP = "CGTGATATTTCACGACTCT"
BRIDGE_BOTTOM = "GAGTCGTGAAATATCACGT"
BRIDGE_BIOTIN_POS = 10   # 1-based position of /iBIOdT/ in BRIDGE_TOP
MULTIOME_SPACER = "CGCGTCTG"   # CG000338 Rev G appendix p. 80, ATAC library


def bridge_linker_ends() -> tuple[int, str, str]:
    """Anneal the two printed strands; return (duplex bp, top 3' overhang, bottom 3' overhang).

    Raises unless they form one duplex with single-base 3'-T overhangs, the ends that
    the Klenow (3'->5' exo-) dA-tailing step makes ligatable (Methods line 179).
    """
    top, bottom = BRIDGE_TOP, revcomp(BRIDGE_BOTTOM)   # bottom in top orientation
    for shift in range(-3, 4):
        a0, b0 = max(0, shift), max(0, -shift)
        n = min(len(top) - a0, len(bottom) - b0)
        if n >= 15 and top[a0:a0 + n] == bottom[b0:b0 + n]:
            top_3 = top[a0 + n:]
            bottom_3 = complement(bottom[:b0])[::-1]
            if a0 or b0 + n != len(bottom) or top_3 != "T" or bottom_3 != "T":
                break
            return n, top_3, bottom_3
    raise ValueError("ChAIR bridge-linker strands no longer form a 3'-T/3'-T duplex")


def contact_junction() -> str:
    """Top strand across fragment A (dA-tailed) + bridge linker + fragment B.

    The top-strand reading is dA + the top oligo; it must equal the reverse complement
    of (dA + bottom oligo), i.e. both dA tails pair with the two 3'-T overhangs. One
    linker orientation is shown; the opposite orientation reads the reverse complement.
    """
    bridge_linker_ends()
    junction = "A" + BRIDGE_TOP
    if revcomp(junction) != "A" + BRIDGE_BOTTOM:
        raise ValueError("bridge linker does not seal to two dA-tailed ends")
    return junction


CELL = role("10x ATAC cell barcode", tag="cbc", length_bp=16,
            molecular_feature=feature("cell_barcode", "cell_barcode", "whitelist",
                                      whitelist="10x Multiome ATAC gel-bead barcodes"))
SAMPLE = role("Sample Index N", tag="cbc", length_bp=8,
              molecular_feature=feature("sample_index", "sample_index", "unknown"))


def dna_library() -> Construct:
    """ChAIR-DNA contact (PET) molecule inside the CG000338 Rev G ATAC-library shell."""
    return Construct([
        Segment("P5", il.P5, "p5"), CELL,
        Segment("10x spacer", MULTIOME_SPACER, None),
        Segment("Read 1N (s5 + mosaic end)", nx.ADAPTOR_S5, "s5"),
        role("Tn5-cut end of restriction fragment A"),
        Segment("dA tail of fragment A", "A", None),
        Segment(f"bridge linker top strand (biotin-dT at {BRIDGE_BIOTIN_POS})",
                BRIDGE_TOP, "w1"),
        role("restriction fragment B to its Tn5-cut end"),
        Segment("mosaic end + s7 (Read 2N side)", nx.ME_RC + nx.S7_RC, "s7"),
        SAMPLE,
        Segment("P7", il.P7_RC, "p7"),
    ], name="ChAIR-DNA library (contact molecule)")


RNA_LIBRARY = Construct([role("10x Multiome GEX P5 end"), role("10x GEX cell barcode"),
                         role("UMI", tag="umi"), role("poly(dT)VN and cDNA insert"),
                         role("10x GEX P7 end with sample index")],
                        name="ChAIR-RNA library roles")

FINAL_LIBRARY = dna_library()
FINAL_CAPTION = ("ChAIR-DNA contact molecule in the Multiome ATAC shell (CG000338 Rev G "
                 "appendix): the 20-bp dA + bridge-linker junction is computed from the two "
                 "printed linker strands; one of two linker orientations is shown. "
                 "Accessibility (ChAIR-ATAC) reads come from the same library without a "
                 "linker; ChAIR-PIPE separates them by linker trimming.")
SEQ_PRIMERS = ()
SEQUENCING_INTRO = FINAL_CAPTION
READ_LENGTHS = None
SEQUENCING_ENDING = ("ChAIR-DNA and ChAIR-RNA libraries were \"prepared and sequenced in "
                     "accordance with the manufacturer's guidelines\" and demultiplexed with "
                     "bcl2fastq (Methods). ChAIR-ATAC and ChAIR-PET are computational "
                     "partitions of the ChAIR-DNA reads (linker trimming in ChAIR-PIPE).")
SEQUENCING_UNAVAILABLE = ("ChAIR states no read lengths and no kit revision; the "
                          "CG000338 cycle recommendations are not assigned to ChAIR here.")


def oligos():
    return (f"bridge linker top — 5′-/5Phos/{BRIDGE_TOP[:BRIDGE_BIOTIN_POS - 1]}/iBIOdT/"
            f"{BRIDGE_TOP[BRIDGE_BIOTIN_POS:]}-3′",
            f"bridge linker bottom — 5′-/5Phos/{BRIDGE_BOTTOM}-3′")


def contact_state():
    return molecule("bridge-linked contact junction", (
        role("restriction fragment A"), role("biotinylated bridge linker"),
        role("restriction fragment B")),
        junctions=(("restriction fragment A", "biotinylated bridge linker", "T/A ligation"),
                   ("biotinylated bridge linker", "restriction fragment B", "T/A ligation")))


def gem_state():
    return molecule("GEM products", (CELL, role("transposed ChAIR DNA"), role("cDNA")))


def _pcr(g, ident, label, template, mix, construct, hifi, out):
    """PCR with a kit primer mix; its two primer species are declared, not sequenced."""
    a = g.input(ident + "_primer_1", mix + ", first primer species (kit-defined)",
                MaterialKind.OLIGO)
    b = g.input(ident + "_primer_2", mix + ", second primer species (kit-defined)",
                MaterialKind.OLIGO)
    g.react(ident, label, StepKind.PCR,
            ((template, InputRole.TEMPLATE), (a, InputRole.STRUCTURAL),
             (b, InputRole.STRUCTURAL), (hifi, InputRole.CATALYST)), (out,),
            geometry=PCRGeometry.unavailable(construct, reason=(
                mix + " is a 10x kit reagent; ChAIR does not print its sequences and the "
                "kit revision is not cited.")))


def protocol_graph():
    source = f"{CITATION}; doi:{DOI}"
    kit = "10x CG000338 (revision not cited by ChAIR); kit-defined reagent"
    g = ProtocolGraph(TITLE, authority=source)
    M, R, C, S, T = (MaterialKind.MOLECULAR, InputRole.TEMPLATE, InputRole.CATALYST,
                     InputRole.STRUCTURAL, InputRole.CONDITION)

    nucleus = g.input("nucleus", "FA-EGS-crosslinked cells or nuclei (50,000-300,000)",
                      MaterialKind.BIOLOGICAL, primary=True)
    sds = g.input("sds", "0.1% SDS FA lysis buffer, then 0.1% SDS", MaterialKind.CONDITION)
    perm = g.product("permeabilized", "lysed, SDS-permeabilized nuclei", M)
    g.react("permeabilize", "Lyse (omit for nuclei) and permeabilize with 0.1% SDS",
            StepKind.PERMEABILIZATION, ((nucleus, R), (sds, T)), (perm,))
    triton = g.input("triton", "Triton X-100 (2.5 µl 20%)", MaterialKind.CONDITION)
    quenched = g.product("quenched", "SDS-quenched nuclei", M)
    g.react("quench", "Quench SDS with Triton X-100 at 37 °C", StepKind.QUENCHING,
            ((perm, R), (triton, T)), (quenched,), changes_molecule=False)
    enzyme = g.input("restriction_enzyme", "AluI, or AluI + HpyCH4V (NEB)", MaterialKind.CATALYST)
    digested = g.product("digested", "in-situ restricted chromatin", M)
    g.react("digest", "Digest in situ (37 °C, ≥2 h or overnight, CutSmart)",
            StepKind.DIGESTION, ((quenched, R), (enzyme, C)), (digested,))
    klenow = g.input("klenow", "Klenow fragment (3′→5′ exo−) with dATP", MaterialKind.CATALYST)
    tailed = g.product("a_tailed", "dA-tailed restriction ends", M)
    g.react("a_tail", "dA-tail restriction ends (37 °C 1 h)", StepKind.TAILING,
            ((digested, R), (klenow, C)), (tailed,))
    stopped = g.product("a_tail_stopped", "heat-stopped dA-tailed chromatin", M)
    heat = g.input("heat_65", "65 °C, 20 min", MaterialKind.CONDITION)
    g.react("heat_stop", "Stop dA-tailing at 65 °C for 20 min", StepKind.THERMAL_TREATMENT,
            ((tailed, R), (heat, T)), (stopped,), changes_molecule=False)

    top = g.input("bridge_top", "bridge linker top oligo (5′Phos, internal biotin-dT)",
                  MaterialKind.OLIGO, state=(Segment("bridge linker top", BRIDGE_TOP, "w1"),))
    bottom = g.input("bridge_bottom", "bridge linker bottom oligo (5′Phos)",
                     MaterialKind.OLIGO,
                     state=(Segment("bridge linker bottom", BRIDGE_BOTTOM, "w1"),))
    linker = g.product("bridge_linker",
                       "annealed bridge linker, 18 bp with 3′-T overhangs (computed)",
                       MaterialKind.ADAPTER)
    g.react("anneal_linker", "Anneal bridge-linker strands (95 °C, stepped ramp to 4 °C)",
            StepKind.HYBRIDIZATION, ((top, R), (bottom, S)), (linker,))
    ligase = g.input("t4_ligase", "T4 DNA ligase, NEB Quick ligation buffer", MaterialKind.CATALYST)
    contacts = g.product("contacts", "bridge-linked chromatin-contact junctions retained with RNA",
                         M, state=contact_state())
    g.react("proximity_ligation", "Ligate the bridge linker in situ (RT ≥4 h or 16 °C overnight)",
            StepKind.PROXIMITY_LIGATION, ((stopped, R), (linker, S), (ligase, C)), (contacts,))

    buffers = g.input("nuclei_buffer", "ATAC-RSB wash; 10x Diluted Nuclei Buffer",
                      MaterialKind.CONDITION)
    strainer = g.input("strainer", "brief sonication and 20 µm strainer", MaterialKind.SELECTION)
    nuclei = g.product("filtered_nuclei", "washed, declumped ligated nuclei", M)
    clumps = g.product("aggregates", "nuclear aggregates retained by the strainer",
                       MaterialKind.WASTE)
    g.react("wash_filter", "Wash, resuspend, sonicate briefly and filter (20 µm)",
            StepKind.SELECTION, ((contacts, R), (buffers, T), (strainer, InputRole.SELECTION)),
            (nuclei, clumps), changes_molecule=False)
    g.discard(clumps, "physically retained on the 20 µm strainer")
    tn5_ends = g.input("atac_transposome_ends", "Tn5 adapter ends loaded in ATAC Enzyme B ("
                       + kit + ")", MaterialKind.ADAPTER)
    tn5 = g.input("atac_enzyme_b", "ATAC Enzyme B in ATAC Buffer B (" + kit + ")",
                  MaterialKind.CATALYST)
    tagged = g.product("tagged", "transposed nuclei carrying contact and accessible fragments", M)
    g.react("tagment", "Tagment 12,000 nuclei, 37 °C 1 h", StepKind.TAGMENTATION,
            ((nuclei, R), (tn5_ends, S), (tn5, C)), (tagged,))
    beads = g.input("gel_beads", "Multiome Gel Beads / GEM master mix, Chip J (" + kit + ")",
                    MaterialKind.ADAPTER)
    gems = g.product("gems", "cell-barcoded DNA and cDNA in GEMs", M,
                     state=gem_state())
    g.react("gem_barcoding", "Partition and barcode in GEMs (37 °C 45 min, 25 °C 30 min)",
            StepKind.MAGIC, ((tagged, R), (beads, S)), (gems,),
            reason=("Chip J partitioning, ATAC-fragment barcoding and reverse transcription "
                    "follow CG000338 'strictly' (SI steps 36-41); kit reactions not re-modelled."))
    recovery = g.input("recovery_agent", "Recovery Agent", MaterialKind.SELECTION)
    broken = g.product("recovered", "aqueous GEM products", M)
    oil = g.product("partitioning_oil", "Recovery Agent/partitioning oil phase", MaterialKind.WASTE)
    g.react("break_emulsion", "Break emulsion; remove oil phase", StepKind.CLEANUP,
            ((gems, R), (recovery, InputRole.SELECTION)), (broken, oil), changes_molecule=False)
    g.discard(oil, "removed from the bottom of the tube")
    pk = g.input("proteinase_k", "Proteinase K", MaterialKind.CATALYST)
    decrosslinked = g.product("decrosslinked", "decrosslinked DNA and cDNA", M)
    g.react("decrosslink", "Proteinase K, 65 °C 2 h (reverse crosslinks)",
            StepKind.DEPROTEINIZATION, ((broken, R), (pk, C)), (decrosslinked,))
    dyna = g.input("dynabeads_cleanup", "Dynabeads Cleanup Mix (" + kit + ")", MaterialKind.SELECTION)
    cleaned = g.product("cleaned", "purified nucleic acids in Elution Solution I", M)
    sup1 = g.product("cleanup_supernatant", "Dynabeads cleanup supernatant and ethanol washes",
                     MaterialKind.WASTE)
    g.react("dynabeads_cleanup_step", "Dynabeads cleanup and elution", StepKind.CLEANUP,
            ((decrosslinked, R), (dyna, InputRole.SELECTION)), (cleaned, sup1),
            changes_molecule=False)
    g.discard(sup1, "removed on the magnet")
    m280 = g.input("m280", "Dynabeads M-280 streptavidin, iBlock- and genomic-DNA-blocked",
                   MaterialKind.SELECTION)
    dna = g.product("dna_component", "biotinylated DNA on M-280 beads (DNA component)", M)
    cdna = g.product("cdna_component", "M-280 supernatant (cDNA component)", M)
    g.react("streptavidin_capture", "Capture biotinylated dsDNA on M-280 beads (RT 45 min)",
            StepKind.AFFINITY_CAPTURE, ((cleaned, R), (m280, InputRole.SELECTION)),
            (dna, cdna), changes_molecule=False)

    hifi = g.input("nebnext_hifi", "NEBNext High-Fidelity 2× PCR Master Mix", MaterialKind.CATALYST)
    dna_amp = g.product("dna_preamp", "pre-amplified ChAIR DNA", M)
    _pcr(g, "dna_preamp_pcr", "Pre-amplify on beads (72 °C 5 min; 7 cycles)", dna,
         "Pre-Amp Primers", FINAL_LIBRARY, hifi, dna_amp)
    spri = g.input("spriselect", "SPRIselect beads", MaterialKind.SELECTION)
    dna_clean = g.product("dna_preamp_clean", "1.6× SPRI-purified DNA pre-amplicons", M)
    m280_waste = g.product("used_m280", "M-280 beads and SPRI-excluded material", MaterialKind.WASTE)
    g.react("dna_preamp_cleanup", "Discard M-280 beads; 1.6× SPRI", StepKind.CLEANUP,
            ((dna_amp, R), (spri, InputRole.SELECTION)), (dna_clean, m280_waste),
            changes_molecule=False)
    g.discard(m280_waste, "beads removed on the magnet; unbound material in SPRI supernatant")
    dna_idx = g.product("dna_indexed", "sample-indexed ChAIR DNA", M)
    _pcr(g, "dna_index_pcr", "Sample-index PCR (7 cycles)", dna_clean,
         "SI-PCR Primer B + Sample Index N", FINAL_LIBRARY, hifi, dna_idx)
    dna_lib = g.product("dna_library", "ChAIR-DNA library", MaterialKind.LIBRARY,
                        state=FINAL_LIBRARY)
    dna_rej = g.product("dna_size_rejects", "fragments outside the 0.55×/0.8× window",
                        MaterialKind.WASTE)
    g.react("dna_size_select", "1× SPRI, then double-sided 0.55×/0.8× SPRI", StepKind.CLEANUP,
            ((dna_idx, R), (spri, InputRole.SELECTION)), (dna_lib, dna_rej),
            changes_molecule=False)
    g.discard(dna_rej, "excluded by double-sided SPRI selection")
    g.sequence(dna_lib, SequencingPlan(limitation=SEQUENCING_ENDING))

    cdna_clean = g.product("cdna_clean", "1.8× SPRI-purified cDNA component", M)
    cdna_rej = g.product("cdna_sup_rejects", "1.8× SPRI supernatant", MaterialKind.WASTE)
    g.react("cdna_cleanup", "1.8× SPRI of the M-280 supernatant", StepKind.CLEANUP,
            ((cdna, R), (spri, InputRole.SELECTION)), (cdna_clean, cdna_rej),
            changes_molecule=False)
    g.discard(cdna_rej, "SPRI supernatant removed; non-biotinylated DNA fate not stated")
    cdna_amp = g.product("cdna_preamp", "pre-amplified cDNA", M)
    _pcr(g, "cdna_preamp_pcr", "Pre-amplify cDNA (72 °C 5 min; 7 cycles), 1.6× SPRI",
         cdna_clean, "Pre-Amp Primers", RNA_LIBRARY, hifi, cdna_amp)
    cdna_full = g.product("cdna_amplified", "amplified cDNA (0.6× SPRI)", M)
    _pcr(g, "cdna_amplification", "cDNA amplification (7 cycles), 0.6× SPRI", cdna_amp,
         "cDNA Primers", RNA_LIBRARY, hifi, cdna_full)
    gex = g.input("gex_reagents", "GEX fragmentation, end-repair/A-tailing, adaptor and "
                  "sample-index reagents (" + kit + ")", MaterialKind.CATALYST)
    rna_lib = g.product("rna_library", "ChAIR-RNA library", MaterialKind.LIBRARY,
                        state=RNA_LIBRARY)
    g.react("gex_library", "Build the GEX library per CG000338", StepKind.MAGIC,
            ((cdna_full, R), (gex, C)), (rna_lib,),
            reason=("ChAIR delegates fragmentation, end repair/A-tailing, adaptor ligation "
                    "and sample-index PCR (with 0.6×/0.8× selections) to CG000338; those "
                    "kit reactions are not re-modelled in this pass."))
    g.sequence(rna_lib, SequencingPlan(limitation=SEQUENCING_ENDING))

    g.annotate(SourceDiscrepancy(
        "permeabilize", "Permeabilization time differs between Methods and protocol",
        (SourceClaim("Methods 'ChAIR' (PMC12955729 JATS line 179)",
                     "permeabilized by 10 µl 0.1% SDS at RT for 2 h"),
         SourceClaim("Supplementary Information steps 4–5",
                     "incubate at RT for 1.5 h (to be optimized per cell type)")),
        "The graph names the SDS permeabilization without a duration.",
        "Both values are in the published article and its Supplementary Information.",
        "The protocol says the time must be optimized per cell type; which value applied to each dataset is not stated."))
    g.annotate(SourceDiscrepancy(
        "quench", "SDS-quench time differs",
        (SourceClaim("Methods 'ChAIR' (JATS line 179)", "Triton X-100, 37 °C for 20 min"),
         SourceClaim("Supplementary Information step 7", "37 °C for 15 min with shaking")),
        "The graph records the quench without a duration.",
        "Both are published originals for the same assay.",
        "Not resolved by the inspected sources."))
    g.annotate(SourceDiscrepancy(
        "wash_filter", "Declumping sonication differs",
        (SourceClaim("Methods 'ChAIR' (JATS line 179)", "sonicated for 1 second"),
         SourceClaim("Supplementary Information step 27", "sonicate with 43% for 1.5 s")),
        "The graph records a brief sonication without a duration or power.",
        "Both are published originals for the same assay.",
        "Not resolved by the inspected sources; immaterial to the drawn molecules."))
    g.declare_reviewed_scope(tuple(g.materials) + tuple(g.steps), source)
    return g.require_publishable()
