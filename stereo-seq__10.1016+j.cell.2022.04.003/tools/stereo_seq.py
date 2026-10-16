"""Stereo-seq DNA-nanoball patterned-array spatial transcriptomics."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"lib"))
from chemdraw import Construct, Segment, feature
from adapter_ligation import AdapterLigationGeometry
from pcr import PCRGeometry
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind)

TITLE="Stereo-seq"
NOTES="01_stereo-seq.html"
SOURCE='Defining source: <a href="https://doi.org/10.1016/j.cell.2022.04.003">Chen et al., <i>Cell</i> (2022)</a>.'
SUMMARY="Spatially indexed DNA nanoballs are arrayed and decoded before receiving a UMI and oligo-dT capture sequence. Tissue mRNA is reverse-transcribed on this surface; released cDNA retains the coordinate identifier and is converted to a DNBSEQ library."
CAVEAT=("The published Cell STAR Methods print the chip, capture, template-switch, "
        "cDNA-PCR, library-PCR and read-primer oligos listed below. Its Key Resources "
        "Table DNB oligo1 (local text lines 1365–1369) has TTGAACT, while bioRxiv v3 "
        "Methods precursor section and the linked MOSTA detailed protocol primer table line 229 "
        "have TTGNAACT (an extra N). This model uses TTGNAACT "
        "because the author preprint and detailed protocol agree on it; which precursor "
        "was synthesized for the published experiment remains uncertain. Cell Methods "
        "(local text line 1684) call the capture tail 22-nt poly(T), but the Key Resources "
        "Table capture oligo (lines 1412–1414) prints 21 T followed by V; the model uses "
        "the printed short oligo. The MOSTA full capture-probe table instead prints 23 T followed by V (local HTML text 1276–1279), while its body prints the same short 21T–V oligo as Cell. Full probe and short ligation input are distinct molecular states; their tail-length difference is unresolved. The in-house Tn5 loading "
        "adapter and exact final junction remain unresolved, so the final construct "
        "is not constructed base by base. Physical circle and primer-binding geometry remain model questions.")
SPATIAL=feature("coordinate_identifier","spatial_barcode","unknown")
UMI_FEATURE=feature("molecular_identifier","umi","random")
SEQ_PRIMERS=()
FINAL_LIBRARY=None
SEQUENCING_ENDING="The article reports DNB generation and prints library-splint, MDA and read-primer bases. It assigns CID to Read 1 cycles 1–25 and MID to cycles 26–35; Read 2 is cDNA. Sample-specific run lengths refer to its Table S1, which remains under review. Exact physical circle, primer-binding sites and the source-matched Tn5 junction are not constructed here."


CID_PRECURSOR=("TGTGAGCCAAGGAGTTGNAACTGCTGACGTACTGAGAGGCATGGCGACCTTATCAG"+
               "N"*25+"TTGTCTTCCTAAGACCG")
FIXED_LIGATION_OLIGO="CTTGGCCTCCGACTTAAGTCGGATCGTAGCCATGTCGTTC"
CID_SPLINT="TCGGAGGCCAAGCGGTCTTAGGAA"
CAPTURE_OLIGO="TTGTCTTCCTAAGAC"+"N"*10+"T"*21+"V"
TSO="CTGCTGACGTACTGAGAGGC/rG//rG//iXNA_G/"
CDNA_PCR_PRIMER="CTGCTGACGTACTGAGAGGC"
LIBRARY_F="CTGCTGACGTACTGAGAGGCA"
LIBRARY_R="GAGACGTTCTCGACTCAGCAGA"
LIBRARY_SPLINT="GTACGTCAGCAGGAGACGTTCTCG"
READ1="CTGCTGACGTACTGAGAGGCATGGCGACCTTATCAG"
READ2="GCCATGTCGTTCTGTGAGCCAAGGAGTT"
MDA_PRIMER="TCTGCTGAGTCGAGAACGTC"

def oligos():
    return [
        f"CID precursor (25-nt random CID; preprint/MOSTA additional N): 5′-{CID_PRECURSOR}-3′",
        f"fixed oligo (5′ phosphate): 5′-{FIXED_LIGATION_OLIGO}-3′",
        f"CID/fixed-oligo splint: 5′-{CID_SPLINT}-3′",
        f"UMI/poly(dT) capture oligo (5′ phosphate): 5′-{CAPTURE_OLIGO}-3′",
        f"Stereo-seq TSO: 5′-{TSO}-3′",
        f"cDNA PCR primer: 5′-{CDNA_PCR_PRIMER}-3′",
        f"library F (5′ phosphate; final C and A phosphorothioated): 5′-/5phos/CTGCTGACGTACTGAGAGG*C*A-3′",
        f"library R: 5′-{LIBRARY_R}-3′",
        f"library splint: 5′-{LIBRARY_SPLINT}-3′",
        f"published read 1 primer: 5′-{READ1}-3′",
        f"published read 2 primer: 5′-{READ2}-3′",
        f"published MDA primer: 5′-{MDA_PRIMER}-3′",
    ]


def roles(name):
    return Construct([Segment(name, "[" + name + "]", placeholder=True)], name=name)


def capture_oligo():
    return Construct([Segment("capture handle", CAPTURE_OLIGO[:15]),
                      Segment("molecular identifier", "N" * 10, "umi", placeholder=True,
                              feature=UMI_FEATURE),
                      Segment("poly(dT)V; short-oligo source", "T" * 21 + "V", placeholder=True)],
                     name="Short phosphorylated capture oligo; not the full CID construct")


def protocol_graph():
    g = ProtocolGraph(TITLE, authority=SOURCE)
    dnb = g.input("decoded_dnb", "Predecoded CID DNA-nanoball array", MaterialKind.MOLECULAR,
                  state=roles("predecoded CID DNB array"), primary=True)
    capture = g.input("capture_extension", "Short phosphorylated capture oligo", MaterialKind.ADAPTER,
                      state=capture_oligo())
    ligase = g.input("surface_ligase", "T4 DNA ligase", MaterialKind.CATALYST)
    surface = g.product("capture_surface", "CID–MID–poly(dT) capture array", MaterialKind.MOLECULAR,
                        state=roles("CID-bearing surface with ligated capture oligo"))
    g.react("extend_surface", "Ligate capture oligo to decoded DNB", StepKind.LIGATION,
            ((dnb, InputRole.TEMPLATE), (capture, InputRole.STRUCTURAL), (ligase, InputRole.CATALYST)),
            (surface,), geometry=AdapterLigationGeometry.unavailable(
                "Both source oligos are known, but the exposed surface-DNB strand and complete nick are not constructed here. This is a model limitation."))
    rna = g.input("tissue_rna", "Polyadenylated tissue RNA", MaterialKind.MOLECULAR)
    captured = g.product("captured", "Surface-hybridized tissue RNA", MaterialKind.MOLECULAR)
    g.react("capture_rna", "Hybridize RNA to surface poly(dT)", StepKind.HYBRIDIZATION,
            ((surface, InputRole.TEMPLATE), (rna, InputRole.STRUCTURAL)), (captured,))
    rt = g.input("reverse_transcriptase", "SuperScript II", MaterialKind.CATALYST)
    cdna = g.product("cdna", "Surface-primed cDNA", MaterialKind.MOLECULAR)
    g.react("reverse_transcribe", "Reverse-transcribe captured RNA", StepKind.REVERSE_TRANSCRIPTION,
            ((captured, InputRole.TEMPLATE), (surface, InputRole.STRUCTURAL), (rt, InputRole.CATALYST)), (cdna,))
    tso = g.input("tso", "TSO; source terminal rG/rG/iXNA_G", MaterialKind.OLIGO,
                  state=Construct([Segment("TSO DNA handle", CDNA_PCR_PRIMER),
                                   Segment("source modified terminal G residues", "/rG//rG//iXNA_G/",
                                           placeholder=True, length_bp=3, note="Literal source modification tokens; not an unmodified DNA GGG tail")] ))
    switched = g.product("switched", "Template-switched cDNA", MaterialKind.MOLECULAR)
    g.react("template_switch", "Template switch with the published TSO", StepKind.TEMPLATE_SWITCHING,
            ((cdna, InputRole.TEMPLATE), (tso, InputRole.STRUCTURAL), (rt, InputRole.CATALYST)), (switched,))
    removal = g.input("tissue_removal", "Source tissue-removal buffer", MaterialKind.CONDITION)
    cleaned = g.product("cleaned", "Washed cDNA-containing chip", MaterialKind.MOLECULAR)
    tissue_waste = g.product("tissue_waste", "Removed tissue and wash material", MaterialKind.WASTE)
    g.react("remove_tissue", "Remove tissue and wash the cDNA-containing chip", StepKind.CLEANUP,
            ((switched, InputRole.TEMPLATE), (removal, InputRole.CONDITION)), (cleaned, tissue_waste))
    g.discard(tissue_waste, "Reported tissue-removal and wash steps.")
    exo = g.input("exo1", "Exonuclease I", MaterialKind.CATALYST)
    processed = g.product("processed", "Exonuclease-I-treated cDNA chip", MaterialKind.MOLECULAR)
    digest = g.product("digest", "Exonuclease-sensitive material", MaterialKind.WASTE)
    g.react("exo_treatment", "Treat the chip with Exonuclease I", StepKind.DIGESTION,
            ((cleaned, InputRole.TEMPLATE), (exo, InputRole.CATALYST)), (processed, digest))
    g.discard(digest, "Source digestion/washes; exact molecular identities are not resolved.")
    primer = g.input("cdna_primer", "Published cDNA PCR primer", MaterialKind.OLIGO,
                     state=Construct([Segment("cDNA PCR primer", CDNA_PCR_PRIMER)]))
    kapa = g.input("kapa", "KAPA HiFi HotStart ReadyMix polymerase", MaterialKind.CATALYST)
    amplified = g.product("amplified", "Amplified spatial cDNA", MaterialKind.MOLECULAR)
    g.react("amplify", "Amplify with the cDNA PCR primer", StepKind.PCR,
            ((processed, InputRole.TEMPLATE), (primer, InputRole.STRUCTURAL),
             (primer, InputRole.STRUCTURAL), (kapa, InputRole.CATALYST)),
            (amplified,), geometry=PCRGeometry.unavailable(roles("spatial cDNA duplex"),
                reason="The source primer is known; its complete duplex binding state remains unconstructed."))
    tn5 = g.input("tn5", "In-house Tn5 transposase", MaterialKind.CATALYST)
    loaded = g.input("loaded_adapter", "Source-matched Tn5 loading adapter; identity unassigned", MaterialKind.ADAPTER,
                     inferred=True)
    fragmented = g.product("fragmented", "Tn5-fragmented spatial cDNA", MaterialKind.MOLECULAR)
    g.react("fragment", "Fragment amplified cDNA with in-house Tn5", StepKind.TAGMENTATION,
            ((amplified, InputRole.TEMPLATE), (loaded, InputRole.STRUCTURAL), (tn5, InputRole.CATALYST)), (fragmented,))
    sds = g.input("sds", "SDS stop condition", MaterialKind.CONDITION)
    stopped = g.product("stopped", "Stopped fragmentation reaction", MaterialKind.MOLECULAR)
    g.react("stop", "Stop fragmentation with SDS", StepKind.QUENCHING,
            ((fragmented, InputRole.TEMPLATE), (sds, InputRole.CONDITION)), (stopped,))
    forward = g.input("library_f", "Library F; 5′P and terminal phosphorothioates", MaterialKind.OLIGO,
                      state=Construct([Segment("library F", LIBRARY_F)]))
    reverse = g.input("library_r", "Library R", MaterialKind.OLIGO,
                      state=Construct([Segment("library R", LIBRARY_R)]))
    linear = g.product("linear", "Amplified platform-library product", MaterialKind.MOLECULAR,
                       state=roles("amplified platform-library product"))
    g.react("library_pcr", "Amplify with the published library primers", StepKind.PCR,
            ((stopped, InputRole.TEMPLATE), (forward, InputRole.STRUCTURAL),
             (reverse, InputRole.STRUCTURAL), (kapa, InputRole.CATALYST)), (linear,),
            geometry=PCRGeometry.unavailable(roles("fragmented spatial cDNA"),
                reason="Loaded-adapter identity and complete primer-binding geometry are unassigned in this model; this is not a confirmed publication absence."))
    beads = g.input("cleanup_beads", "AMPure XP cleanup beads", MaterialKind.SELECTION)
    purified = g.product("purified", "Purified amplified library", MaterialKind.MOLECULAR)
    cleanup_waste = g.product("cleanup_waste", "Material excluded by bead cleanup", MaterialKind.WASTE)
    g.react("cleanup_library", "Purify the amplified library with beads", StepKind.CLEANUP,
            ((linear, InputRole.TEMPLATE), (beads, InputRole.SELECTION)), (purified, cleanup_waste))
    g.discard(cleanup_waste, "Source-reported AMPure cleanup; no exact excluded molecular species assigned.")
    library = g.product("library", "Reported sequencing DNB endpoint", MaterialKind.LIBRARY,
                        state=roles("sequencing DNB; physical circle unconstructed"))
    g.react("form_dnbs", "Reported DNB generation", StepKind.MAGIC,
            ((purified, InputRole.TEMPLATE),), (library,),
            reason="The published article reports DNB generation and prints a library splint/MDA primer. Complete strand selection, circularization and RCA geometry are not reconstructed; no generic platform adapter or polymerase is assigned by analogy.")
    g.sequence(library, SequencingPlan(limitation=SEQUENCING_ENDING))
    g.declare_reviewed_scope((
        "decoded_dnb", "capture_extension", "surface_ligase", "extend_surface",
        "tissue_rna", "capture_rna", "reverse_transcriptase", "reverse_transcribe",
        "tso", "template_switch", "tissue_removal", "remove_tissue", "exo1", "exo_treatment",
        "cdna_primer", "kapa", "amplify", "tn5", "loaded_adapter", "fragment", "sds", "stop",
        "library_f", "library_r", "library_pcr", "cleanup_beads", "cleanup_library", "form_dnbs", "library"), SOURCE)
    return g.require_publishable()
