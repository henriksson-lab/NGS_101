"""Slide-seqV2 bead-puck capture and library architecture."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2]/"lib"))
from batch_ngs import seg, truseq_library
from chemdraw import Row, feature
from multimodal_spatial import surface_capture_rows
from pcr import PCRGeometry
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind)

TITLE="Slide-seqV2"
NOTES="01_slide-seqv2.html"
SOURCE='Defining source: <a href="https://doi.org/10.1038/s41587-020-0739-1">Stickels et al., <i>Nature Biotechnology</i> (2021)</a>.'
SUMMARY="A lawn of randomly positioned 10-µm beads captures tissue mRNA. Each bead carries a 15-base split-pool barcode that is decoded in situ, an 8-base UMI and oligo-dT; template switching plus second-strand synthesis increases recovery before Nextera library preparation."
CAVEAT="The Methods print Chemgenes and two custom bead-oligo versions. This schematic represents the Chemgenes barcode–UMI–oligo-dT core; the full bead oligo is listed below. Its final duplex uses normalized Illumina arms and a variable transcript insert, so it is not a literal transcription of one sequenced molecule."
BEAD_CHEMGENES="TTTTTTTTCTACACGACGCTCTTCCGATCTJJJJJJJJTCTTCAGCGTTCCCGAGAJJJJJJJNNNNNNNNT30"
BEAD_CUSTOM_V1="TTT_PC_GCCGGTAATACGACTCACTATAGGGCTACACGACGCTCTTCCGATCTJJJJJJJJTCTTCAGCGTTCCCGAGAJJJJJJJTCNNNNNNNNT25"
BEAD_CUSTOM_V2="TTT_PC_GCCGGTAATACGACTCACTATAGGGCTACACGACGCTCTTCCGATCTJJJJJJJJTCTTCAGCGTTCCCGAGAJJJJJJNNNNNNNVVT30"
TSO="AAGCAGTGGTATCAACGCAGAGTGAATrG+GrG"
SECOND_STRAND="AAGCAGTGGTATCAACGCAGAGTGANNNGGNNNB"
TRUSEQ_PCR="CTACACGACGCTCTTCCGATCT"
SMART_PCR="AAGCAGTGGTATCAACGCAGAGT"
TRUSEQ_P5="AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT"
SPATIAL=feature("spatial_barcode","spatial_barcode","combinatorial")
UMI_FEATURE=feature("umi","umi","random")

LIB, PRIMERS=truseq_library([
    seg("bead barcode", "J"*15, "cbc", placeholder=True, feature=SPATIAL),
    seg("bead constant", "TCTTCAGCGTTCCCGAGA"),
    seg("UMI", "N"*8, "umi", placeholder=True, feature=UMI_FEATURE),
    seg("poly(dT) junction", "T"*30, placeholder=True),
    seg("captured cDNA", "X"*36, placeholder=True)],
    "Slide-seqV2 library", dual_index=False)
FINAL_LIBRARIES=(("Final Slide-seqV2 library",LIB,PRIMERS,
                  "The barcode-bearing read contains the bead barcode and UMI; the opposite read enters transcript-derived sequence.",
                  "The final duplex shows the Illumina primer sites used to read bead identity and transcript sequence."),)

def oligos():
    return [
        "<b>Chemgenes bead capture oligo</b> — 5′-"+BEAD_CHEMGENES+"-3′ (J: split-pool barcode; N: UMI)",
        "<b>Custom bead capture oligo vs1</b> — 5′-"+BEAD_CUSTOM_V1+"-3′ (PC: photocleavable linker)",
        "<b>Custom bead capture oligo vs2</b> — 5′-"+BEAD_CUSTOM_V2+"-3′ (V: A/C/G)",
        "<b>Template-switch oligo</b> — 5′-"+TSO+"-3′",
        "<b>Second-strand dN-SMRT oligo</b> — 5′-"+SECOND_STRAND+"-3′",
        "<b>TruSeq PCR handle</b> — 5′-"+TRUSEQ_PCR+"-3′",
        "<b>SMART PCR primer</b> — 5′-"+SMART_PCR+"-3′",
        "<b>TruSeq P5 hybrid primer</b> — 5′-"+TRUSEQ_P5+"-3′",
    ]

def sections():
    return [
        ("Decode the randomly packed bead puck",
         [Row(chunks=[("10-µm bead lawn:  ● ● ● ● ●", "cbc", False)]),
          Row(chunks=[("in-situ ligation sequencing → each 15-base bead barcode gets an (x,y) coordinate", None, False)])],
         "Fourteen of the fifteen split-pool barcode bases are decoded on the puck and mapped to bead coordinates."),
        ("Capture tissue RNA on spatially decoded beads",
         surface_capture_rows("bead", "15-base bead barcode", "tissue mRNA"),
         "Oligo-dT captures polyadenylated RNA at the bead beneath the tissue; the same oligo contributes bead barcode and UMI."),
        ("Copy, add the distal handle and amplify",
         [Row(chunks=[("first-strand RT + template switch → tissue digestion → bead recovery", None, False)]),
          Row(chunks=[("→ second-strand synthesis → PCR → Nextera tagmentation", "me", False)])],
         "Slide-seqV2 adds bead-bound second-strand synthesis after RT, then amplifies and converts cDNA to an Illumina library."),
    ]


def protocol_graph():
    final,primers=LIB,PRIMERS
    graph=ProtocolGraph(TITLE,authority=SOURCE)
    rna=graph.input("tissue_rna","polyadenylated tissue RNA",MaterialKind.MOLECULAR,
                    primary=True)
    bead=graph.input("bead_oligo","decoded bead barcode-UMI-oligo-dT capture oligo",
                     MaterialKind.OLIGO)
    captured=graph.product("captured","bead-captured tissue mRNA",MaterialKind.MOLECULAR)
    graph.react("capture","Hybridize tissue mRNA to spatially decoded bead oligo-dT",
                StepKind.HYBRIDIZATION,
                ((rna,InputRole.TEMPLATE),(bead,InputRole.STRUCTURAL)),(captured,))
    rt=graph.input("reverse_transcriptase","template-switching reverse transcriptase",
                   MaterialKind.CATALYST)
    cdna=graph.product("cdna","bead-barcode- and UMI-bearing first-strand cDNA",
                       MaterialKind.MOLECULAR)
    graph.react("reverse_transcribe","Copy bead-captured tissue RNA",
                StepKind.REVERSE_TRANSCRIPTION,
                ((captured,InputRole.TEMPLATE),(bead,InputRole.STRUCTURAL),
                 (rt,InputRole.CATALYST)),(cdna,))
    tso=graph.input("tso","Slide-seqV2 template-switch oligo "+TSO,MaterialKind.OLIGO,
                    source=SOURCE+" Supplementary Table 6")
    switched=graph.product("switched","template-switched spatial cDNA",
                           MaterialKind.MOLECULAR)
    graph.react("template_switch","Install the distal cDNA handle",
                StepKind.TEMPLATE_SWITCHING,
                ((cdna,InputRole.TEMPLATE),(tso,InputRole.STRUCTURAL),
                 (rt,InputRole.CATALYST)),(switched,))
    second_primer=graph.input("second_primer","dN-SMRT second-strand oligo "+SECOND_STRAND,
                              MaterialKind.OLIGO,source=SOURCE+" Supplementary Table 6")
    second_pol=graph.input("second_polymerase","second-strand polymerase",
                           MaterialKind.CATALYST)
    duplex=graph.product("duplex","double-stranded spatial cDNA",MaterialKind.MOLECULAR)
    graph.react("second_strand","Synthesize the bead-bound second strand",
                StepKind.EXTENSION,
                ((switched,InputRole.TEMPLATE),(second_primer,InputRole.STRUCTURAL),
                 (second_pol,InputRole.CATALYST)),(duplex,))
    fwd=graph.input("preamp_forward","TruSeq PCR handle "+TRUSEQ_PCR,MaterialKind.OLIGO,
                    source=SOURCE+" Supplementary Table 6")
    rev=graph.input("preamp_reverse","SMART PCR primer "+SMART_PCR,MaterialKind.OLIGO,
                    source=SOURCE+" Supplementary Table 6")
    preamp_pol=graph.input("preamp_polymerase","cDNA PCR polymerase",MaterialKind.CATALYST)
    amplified=graph.product("amplified","amplified spatial cDNA",MaterialKind.MOLECULAR)
    graph.react("preamp","Amplify spatial cDNA",StepKind.PCR,
                ((duplex,InputRole.TEMPLATE),(fwd,InputRole.STRUCTURAL),
                 (rev,InputRole.STRUCTURAL),(preamp_pol,InputRole.CATALYST)),
                (amplified,),geometry=PCRGeometry.unavailable(
                    final,reason="Both preamplification primers are published, but the variable tissue cDNA insert has no fixed sequence for an exact amplicon construct."))
    transposome=graph.input("transposome","Nextera loaded transposomes",
                            MaterialKind.ADAPTER)
    tn5=graph.input("tn5","Tn5 transposase",MaterialKind.CATALYST)
    tagged=graph.product("tagged","tagmented spatial cDNA",MaterialKind.MOLECULAR)
    waste=graph.product("tagment_waste","non-amplifiable tagmentation products",
                        MaterialKind.WASTE)
    graph.react("tagment","Tagment amplified spatial cDNA",StepKind.TAGMENTATION,
                ((amplified,InputRole.TEMPLATE),(transposome,InputRole.STRUCTURAL),
                 (tn5,InputRole.CATALYST)),(tagged,waste))
    graph.discard(waste,"non-amplifiable Tn5 outcomes are not sequenced")
    i5=graph.input("i5_primer","Nextera i5 indexing primer",MaterialKind.OLIGO)
    i7=graph.input("i7_primer","Nextera i7 indexing primer",MaterialKind.OLIGO)
    index_pol=graph.input("index_polymerase","index PCR polymerase",MaterialKind.CATALYST)
    library=graph.product("library",final.name,MaterialKind.LIBRARY,state=final)
    graph.react("index_pcr","Complete the Illumina library",StepKind.PCR,
                ((tagged,InputRole.TEMPLATE),(i5,InputRole.STRUCTURAL),
                 (i7,InputRole.STRUCTURAL),(index_pol,InputRole.CATALYST)),
                (library,),geometry=PCRGeometry.unavailable(
                    final,reason="The final sites are exact but index identities vary."))
    graph.sequence(library,SequencingPlan(library=final,primers=tuple(primers)))
    graph.declare_reviewed_scope((
        "tissue_rna","bead_oligo","capture","reverse_transcriptase",
        "reverse_transcribe","tso","template_switch","second_primer",
        "second_polymerase","second_strand","preamp_forward","preamp_reverse",
        "preamp_polymerase","preamp","transposome","tn5","tagment",
        "i5_primer","i7_primer","index_polymerase","index_pcr","library"),SOURCE)
    return graph.require_publishable()
