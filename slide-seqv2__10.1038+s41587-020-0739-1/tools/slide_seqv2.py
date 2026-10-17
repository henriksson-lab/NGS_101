"""Slide-seqV2 bead-puck capture and library architecture."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2]/"lib"))
import illumina as il
import nextera as nx
import seqprimers as sp
from batch_ngs import seg
from chemdraw import Construct, Row, feature
from multimodal_spatial import surface_capture_rows
from pcr import PCRGeometry
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            SourceClaim, SourceDiscrepancy, StepKind)

TITLE="Slide-seqV2"
NOTES="01_slide-seqv2.html"
SOURCE='Defining source: <a href="https://doi.org/10.1038/s41587-020-0739-1">Stickels et al., <i>Nature Biotechnology</i> (2021)</a>.'
SUMMARY="A lawn of randomly positioned 10-µm beads captures tissue mRNA. Each bead carries a 15-base split-pool barcode that is decoded in situ, an 8-base UMI and oligo-dT; template switching plus second-strand synthesis increases recovery before Nextera library preparation."
CAVEAT="The Methods print Chemgenes and two custom bead-oligo versions. This schematic represents the Chemgenes J8–bridge–J7–UMI–oligo-dT core; the full bead oligos are listed below. The final duplex is the source-defined mixed endpoint (Truseq5 P5/TruSeq Read 1 arm on the bead end, Nextera XT mosaic end with an N7xx index on the tagmented end) with a variable transcript insert, so it is not a literal transcription of one sequenced molecule."
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

BEAD_BRIDGE="TCTTCAGCGTTCCCGAGA"


def _final_library():
    """Truseq5 x Nextera N7xx endpoint (Online Methods 'PCR cleanup and Nextera Tagmentation')."""
    lib=Construct([
        seg("P5 (Truseq5 primer 5' part)", il.P5[:-4], "p5"),
        seg("TruSeq Read 1 arm (bead PCR handle at its 3' end)", il.TRUSEQ_READ1, "r1"),
        seg("bead barcode block 1", "J"*8, "cbc", placeholder=True, feature=SPATIAL),
        seg("bead constant (UP site)", BEAD_BRIDGE),
        seg("bead barcode block 2", "J"*7, "cbc", placeholder=True, feature=SPATIAL),
        seg("UMI", "N"*8, "umi", placeholder=True, feature=UMI_FEATURE),
        seg("poly(dT)", "T"*30),
        seg("captured cDNA", "X"*36, placeholder=True),
        seg("Nextera mosaic end reverse complement", nx.ME_RC, "me"),
        seg("S7 reverse complement", nx.S7_RC, "s7"),
        seg("i7 reverse complement (N7xx)", "I"*8, "cbc", placeholder=True,
            feature=feature("slide_seqv2.i7", "sample_index", "whitelist",
                            whitelist="Nextera N700 index primers")),
        seg("P7 reverse complement", il.P7_RC, "p7"),
    ], name="Slide-seqV2 library")
    if il.P5[:-4]+il.TRUSEQ_READ1!=TRUSEQ_P5:
        raise AssertionError("bead-end arm differs from the published Truseq P5 primer")
    primers=(sp.TRUSEQ["R1"], sp.NEXTERA["I1"], sp.NEXTERA["R2"])
    problems=sp.verify(lib, primers, required_roles=("Read 1", "Index 1 (i7)", "Read 2"))
    if problems:
        raise ValueError("invalid Slide-seqV2 endpoint: "+"; ".join(problems))
    return lib, primers


LIB, PRIMERS=_final_library()
READ_LENGTHS={"Read 1": 42, "Index 1 (i7)": 8, "Read 2": 50}
READ_USES={"Read 1": "extract", "Index 1 (i7)": "demultiplex", "Read 2": "include"}
FINAL_LIBRARIES=(("Final Slide-seqV2 library",LIB,PRIMERS,
                  "TruSeq Read 1 (42 cycles) reads the bead barcode blocks, constant and UMI; Nextera Read 2 (50 cycles) enters transcript sequence from the tagmented end; i7 (8 cycles) reads the N7xx sample index.",
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
        "<b>Truseq P5 (Truseq5) library primer</b> — 5′-"+TRUSEQ_P5+"-3′ (with Nextera XT N700 index primers)",
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
    puck=graph.input("puck","randomly deposited 10-µm barcoded bead puck",
                     MaterialKind.BIOLOGICAL,primary=True)
    decoders=graph.input("decoders","TruSeq/UP hybridization primer series and 5(base)/3(base) "
                         "fluorescent ligation oligos (Supplementary Table 6); T4 DNA ligase",
                         MaterialKind.OLIGO,source=SOURCE+" Online Methods 'Puck Sequencing'; Supplementary Table 6")
    mapped=graph.product("mapped_puck","bead puck with 14 of 15 barcode bases mapped to x–y",
                         MaterialKind.MOLECULAR)
    graph.react("decode","Decode bead barcodes by 8 five-prime, 3 three-prime and 3 SEDAL "
                "ligations, formamide-stripping between primers",StepKind.HYBRIDIZATION,
                ((puck,InputRole.TEMPLATE),(decoders,InputRole.STRUCTURAL)),(mapped,),
                changes_molecule=False)
    rna=graph.input("tissue_rna","polyadenylated RNA from the 10-µm tissue section",
                    MaterialKind.MOLECULAR)
    hyb=graph.input("hyb_buffer","6× SSC, 2 U/µl RNase inhibitor, 30 min room temperature",
                    MaterialKind.CONDITION)
    captured=graph.product("captured","bead-captured tissue mRNA",MaterialKind.MOLECULAR)
    graph.react("capture","Hybridize tissue mRNA to spatially decoded bead oligo-dT",
                StepKind.HYBRIDIZATION,
                ((mapped,InputRole.TEMPLATE),(rna,InputRole.STRUCTURAL),
                 (hyb,InputRole.CONDITION)),(captured,))
    rt=graph.input("reverse_transcriptase","Maxima H− RT in Maxima RT buffer, dNTPs, RNase "
                   "inhibitor; 1.5 h at 52 °C",MaterialKind.CATALYST)
    bead=graph.input("bead_oligo","bead-tethered barcode–UMI–oligo-dT capture oligo (RT primer)",
                     MaterialKind.OLIGO)
    cdna=graph.product("cdna","bead-barcode- and UMI-bearing first-strand cDNA",
                       MaterialKind.MOLECULAR)
    graph.react("reverse_transcribe","Copy bead-captured tissue RNA from the bead oligo-dT",
                StepKind.REVERSE_TRANSCRIPTION,
                ((captured,InputRole.TEMPLATE),(bead,InputRole.STRUCTURAL),
                 (rt,InputRole.CATALYST)),(cdna,))
    tso=graph.input("tso","Slide-seqV2 template-switch oligo "+TSO+" (50 µM, 10 µl in 200 µl RT)",
                    MaterialKind.OLIGO,source=SOURCE+" Supplementary Table 6")
    switched=graph.product("switched","template-switched spatial cDNA (where switching succeeded)",
                           MaterialKind.MOLECULAR)
    graph.react("template_switch","Install the distal SMART handle in the same RT reaction",
                StepKind.TEMPLATE_SWITCHING,
                ((cdna,InputRole.TEMPLATE),(tso,InputRole.STRUCTURAL),
                 (rt,InputRole.CATALYST)),(switched,))
    digest=graph.input("tissue_digestion","2× tissue digestion buffer (SDS, EDTA, proteinase K), "
                       "37 °C 30 min",MaterialKind.CATALYST)
    digested=graph.product("digested","bead cDNA in digested tissue/RT mixture",
                           MaterialKind.MOLECULAR)
    graph.react("digest_tissue","Digest the tissue in the RT mixture",StepKind.DEPROTEINIZATION,
                ((switched,InputRole.TEMPLATE),(digest,InputRole.CATALYST)),(digested,),
                changes_molecule=False)
    wash1=graph.input("bead_washes","beads pipetted off, glass discarded; three wash-buffer "
                      "spins at 3000 RCF",MaterialKind.SELECTION)
    washed=graph.product("washed","washed cDNA-bearing beads",MaterialKind.MOLECULAR)
    wash_waste=graph.product("wash_waste","glass substrate, digested tissue and RT reagents",
                             MaterialKind.WASTE)
    graph.react("wash_beads","Release beads and wash by centrifugation",StepKind.CLEANUP,
                ((digested,InputRole.TEMPLATE),(wash1,InputRole.SELECTION)),(washed,wash_waste))
    graph.discard(wash_waste,"glass discarded and supernatants removed")
    exo=graph.input("exo_i","Exonuclease I (NEB M0568) in ExoI buffer, 37 °C 50 min, then three washes",
                    MaterialKind.CATALYST)
    exo_beads=graph.product("exo_beads","ExoI-treated cDNA beads",MaterialKind.MOLECULAR)
    exo_waste=graph.product("exo_waste","ExoI digestion products and washes",MaterialKind.WASTE)
    graph.react("exo_digest","Treat beads with Exonuclease I",StepKind.DIGESTION,
                ((washed,InputRole.TEMPLATE),(exo,InputRole.CATALYST)),(exo_beads,exo_waste))
    graph.discard(exo_waste,"removed in the three post-ExoI wash spins")
    naoh=graph.input("naoh","0.1 N NaOH, 5 min room temperature; quenched with wash buffer, three washes",
                     MaterialKind.CONDITION)
    ss_beads=graph.product("ss_beads","bead-bound single-stranded first-strand cDNA",
                           MaterialKind.MOLECULAR)
    released=graph.product("released","released RNA template and non-covalently paired strands",
                           MaterialKind.WASTE)
    graph.react("denature","Denature with NaOH and wash",StepKind.STRAND_SEPARATION,
                ((exo_beads,InputRole.TEMPLATE),(naoh,InputRole.CONDITION)),(ss_beads,released))
    graph.discard(released,"removed with the quench and wash supernatants")
    second_primer=graph.input("second_primer","dN-SMRT second-strand oligo "+SECOND_STRAND+
                              " (2 µl of 1 mM in 200 µl)",
                              MaterialKind.OLIGO,source=SOURCE+" Supplementary Table 6")
    second_pol=graph.input("second_polymerase","Klenow fragment (NEB M0210) in Maxima RT buffer "
                           "with dNTPs, 37 °C 1 h, then three washes",MaterialKind.CATALYST)
    duplex=graph.product("duplex","bead-bound double-stranded spatial cDNA with a SMART end",
                         MaterialKind.MOLECULAR)
    graph.react("second_strand","Randomly prime and extend the bead-bound second strand",
                StepKind.EXTENSION,
                ((ss_beads,InputRole.TEMPLATE),(second_primer,InputRole.STRUCTURAL),
                 (second_pol,InputRole.CATALYST)),(duplex,))
    fwd=graph.input("preamp_forward","TruSeq PCR handle "+TRUSEQ_PCR,MaterialKind.OLIGO,
                    source=SOURCE+" Supplementary Table 6")
    rev=graph.input("preamp_reverse","SMART PCR primer "+SMART_PCR,MaterialKind.OLIGO,
                    source=SOURCE+" Supplementary Table 6")
    preamp_pol=graph.input("preamp_polymerase","Terra PCR Direct polymerase mix (Takara 639270); "
                           "95 °C 3 min, 4 cycles (65 °C anneal), 9 cycles (67 °C anneal), 72 °C 5 min",
                           MaterialKind.CATALYST)
    amplified=graph.product("amplified","amplified spatial cDNA",MaterialKind.MOLECULAR)
    graph.react("preamp","Amplify spatial cDNA on beads",StepKind.PCR,
                ((duplex,InputRole.TEMPLATE),(fwd,InputRole.STRUCTURAL),
                 (rev,InputRole.STRUCTURAL),(preamp_pol,InputRole.CATALYST)),
                (amplified,),geometry=PCRGeometry.unavailable(
                    final,reason="Both preamplification primers are published, but the variable tissue cDNA insert has no fixed sequence for an exact amplicon construct."))
    ampure=graph.input("ampure_cdna","AMPure XP 30 µl per 50 µl PCR (0.6×), performed twice; "
                       "final 10 µl",MaterialKind.SELECTION)
    clean=graph.product("clean_cdna","purified amplified cDNA (600 pg taken forward)",
                        MaterialKind.MOLECULAR)
    ampure_waste=graph.product("ampure_cdna_waste","0.6× AMPure supernatants and washes",
                               MaterialKind.WASTE)
    graph.react("cdna_cleanup","Purify amplified cDNA twice with 0.6× AMPure XP",StepKind.CLEANUP,
                ((amplified,InputRole.TEMPLATE),(ampure,InputRole.SELECTION)),(clean,ampure_waste))
    graph.discard(ampure_waste,"removed with the bead supernatants")
    transposome=graph.input("transposome","Nextera XT tagmentation reagent (FC-131-1096)",
                            MaterialKind.ADAPTER)
    tn5=graph.input("tn5","Tn5 transposase (in Nextera XT)",MaterialKind.CATALYST)
    tagged=graph.product("tagged","tagmented spatial cDNA",MaterialKind.MOLECULAR)
    waste=graph.product("tagment_waste","internal and SMART-end fragments without the bead TruSeq handle",
                        MaterialKind.WASTE)
    graph.react("tagment","Tagment 600 pg amplified cDNA (Nextera XT, manufacturer's instructions)",
                StepKind.TAGMENTATION,
                ((clean,InputRole.TEMPLATE),(transposome,InputRole.STRUCTURAL),
                 (tn5,InputRole.CATALYST)),(tagged,waste))
    graph.discard(waste,"not exponentially amplified by Truseq P5 + N7xx; carried into the 0.6× "
                  "AMPure cleanup rather than physically separated by sequence")
    p5=graph.input("p5_primer","Truseq P5 primer "+TRUSEQ_P5,MaterialKind.OLIGO,
                   source=SOURCE+" Supplementary Table 6")
    i7=graph.input("i7_primer","Nextera XT N700-series i7 index primer",MaterialKind.OLIGO)
    index_pol=graph.input("index_polymerase","Nextera XT PCR mix; 72 °C 3 min, 95 °C 30 s, "
                          "12 cycles (95/55/72 °C), 72 °C 5 min",MaterialKind.CATALYST)
    pcr_lib=graph.product("pcr_library","amplified Slide-seqV2 library before cleanup",
                          MaterialKind.MOLECULAR,state=final)
    graph.react("index_pcr","Amplify with Truseq P5 and an N7xx index primer",StepKind.PCR,
                ((tagged,InputRole.TEMPLATE),(p5,InputRole.STRUCTURAL),
                 (i7,InputRole.STRUCTURAL),(index_pol,InputRole.CATALYST)),
                (pcr_lib,),geometry=PCRGeometry.unavailable(
                    final,reason="The final sites are exact but N7xx index identities vary."))
    ampure2=graph.input("ampure_library","AMPure XP 0.6× (30 µl to 50 µl); 10 µl water",
                        MaterialKind.SELECTION)
    library=graph.product("library",final.name,MaterialKind.LIBRARY,state=final)
    lib_waste=graph.product("library_ampure_waste","0.6× AMPure supernatant, primers and washes",
                            MaterialKind.WASTE)
    graph.react("library_cleanup","Purify the library with 0.6× AMPure XP",StepKind.CLEANUP,
                ((pcr_lib,InputRole.TEMPLATE),(ampure2,InputRole.SELECTION)),(library,lib_waste))
    graph.discard(lib_waste,"removed with the bead supernatant")
    graph.sequence(library,SequencingPlan(library=final,primers=tuple(primers),
                                          read_lengths=READ_LENGTHS,read_uses=READ_USES))
    graph.annotate(SourceDiscrepancy("decode","Custom bead version 2 barcode length",
        (SourceClaim("Online Methods 'Barcoded beads', JATS text line 105 (vs2 oligo)",
                     "JJJJJJJJ…JJJJJJ: 8 + 6 split-pool J bases, then 7 N and VV before T30"),
         SourceClaim("Online Methods 'Bead synthesis', JATS text line 111",
                     "split-pool procedure repeated 15 times (two blocks of 8 and 7 cycles)")),
        "Model only the Chemgenes oligo (J8 + J7); vs1/vs2 are listed as printed.",
        "Both passages were read in the PMC author-manuscript JATS; the Chemgenes and vs1 oligos print 8 + 7 J.",
        "Whether vs2 has 14 split-pool bases or the printed oligo omits one J is unresolved; "
        "the publisher PDF was not available for rendered comparison."))
    graph.declare_reviewed_scope((
        "puck","decoders","decode","tissue_rna","hyb_buffer","capture","bead_oligo","reverse_transcriptase",
        "reverse_transcribe","tso","template_switch","tissue_digestion","digest_tissue",
        "bead_washes","wash_beads","exo_i","exo_digest","naoh","denature","second_primer",
        "second_polymerase","second_strand","preamp_forward","preamp_reverse",
        "preamp_polymerase","preamp","ampure_cdna","cdna_cleanup","transposome","tn5","tagment",
        "p5_primer","i7_primer","index_polymerase","index_pcr","ampure_library",
        "library_cleanup","library"),SOURCE)
    return graph.require_publishable()
