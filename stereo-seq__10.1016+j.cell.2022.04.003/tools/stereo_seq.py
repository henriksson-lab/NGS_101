"""Stereo-seq DNA-nanoball patterned-array spatial transcriptomics."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"lib"))
from batch_ngs import seg
from adapter_ligation import AdapterLigationGeometry
from chemdraw import Construct, Row, feature
from circular import circularize_ssdna, rolling_circle
from multimodal_spatial import surface_capture_rows
from circular import CircularizationGeometry
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
        "Methods line 904 and the linked MOSTA detailed protocol primer table line 229 "
        "have TTGNAACT (an extra N). This model uses TTGNAACT "
        "because the author preprint and detailed protocol agree on it; which precursor "
        "was synthesized for the published experiment remains uncertain. Cell Methods "
        "(local text line 1684) call the capture tail 22-nt poly(T), but the Key Resources "
        "Table capture oligo (lines 1412–1414) prints 21 T followed by V; the model uses "
        "the printed oligo. The in-house Tn5 loading "
        "adapter and exact final junction remain unresolved, so the final construct "
        "retains inferred placeholders only for those regions.")
SPATIAL=feature("coordinate_identifier","spatial_barcode","unknown")
UMI_FEATURE=feature("molecular_identifier","umi","random")
SEQ_PRIMERS=()
FINAL_LIBRARY=None
SEQUENCING_ENDING="The completed library is circularized and amplified into sequencing DNA nanoballs. The published article prints Stereo-seq-read1 and -read2 primer bases and assigns CID to Read 1 cycles 1–25 and MID to cycles 26–35; Read 2 is cDNA. Their exact sites on a final library strand cannot be drawn base by base while the in-house Tn5 loading adapter and complete junction remain unidentified."

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

if CID_SPLINT != "TCGGAGGCCAAG"+"CGGTCTTAGGAA":
    raise ValueError("Stereo-seq splint transcription changed")

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
    ]

LINEAR=Construct([seg("DNBSEQ primer site","P"*18,"r1",placeholder=True,inferred=True),seg("coordinate identifier","C"*25,"cbc",placeholder=True,feature=SPATIAL),seg("molecular identifier","M"*10,"umi",placeholder=True,feature=UMI_FEATURE),seg("captured cDNA","X"*38,placeholder=True),seg("distal DNBSEQ adapter","Q"*18,placeholder=True,inferred=True)],name="Stereo-seq ssDNA library")
CIRCLE=circularize_ssdna(LINEAR,five_prime_phosphate=True)
RCA=rolling_circle(CIRCLE,"P"*18,copies=3)

def sections():
    return [
        ("Pattern and decode coordinate-identifier DNA nanoballs",
         [Row(chunks=[("DNB array:  ●CID₁  ●CID₂  ●CID₃  …", "cbc", False)]),
          Row(chunks=[("sequencing-by-synthesis before tissue placement → CID-to-(x,y) map", None, False)])],
         "Each surface DNB carries many copies of one coordinate identifier; pre-decoding maps molecular barcode to a physical location."),
        ("Add UMI/oligo-dT and capture tissue RNA",
         surface_capture_rows("decoded DNB", "coordinate identifier", "tissue mRNA"),
         "A UMI-bearing poly(dT) capture sequence is ligated onto the decoded chip oligo before tissue RNA capture."),
        ("Reverse-transcribe, release and amplify spatial cDNA",
         [Row(chunks=[("surface CID—MID—poly(dT) || RNA → reverse transcription + template switch", "umi", False)]),
          Row(chunks=[("→ release cDNA → PCR → fragmentation / platform-library conversion", None, False)])],
         "The coordinate identifier (CID) and molecular identifier (MID) remain attached to the cDNA through bulk library preparation."),
        ("Circularize the sequencing strand and form library DNBs",
         CIRCLE.rows(closure_label="splint ligation"),
         "INFERRED — A platform library strand is splint-ligated into a circle; the printed library primers and splint are listed above, but the in-house Tn5 adapter and final junction remain unresolved."),
        ("Bind the sequencing primer at every repeat",
         [Row(chunks=[("…[primer site]—CID—MID—cDNA—[adapter][primer site]—CID—MID—cDNA…", "r1", False)]),
          Row(chunks=[("   sequencing primer →                    sequencing primer →", "r2", False)])],
         "INFERRED — The published read primers are listed above; their exact placement on this unresolved final-junction drawing is inferred."),]


def protocol_graph():
    graph=ProtocolGraph(TITLE,authority=SOURCE)
    dnb=graph.input("decoded_dnb","decoded coordinate-identifier DNA nanoball array",
                    MaterialKind.MOLECULAR,
                    state=(seg("CID precursor",CID_PRECURSOR,"cbc",placeholder=True,
                               feature=SPATIAL),
                           seg("fixed ligation oligo",FIXED_LIGATION_OLIGO)),primary=True)
    capture_extension=graph.input("capture_extension","UMI and oligo-dT capture extension",
                                  MaterialKind.ADAPTER,
                                  state=(seg("capture handle",CAPTURE_OLIGO[:15]),
                                         seg("molecular identifier","N"*10,"umi",placeholder=True,
                                             feature=UMI_FEATURE),
                                         seg("poly(dT)V","T"*21+"V",placeholder=True)))
    ligase=graph.input("surface_ligase","surface capture-oligo ligase",
                       MaterialKind.CATALYST)
    capture_surface=graph.product("capture_surface","CID-UMI-oligo-dT capture array",
                                  MaterialKind.MOLECULAR)
    graph.react("extend_surface","Ligate UMI/oligo-dT onto decoded chip oligos",
                StepKind.LIGATION,
                ((dnb,InputRole.TEMPLATE),(capture_extension,InputRole.STRUCTURAL),
                 (ligase,InputRole.CATALYST)),(capture_surface,),
                geometry=AdapterLigationGeometry.unavailable(
                    "The source prints the capture oligo, but the decoded surface DNB's exposed ligation end is not represented as one complete linear strand."))
    rna=graph.input("tissue_rna","polyadenylated tissue RNA",MaterialKind.MOLECULAR)
    captured=graph.product("captured","surface-captured tissue mRNA",MaterialKind.MOLECULAR)
    graph.react("capture_rna","Hybridize tissue mRNA to surface oligo-dT",
                StepKind.HYBRIDIZATION,
                ((capture_surface,InputRole.TEMPLATE),(rna,InputRole.STRUCTURAL)),
                (captured,))
    rt_primer=graph.input("rt_primer","surface UMI/oligo-dT RT primer",MaterialKind.OLIGO,
                          state=(seg("capture handle",CAPTURE_OLIGO[:15]),
                                 seg("molecular identifier","N"*10,"umi",placeholder=True,
                                     feature=UMI_FEATURE),
                                 seg("poly(dT)V","T"*21+"V",placeholder=True)))
    rt=graph.input("reverse_transcriptase","template-switching reverse transcriptase",
                   MaterialKind.CATALYST)
    cdna=graph.product("cdna","CID-MID-bearing spatial cDNA",MaterialKind.MOLECULAR)
    graph.react("reverse_transcribe","Copy captured tissue RNA",StepKind.REVERSE_TRANSCRIPTION,
                ((captured,InputRole.TEMPLATE),(rt_primer,InputRole.STRUCTURAL),
                 (rt,InputRole.CATALYST)),(cdna,))
    tso=graph.input("tso","published Stereo-seq template-switch oligo",MaterialKind.OLIGO,
                    state=(seg("TSO DNA handle",CDNA_PCR_PRIMER),),
                    source=SOURCE+" (STAR Methods reagent table; modified terminal rG/rG/iXNA_G tail)")
    switched=graph.product("switched","template-switched spatial cDNA",
                           MaterialKind.MOLECULAR)
    graph.react("template_switch","Install the distal cDNA handle",
                StepKind.TEMPLATE_SWITCHING,
                ((cdna,InputRole.TEMPLATE),(tso,InputRole.STRUCTURAL),
                 (rt,InputRole.CATALYST)),(switched,))
    pcr_f=graph.input("pcr_f","published spatial cDNA PCR primer at one end",MaterialKind.OLIGO,
                      state=(seg("cDNA PCR primer",CDNA_PCR_PRIMER),),source=SOURCE)
    pcr_r=graph.input("pcr_r","same published cDNA PCR primer at opposite end",MaterialKind.OLIGO,
                      state=(seg("cDNA PCR primer",CDNA_PCR_PRIMER),),source=SOURCE)
    polymerase=graph.input("pcr_polymerase","cDNA amplification polymerase",
                           MaterialKind.CATALYST)
    amplified=graph.product("amplified","amplified CID-MID-bearing cDNA",
                            MaterialKind.MOLECULAR)
    graph.react("amplify","Amplify released spatial cDNA",StepKind.PCR,
                ((switched,InputRole.TEMPLATE),(pcr_f,InputRole.STRUCTURAL),
                 (pcr_r,InputRole.STRUCTURAL),
                 (polymerase,InputRole.CATALYST)),
                (amplified,),geometry=PCRGeometry.unavailable(
                    LINEAR,reason="The cDNA PCR primer sequence is printed, but the full amplified duplex is not represented base by base."))
    platform_adapter=graph.input("platform_adapter","in-house Tn5-loaded platform adapter",
                                 MaterialKind.ADAPTER,inferred=True,
                                 source=SOURCE+" (in-house Tn5 loading adapter and junction unresolved)")
    library_f=graph.input("library_f","published Stereo-seq library forward primer",
                          MaterialKind.OLIGO,state=(seg("library F",LIBRARY_F),),source=SOURCE)
    library_r=graph.input("library_r","published Stereo-seq library reverse primer",
                          MaterialKind.OLIGO,state=(seg("library R",LIBRARY_R),),source=SOURCE)
    linear=graph.product("linear","linear Stereo-seq sequencing strand",
                         MaterialKind.MOLECULAR,state=LINEAR,inferred=True,
                         source=SOURCE+" (platform adapter bases inferred)")
    graph.react("platform_convert","Tn5-fragment and amplify with published library primers",
                StepKind.ADAPTER_ATTACHMENT,
                ((amplified,InputRole.TEMPLATE),(platform_adapter,InputRole.STRUCTURAL),
                 (library_f,InputRole.STRUCTURAL),(library_r,InputRole.STRUCTURAL)),
                (linear,))
    circle_ligase=graph.input("circle_ligase","splint-assisted circle ligase",
                              MaterialKind.CATALYST)
    circle=graph.product("circle","circular Stereo-seq sequencing strand",
                         MaterialKind.MOLECULAR,state=CIRCLE)
    graph.react("circularize","Splint-ligate the sequencing strand into a circle",
                StepKind.CIRCULARIZATION,
                ((linear,InputRole.TEMPLATE),(circle_ligase,InputRole.CATALYST)),
                (circle,),geometry=CircularizationGeometry.known(
                    "splinted single-strand ligation",
                    statement="The platform library strand is splint-ligated into a closed circle before DNB formation."))
    phi29=graph.input("phi29","rolling-circle amplification polymerase",
                      MaterialKind.CATALYST)
    library=graph.product("library","Stereo-seq sequencing DNA nanoballs",
                          MaterialKind.LIBRARY,state=RCA,inferred=True,
                          source=SOURCE+" (platform adapter and priming-site bases inferred)")
    graph.react("form_dnbs","Amplify the circle into sequencing DNA nanoballs",
                StepKind.EXTENSION,
                ((circle,InputRole.TEMPLATE),(phi29,InputRole.CATALYST)),(library,))
    graph.sequence(library,SequencingPlan(
        limitation="Published read-primer sequences and R1/R2 roles are known, but the in-house Tn5 loading adapter and exact completed-library junction remain unresolved."))
    graph.declare_reviewed_scope((
        "decoded_dnb","capture_extension","surface_ligase","extend_surface",
        "tissue_rna","capture_rna","rt_primer","reverse_transcriptase",
        "reverse_transcribe","tso","template_switch","pcr_f","pcr_r",
        "pcr_polymerase","amplify","platform_adapter","library_f","library_r","platform_convert",
        "circle_ligase","circularize","phi29","form_dnbs","library"),SOURCE)
    return graph.require_publishable()
