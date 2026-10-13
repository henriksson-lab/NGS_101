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
CAVEAT="The paper establishes coordinate-identifier and molecular-identifier roles, but complete STOmics chip and sequencing-primer bases are commercial. They are shown as inferred placeholders and no proprietary sequence is invented."
SPATIAL=feature("coordinate_identifier","spatial_barcode","unknown")
UMI_FEATURE=feature("molecular_identifier","umi","random")
SEQ_PRIMERS=()
FINAL_LIBRARY=None
SEQUENCING_ENDING="The completed library is circularized and amplified into sequencing DNA nanoballs. A DNBSEQ/cPAS primer anneals to the repeated platform-adapter site immediately before the coordinate-identifier/UMI-bearing insert. The current primer and adapter bases are proprietary, so the final panel shows the binding geometry and extension direction rather than a fabricated sequence."

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
         "INFERRED — A platform library strand is splint-ligated into a circle; undisclosed platform-adapter bases are visibly inferred."),
        ("Bind the sequencing primer at every repeat",
         [Row(chunks=[("…[primer site]—CID—MID—cDNA—[adapter][primer site]—CID—MID—cDNA…", "r1", False)]),
          Row(chunks=[("   sequencing primer →                    sequencing primer →", "r2", False)])],
         "INFERRED — The undisclosed cPAS primer binds the repeated adapter role and extends toward the coordinate/UMI-bearing insert."),]


def protocol_graph():
    graph=ProtocolGraph(TITLE,authority=SOURCE)
    dnb=graph.input("decoded_dnb","decoded coordinate-identifier DNA nanoball array",
                    MaterialKind.MOLECULAR,primary=True)
    capture_extension=graph.input("capture_extension","UMI and oligo-dT capture extension",
                                  MaterialKind.ADAPTER,inferred=True,
                                  source=SOURCE+" (commercial bases undisclosed)")
    ligase=graph.input("surface_ligase","surface capture-oligo ligase",
                       MaterialKind.CATALYST)
    capture_surface=graph.product("capture_surface","CID-UMI-oligo-dT capture array",
                                  MaterialKind.MOLECULAR,inferred=True,
                                  source=SOURCE+" (commercial bases undisclosed)")
    graph.react("extend_surface","Ligate UMI/oligo-dT onto decoded chip oligos",
                StepKind.LIGATION,
                ((dnb,InputRole.TEMPLATE),(capture_extension,InputRole.STRUCTURAL),
                 (ligase,InputRole.CATALYST)),(capture_surface,),
                geometry=AdapterLigationGeometry.unavailable(
                    "The paper establishes surface ligation but complete commercial oligo ends are undisclosed."))
    rna=graph.input("tissue_rna","polyadenylated tissue RNA",MaterialKind.MOLECULAR)
    captured=graph.product("captured","surface-captured tissue mRNA",MaterialKind.MOLECULAR)
    graph.react("capture_rna","Hybridize tissue mRNA to surface oligo-dT",
                StepKind.HYBRIDIZATION,
                ((capture_surface,InputRole.TEMPLATE),(rna,InputRole.STRUCTURAL)),
                (captured,))
    rt_primer=graph.input("rt_primer","surface oligo-dT RT primer",MaterialKind.OLIGO,
                          inferred=True,source=SOURCE+" (commercial bases undisclosed)")
    rt=graph.input("reverse_transcriptase","template-switching reverse transcriptase",
                   MaterialKind.CATALYST)
    cdna=graph.product("cdna","CID-MID-bearing spatial cDNA",MaterialKind.MOLECULAR)
    graph.react("reverse_transcribe","Copy captured tissue RNA",StepKind.REVERSE_TRANSCRIPTION,
                ((captured,InputRole.TEMPLATE),(rt_primer,InputRole.STRUCTURAL),
                 (rt,InputRole.CATALYST)),(cdna,))
    tso=graph.input("tso","template-switch oligo",MaterialKind.OLIGO,inferred=True,
                    source=SOURCE+" (commercial bases undisclosed)")
    switched=graph.product("switched","template-switched spatial cDNA",
                           MaterialKind.MOLECULAR)
    graph.react("template_switch","Install the distal cDNA handle",
                StepKind.TEMPLATE_SWITCHING,
                ((cdna,InputRole.TEMPLATE),(tso,InputRole.STRUCTURAL),
                 (rt,InputRole.CATALYST)),(switched,))
    fwd=graph.input("pcr_forward","spatial cDNA PCR forward primer",MaterialKind.OLIGO,
                    inferred=True,source=SOURCE+" (commercial bases undisclosed)")
    rev=graph.input("pcr_reverse","spatial cDNA PCR reverse primer",MaterialKind.OLIGO,
                    inferred=True,source=SOURCE+" (commercial bases undisclosed)")
    polymerase=graph.input("pcr_polymerase","cDNA amplification polymerase",
                           MaterialKind.CATALYST)
    amplified=graph.product("amplified","amplified CID-MID-bearing cDNA",
                            MaterialKind.MOLECULAR)
    graph.react("amplify","Amplify released spatial cDNA",StepKind.PCR,
                ((switched,InputRole.TEMPLATE),(fwd,InputRole.STRUCTURAL),
                 (rev,InputRole.STRUCTURAL),(polymerase,InputRole.CATALYST)),
                (amplified,),geometry=PCRGeometry.unavailable(
                    LINEAR,reason="Commercial amplification-primer bases are not disclosed."))
    platform_adapter=graph.input("platform_adapter","DNBSEQ platform-library adapter",
                                 MaterialKind.ADAPTER,inferred=True,
                                 source=SOURCE+" (commercial bases undisclosed)")
    linear=graph.product("linear","linear Stereo-seq sequencing strand",
                         MaterialKind.MOLECULAR,state=LINEAR,inferred=True,
                         source=SOURCE+" (platform adapter bases inferred)")
    graph.react("platform_convert","Fragment and install DNBSEQ platform adapters",
                StepKind.ADAPTER_ATTACHMENT,
                ((amplified,InputRole.TEMPLATE),(platform_adapter,InputRole.STRUCTURAL)),
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
        limitation="DNBSEQ/cPAS uses a proprietary primer at repeated platform-adapter sites; complete primer and adapter bases are unavailable."))
    graph.declare_reviewed_scope((
        "decoded_dnb","capture_extension","surface_ligase","extend_surface",
        "tissue_rna","capture_rna","rt_primer","reverse_transcriptase",
        "reverse_transcribe","tso","template_switch","pcr_forward","pcr_reverse",
        "pcr_polymerase","amplify","platform_adapter","platform_convert",
        "circle_ligase","circularize","phi29","form_dnbs","library"),SOURCE)
    return graph.require_publishable()
