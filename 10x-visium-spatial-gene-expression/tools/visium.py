"""Visium Spatial Gene Expression (fresh-frozen workflow) architecture."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
from batch_ngs import seg, spatial_rt_scene, truseq_library
from chemdraw import Row, feature
from adapter_ligation import AdapterLigationGeometry
from multimodal_spatial import finish_pcr_branch
from pcr import PCRGeometry
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, StepKind
from protocol_model import molecule, role

TITLE = "10x Visium Spatial Gene Expression — array capture"
NOTES = "01_visium.html"
SOURCE = 'Commercial protocol: 10x Genomics <a href="https://cdn.10xgenomics.com/image/upload/v1660261286/support-documents/CG000239_Visium_Spatial_Gene_Expression_User_Guide_Rev_F.pdf">CG000239 Rev F</a>.'
SUMMARY = "A tissue section is permeabilized over an array of spatially indexed oligo-dT spots. Released polyadenylated RNA is copied on the slide so each cDNA receives the spot barcode and a UMI before standard short-read library construction."
CAVEAT = "The spot barcode whitelist and complete surface oligo are commercial. Their published roles and lengths are drawn as placeholders; no undisclosed bases are invented."
SPATIAL=feature("spatial_barcode","spatial_barcode","whitelist",whitelist="10x-visium-spot-barcode-set")
UMI_FEATURE=feature("umi","umi","random")

FINAL_LIBRARY, SEQ_PRIMERS = truseq_library([
    seg("spatial barcode", "B" * 16, "cbc", placeholder=True, feature=SPATIAL),
    seg("UMI", "U" * 12, "umi", placeholder=True, feature=UMI_FEATURE),
    seg("poly(dT) junction", "T" * 12, placeholder=True),
    seg("captured cDNA", "X" * 38, placeholder=True)],
    "Visium gene-expression library", dual_index=False, inferred_adapters=True)
FINAL_CAPTION = "Read 1 contains the 16-nt spatial barcode followed by the 12-nt UMI; Read 2 enters the captured transcript. Adapter bases are canonical where the commercial guide specifies standard Illumina functions."
SEQUENCING_INTRO = "The 28-cycle Read 1 reports spot barcode plus UMI. Read 2 reports transcript sequence, and the i7 read identifies the sample library."

def sections():
    return [
        ("Place, image and permeabilize tissue over the capture array",
         [Row(chunks=[("tissue section:       [cell][cell][cell]", None, False)]),
          Row(chunks=[("                         ↓ released poly(A) RNA", None, False)]),
          Row(chunks=[("capture area:    ● ● ● ● ● spatially barcoded spots", "cbc", False)])],
         "Each printed spot carries many copies of one spatial-barcode species. Imaging preserves the mapping from spot coordinates to tissue morphology."),
        ("Capture RNA and reverse-transcribe on the slide",
         spatial_rt_scene(barcode_parts=(("16-nt spatial barcode",16),), umi=12,
                          surface="slide attachment").rows(),
         "The surface oligo contributes spatial barcode, UMI and oligo-dT. Reverse transcription copies the locally captured RNA while it remains associated with its spot."),
        ("Release cDNA and construct the sequencing library",
         [Row(chunks=[("barcoded cDNA → second strand → amplification → fragmentation", None, False)]),
          Row(chunks=[("→ end repair / dA-tail → adapter ligation → sample-index PCR", "me", False)])],
         "The spatial identity is already inside the cDNA before bulk recovery; downstream fragmentation must retain the barcode-bearing end used for Read 1."),
    ]

def protocol_graph():
    g=ProtocolGraph(TITLE,authority=SOURCE)
    tissue=g.input("tissue","imaged tissue section containing polyadenylated RNA",MaterialKind.BIOLOGICAL,state=molecule("tissue section",(role("cells"),role("poly(A) RNA"))),primary=True)
    permeabilization=g.input("permeabilization","tissue permeabilization condition",MaterialKind.CONDITION)
    released=g.product("released_rna","locally released poly(A) RNA",MaterialKind.MOLECULAR,state=molecule("released RNA",(role("transcript"),role("poly(A)"))))
    g.react("permeabilize","Permeabilize tissue over the capture array",StepKind.EXTRACTION,((tissue,InputRole.TEMPLATE),(permeabilization,InputRole.CONDITION)),(released,))
    surface=g.input("surface_oligo","surface-bound spatial-barcode/UMI oligo(dT)",MaterialKind.OLIGO,inferred=True,source=SOURCE+" (complete commercial surface oligo not disclosed)")
    rt=g.input("rt","reverse transcriptase",MaterialKind.CATALYST)
    surface_cdna=g.product("surface_cdna","spatial-barcode/UMI-tagged surface cDNA",MaterialKind.MOLECULAR,state=molecule("surface cDNA",(role("spatial barcode"),role("UMI"),role("captured cDNA"))))
    g.react("reverse_transcribe","Capture RNA and reverse-transcribe on the slide",StepKind.REVERSE_TRANSCRIPTION,((released,InputRole.TEMPLATE),(surface,InputRole.STRUCTURAL),(rt,InputRole.CATALYST)),(surface_cdna,))
    release=g.input("release_condition","surface-cDNA release condition",MaterialKind.CONDITION)
    cdna=g.product("released_cdna","released spatially barcoded cDNA",MaterialKind.MOLECULAR,state=surface_cdna.state)
    g.react("release","Release spatially barcoded cDNA from the array",StepKind.EXTRACTION,((surface_cdna,InputRole.TEMPLATE),(release,InputRole.CONDITION)),(cdna,))
    second_primer=g.input("second_primer","second-strand synthesis primer",MaterialKind.OLIGO,inferred=True,source=SOURCE+" (complete commercial sequence not disclosed)")
    second_pol=g.input("second_polymerase","second-strand polymerase",MaterialKind.CATALYST)
    duplex=g.product("duplex","double-stranded spatial cDNA",MaterialKind.MOLECULAR,state=cdna.state)
    g.react("second_strand","Synthesize the second cDNA strand",StepKind.EXTENSION,((cdna,InputRole.TEMPLATE),(second_primer,InputRole.STRUCTURAL),(second_pol,InputRole.CATALYST)),(duplex,))
    wf=g.input("wta_forward","cDNA preamplification forward primer",MaterialKind.OLIGO,inferred=True,source=SOURCE+" (complete commercial sequence not disclosed)")
    wr=g.input("wta_reverse","cDNA preamplification reverse primer",MaterialKind.OLIGO,inferred=True,source=SOURCE+" (complete commercial sequence not disclosed)")
    wp=g.input("wta_polymerase","cDNA preamplification polymerase",MaterialKind.CATALYST)
    amplified=g.product("amplified","amplified spatial cDNA",MaterialKind.MOLECULAR,state=duplex.state)
    g.react("wta","Amplify spatially barcoded cDNA",StepKind.PCR,((duplex,InputRole.TEMPLATE),(wf,InputRole.STRUCTURAL),(wr,InputRole.STRUCTURAL),(wp,InputRole.CATALYST)),(amplified,),geometry=PCRGeometry.unavailable(amplified.state,reason="The complete commercial preamplification primers are not disclosed."))
    frag=g.input("fragmentation","cDNA fragmentation condition",MaterialKind.CONDITION)
    fragments=g.product("fragments","fragmented barcode-bearing cDNA",MaterialKind.MOLECULAR,state=amplified.state)
    g.react("fragment","Fragment amplified cDNA while retaining barcode-bearing ends",StepKind.FRAGMENTATION,((amplified,InputRole.TEMPLATE),(frag,InputRole.CONDITION)),(fragments,))
    endmix=g.input("endprep","end-repair and dA-tailing enzyme mix",MaterialKind.CATALYST)
    repaired=g.product("repaired","end-repaired cDNA fragments",MaterialKind.MOLECULAR,state=fragments.state)
    g.react("end_repair","Repair fragmented cDNA ends",StepKind.END_REPAIR,((fragments,InputRole.TEMPLATE),(endmix,InputRole.CATALYST)),(repaired,))
    tailed=g.product("tailed","dA-tailed cDNA fragments",MaterialKind.MOLECULAR,state=repaired.state)
    g.react("da_tail","Add one 3-prime dA to each repaired end",StepKind.TAILING,((repaired,InputRole.TEMPLATE),(endmix,InputRole.CATALYST)),(tailed,))
    adapter=g.input("adapter","Visium library adapter",MaterialKind.ADAPTER,inferred=True,source=SOURCE+" (complete commercial sequence not disclosed)")
    ligase=g.input("ligase","library adapter ligase",MaterialKind.CATALYST)
    adapted=g.product("adapted","adapter-ligated spatial cDNA",MaterialKind.MOLECULAR,state=FINAL_LIBRARY,inferred=True,source=SOURCE+" (complete adapter boundary not disclosed)")
    g.react("ligate","Ligate sequencing adapters",StepKind.LIGATION,((tailed,InputRole.TEMPLATE),(adapter,InputRole.STRUCTURAL),(ligase,InputRole.CATALYST)),(adapted,),geometry=AdapterLigationGeometry.unavailable("The guide does not disclose the complete adapter-terminal geometry."))
    reviewed=["tissue","permeabilization","permeabilize","surface_oligo","rt","reverse_transcribe","release_condition","release","second_primer","second_polymerase","second_strand","wta_forward","wta_reverse","wta_polymerase","wta","fragmentation","fragment","endprep","end_repair","da_tail","adapter","ligase","ligate"]
    finish_pcr_branch(g,"visium",adapted,FINAL_LIBRARY,SEQ_PRIMERS,reviewed,label="Visium gene-expression",read_lengths={"Read 1":28},inferred_primers=True)
    g.declare_reviewed_scope(reviewed,SOURCE)
    return g.require_publishable()
