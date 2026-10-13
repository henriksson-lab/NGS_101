"""SureSelect XT HS2 target-capture model."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
from adapter_ligation import AdapterLigationGeometry
from batch_ngs import seg
from chemdraw import Construct, MolecularState, Row, Scene, Workflow
from endprep import dA_tailed_scene, repair_and_dA_tail
from pcr import PCRGeometry
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind)
from targeted_ngs import hybrid_capture_scene, sureselect_library

TITLE = "SureSelect XT HS2 DNA target enrichment"
NOTES = "01_sureselect-xt-hs2.html"
SOURCE = 'Agilent <a href="https://www.agilent.com/cs/library/usermanuals/public/G9957-90000.pdf">SureSelect XT HS2 DNA Kits, G9957-90000 B0 (2026)</a>.'
SUMMARY = "A dual-indexed, optionally duplex-MBC-tagged Illumina library is pooled, hybridized to biotinylated RNA baits, recovered on streptavidin beads and amplified after capture."
CAVEAT = "The current manual publishes the molecular-barcode length and index sequences but not every adaptor-oligo base. Kit-specific hidden sequence is represented by its role; no sequence is fabricated."
FINAL_LIBRARY, SEQ_PRIMERS = sureselect_library()
FINAL_CAPTION = "Captured HS2 library in the MBC configuration. Each read begins with its five-base end tag; the two tags identify the original duplex family."
SEQUENCING_INTRO = "The manual specifies standard Illumina paired-end and dual-index primers with no custom sequencing primers."

def workflow():
    insert=Construct([seg("sheared genomic insert","X"*40,placeholder=True)])
    w=Workflow(MolecularState("sheared genomic DNA",tuple(Scene.duplex(list(insert)).rows())))
    tailed=repair_and_dA_tail(insert)
    w.react("End-repair and add 3′ dA",dA_tailed_scene(tailed).rows(),name="dA-tailed genomic DNA",note="The manual explicitly specifies the end-repair/dA-tailing step.")
    geometry=AdapterLigationGeometry.unavailable("Agilent's current SureSelect XT HS2 manual specifies HS2 adaptor ligation after dA-tailing but does not publish the adaptor oligo or its terminal bases; a complementary T is not fabricated here.")
    product=Scene.duplex(list(FINAL_LIBRARY),label="indexed HS2 library",unpaired=tuple(s.name+"'" for s in FINAL_LIBRARY if s.is_role_token()))
    w.ligate_adapter("Ligate the HS2 adaptor",geometry,product.rows(),name="HS2 adapter-ligated library",note="Optional five-base MBCs sit beside both insert ends; pre-capture PCR completes the flow-cell arms.")
    w.react("Hybridize the target interval to a biotinylated bait",hybrid_capture_scene().rows(),name="bait-bound target library",note="Adaptor blockers suppress arm-to-arm hybridization.")
    w.react("Capture, wash and amplify",[
            Row(chunks=[("bait—library duplex + streptavidin bead → retained", None, False)]),
            Row(chunks=[("unbound library → stringent wash → discarded", None, False)]),
            Row(chunks=[("bead-bound target library → post-capture PCR", None, False)])],name="captured HS2 library",note="Affinity selection changes abundance, not adapter orientation.")
    return w


def protocol_graph():
    graph=ProtocolGraph(TITLE,authority=SOURCE)
    dna=graph.input("genomic_dna","genomic DNA",MaterialKind.MOLECULAR,primary=True)
    fragmentation=graph.input("fragmentation_condition","mechanical or enzymatic fragmentation",
                              MaterialKind.CONDITION)
    fragmented=graph.product("fragmented","fragmented genomic DNA",MaterialKind.MOLECULAR)
    graph.react("fragment","Fragment genomic DNA",StepKind.FRAGMENTATION,
                ((dna,InputRole.TEMPLATE),(fragmentation,InputRole.CONDITION)),
                (fragmented,))
    endprep=graph.input("endprep_mix","end-repair and dA-tailing enzyme mix",
                        MaterialKind.CATALYST)
    tailed=graph.product("tailed","end-repaired 3-prime-dA-tailed DNA",
                         MaterialKind.MOLECULAR)
    graph.react("endprep","Repair fragment ends and add 3-prime dA",
                StepKind.END_REPAIR,
                ((fragmented,InputRole.TEMPLATE),(endprep,InputRole.CATALYST)),
                (tailed,))
    adapter=graph.input("hs2_adapter","SureSelect XT HS2 adapter mix",
                        MaterialKind.ADAPTER,inferred=True,
                        source=SOURCE+" (complete oligo bases undisclosed)")
    ligase=graph.input("adapter_ligase","HS2 adapter ligase",MaterialKind.CATALYST)
    ligated=graph.product("ligated","HS2 adapter-ligated molecules",
                          MaterialKind.MOLECULAR,inferred=True,
                          source=SOURCE+" (hidden adapter sequence)")
    geometry=AdapterLigationGeometry.unavailable(
        "Agilent specifies ligation after dA-tailing but does not publish the HS2 adapter terminal bases.")
    graph.react("ligate_adapter","Ligate the HS2 adapter mix",StepKind.LIGATION,
                ((tailed,InputRole.TEMPLATE),(adapter,InputRole.STRUCTURAL),
                 (ligase,InputRole.CATALYST)),(ligated,),geometry=geometry)
    i5=graph.input("i5_primer","HS2 i5 indexing PCR primer",MaterialKind.OLIGO)
    i7=graph.input("i7_primer","HS2 i7 indexing PCR primer",MaterialKind.OLIGO)
    polymerase=graph.input("pcr_polymerase","library PCR polymerase",MaterialKind.CATALYST)
    precapture=graph.product("precapture","indexed pre-capture HS2 library",
                             MaterialKind.MOLECULAR,state=FINAL_LIBRARY)
    graph.react("precapture_pcr","Complete indexed arms before capture",StepKind.PCR,
                ((ligated,InputRole.TEMPLATE),(i5,InputRole.STRUCTURAL),
                 (i7,InputRole.STRUCTURAL),(polymerase,InputRole.CATALYST)),
                (precapture,),geometry=PCRGeometry.unavailable(
                    FINAL_LIBRARY,reason="Index sequences are published, but the hidden adapter oligos prevent an exact precursor geometry."))
    bait=graph.input("rna_baits","biotinylated target-specific RNA baits",
                     MaterialKind.MOLECULAR)
    hybridized=graph.product("hybridized","bait-bound target library molecules",
                             MaterialKind.MOLECULAR,
                             state=tuple(hybrid_capture_scene().rows()))
    graph.react("hybridize","Hybridize target intervals to biotinylated RNA baits",
                StepKind.HYBRIDIZATION,
                ((precapture,InputRole.TEMPLATE),(bait,InputRole.STRUCTURAL)),
                (hybridized,))
    beads=graph.input("streptavidin","streptavidin capture beads",
                      MaterialKind.SELECTION)
    captured=graph.product("captured","captured target-enriched library",
                           MaterialKind.MOLECULAR)
    unbound=graph.product("unbound","unbound off-target library molecules",
                          MaterialKind.WASTE)
    graph.react("capture","Capture bait-bound molecules and wash off unbound library",
                StepKind.AFFINITY_CAPTURE,
                ((hybridized,InputRole.TEMPLATE),(beads,InputRole.SELECTION)),
                (captured,unbound))
    graph.discard(unbound,"unbound molecules are removed by stringent washes")
    post_f=graph.input("post_forward","post-capture PCR forward primer",MaterialKind.OLIGO)
    post_r=graph.input("post_reverse","post-capture PCR reverse primer",MaterialKind.OLIGO)
    library=graph.product("library",FINAL_LIBRARY.name,MaterialKind.LIBRARY,
                          state=FINAL_LIBRARY)
    graph.react("postcapture_pcr","Amplify the captured target library",StepKind.PCR,
                ((captured,InputRole.TEMPLATE),(post_f,InputRole.STRUCTURAL),
                 (post_r,InputRole.STRUCTURAL),(polymerase,InputRole.CATALYST)),
                (library,),geometry=PCRGeometry.unavailable(
                    FINAL_LIBRARY,reason="Post-capture PCR preserves the final construct, but the manual does not disclose every kit primer base."))
    graph.sequence(library,SequencingPlan(library=FINAL_LIBRARY,
                                          primers=tuple(SEQ_PRIMERS)))
    graph.declare_reviewed_scope((
        "genomic_dna","fragmentation_condition","fragment","endprep_mix","endprep",
        "hs2_adapter","adapter_ligase","ligate_adapter","i5_primer","i7_primer",
        "pcr_polymerase","precapture_pcr","rna_baits","hybridize","streptavidin",
        "capture","post_forward","post_reverse","postcapture_pcr","library"),SOURCE)
    return graph.require_publishable()
