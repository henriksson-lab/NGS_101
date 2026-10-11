"""SureSelect XT HS2 target-capture model."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
from adapter_ligation import AdapterLigationGeometry
from batch_ngs import seg
from chemdraw import Construct, MolecularState, Row, Scene, Workflow
from endprep import dA_tailed_scene, repair_and_dA_tail
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
