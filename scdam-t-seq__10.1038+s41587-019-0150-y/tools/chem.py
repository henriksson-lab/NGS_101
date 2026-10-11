"""scDam&T-seq DamID and transcriptome branches."""
from adapter_ligation import AdapterLigationGeometry
from batch_ngs import seg
from chemdraw import MolecularState, Row, Scene, Workflow
from multibranch import illumina_branch
from protocol_extensions import analyte_split
TITLE="scDam&T-seq — protein–DNA contacts and transcriptome"
NOTES="01_scdam-t-seq.html"
SOURCE='Defining source: <a href="https://doi.org/10.1038/s41587-019-0150-y">Rooijers et al., <i>Nature Biotechnology</i> (2019)</a>; detailed protocol <a href="https://doi.org/10.1038/s41596-020-0314-8">Nat. Protocols (2020)</a>.'
SUMMARY="Dam fusion proteins methylate nearby GATC sites in vivo; mRNA is removed on oligo(dT) beads, leaving genomic DNA for DpnI-selective DamID amplification."
RNA,RP=illumina_branch("scDam&T-seq RNA","CEL-seq2 transcript insert",umi=6)
DNA,DP=illumina_branch("scDam&T-seq DamID","methylated-GATC genomic insert")
FINAL_LIBRARIES=(("Transcriptome library",RNA,RP,"Bead-bound RNA is processed by CEL-seq2.","Illumina primers on the RNA branch."),("DamID library",DNA,DP,"DpnI-selective fragments report protein-proximal adenine methylation.","Illumina primers on the DamID branch."))
def _damid_geometry():
    sc=Scene(); sc.strand("AdRt",[seg("unpaired AdRt tail","X"*20,placeholder=True),seg("paired blunt end","CGGCCGAGGA")],label="DamID AdRt")
    sc.anneal("AdRb",[seg("paired blunt end reverse complement","TCCTCGGCCG")],to="AdRt",pair=("paired blunt end reverse complement","paired blunt end"),label="DamID AdRb")
    sc.ligation_end("AdRt","paired blunt end",kind="blunt",owner="adapter",text="blunt DamID adapter end")
    return AdapterLigationGeometry.known("blunt",sc,statement="The detailed DamID protocol prints AdRt/AdRb; their ten-base paired terminus is blunt and DpnI supplies a blunt methylated-GATC substrate.")

def _workflow():
    w=Workflow(MolecularState("one cell",(Row(chunks=[("one cell expressing a Dam fusion",None,False)]),)))
    w.react("Mark protein-proximal GATC sites in vivo",[Row(chunks=[("DNA — G^mATC — DNA   ← Dam fusion at the chromatin target","w1",False)])],name="Dam-marked cell",note="Dam deposition records protein–DNA proximity.")
    w.react("Capture poly(A) RNA on magnetic beads",analyte_split("bead-bound RNA → CEL-seq2","supernatant DNA → DamID"),name="RNA and DNA branches",note="Plate position preserves cell identity.")
    cut=Scene.duplex([seg("DpnI-selected GATC fragment","GATC"+"X"*28+"GATC",placeholder=True)],label="blunt DpnI fragment")
    w.react("Cut methylated GATC with DpnI",cut.rows(),name="blunt DpnI fragments",note="DpnI creates amplifiable ends only at methylated GATC sites.")
    joined=Scene.duplex([seg("DamID adapter","X"*18,placeholder=True),seg("DpnI-selected genomic DNA","X"*28,placeholder=True)],label="adapter-ligated DamID fragment"); joined.junction("top","DamID adapter","DpnI-selected genomic DNA","blunt ligation")
    w.ligate_adapter("Ligate DamID adapters",_damid_geometry(),joined.rows(),name="adapter-ligated DamID DNA",note="Selective PCR follows from the restored adapter/GATC junctions.")
    return w

def sections():
    return [(r.action,list(r.after.rows),r.note) for r in _workflow().reactions]
