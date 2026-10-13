"""Quantum Barcoding (QBC) molecular model from the defining supplement."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
import illumina as il
import seqprimers as sp
from adapter_ligation import reported_adapter_end
from chemdraw import Construct, Scene, Segment, feature
from pcr import PCRGeometry, PCRPrimer
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind
from topology import InterruptedSplintTopology, SplintProductSegment

TITLE="Quantum Barcoding (QBC) — targeted RNA and protein indexing"
SOURCE="Roy et al., Communications Biology (2020), doi:10.1038/s42003-020-0896-2"
def seg(name,seq,tag=None,**kw): return Segment(name,seq,tag,**kw)

LEFT_PARTIAL_R1="CCCTACACGACGCTCTTCCGATCT"
SC3_RIGHT="GCTCACCAAGTCCTAGT"; SC2_LEFT="GCGTTACCC"; SC2_RIGHT="CCCATTCTC"
SCN_LEFT="GTCTTTGCC"; SCN_RIGHT="CCGTTGTGG"; SC1_LEFT="GCGACTTCC"; SC1_RIGHT="CCAACCACC"
ANCHOR="GCTCCCTGTCTGACG"
RIGHT_ARM="AGATCGGAAGAGCACACGTCTGAACTCCAGTCACGAACTCCATCAGC"
SPLINT6="CGTCAGACAGGGAGCGGTGGTTGGGGAAGTCGCCCACAACGGGGCAAAGACGAGAATGGGGGGTAACGCACTAGGACTTGGTGAGC"

UMI_F=feature("qbc_umi","umi","random")
def bc_feature(n): return feature(f"qbc_cell_bc_{n}","cell_barcode","whitelist",group="cell_barcode",part=f"round {n}",whitelist="QBC 60-well subcode set")
AHCA_F=feature("qbc_antibody_code","feature_barcode","whitelist",whitelist="QBC AHCA set")

def indexed_core(protein=True) -> Construct:
    parts=[seg("10-nt UMI","N"*10,"umi",placeholder=True,feature=UMI_F),
      seg("SC3 cell subcode","B"*7,"cbc",placeholder=True,feature=bc_feature(4)),seg("SC3 right arm",SC3_RIGHT),
      seg("SC2 left arm",SC2_LEFT),seg("SC2 cell subcode","C"*7,"cbc",placeholder=True,feature=bc_feature(3)),seg("SC2 right arm",SC2_RIGHT),
      seg("SCN left arm",SCN_LEFT),seg("SCN cell subcode","D"*7,"cbc",placeholder=True,feature=bc_feature(2)),seg("SCN right arm",SCN_RIGHT),
      seg("SC1 left arm",SC1_LEFT),seg("SC1 cell subcode","E"*7,"cbc",placeholder=True,feature=bc_feature(1)),seg("SC1 right arm",SC1_RIGHT),
      seg("QBC anchor",ANCHOR)]
    if protein: parts.append(seg("9-nt antibody-handle code","H"*9,"cbc",placeholder=True,feature=AHCA_F))
    else: parts.append(seg("target-specific cDNA","X"*40,placeholder=True))
    return Construct(parts, "QBC antibody reporter" if protein else "QBC targeted-RNA reporter")

def final_library(protein=True) -> Construct:
    payload=list(indexed_core(protein))
    return Construct([seg("P5",il.P5,"p5"),seg("i5 sample index","I"*8,"cbc",placeholder=True,
      feature=feature("qbc_i5","sample_index","random")),seg("Read 1 site",il.TRUSEQ_READ1,"r1"),
      *payload,seg("Read 2 site reverse complement",il.TRUSEQ_READ2[::-1].translate(str.maketrans("ACGT","TGCA")),"r2"),
      seg("i7 sample index reverse complement","J"*8,"cbc",placeholder=True,
          feature=feature("qbc_i7","sample_index","random")),seg("P7 reverse complement",il.P7_RC,"p7")],
      "QBC antibody Illumina library" if protein else "QBC targeted-RNA Illumina library")

def pcr_template(protein=True) -> Construct:
    return Construct([seg("partial Read 1 arm",LEFT_PARTIAL_R1,"r1"),
                      *list(indexed_core(protein)),
                      seg("Read 2-side arm",RIGHT_ARM,"r2")],
                     "QBC antibody pre-PCR template" if protein else "QBC RNA pre-PCR template")

PROTEIN_LIBRARY=final_library(True); RNA_LIBRARY=final_library(False)
PRIMERS=(sp.TRUSEQ["R1"],sp.TRUSEQ["I1"],sp.TRUSEQ["I2"],sp.TRUSEQ["R2"])
for _lib in (PROTEIN_LIBRARY,RNA_LIBRARY):
    _problems=sp.verify(_lib,PRIMERS)
    if _problems: raise ValueError(_lib.name+": "+"; ".join(_problems))

def splint_scene() -> Scene:
    sc=Scene(); sc.strand("Splint6",[seg("SC3 docking",SPLINT6[:24]),seg("C3 spacer 1","[C3]",placeholder=True),
      seg("SC2 docking",SPLINT6[24:42]),seg("C3 spacer 2","[C3]",placeholder=True),
      seg("SCN docking",SPLINT6[42:60]),seg("C3 spacer 3","[C3]",placeholder=True),
      seg("SC1 docking",SPLINT6[60:])],label="Splint6, written 5′→3′")
    sc.note("Splint6","one scaffold co-localizes the four successive barcode junctions")
    return sc


SPLINT6_CHUNKS = (
    ("anchor / SC1 chunk", SPLINT6[:24]),
    ("SC1 / SCN chunk", SPLINT6[24:42]),
    ("SCN / SC2 chunk", SPLINT6[42:60]),
    ("SC2 / SC3 chunk", SPLINT6[60:]),
)
SPLINT6_TOPOLOGY = InterruptedSplintTopology(
    "QBC four-round assembly on Splint6",
    SPLINT6_CHUNKS,
    (
        SplintProductSegment("SC3 PCR handle, UMI and cell subcode",
                             LEFT_PARTIAL_R1 + "N" * 10 + "N" * 7),
        SplintProductSegment("SC3 right arm", SC3_RIGHT, docking=True),
        SplintProductSegment("SC2 left arm", SC2_LEFT, docking=True,
                             junction_before=True),
        SplintProductSegment("SC2 cell subcode", "N" * 7),
        SplintProductSegment("SC2 right arm", SC2_RIGHT, docking=True),
        SplintProductSegment("SCN left arm", SCN_LEFT, docking=True,
                             junction_before=True),
        SplintProductSegment("SCN cell subcode", "N" * 7),
        SplintProductSegment("SCN right arm", SCN_RIGHT, docking=True),
        SplintProductSegment("SC1 left arm", SC1_LEFT, docking=True,
                             junction_before=True),
        SplintProductSegment("SC1 cell subcode", "N" * 7),
        SplintProductSegment("SC1 right arm", SC1_RIGHT, docking=True),
        SplintProductSegment("reporter anchor", ANCHOR, docking=True,
                             junction_before=True),
    ),
)

def ligation_geometry(round_no):
    return reported_adapter_end(f"QBC round-{round_no} barcode",kind="splinted",
      polarity="5-prime",terminal_sequence="N"*7,
      statement="Supplementary Data 2 prints each phosphorylated subcode oligo and Splint6-directed junction.",
      substrate_scenes=(splint_scene(),))

def pcr_geometry(template):
    f=PCRPrimer("NEBNext universal primer",il.NEBNEXT_UNIVERSAL_PRIMER,anneal_nt=24)
    r=PCRPrimer("NEBNext indexed primer",il.nebnext_i7_primer("N"*8),anneal_nt=len(il.TRUSEQ_READ2))
    return PCRGeometry.exact(template,f,r,statement="The paper specifies NEBNext E7600/E7335 indexed PCR; the supplement prints both internal primer-binding arms.")

def protocol_graph() -> ProtocolGraph:
    g=ProtocolGraph(TITLE,authority=SOURCE)
    cells=g.input("fixed_cells","fixed cells containing RNA and antibody-bound proteins",MaterialKind.BIOLOGICAL,primary=True)
    rtprimer=g.input("target_rt_primers","phosphorylated anchor–target-specific RT primers",MaterialKind.OLIGO)
    antibodies=g.input("dna_antibodies","AHCA-tagged antibodies",MaterialKind.OLIGO)
    rt=g.input("reverse_transcriptase","reverse transcriptase",MaterialKind.CATALYST)
    reporters=g.product("reporters","cells carrying targeted cDNA and antibody DNA reporters",MaterialKind.MOLECULAR,state=(indexed_core(False),indexed_core(True)))
    g.react("capture","Bind antibody reporters and reverse-transcribe selected RNAs",StepKind.REVERSE_TRANSCRIPTION,
      ((cells,InputRole.TEMPLATE),(rtprimer,InputRole.STRUCTURAL),(antibodies,InputRole.STRUCTURAL),(rt,InputRole.CATALYST)),(reporters,))
    current=reporters; reviewed=["fixed_cells","target_rt_primers","dna_antibodies","reverse_transcriptase","capture"]
    ligase=g.input("ligase","T4 DNA ligase",MaterialKind.CATALYST); splint=g.input("splint6","Splint6 scaffold",MaterialKind.OLIGO,state=splint_scene())
    reviewed += ["ligase","splint6"]
    for n,name in enumerate(("SC1","SCN","SC2","SC3"),1):
      ol=g.input(f"round_{n}_subcodes",f"{name} 60-member phosphorylated subcode plate",MaterialKind.ADAPTER)
      out=g.product(f"after_round_{n}",f"reporters after {name} split-pool ligation",MaterialKind.MOLECULAR,state=(indexed_core(False),indexed_core(True)))
      g.react(f"ligate_round_{n}",f"Split into 60 wells, ligate {name}, wash and pool",StepKind.LIGATION,
        ((current,InputRole.TEMPLATE),(ol,InputRole.STRUCTURAL),(splint,InputRole.STRUCTURAL),(ligase,InputRole.CATALYST)),(out,),geometry=ligation_geometry(n))
      current=out; reviewed += [f"round_{n}_subcodes",f"ligate_round_{n}"]
    rna_t=g.product("rna_template","barcoded targeted-RNA PCR template",MaterialKind.MOLECULAR,state=pcr_template(False))
    ab_t=g.product("protein_template","barcoded antibody-reporter PCR template",MaterialKind.MOLECULAR,state=pcr_template(True))
    g.react("separate_branches","Release reporters into targeted-RNA and antibody branches",StepKind.SPLIT,((current,InputRole.TEMPLATE),),(rna_t,ab_t),changes_molecule=False)
    reviewed += ["separate_branches"]
    fp=g.input("universal_pcr_primer","NEBNext universal PCR primer",MaterialKind.OLIGO)
    rp=g.input("indexed_pcr_primer","NEBNext indexed PCR primer",MaterialKind.OLIGO)
    reviewed += ["universal_pcr_primer","indexed_pcr_primer"]
    for key,temp,lib in (("rna",rna_t,RNA_LIBRARY),("protein",ab_t,PROTEIN_LIBRARY)):
      out=g.product(f"{key}_library",lib.name,MaterialKind.LIBRARY,state=lib)
      g.react(f"{key}_pcr",f"Indexed library PCR for the {key} branch",StepKind.PCR,
        ((temp,InputRole.TEMPLATE),(fp,InputRole.STRUCTURAL),(rp,InputRole.STRUCTURAL)),(out,),geometry=pcr_geometry(pcr_template(key=="protein")))
      roles=tuple(p.role for p in PRIMERS)
      g.sequence(out,SequencingPlan(library=lib,primers=PRIMERS,acquired_roles=roles,
        read_uses={"Read 1":"include","Index 1 (i7)":"extract","Index 2 (i5)":"extract","Read 2":"include"}))
      reviewed.extend((f"{key}_pcr",f"{key}_library"))
    g.declare_reviewed_scope(tuple(reviewed),SOURCE)
    return g.require_publishable()
