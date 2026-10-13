"""PLAYR proximity-ligation/RCA assay model."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"lib"))
from chemdraw import Construct, MolecularState, Row, Scene, Segment, Workflow, circle_rows, duplex_rows, revcomp
from circular import circularize_ssdna
from topology import ProximityCircleTopology, ProximityProbe

TITLE="PLAYR — proximity ligation assay for RNA"
SOURCE=('Defining source: <a href="https://doi.org/10.1038/nmeth.3742">Frei et al., '
        '<i>Nature Methods</i> (2016)</a>; reaction detail from '
        '<a href="https://doi.org/10.1038/s41596-018-0120-8">Frei et al., '
        '<i>Nature Protocols</i> (2019)</a>.')
LINKER="A"*10
PROBE1_INSERT_SITE="CTCAGTCGT"
PROBE1_BACKBONE_SITE="GACACTCTT"
PROBE2_BACKBONE_SITE="GACGCTAAT"
PROBE2_INSERT_SITE="ATCGTGACC"
PROBE1_TAIL=LINKER+PROBE1_INSERT_SITE+PROBE1_BACKBONE_SITE
PROBE2_TAIL=LINKER+PROBE2_BACKBONE_SITE+PROBE2_INSERT_SITE
INSERT="ACGACTGAGTTTGGTCACGAT"
BACKBONE="ATTAGCGTCCAGTGAATGCGAGTCCGTCTAGGAGAGTAGTACAGCAGCCGTCAAGAGTGTC"
DETECTION=INSERT

def seg(name,seq,tag=None,**kw): return Segment(name,seq,tag,**kw)

TOPOLOGY=ProximityCircleTopology(
  upstream=ProximityProbe("Probe 2","N"*24,LINKER,"outer",PROBE2_BACKBONE_SITE,
                          "inner",PROBE2_INSERT_SITE),
  downstream=ProximityProbe("Probe 1","N"*24,LINKER,"inner",PROBE1_INSERT_SITE,
                            "outer",PROBE1_BACKBONE_SITE),
  inner=INSERT,outer=BACKBONE,inner_spacer="TTT",outer_body=BACKBONE[9:-9],
  name="PLAYR adjacent probes and ligation-ready reporter circle")

def reporter_circle_template() -> Construct:
  return Construct([seg("target-specific insert",INSERT,"cbc"),seg("universal backbone",BACKBONE,"r3")],"PLAYR system-1 reporter circle")

CIRCLE=circularize_ssdna(reporter_circle_template(),five_prime_phosphate=True)

def rolony() -> Construct:
  parts=[]
  for i in range(1,4):
    parts.extend([seg(f"insert complement {i}",revcomp(INSERT),"cbc"),
                  seg(f"backbone complement {i}",revcomp(BACKBONE),"r3")])
  return Construct(parts,"PLAYR RCA concatemer (three repeats drawn)")

def detection_scene() -> Scene:
  product=rolony(); target="insert complement 2"
  sc=Scene(); sc.strand("rolony",list(product),label="phi29 RCA concatemer")
  # Detection oligo is identical to the published insert; its complement occurs in the RCA product.
  sc.anneal("detection oligo",[seg("system-1 detection sequence",DETECTION,"umi")],to="rolony",pair=("system-1 detection sequence",target),label="fluorescent or metal-tagged detection oligo")
  sc.arrow("detection oligo","measure reporter abundance by flow or mass cytometry")
  return sc

def workflow() -> Workflow:
  w=Workflow(MolecularState("adjacent PLAYR probes on target RNA",
    (Row(chunks=[("Topology above: Probe 2 and Probe 1 expose the complete split docking platform",None,False)]),)),track_graph=False)
  w.react("Hybridize phosphorylated insert and universal backbone; T4 ligase closes the reporter circle",circle_rows(reporter_circle_template(),"two insert/backbone ligation junctions"),name="RNA-templated reporter circle",note="Both circle oligos and their probe-tail docking sequences are published.")
  w.react("Phi29 rolling-circle amplification",duplex_rows(rolony()),name="reporter RCA concatemer",note="The circle is copied into tandem repeats in the fixed cell.")
  w.react("Hybridize the labelled detection oligo",detection_scene().rows(),name="cytometry-ready PLAYR signal",note="The detection oligo reports the insert identity; no nucleotide sequencing follows.")
  return w
