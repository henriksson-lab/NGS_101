"""PLAYR proximity-ligation/RCA assay model."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"lib"))
from chemdraw import Construct, MolecularState, Row, Scene, Segment, Workflow, revcomp, strand_row
from circular import circularize_ssdna, rolling_circle
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
INSERT_SPACER=INSERT[len(PROBE1_INSERT_SITE):-len(PROBE2_INSERT_SITE)]
DETECTION=INSERT
# Nature Methods Supplementary Table 1 (MOESM647 "PLAYR sequences") and Nature
# Protocols Table 1 write every detection oligo as "Z-" + the insert sequence; the
# Nature Protocols footnote defines Z as a 5'-C6-thiol (mass cytometry, metal X8
# polymer conjugate) or a fluorophore (flow cytometry).
DETECTION_5PRIME_MOD="Z (5′-C6-thiol for metal conjugation, or fluorophore)"
# MOESM647 also lists a universal "Backbone Detection" oligo.  It is not used in the
# drawn system-1 workflow.
BACKBONE_DETECTION="CAGTGAATGCGAGTCCGTCT"
if BACKBONE_DETECTION not in BACKBONE:
  raise AssertionError("backbone detection oligo must share the circle's sense to bind the RCA product")

def seg(name,seq,tag=None,**kw): return Segment(name,seq,tag,**kw)

TOPOLOGY=ProximityCircleTopology(
  upstream=ProximityProbe("Probe 2","N"*24,LINKER,"outer",PROBE2_BACKBONE_SITE,
                          "inner",PROBE2_INSERT_SITE),
  downstream=ProximityProbe("Probe 1","N"*24,LINKER,"inner",PROBE1_INSERT_SITE,
                            "outer",PROBE1_BACKBONE_SITE),
  inner=INSERT,outer=BACKBONE,inner_spacer=INSERT_SPACER,
  outer_body=BACKBONE[len(PROBE2_BACKBONE_SITE):-len(PROBE1_BACKBONE_SITE)],
  name="PLAYR adjacent probes and ligation-ready reporter circle")

def reporter_circle_template() -> Construct:
  return Construct([seg("target-specific insert",INSERT,"cbc"),seg("universal backbone",BACKBONE,"r3")],"PLAYR system-1 reporter circle")

CIRCLE=circularize_ssdna(reporter_circle_template(),five_prime_phosphate=True)

# PLAYR does not add a separate RCA primer: Phi29 extends the 3'-OH of a probe-pair
# member.  In the published orientation, Probe 1's complete 3' docking region spans the
# circle closure (insert then backbone in the probe, hence backbone then insert on the
# circular template).  rolling_circle() makes the unique landing site and product phase
# construction invariants rather than labels added by the page.
RCA_PRIMER=Construct([
  seg("Probe 1 insert-binding region",PROBE1_INSERT_SITE,"cbc"),
  seg("Probe 1 backbone-binding 3'-OH",PROBE1_BACKBONE_SITE,"r3"),
],"Probe 1 RCA-priming 3' tail")
RCA=rolling_circle(CIRCLE,RCA_PRIMER.top(),copies=3)

def hybridized_reporter_rows() -> list[Row]:
  """The two phosphorylated oligos docked on the probe tails, before T4 ligation."""
  return [
    strand_row(Construct([seg("target-specific insert",INSERT,"cbc")]),prefix="Insert 1  5′-p-",suffix="-3′"),
    strand_row(Construct([seg("universal backbone",BACKBONE,"r3")]),prefix="Backbone  5′-p-",suffix="-3′"),
    Row(chunks=[("Both oligos are docked on the probe tails as in the topology panel above; two nicks remain (one insert/backbone junction on each probe).",None,False)]),
  ]

def rca_priming_scene() -> Scene:
  """Draw the circular template opened away from Probe 1's exact landing site."""
  template=Construct([
    seg("backbone remainder",BACKBONE[:-len(PROBE1_BACKBONE_SITE)],"r3"),
    seg("Probe 1 backbone template site",BACKBONE[-len(PROBE1_BACKBONE_SITE):],"r3"),
    seg("Probe 1 insert template site",INSERT[:len(PROBE1_INSERT_SITE)],"cbc"),
    seg("insert remainder",INSERT[len(PROBE1_INSERT_SITE):],"cbc"),
  ],"PLAYR reporter circle opened for primer placement")
  sc=Scene()
  sc.strand("circle",list(template),label="primer-binding site expanded from the circle above")
  sc.anneal("probe1",list(RCA_PRIMER),to="circle",
            pair=("Probe 1 backbone-binding 3'-OH","Probe 1 backbone template site"),
            label="Probe 1 3' tail")
  sc.arrow("probe1","Phi29 extends the probe's 3'-OH")
  sc.mark("circle","Probe 1 backbone template site","primer-binding site",
          through="Probe 1 insert template site")
  sc.footer("The covalently closed template is opened at a non-priming position only for sequence display.","circle")
  return sc

def rca_priming_rows() -> list[Row]:
  """Keep the closed topology visible above the sequence-resolved binding detail."""
  return [*CIRCLE.rows("closed reporter circle"),
          Row(chunks=[("Primer-binding junction expanded below",None,False)]),
          *rca_priming_scene().rows()]

def rolony() -> Construct:
  parts=[]
  for i in range(1,4):
    p1_insert=("Probe 1 insert-binding primer" if i == 1
               else f"Probe 1 insert-site copy {i}")
    p1_backbone=("Probe 1 backbone-binding 3'-OH" if i == 1
                 else f"Probe 1 backbone-site copy {i}")
    parts.extend([
      seg(p1_insert,PROBE1_INSERT_SITE,"cbc"),
      seg(p1_backbone,PROBE1_BACKBONE_SITE,"r3"),
      seg(f"backbone body complement {i}",revcomp(BACKBONE[len(PROBE2_BACKBONE_SITE):-len(PROBE1_BACKBONE_SITE)]),"r3"),
      seg(f"Probe 2 backbone-binding copy {i}",PROBE2_BACKBONE_SITE,"r3"),
      seg(f"Probe 2 insert-binding copy {i}",PROBE2_INSERT_SITE,"cbc"),
      seg(f"insert spacer complement {i}",revcomp(INSERT_SPACER),"cbc"),
    ])
  product=Construct(parts,"PLAYR RCA concatemer (three repeats drawn)")
  if product.top() != RCA.sequence:
    raise AssertionError("drawn PLAYR concatemer drifted from rolling-circle construction")
  return product

def rolony_scene() -> Scene:
  sc=Scene()
  sc.strand("rolony",list(rolony()),label="single-stranded RCA concatemer")
  sc.note("rolony","three circle traversals shown; the physical product can contain hundreds")
  return sc

def detection_scene() -> Scene:
  product=rolony(); target="Probe 2 insert-binding copy 2"
  sc=Scene(); sc.strand("rolony",list(product),label="phi29 RCA concatemer")
  # Detection oligo is identical to the published insert; its complement occurs in the RCA product.
  sc.anneal("detection oligo",[seg("system-1 detection sequence",DETECTION,"umi")],to="rolony",pair=("system-1 detection sequence",target),label="fluorescent or metal-tagged detection oligo")
  sc.arrow("detection oligo","measure reporter abundance by flow or mass cytometry")
  return sc

def workflow() -> Workflow:
  w=Workflow(MolecularState("adjacent PLAYR probes on target RNA",
    (Row(chunks=[("Topology above: Probe 2 and Probe 1 expose the complete split docking platform",None,False)]),)),track_graph=False)
  w.react("Hybridize the phosphorylated insert and universal backbone to the probe tails (100 nM each, 1× SSC, 37 °C, 30 min; wash twice)",
    hybridized_reporter_rows(),name="ligation-ready insert/backbone on the probe pair",
    note="Frei 2016 Online Methods 'PLAYR protocol'; Nature Protocols reagent setup and Steps 23–24. Insert/backbone hybridization is a separate incubation before ligation.")
  w.react("T4 DNA ligase seals both nicks (room temperature/21 °C, 30 min; wash twice)",CIRCLE.rows("closed reporter circle"),name="RNA-templated reporter circle",note="Both circle oligos and their probe-tail docking sequences are published.")
  w.react("Position a probe-pair member's 3'-OH as the RCA primer",rca_priming_rows(),name="Probe-primed reporter circle (Probe 1 shown)",note="The return path above is the covalently closed circle; the lower expansion shows its exact primer-binding junction. No separate RCA primer is added.")
  w.react("Phi29 rolling-circle amplification (30 °C, ≥2 h flow cytometry, ≥4 h mass cytometry)",rolony_scene().rows(),name="single-stranded reporter RCA concatemer",note="Phi29 extends the bound probe around the circle into complementary tandem repeats in the fixed cell.")
  w.react("Hybridize the 5′-Z-labelled detection oligo (37 °C, 30 min)",detection_scene().rows(),name="cytometry-ready PLAYR signal",note="Z is a 5′-C6-thiol (metal-polymer conjugate, mass cytometry) or a fluorophore (flow cytometry). The detection oligo reports the insert identity; no nucleotide sequencing follows.")
  return w
