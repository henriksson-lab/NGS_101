"""10x Chromium 3' v3 Cell Surface Protein chemistry (CG000185 Rev B)."""
from __future__ import annotations
import illumina as il
import nextera as nx
import seqprimers as sp
from chemdraw import Construct,Scene,Segment,feature,revcomp
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind
from protocol_model import molecule, role
import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[2]))
from protocol_families.tenx_graph import chromium_pcr, partition_material

CS1="TTGCTAGGACCGGCCTTAAAGC"; CS1_RC=revcomp(CS1)
PARTIAL_R1N=nx.READ1_PRIMER[-22:]
CELL_BARCODE=feature("cell_barcode","cell_barcode","whitelist",whitelist="10x-chromium-3prime-v3")
UMI_FEATURE=feature("umi","umi","random")
ANTIBODY_FEATURE=feature("antibody_feature_barcode","feature_barcode","unknown")
I7_FEATURE=feature("sample_index_i7","sample_index","unknown")
FEATURES={"cell barcode":CELL_BARCODE,"UMI":UMI_FEATURE,"feature barcode":ANTIBODY_FEATURE,
          "i7 reverse complement":I7_FEATURE}
def seg(name,top,tag=None,**kw):
    kw.setdefault("feature",FEATURES.get(name)); return Segment(name=name,top=top,tag=tag,**kw)
def bead_primer(): return [seg("Partial Read 1N",PARTIAL_R1N,"r1"),seg("cell barcode","B"*16,"cbc",placeholder=True),seg("UMI","U"*12,"umi",placeholder=True),seg("Capture Sequence 1",CS1,"tso")]
def antibody_oligo(): return [seg("TruSeq Read 2",il.TRUSEQ_READ2,"r2"),seg("diversity N10","N"*10,placeholder=True),seg("feature barcode","F"*15,"cbc",placeholder=True),seg("diversity N9","N"*9,placeholder=True),seg("Capture Sequence 1 reverse complement",CS1_RC,"tso")]
def capture_scene():
    sc=Scene(); sc.strand("antibody oligo",antibody_oligo(),label="antibody-conjugated DNA",mod5="antibody")
    sc.anneal("bead primer",bead_primer(),to="antibody oligo",pair=("Capture Sequence 1","Capture Sequence 1 reverse complement"),label="gel-bead primer",mod5="bead")
    sc.arrow("bead primer","extension"); return sc
def feature_cdna(): return Construct([*bead_primer(),seg("diversity N9","N"*9,placeholder=True),seg("feature barcode","F"*15,"cbc",placeholder=True),seg("diversity N10","N"*10,placeholder=True),seg("Read 2 reverse complement",revcomp(il.TRUSEQ_READ2),"r2")],name="Cell Surface Protein feature cDNA")
SEQ_PRIMERS=(sp.NEXTERA["R1"],sp.TRUSEQ["I1"],sp.TRUSEQ["R2"])
def final_library():
    lib=Construct([seg("P5",il.P5,"p5"),seg("Nextera Read 1",nx.READ1_PRIMER,"r1"),seg("cell barcode","B"*16,"cbc",placeholder=True),seg("UMI","U"*12,"umi",placeholder=True),seg("Capture Sequence 1",CS1,"tso"),seg("diversity N9","N"*9,placeholder=True),seg("feature barcode","F"*15,"cbc",placeholder=True),seg("diversity N10","N"*10,placeholder=True),seg("Read 2 leading A","A","r2"),seg("Index 1 / Read 2 arm",il.INDEX1_PRIMER,"r2"),seg("i7 reverse complement","I"*8,"cbc",placeholder=True),seg("P7 reverse complement",il.P7_RC,"p7")],name="10x 3' v3 Cell Surface Protein library")
    problems=sp.verify(lib,SEQ_PRIMERS,required_roles=("Read 1","Index 1 (i7)","Read 2"))
    if problems: raise ValueError("invalid Cell Surface Protein library: "+"; ".join(problems))
    return lib
def read_layout():
    spans=sp.feature_spans(final_library(),SEQ_PRIMERS,{"Read 1":28,"Read 2":25})
    rows=lambda role:[(f"{x.read}, {x.cycles}",x.feature.label) for x in spans if x.read==role]
    return [*rows("Read 1"),("Read 2, 1–10","diversity bases"),*rows("Read 2")]


def protocol_graph():
    authority=('Primary chemistry source: <a href="https://teichlab.github.io/scg_lib_structs/data/10X-Genomics/CG000185_ChromiumSingleCell3__FeatureBarcode_CellSurfaceProtein_Rev_B.pdf">'
               '10x Genomics User Guide CG000185 Rev B</a>.')
    g=ProtocolGraph("10x Chromium 3-prime v3 Cell Surface Protein",authority=authority)
    cells=g.input("cells","cells with surface epitopes",MaterialKind.BIOLOGICAL,
                  state=molecule("cells",(role("surface epitopes"),)),primary=True)
    antibodies=g.input("antibodies","antibody-conjugated Feature Barcode oligos",
                       MaterialKind.SELECTION)
    stained=g.product("stained","cells retaining epitope-bound Feature Barcode oligos",
                      MaterialKind.MOLECULAR,state=molecule("stained cells",(role("surface epitopes"),role("antibody feature tags"))))
    unbound=g.product("unbound","unbound antibody tags",MaterialKind.WASTE)
    g.react("stain","Stain cells with antibody-conjugated Feature Barcode oligos",
            StepKind.BINDING,((cells,InputRole.TEMPLATE),(antibodies,InputRole.SELECTION)),
            (stained,unbound)); g.discard(unbound,"does not report a cell-bound epitope")
    reviewed=["cells","antibodies","stain"]
    occupied=partition_material(g,stained,reviewed,label="stained cells")
    bead=g.input("bead_primer","Capture Sequence 1 cell-barcode/UMI bead primer",
                 MaterialKind.OLIGO)
    polymerase=g.input("extension_polymerase","Feature Barcode extension polymerase",
                       MaterialKind.CATALYST)
    copied=g.product("copied","cell-barcoded antibody feature cDNA",MaterialKind.MOLECULAR,
                     state=feature_cdna())
    g.react("copy","Copy antibody tags from the Capture Sequence 1 bead primer",
            StepKind.EXTENSION,((occupied,InputRole.TEMPLATE),(bead,InputRole.STRUCTURAL),
            (polymerase,InputRole.CATALYST)),(copied,))
    reviewed.extend(("bead_primer","extension_polymerase","copy"))
    amplified=chromium_pcr(g,"feature_pcr",copied,feature_cdna(),reviewed,
                           label="Co-amplify short feature cDNA",inferred_primers=True)
    beads=g.input("size_selection","feature-cDNA size-selection beads",MaterialKind.SELECTION)
    selected=g.product("selected","size-selected feature cDNA",MaterialKind.MOLECULAR,
                       state=feature_cdna())
    rejected=g.product("rejected","non-feature cDNA size fraction",MaterialKind.WASTE)
    g.react("size_select","Separate short feature cDNA from gene-expression cDNA",
            StepKind.SELECTION,((amplified,InputRole.TEMPLATE),(beads,InputRole.SELECTION)),
            (selected,rejected),changes_molecule=False)
    g.discard(rejected,"belongs to a different library branch")
    reviewed.extend(("size_selection","size_select"))
    library=chromium_pcr(g,"index_pcr",selected,final_library(),reviewed,
                         label="Add P5, P7 and sample index",library=True,
                         inferred_primers=True)
    g.sequence(library,SequencingPlan(library=final_library(),primers=tuple(SEQ_PRIMERS),
        acquired_roles=tuple(p.role for p in SEQ_PRIMERS),
        read_lengths={"Read 1":28,"Index 1 (i7)":8,"Read 2":25}))
    reviewed.append(library.id)
    g.declare_reviewed_scope(reviewed,authority)
    return g.require_publishable()
