"""10x Chromium 3' v3 Cell Surface Protein chemistry (CG000185 Rev B)."""
from __future__ import annotations
import illumina as il
import nextera as nx
import seqprimers as sp
from chemdraw import Construct,Scene,Segment,feature,revcomp
from pcr import PCRGeometry, PCRPrimer
from protocol_graph import Disposition, InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind, Terminal
from protocol_model import molecule, role
import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[2]))
from protocol_families.tenx_graph import partition_material

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
# CG000185 Rev B Step 5 run parameters (txt 5310-5330): recommended cycles; Read 2 minimum is 25.
READ_LENGTHS={"Read 1":28,"Index 1 (i7)":8,"Read 2":91}
READ2_MINIMUM=25
def final_library():
    lib=Construct([seg("P5",il.P5,"p5"),seg("Nextera Read 1",nx.READ1_PRIMER,"r1"),seg("cell barcode","B"*16,"cbc",placeholder=True),seg("UMI","U"*12,"umi",placeholder=True),seg("Capture Sequence 1",CS1,"tso"),seg("diversity N9","N"*9,placeholder=True),seg("feature barcode","F"*15,"cbc",placeholder=True),seg("diversity N10","N"*10,placeholder=True),seg("Read 2 leading A","A","r2"),seg("Index 1 / Read 2 arm",il.INDEX1_PRIMER,"r2"),seg("i7 reverse complement","I"*8,"cbc",placeholder=True),seg("P7 reverse complement",il.P7_RC,"p7")],name="10x 3' v3 Cell Surface Protein library")
    problems=sp.verify(lib,SEQ_PRIMERS,required_roles=("Read 1","Index 1 (i7)","Read 2"))
    if problems: raise ValueError("invalid Cell Surface Protein library: "+"; ".join(problems))
    return lib
def read_layout():
    spans=sp.feature_spans(final_library(),SEQ_PRIMERS,READ_LENGTHS)
    rows=lambda role:[(f"{x.read}, {x.cycles}",x.feature.label) for x in spans if x.read==role]
    return [*rows("Read 1"),("Read 2, 1–10","diversity bases"),*rows("Read 2")]


SOURCE="10x Genomics User Guide CG000185 Rev B"
# CG000185 Rev B appendix p.65 (Step 2.2) and p.66 (Step 4.1); expressed via canonical parts.
ADT_FORWARD=nx.READ1_PRIMER[-27:]            # Feature cDNA Primers 2 PN-2000097, Partial Read 1N
ADT_REVERSE=il.TRUSEQ_READ2[:19]             # Feature cDNA Primers 2 PN-2000097, Partial Read 2
SI_FORWARD=il.P5+nx.READ1_PRIMER             # Feature SI Primers 2 PN-2000099
I7_PRIMER=il.P7+"N"*8+il.TRUSEQ_READ2[:21]   # Chromium i7 Sample Index PN-220103 (4-oligo well mix)
def oligo_state(name,parts): return Construct([seg(*x) if isinstance(x,tuple) else x for x in parts],name=name)
def amplified_feature():
    """Step 2.2 product: the 27-nt Partial Read 1N primer extends the 22-nt bead handle."""
    return Construct([seg("Partial Read 1N",ADT_FORWARD,"r1"),*list(feature_cdna())[1:]],name="amplified Cell Surface Protein feature DNA")


def protocol_graph():
    authority=('Primary chemistry source: <a href="https://teichlab.github.io/scg_lib_structs/data/10X-Genomics/CG000185_ChromiumSingleCell3__FeatureBarcode_CellSurfaceProtein_Rev_B.pdf">'
               '10x Genomics User Guide CG000185 Rev B</a> (vendor PDF mirrored by scg_lib_structs); '
               'labeling per linked Demonstrated Protocol CG000149 Rev D.')
    g=ProtocolGraph("10x Chromium 3-prime v3 Cell Surface Protein",authority=authority)
    cells=g.input("cells","cells with surface epitopes",MaterialKind.BIOLOGICAL,
                  state=molecule("cells",(role("surface epitopes"),)),primary=True)
    antibodies=g.input("antibodies","antibody-conjugated Feature Barcode oligos (e.g. TotalSeq-B), then washes",
                       MaterialKind.SELECTION)
    stained=g.product("stained","cells retaining epitope-bound Feature Barcode oligos",
                      MaterialKind.MOLECULAR,state=molecule("stained cells",(role("surface epitopes"),role("antibody feature tags"))))
    unbound=g.product("unbound","unbound antibody tags",MaterialKind.WASTE)
    g.react("stain","Stain cells (Fc block, 30 min 4 °C) and wash away unbound conjugates",
            StepKind.BINDING,((cells,InputRole.TEMPLATE),(antibodies,InputRole.SELECTION)),
            (stained,unbound)); g.discard(unbound,"removed in the CG000149 cell washes")
    reviewed=["cells","antibodies","stain"]
    occupied=partition_material(g,stained,reviewed,label="stained cells")
    bead=g.input("bead_primer","Capture Sequence 1 cell-barcode/UMI gel-bead primer",
                 MaterialKind.OLIGO,state=Construct(bead_primer(),name="Capture Sequence 1 gel-bead primer"))
    polymerase=g.input("extension_polymerase","RT Enzyme C (GEM-RT, 53 °C 45 min, 85 °C 5 min)",
                       MaterialKind.CATALYST)
    copied=g.product("copied","cell-barcoded antibody feature DNA",MaterialKind.MOLECULAR,
                     state=feature_cdna())
    g.react("copy","Extend the Capture Sequence 1 bead primer across the antibody tag",
            StepKind.EXTENSION,((occupied,InputRole.TEMPLATE),(bead,InputRole.STRUCTURAL),
            (polymerase,InputRole.CATALYST)),(copied,))
    recovery=g.input("gem_recovery","Recovery Agent GEM break and Dynabeads MyOne SILANE cleanup",
                     MaterialKind.SELECTION)
    recovered=g.product("recovered","recovered GEM-RT DNA",MaterialKind.MOLECULAR,state=feature_cdna())
    recovery_waste=g.product("recovery_waste","Recovery Agent/partitioning oil phase and bead-cleanup supernatant",
                             MaterialKind.WASTE)
    g.react("recover","Break GEMs and purify GEM-RT products (step 2.1)",StepKind.SELECTION,
            ((copied,InputRole.TEMPLATE),(recovery,InputRole.SELECTION)),(recovered,recovery_waste),
            changes_molecule=False)
    g.discard(recovery_waste,"removed during post GEM-RT cleanup")
    reviewed.extend(("bead_primer","extension_polymerase","copy","gem_recovery","recover"))
    fwd=g.input("feature_pcr_forward","Feature cDNA Primers 2 Partial Read 1N primer (PN-2000097)",
                MaterialKind.OLIGO,state=oligo_state("Partial Read 1N primer",[("Partial Read 1N",ADT_FORWARD,"r1")]))
    rev=g.input("feature_pcr_reverse","Feature cDNA Primers 2 Partial Read 2 primer (PN-2000097)",
                MaterialKind.OLIGO,state=oligo_state("Partial Read 2 primer",[("Partial Read 2",ADT_REVERSE,"r2")]))
    amp=g.input("feature_pcr_polymerase","Amp Mix (PN-2000047)",MaterialKind.CATALYST)
    amplified=g.product("feature_pcr_product","co-amplified feature DNA and cDNA",MaterialKind.MOLECULAR,
                        state=amplified_feature())
    g.react("feature_pcr","Co-amplify feature DNA with Feature cDNA Primers 2 (11–13 cycles)",StepKind.PCR,
            ((recovered,InputRole.TEMPLATE),(fwd,InputRole.STRUCTURAL),(rev,InputRole.STRUCTURAL),
             (amp,InputRole.CATALYST)),(amplified,),
            geometry=PCRGeometry.exact(feature_cdna(),PCRPrimer("Partial Read 1N",ADT_FORWARD),
                                       PCRPrimer("Partial Read 2",ADT_REVERSE),
                                       statement="CG000185 Rev B appendix Step 2.2 prints both feature primers; the poly(dT) cDNA pair in the same mix is not drawn."))
    reviewed.extend(("feature_pcr_forward","feature_pcr_reverse","feature_pcr_polymerase","feature_pcr"))
    beads=g.input("size_selection","0.6X SPRIselect",MaterialKind.SELECTION)
    selected=g.product("selected","transferred 80 µl supernatant with feature DNA",MaterialKind.MOLECULAR,
                       state=amplified_feature())
    gex=g.product("gex_pellet","bead-bound gene-expression cDNA (pellet)",MaterialKind.MOLECULAR)
    g.react("size_select","Separate short feature DNA (supernatant) from gene-expression cDNA (pellet)",
            StepKind.SELECTION,((amplified,InputRole.TEMPLATE),(beads,InputRole.SELECTION)),
            (selected,gex),changes_molecule=False)
    g.mark_terminal(Terminal(gex.id,Disposition.STORED,
        "pellet cleanup (2.3A) feeds the separate 3-prime Gene Expression library, step 3; not drawn"))
    clean=g.input("supernatant_cleanup","2.0X SPRIselect on the transferred supernatant",MaterialKind.SELECTION)
    cleaned=g.product("cleaned","purified feature DNA (5 µl of 40 µl used)",MaterialKind.MOLECULAR,
                      state=amplified_feature())
    clean_waste=g.product("cleanup_waste","2.0X SPRI supernatant and ethanol washes",MaterialKind.WASTE)
    g.react("supernatant_cleanup_step","Bind feature DNA at 2.0X SPRIselect and elute (2.3B)",
            StepKind.SELECTION,((selected,InputRole.TEMPLATE),(clean,InputRole.SELECTION)),
            (cleaned,clean_waste),changes_molecule=False)
    g.discard(clean_waste,"removed by the 2.3B bead cleanup")
    reviewed.extend(("size_selection","size_select","supernatant_cleanup","supernatant_cleanup_step"))
    lib=final_library()
    sif=g.input("index_pcr_forward","Feature SI Primers 2 P5/Read 1N primer (PN-2000099)",MaterialKind.OLIGO,
                state=oligo_state("Feature SI Primers 2",[("P5",il.P5,"p5"),("Nextera Read 1",nx.READ1_PRIMER,"r1")]))
    sir=g.input("index_pcr_reverse","Chromium i7 Sample Index (PN-220103)",MaterialKind.OLIGO,
                state=oligo_state("Chromium i7 Sample Index",[("P7",il.P7,"p7"),seg("i7","I"*8,"cbc",placeholder=True,feature=I7_FEATURE),("Partial Read 2",il.TRUSEQ_READ2[:21],"r2")]))
    sipol=g.input("index_pcr_polymerase","Amp Mix (PN-2000047)",MaterialKind.CATALYST)
    indexed=g.product("index_pcr_product","indexed Cell Surface Protein library",MaterialKind.MOLECULAR,state=lib)
    g.react("index_pcr","Add P5, P7 and sample index by Sample Index PCR (10 cycles)",StepKind.PCR,
            ((cleaned,InputRole.TEMPLATE),(sif,InputRole.STRUCTURAL),(sir,InputRole.STRUCTURAL),
             (sipol,InputRole.CATALYST)),(indexed,),
            geometry=PCRGeometry.exact(amplified_feature(),PCRPrimer("Feature SI Primers 2",SI_FORWARD),
                                       PCRPrimer("Chromium i7 Sample Index",I7_PRIMER),
                                       statement="CG000185 Rev B appendix Step 4.1 prints both primers; each i7 well is a mix of four 8-nt index oligos."))
    final_beads=g.input("final_size_selection","1.2X SPRIselect",MaterialKind.SELECTION)
    library=g.product("library","10x 3' v3 Cell Surface Protein library",MaterialKind.LIBRARY,state=lib)
    final_waste=g.product("final_waste","1.2X SPRI supernatant (primers, short products)",MaterialKind.WASTE)
    g.react("final_cleanup","Post Sample Index PCR 1.2X SPRIselect (4.2)",StepKind.SELECTION,
            ((indexed,InputRole.TEMPLATE),(final_beads,InputRole.SELECTION)),(library,final_waste),
            changes_molecule=False)
    g.discard(final_waste,"removed by the 4.2 bead cleanup")
    reviewed.extend(("index_pcr_forward","index_pcr_reverse","index_pcr_polymerase","index_pcr",
                     "final_size_selection","final_cleanup"))
    g.sequence(library,SequencingPlan(library=lib,primers=tuple(SEQ_PRIMERS),
        acquired_roles=tuple(p.role for p in SEQ_PRIMERS),
        read_lengths=dict(READ_LENGTHS)))
    reviewed.append(library.id)
    g.declare_reviewed_scope(reviewed,authority)
    return g.require_publishable()
