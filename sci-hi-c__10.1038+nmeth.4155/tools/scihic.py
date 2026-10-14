"""Ramani et al. 2017 sci-Hi-C chemistry."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
import illumina as il
import seqprimers as sp
from adapter_ligation import AdapterLigationGeometry
from chemdraw import Construct, MolecularState, Row, Scene, Segment, Workflow, feature, revcomp
from restriction import DPNII
from single_cell_hic import ContactWorkflow
WORKFLOW = ContactWorkflow("sci-Hi-C", DPNII, "combinatorial indexing", "bridge-adaptor", "PCR", barcoding_rounds=2)
TITLE = "sci-Hi-C — two-level combinatorial indexing of contacts"
NOTES = "01_sci-hi-c.html"
CITATION = 'Defining source: <a href="https://doi.org/10.1038/nmeth.4155">Ramani et al., <i>Nature Methods</i> (2017)</a>.'
SUMMARY = "A first well barcode enters through indexed bridge adaptors before proximity ligation; a second well barcode enters through indexed library adaptors after redistribution."
STEPS = (("Digest in 96 wells", "Distribute fixed nuclei and digest chromatin with DpnII."),
         ("Barcode round 1", "Ligate biotinylated indexed bridge adaptors to restriction ends in each first-round well."),
         ("Pool, ligate, redistribute", "Pool nuclei, proximity-ligate bridged ends, then sort a small number of nuclei into second-round wells."),
         ("Barcode round 2 and PCR", "Add indexed Illumina Y adaptors, enrich biotin-bearing molecules and PCR-amplify; the barcode pair identifies one nucleus."))
READOUT = "The first eight bases of each genomic read carry BC2; the two BC1 copies flank the internal bridge. The barcode pair assigns the contact to a nucleus."
JUNCTION_CAPTION = "The contact junction contains the first-round indexed bridge (BC1); ** marks each ligation boundary."
LIBRARY_CAPTION = "BC1 remains internal. BC2 is an 8-bp inline barcode at both read starts, followed by the source-printed CGT fixed bases; sci-Hi-C does not use an index read for BC2."
JUNCTIONS = (("locus A", "DpnII end A", "restriction end"),
             ("DpnII end A", "round-1 barcode A", "adapter ligation"),
             ("round-1 barcode B", "DpnII end B", "adapter ligation"))
BC1_EXAMPLE = "GAGCTCGA"  # A1 bridge-adaptor barcode, Supplementary Spreadsheet 1
BRIDGE_CORE = BC1_EXAMPLE + "G"  # pairs to MKVII_R_1; final G begins the EcoRI site
BRIDGE_ECORI_OVERHANG = "AATT"

def contact_product():
    return Construct([Segment("locus A", "X"*22, placeholder=True), Segment("DpnII end A", "GATC", "me"),
        Segment("round-1 barcode A", BC1_EXAMPLE, "cbc",
                feature=feature("cell_bc1_a", "cell_barcode", "combinatorial", group="cell_id", part="round 1, end A")), Segment("bridge EcoRI site", "GAATTC", "me"),
        Segment("round-1 barcode B", revcomp(BC1_EXAMPLE), "cbc",
                feature=feature("cell_bc1_b", "cell_barcode", "combinatorial", group="cell_id", part="round 1, end B")), Segment("DpnII end B", "GATC", "me"),
        Segment("locus B", "X"*22, placeholder=True)], name="sci-Hi-C bridged contact")


def bridge_adapter_geometry():
    """MKVII_F_1/MKVII_R_1: DpnII-compatible end, BC1 and EcoRI overhang."""
    sc = Scene()
    sc.strand("MKVII_F_1", [Segment("5′ GATC", "GATC", "me"),
                             Segment("BC1 + EcoRI G", BRIDGE_CORE, "cbc"),
                             Segment("3′ AATT", BRIDGE_ECORI_OVERHANG, "me")],
              mod3="internal biotin-dT in AATT")
    sc.anneal("MKVII_R_1", [Segment("paired reverse", revcomp(BRIDGE_CORE), "cbc")],
              to="MKVII_F_1", pair=("paired reverse", "BC1 + EcoRI G"))
    sc.ligation_end("MKVII_F_1", "5′ GATC", kind="cohesive", polarity="5-prime",
                    owner="adapter", text="5′-GATC DpnII-compatible end")
    sc.overhang("MKVII_F_1", "3′ AATT", polarity="3-prime",
                text="AATT proximity-ligation overhang")
    return AdapterLigationGeometry.known(
        "cohesive", sc,
        statement="Ramani et al. Supplementary Spreadsheet 1 prints MKVII_F_1 and MKVII_R_1, including the 5′ GATC end, eight-base barcode and internal biotin-dT.")


def y_adapter_geometry():
    """SeqShort_BC1 HI/LO: forked TruSeq adapter with a one-base 3′ dT end."""
    hi = il.TRUSEQ_READ1 + BC2_EXAMPLE + "CGT"
    # The printed LO arm joins ``...AAG`` directly to ``ATCG...``: it omits the first G
    # of the standalone 34-nt Index-1 primer constant.
    lo = "CGCAAGGTCAAG" + il.INDEX1_PRIMER[1:]
    paired_hi = hi[-24:-1]
    paired_lo = lo[:23]
    if paired_hi != revcomp(paired_lo):
        raise ValueError("sci-Hi-C Y-adapter HI/LO core does not pair")
    sc = Scene()
    sc.strand("SeqShort_BC1_HI", [Segment("unpaired Read-1 arm", hi[:-24], "r1"),
                                  Segment("paired HI core", paired_hi, "cbc"),
                                  Segment("3′ dT", hi[-1], "me")], mod3="phosphorothioate")
    sc.anneal("SeqShort_BC1_LO paired end", [Segment("paired LO core", paired_lo, "cbc")],
              to="SeqShort_BC1_HI", pair=("paired LO core", "paired HI core"), mod5="phosphate")
    sc.ligation_end("SeqShort_BC1_HI", "3′ dT", kind="cohesive", polarity="3-prime",
                    owner="adapter", text="3′-dT pairs with a dA-tailed insert")
    lo_sc = Scene()
    lo_sc.strand("complete SeqShort_BC1_LO",
                 [Segment("paired LO core", paired_lo, "cbc"),
                  Segment("unpaired Index-1 arm", lo[23:], "r2")], mod5="phosphate")
    lo_sc.note("complete SeqShort_BC1_LO",
               "The forked unpaired arm is shown separately so its bases remain readable.")
    return AdapterLigationGeometry.known(
        "cohesive", sc, lo_sc,
        statement="Ramani et al. Supplementary Spreadsheet 1 prints the complete SeqShort_BC1 HI/LO pair and its modifications; BC1 is representative of 96 Y adapters.")

def adapter_workflow():
    start=Scene.duplex([Segment("DpnII fragment", "X"*20+"GATC", placeholder=True)],label="DpnII-digested chromatin")
    w=Workflow(MolecularState("DpnII restriction ends",tuple(start.rows())))
    bridge=bridge_adapter_geometry()
    w.ligate_adapter("Ligate first-round bridge adapters",bridge,Scene.duplex(list(contact_product()),label="bridged contact").rows(),name="bridge-adapter contact",note="BC1 remains internal after proximity ligation.")
    yend=y_adapter_geometry()
    w.ligate_adapter("Ligate second-round indexed Y adapters",yend,Scene.duplex(list(final_library()),label="sci-Hi-C library").rows(),name="BC2 library",note="BC2 begins each genomic read.")
    return w
BC2_EXAMPLE = "TGACCTTG"  # A1 adapter, Supplementary Data oligo table
SEQ_PRIMERS = (sp.TRUSEQ["R1"], sp.TRUSEQ["R2"])
REQUIRED_ROLES = ("Read 1", "Read 2")
def final_library():
    read_prefix = BC2_EXAMPLE + "CGT"
    lib = Construct([
        Segment("P5", il.P5, "p5"),
        Segment("Read 1 arm", il.TRUSEQ_READ1[4:], "r1"),
        Segment("BC2 A1", BC2_EXAMPLE, "cbc",
                feature=feature("cell_bc2_a", "cell_barcode", "combinatorial", group="cell_id", part="round 2, end A")), Segment("fixed read prefix", "CGT", "me"),
        *list(contact_product()),
        Segment("opposite fixed prefix and BC2", revcomp(read_prefix), "cbc",
                feature=feature("cell_bc2_b_with_prefix", "cell_barcode", "combinatorial", group="cell_id", part="round 2, end B", note="segment also contains the fixed CGT prefix")),
        Segment("Read 2 site", revcomp(il.TRUSEQ_READ2), "r2"),
        Segment("P7 reverse complement", il.P7_RC, "p7"),
    ], name="sci-Hi-C two-barcode library")
    errors = sp.verify(lib, SEQ_PRIMERS, required_roles=REQUIRED_ROLES)
    if errors: raise ValueError("invalid sci-Hi-C library: " + "; ".join(errors))
    return lib
