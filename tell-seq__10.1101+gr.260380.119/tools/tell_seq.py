"""Molecular architecture of the original TELL-seq linked-read protocol."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
from chemdraw import (Construct, MolecularState, Scene, Segment, Workflow,
                      annotation_rows, feature)
import illumina as il
import seqprimers as sp
from pcr import PCRGeometry
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind)

TITLE = "TELL-seq — single-tube bead-linked reads"
NOTES = "01_tell-seq.html"
SOURCE = 'Defining source: <a href="https://doi.org/10.1101/gr.260380.119">Chen et al., <i>Genome Research</i> (2020)</a>.'
SUMMARY = "Stable MuA strand-transfer complexes and hybrid capture keep fragments from one long DNA molecule associated with one clonally barcoded bead in an unpartitioned tube."

def seg(name, top, tag=None, **kw): return Segment(name, top, tag, **kw)

LINKED = feature("tell.linked-molecule", "linked_read_barcode", "random",
                 note="18-base TELL bead molecular barcode read as Index 1")
SAMPLE = feature("tell.sample", "sample_index", "combinatorial",
                 note="8-base library multiplex barcode read as Index 2")

LONG_DNA = Construct([seg("high-molecular-weight genomic DNA", "X"*58, placeholder=True,
                          length_bp=50000)], name="long genomic molecule")
CAPTURED = Construct([
    seg("bead-tethered barcode adaptor", "[bead adaptor]", "p7", placeholder=True),
    seg("18-nt linked-read barcode", "N"*18, "cbc", placeholder=True, feature=LINKED),
    seg("captured MuA strand-transfer site", "[MuA end]", "r1", placeholder=True),
    seg("captured long genomic DNA", "X"*58, placeholder=True, length_bp=50000),
], name="bead-captured strand-transfer complex")
SECOND_TAGGED = Construct([
    seg("bead-tethered barcode adaptor", "[bead adaptor]", "p7", placeholder=True),
    seg("18-nt linked-read barcode", "N"*18, "cbc", placeholder=True, feature=LINKED),
    seg("captured MuA strand-transfer site", "[MuA end]", "r1", placeholder=True),
    seg("genomic subfragment", "X"*34, placeholder=True),
    seg("second MuA strand-transfer site", "[MuA end 2]", "r2", placeholder=True),
], name="two-ended bead-captured linked fragment")
RELEASED = Construct([
    seg("18-nt linked-read barcode", "N"*18, "cbc", placeholder=True, feature=LINKED),
    seg("captured MuA strand-transfer site", "[MuA end]", "r1", placeholder=True),
    seg("genomic subfragment", "X"*34, placeholder=True),
    seg("second MuA strand-transfer site", "[MuA end 2]", "r2", placeholder=True),
], name="released barcoded subfragment")
FINAL = Construct([
    seg("P5", il.P5, "p5"),
    seg("Read 1 primer site", "[custom R1 site]", "r1", placeholder=True),
    seg("genomic insert", "X"*34, placeholder=True),
    seg("Read 2 primer site", "[custom R2 site]", "r2", placeholder=True),
    seg("18-nt linked-read barcode", "N"*18, "cbc", placeholder=True, feature=LINKED),
    seg("Index 1 primer site", "[custom I1 site]", "r1", placeholder=True),
    seg("8-nt sample index", "N"*8, "cbc", placeholder=True, feature=SAMPLE),
    seg("Index 2 primer site", "[custom I2 site]", "r2", placeholder=True),
    seg("P7 reverse complement", il.P7_RC, "p7"),
], name="TELL-seq Illumina library")
FINAL_LIBRARY = FINAL
SEQ_PRIMERS = ()
FINAL_CAPTION = "Index 1 is the 18-base linked-molecule barcode; Index 2 is the 8-base sample index. Adapter ordering is published; proprietary bases remain role tokens."
SEQUENCING_ENDING = "TELL-seq uses vendor custom Read 1, Read 2, Index 1 and Index 2 primers. The guide specifies 2 × 146 reads, an 18-cycle linked-molecule-barcode read and an 8-cycle sample-index read."
SEQUENCING_UNAVAILABLE = "The custom sequencing-primer sequences are proprietary, so exact annealing cannot be computed without fabricating bases."
READ_SPECS = (
    sp.read_spec("Read 1",146,"genomic insert","include"),
    sp.read_spec("Index 1 (i7)",18,"linked-molecule barcode","extract"),
    sp.read_spec("Index 2 (i5)",8,"sample index","demultiplex"),
    sp.read_spec("Read 2",146,"genomic insert","include"),
)

def rows(con, label=""):
    unknown = tuple(s.name + "'" for s in con if s.is_role_token())
    return [*Scene.duplex(list(con), label=label, unpaired=unknown).rows(),
            *annotation_rows(con)]

INITIAL_NAME = "High-molecular-weight genomic DNA"
INITIAL_ROWS = tuple(rows(LONG_DNA,"HMW DNA"))

def workflow():
    wf=Workflow(MolecularState(INITIAL_NAME, INITIAL_ROWS))
    wf.react("Simultaneous MuA strand transfer and bead hybrid capture", rows(CAPTURED,"on bead"),
             note="Stable strand-transfer complexes keep subfragments of one long molecule near one clonal bead barcode.")
    wf.react("Second transpososome tagging between capture complexes", rows(SECOND_TAGGED,"on bead"),
             note="A second priming site is installed before complexes are disrupted.")
    wf.react("Break strand-transfer complexes and wash beads", rows(RELEASED,"released fragment"),
             note="Protein-complex disruption fragments the target while retaining bead-linked barcode identity.")
    wf.react("PCR off beads; add P5, P7 and sample index", rows(FINAL,"library"),
             note="Eight cycles were used for human DNA and 13–14 for microbial DNA in the defining paper.")
    return wf


def protocol_graph():
    graph=ProtocolGraph(TITLE,authority=SOURCE)
    dna=graph.input("hmw_dna","high-molecular-weight genomic DNA",
                    MaterialKind.MOLECULAR,state=LONG_DNA,primary=True)
    bead_adapter=graph.input("bead_adapter","clonally barcoded bead capture adapter",
                             MaterialKind.ADAPTER,inferred=True,
                             source=SOURCE+" (proprietary adapter bases)")
    mua=graph.input("mua","MuA transposase",MaterialKind.CATALYST)
    captured=graph.product("captured","bead-captured strand-transfer complexes",
                           MaterialKind.MOLECULAR,state=CAPTURED,inferred=True,
                           source=SOURCE+" (proprietary adapter bases)")
    first_waste=graph.product("first_waste","uncaptured first-transfer products",
                              MaterialKind.WASTE)
    graph.react("capture_transfer","Combine MuA transfer with bead hybrid capture",
                StepKind.TAGMENTATION,
                ((dna,InputRole.TEMPLATE),(bead_adapter,InputRole.STRUCTURAL),
                 (mua,InputRole.CATALYST)),(captured,first_waste))
    graph.discard(first_waste,"fragments not retained near a clonal bead barcode are washed away")
    second_adapter=graph.input("second_adapter","second MuA strand-transfer adapter",
                               MaterialKind.ADAPTER,inferred=True,
                               source=SOURCE+" (proprietary adapter bases)")
    double_tagged=graph.product("double_tagged","two-ended bead-captured linked fragments",
                                MaterialKind.MOLECULAR,state=SECOND_TAGGED,inferred=True,
                                source=SOURCE+" (proprietary adapter bases)")
    second_waste=graph.product("second_waste","non-amplifiable second-transfer products",
                               MaterialKind.WASTE)
    graph.react("second_transfer","Install the second transpososome priming site",
                StepKind.TAGMENTATION,
                ((captured,InputRole.TEMPLATE),(second_adapter,InputRole.STRUCTURAL),
                 (mua,InputRole.CATALYST)),(double_tagged,second_waste))
    graph.discard(second_waste,"products without the required unlike transfer ends are not amplified")
    disruption=graph.input("disruption_condition","strand-transfer-complex disruption condition",
                           MaterialKind.CONDITION)
    released=graph.product("released","released barcoded genomic subfragments",
                           MaterialKind.MOLECULAR,state=RELEASED,inferred=True,
                           source=SOURCE+" (proprietary adapter bases)")
    graph.react("release","Break strand-transfer complexes and wash beads",
                StepKind.DEPROTEINIZATION,
                ((double_tagged,InputRole.TEMPLATE),(disruption,InputRole.CONDITION)),
                (released,))
    forward=graph.input("pcr_forward","TELL-seq library PCR forward primer",
                        MaterialKind.OLIGO,inferred=True,
                        source=SOURCE+" (primer sequence proprietary)")
    reverse=graph.input("pcr_reverse","TELL-seq indexed PCR reverse primer",
                        MaterialKind.OLIGO,inferred=True,
                        source=SOURCE+" (primer sequence proprietary)")
    polymerase=graph.input("pcr_polymerase","library PCR polymerase",MaterialKind.CATALYST)
    library=graph.product("library",FINAL.name,MaterialKind.LIBRARY,state=FINAL,
                          inferred=True,source=SOURCE+" (proprietary primer sites)")
    graph.react("library_pcr","Add P5, P7 and the sample index",StepKind.PCR,
                ((released,InputRole.TEMPLATE),(forward,InputRole.STRUCTURAL),
                 (reverse,InputRole.STRUCTURAL),(polymerase,InputRole.CATALYST)),
                (library,),geometry=PCRGeometry.unavailable(
                    FINAL,reason="The architecture is published but custom primer bases are proprietary."))
    graph.sequence(library,SequencingPlan(limitation=SEQUENCING_UNAVAILABLE,
                                          read_lengths={"Read 1":146,"Index 1 (i7)":18,
                                                        "Index 2 (i5)":8,"Read 2":146},
                                          read_uses={"Read 1":"include","Index 1 (i7)":"extract",
                                                     "Index 2 (i5)":"demultiplex","Read 2":"include"}))
    graph.declare_reviewed_scope((
        "hmw_dna","bead_adapter","mua","capture_transfer","second_adapter",
        "second_transfer","disruption_condition","release","pcr_forward",
        "pcr_reverse","pcr_polymerase","library_pcr","library"),SOURCE)
    return graph.require_publishable()
