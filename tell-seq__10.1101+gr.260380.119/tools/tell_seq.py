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
from protocol_graph import (Disposition, InputRole, MaterialKind, ProtocolGraph,
                            SequencingPlan, StepKind, Terminal)

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
# Figure 1B (published, rendered p. 899): P5 | 8-base Index 2 | Read 1 region | insert
# (100-400 bp) | Read 2 region | 18-base Index 1 | P7.  Figure 1A shows the bead barcode
# oligo joined to the first (blue) transposon end, which is therefore the Read 2 region;
# the second transpososome (red) supplies the Read 1 region.  Bead constructs are drawn
# from the bead-tethered 5' end, i.e. the strand opposite to the final top strand.
CAPTURED = Construct([
    seg("bead-tethered barcode adaptor", "[bead adaptor]", "p7", placeholder=True),
    seg("18-nt linked-read barcode", "N"*18, "cbc", placeholder=True, feature=LINKED),
    seg("first MuA transposon end (Read 2 region)", "[MuA end, Read 2 region]", "r2", placeholder=True),
    seg("captured long genomic DNA", "X"*58, placeholder=True, length_bp=50000),
], name="bead-captured strand-transfer complex")
SECOND_TAGGED = Construct([
    seg("bead-tethered barcode adaptor", "[bead adaptor]", "p7", placeholder=True),
    seg("18-nt linked-read barcode", "N"*18, "cbc", placeholder=True, feature=LINKED),
    seg("first MuA transposon end (Read 2 region)", "[MuA end, Read 2 region]", "r2", placeholder=True),
    seg("genomic subfragment", "X"*34, placeholder=True),
    seg("second MuA transposon end (Read 1 region)", "[MuA end, Read 1 region]", "r1", placeholder=True),
], name="two-ended bead-captured linked fragment")
RELEASED = Construct(list(SECOND_TAGGED), name="bead-retained barcoded subfragment after STC disruption")
FINAL = Construct([
    seg("P5", il.P5, "p5"),
    seg("8-nt sample index", "N"*8, "cbc", placeholder=True, feature=SAMPLE),
    seg("Read 1 region", "[Read 1 region]", "r1", placeholder=True),
    seg("genomic insert (100-400 bp)", "X"*34, placeholder=True),
    seg("Read 2 region", "[Read 2 region]", "r2", placeholder=True),
    seg("18-nt linked-read barcode", "N"*18, "cbc", placeholder=True, feature=LINKED),
    seg("P7 reverse complement", il.P7_RC, "p7"),
], name="TELL-seq Illumina library")
FINAL_LIBRARY = FINAL
SEQ_PRIMERS = ()
FINAL_CAPTION = "Figure 1B order: P5, 8-base sample index (Index 2), Read 1 region, 100–400 bp insert, Read 2 region, 18-base molecular barcode (Index 1), P7. Read-region bases are not published and remain role tokens."
SEQUENCING_ENDING = "The published Methods report 2 × 146 paired-end reads, an 18-cycle Index 1 read and an 8-cycle Index 2 read, run according to the manufacturer's protocols. Figure 1B places Index 1 primed from the Read 2 region and Index 2 primed from the Read 1 region. Whether each read needs a custom primer comes from the vendor sequencing guide, which this pass did not re-read."
SEQUENCING_UNAVAILABLE = "The paper and its supplements do not print the Read 1/Read 2 region or primer sequences, so exact annealing cannot be computed without fabricating bases."
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
    wf.react("Simultaneous MuA strand transfer and bead capture", rows(CAPTURED,"on bead"),
             note="0.1–5 ng DNA, 3–10 million barcoded TELL beads, non-barcoded spacer beads and transpososomes in one 0.2-mL tube; "
                  "STCs are captured on bead barcode oligos by hybridization (Supplemental Fig. S1A also lists a ligase in this mix).")
    wf.react("Second transpososome tags between STCs", rows(SECOND_TAGGED,"on bead"),
             note="The second transposon end supplies the other library-amplification priming site (Read 1 region in Figure 1B).")
    wf.react("Break STCs and wash the magnetic beads", rows(RELEASED,"on bead"),
             note="Disrupting the protein complexes fragments the target; captured barcoded DNA stays on the TELL beads.")
    wf.react("Amplify off beads; add P5 with sample index and P7", rows(FINAL,"library"),
             note="Eight cycles for human DNA and 13–14 for microbial DNA in the defining paper.")
    return wf


def protocol_graph():
    graph=ProtocolGraph(TITLE,authority=SOURCE)
    dna=graph.input("hmw_dna","high-molecular-weight genomic DNA",
                    MaterialKind.MOLECULAR,state=LONG_DNA,primary=True)
    bead_adapter=graph.input("bead_adapter","clonally barcoded TELL bead oligos (3-µm magnetic beads)",
                             MaterialKind.ADAPTER,inferred=True,
                             source=SOURCE+" (adapter bases not published)")
    spacer=graph.input("spacer_beads","non-barcoded spacer beads and viscous reaction solution",
                       MaterialKind.CONDITION)
    mua=graph.input("mua","MuA transpososomes",MaterialKind.CATALYST)
    ligase=graph.input("capture_ligase","ligase in the simultaneous barcoding mix",MaterialKind.CATALYST,
                       source=SOURCE+" (Supplemental Fig. S1A Method 1 lists transpososomes, barcoded beads, ligase and genomic DNA; "
                              "mapped to TELL-seq by the Results description of the simultaneous reaction; ligated junction not described)")
    captured=graph.product("captured","bead-captured strand-transfer complexes",
                           MaterialKind.MOLECULAR,state=CAPTURED,inferred=True,
                           source=SOURCE+" (adapter bases not published)")
    graph.react("capture_transfer","Combine MuA strand transfer with bead capture in one tube",
                StepKind.TAGMENTATION,
                ((dna,InputRole.TEMPLATE),(bead_adapter,InputRole.STRUCTURAL),
                 (mua,InputRole.CATALYST),(ligase,InputRole.CATALYST),
                 (spacer,InputRole.CONDITION)),(captured,))
    second_adapter=graph.input("second_adapter","second MuA transpososome (Read 1 region end)",
                               MaterialKind.ADAPTER,inferred=True,
                               source=SOURCE+" (adapter bases not published)")
    double_tagged=graph.product("double_tagged","two-ended bead-captured linked fragments",
                                MaterialKind.MOLECULAR,state=SECOND_TAGGED,inferred=True,
                                source=SOURCE+" (adapter bases not published)")
    like_ended=graph.product("like_ended_candidates",
                             "candidate fragments lacking one of the two unlike priming ends",
                             MaterialKind.MOLECULAR,inferred=True,
                             source=SOURCE+" (candidate from end geometry; not described in the paper)")
    graph.react("second_transfer","Second transpososome tags DNA between STCs",
                StepKind.TAGMENTATION,
                ((captured,InputRole.TEMPLATE),(second_adapter,InputRole.STRUCTURAL),
                 (mua,InputRole.CATALYST)),(double_tagged,like_ended))
    graph.mark_terminal(Terminal(like_ended.id,Disposition.UNRESOLVED,
        "Candidate only: such fragments would not be amplified exponentially, but the paper describes no step that removes them."))
    disruption=graph.input("disruption_condition","STC disruption condition (kit-specific; paper names heat, protease or SDS as STC-disassembling treatments)",
                           MaterialKind.CONDITION)
    released=graph.product("released","barcoded subfragments retained on TELL beads",
                           MaterialKind.MOLECULAR,state=RELEASED,inferred=True,
                           source=SOURCE+" (adapter bases not published)")
    graph.react("release","Break strand-transfer complexes",
                StepKind.DEPROTEINIZATION,
                ((double_tagged,InputRole.TEMPLATE),(disruption,InputRole.CONDITION)),
                (released,))
    magnet=graph.input("magnetic_wash","magnetic separation and bead washes",MaterialKind.SELECTION)
    washed=graph.product("washed_beads","washed TELL beads with captured barcoded DNA",
                         MaterialKind.MOLECULAR,state=RELEASED,inferred=True,
                         source=SOURCE+" (adapter bases not published)")
    wash_waste=graph.product("wash_supernatant","bead wash supernatant: unbound DNA, released proteins and reagents",
                             MaterialKind.WASTE)
    graph.react("wash_beads","Wash the magnetic TELL beads",StepKind.CLEANUP,
                ((released,InputRole.TEMPLATE),(magnet,InputRole.SELECTION)),(washed,wash_waste))
    graph.discard(wash_waste,"removed with the bead wash supernatant")
    forward=graph.input("pcr_forward","TELL-seq P5-side library PCR primer with 8-nt sample index",
                        MaterialKind.OLIGO,inferred=True,
                        source=SOURCE+" (primer sequence not published)")
    reverse=graph.input("pcr_reverse","TELL-seq P7-side library PCR primer",
                        MaterialKind.OLIGO,inferred=True,
                        source=SOURCE+" (primer sequence not published)")
    polymerase=graph.input("pcr_polymerase","library PCR polymerase (kit reagent, not named)",MaterialKind.CATALYST)
    library=graph.product("library",FINAL.name,MaterialKind.LIBRARY,state=FINAL,
                          inferred=True,source=SOURCE+" (Figure 1B architecture; primer-site bases not published)")
    graph.react("library_pcr","Amplify off beads (8 cycles human, 13–14 microbial); add P5, sample index and P7",StepKind.PCR,
                ((washed,InputRole.TEMPLATE),(forward,InputRole.STRUCTURAL),
                 (reverse,InputRole.STRUCTURAL),(polymerase,InputRole.CATALYST)),
                (library,),geometry=PCRGeometry.unavailable(
                    FINAL,reason="Figure 1B publishes the architecture but no primer or read-region bases."))
    graph.sequence(library,SequencingPlan(limitation=SEQUENCING_UNAVAILABLE,
                                          read_lengths={"Read 1":146,"Index 1 (i7)":18,
                                                        "Index 2 (i5)":8,"Read 2":146},
                                          read_uses={"Read 1":"include","Index 1 (i7)":"extract",
                                                     "Index 2 (i5)":"demultiplex","Read 2":"include"}))
    graph.declare_reviewed_scope((
        "hmw_dna","bead_adapter","spacer_beads","mua","capture_ligase","capture_transfer",
        "second_adapter","second_transfer","like_ended_candidates","disruption_condition",
        "release","magnetic_wash","wash_beads","pcr_forward",
        "pcr_reverse","pcr_polymerase","library_pcr","library"),SOURCE)
    return graph.require_publishable()
