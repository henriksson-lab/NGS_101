"""Molecular architecture of the original TELL-seq linked-read protocol."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
from chemdraw import (Construct, MolecularState, Scene, Segment, Workflow,
                      annotation_rows, feature)
import illumina as il

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
