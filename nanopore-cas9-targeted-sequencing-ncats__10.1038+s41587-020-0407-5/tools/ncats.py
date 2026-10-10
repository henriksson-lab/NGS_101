"""Molecular model for nanopore Cas9-targeted sequencing (nCATS)."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
from chemdraw import (Construct, MolecularState, Scene, Segment, Workflow,
                      annotation_rows, complement_segments)
from endprep import LigationJunction, dA_tailed_scene, repair_and_dA_tail
from selective_ligation import dephosphorylate_then_cut

def seg(name, top, tag=None, **kw): return Segment(name, top, tag, **kw)

GENOME = Construct([
    seg("left genomic flank", "X"*22, placeholder=True),
    seg("left Cas9 target", "N"*20, "me", placeholder=True),
    seg("region of interest", "X"*44, placeholder=True, length_bp=20000),
    seg("right Cas9 target", "N"*20, "me", placeholder=True),
    seg("right genomic flank", "X"*22, placeholder=True),
], name="high-molecular-weight genomic DNA")

DEPHOSPHORYLATED = GENOME

CUT_TARGET = Construct([
    seg("left Cas9-cleavage flank", "N"*10, "r1", placeholder=True),
    seg("region of interest", "X"*44, placeholder=True, length_bp=20000),
    seg("right Cas9-cleavage flank", "N"*10, "r2", placeholder=True),
], name="Cas9-flanked target")

SELECTIVE_ENDS = dephosphorylate_then_cut(GENOME, CUT_TARGET)

# Taq adds dA only after Cas9 has exposed ligatable, phosphorylated ends.  The shared
# end-prep object derives opposed 3' A overhangs from the target rather than drawing them.
DA_TARGET = repair_and_dA_tail(CUT_TARGET)

TA_JUNCTION = LigationJunction("T", "A")
FINAL_LIBRARY = Construct([
    seg("left motor-loaded ONT adapter", "[LSK109 motor adapter]", "r1", placeholder=True),
    seg("left T:A junction", TA_JUNCTION.adapter_overhang,
        bottom=TA_JUNCTION.insert_overhang),
    *CUT_TARGET.segments,
    seg("right A:T junction", TA_JUNCTION.insert_overhang,
        bottom=TA_JUNCTION.adapter_overhang),
    seg("right ONT adapter", "[LSK109 adapter]", "r2", placeholder=True),
], name="nCATS target library")

def rows(con, label=""):
    unknown=tuple(s.name+"'" for s in con if s.is_role_token())
    return [*Scene.duplex(list(con),label=label,unpaired=unknown).rows(),*annotation_rows(con)]

def end_state_rows(con, label, mod5):
    top=list(con); bottom=complement_segments(top)
    sc=Scene(); sc.strand("top",top,label=label,mod5=mod5)
    sc.anneal("bottom",bottom,to="top",pair=(top[0].name+"'",top[0].name),
              label=label,mod5=mod5)
    return [*sc.rows(),*annotation_rows(con)]

def final_rows():
    unknown=tuple(s.name+"'" for s in FINAL_LIBRARY if s.is_role_token())
    sc=Scene.duplex(list(FINAL_LIBRARY),label="nCATS target",unpaired=unknown)
    sc.junction("top","left T:A junction","left Cas9-cleavage flank","adapter ligation")
    sc.junction("top","right Cas9-cleavage flank","right A:T junction","adapter ligation")
    return [*sc.rows(),*annotation_rows(FINAL_LIBRARY)]

INITIAL_NAME = "High-molecular-weight genomic DNA"
INITIAL_ROWS = tuple(end_state_rows(GENOME,"gDNA","p"))

def workflow():
    wf=Workflow(MolecularState(INITIAL_NAME,INITIAL_ROWS))
    wf.react("Quick CIP dephosphorylation",end_state_rows(DEPHOSPHORYLATED,"old ends","OH"),
             note="Pre-existing DNA ends lose the 5-prime phosphate required for ligation.")
    wf.react("Cas9–guide RNP cleavage around the target",end_state_rows(SELECTIVE_ENDS.selected_fragment,"fresh cuts","p"),
             note="Fresh cuts retain ligatable 5-prime phosphates; the old ends remain blocked.")
    wf.react("Taq plus dATP at 72 °C",dA_tailed_scene(DA_TARGET).rows(),
             note="A-tailing prepares the newly exposed target ends for ONT adapter ligation.")
    wf.react("Quick Ligase plus LSK109 adapters",final_rows(),
             note="Adapters preferentially ligate at Cas9-created, phosphorylated ends.")
    return wf
