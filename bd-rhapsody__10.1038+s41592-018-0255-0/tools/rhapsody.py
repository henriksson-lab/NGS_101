"""Construct model for the documented BD Rhapsody WTA alpha workflow."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"lib"))
from chemdraw import Construct,Scene,Segment,feature
from pcr import PCRGeometry
from protocol_graph import (InputRole,MaterialKind,ProtocolGraph,SequencingPlan,StepKind)

CELL_LABEL=feature("cell_label","cell_barcode","whitelist",whitelist="bd-rhapsody-cell-label-set")
UMI_FEATURE=feature("umi","umi","random")
FEATURES={"cell label":CELL_LABEL,"cell label'":CELL_LABEL,"UMI":UMI_FEATURE,"UMI'":UMI_FEATURE}

# Printed by BD in Whole Transcriptome Analysis Alpha Protocol 23-21179-00 (12/2018).
RANDOMER_HANDLE="TCAGACGTGTGCTCTTCCGATCT"
RANDOM_NT=9
DRAWN_POLYT=15  # visual tract only; BD does not publish the bead oligo's exact sequence

def _s(n,t,g=None,**kw):
    kw.setdefault("feature",FEATURES.get(n)); return Segment(n,t,g,**kw)

def bead_capture_segments():
    """Vendor-published architecture, with proprietary bases represented symbolically."""
    return [_s("bead universal region","XXXXXXXXXXXX",inferred=True,placeholder=True),
            _s("cell label","NNNNNNNNN","cbc",inferred=True,placeholder=True),
            _s("UMI","NNNNNNNN","umi",inferred=True,placeholder=True),
            _s("poly(T)","T"*DRAWN_POLYT,inferred=True)]

def randomer_segments():
    return [_s("published R2-side handle",RANDOMER_HANDLE,"r2"),
            _s("random N9","N"*RANDOM_NT,placeholder=True)]

def capture_scene():
    mrna=[_s("mRNA body","XXXXXXXXXXXX",placeholder=True),
          _s("poly(A)","A"*DRAWN_POLYT)]
    sc=Scene();sc.strand("mRNA",mrna,label="mRNA")
    sc.anneal("bead oligo",bead_capture_segments(),to="mRNA",
              pair=("poly(T)","poly(A)"),label="bead oligo")
    sc.arrow("bead oligo","reverse transcriptase")
    return sc

def bead_cdna():
    return Construct([*bead_capture_segments(),
                      _s("first-strand cDNA","XXXXXXXX...XXXXXXXX",placeholder=True)],
                     name="bead-bound first-strand cDNA")

def random_priming_scene():
    template=[_s("bead-proximal cDNA","XXXXXXXX...XXXXXXXX",placeholder=True),
              _s("random site","XXXXXXXXX",placeholder=True)]
    sc=Scene();sc.strand("first strand",template,label="bead-bound cDNA")
    sc.anneal("Randomer",randomer_segments(),to="first strand",
              pair=("random N9","random site"),label="Randomer")
    sc.arrow("Randomer","Klenow exo-")
    return sc

def rpe_product():
    return Construct([*randomer_segments()[:-1],
                      _s("random-primed cDNA","XXXXXXXX...XXXXXXXX",placeholder=True),
                      _s("poly(A) copy","AAA",placeholder=True),
                      _s("UMI'","N"*8,"umi",inferred=True,placeholder=True),
                      _s("cell label'","N"*9,"cbc",inferred=True,placeholder=True),
                      _s("bead universal complement","XXXXXXXXXXXX",inferred=True,
                         placeholder=True)],name="released RPE product")

def final_library():
    return Construct([
        _s("unpublished P5/Read-1 arm","[P5 / READ 1 ARM]",inferred=True,placeholder=True),
        _s("cell label","[CELL LABEL]","cbc",inferred=True,placeholder=True),
        _s("UMI","[UMI]","umi",inferred=True,placeholder=True),
        _s("poly(T)/cDNA","TTTXXXXXXXX...XXXXXXXX",placeholder=True),
        _s("published Randomer handle",RANDOMER_HANDLE,"r2"),
        _s("unpublished index/P7 arm","[i7 / P7 ARM]",inferred=True,placeholder=True),
    ],name="BD Rhapsody WTA library boundary")

def _validate():
    if "".join(s.top for s in randomer_segments()) != RANDOMER_HANDLE+"N"*RANDOM_NT:
        raise ValueError("Randomer no longer matches the BD alpha protocol")
    capture_scene().rows();random_priming_scene().rows()
_validate()

def protocol_graph():
    source=("Chemistry authority: BD Whole Transcriptome Analysis Alpha Protocol "
            "23-21179-00 (December 2018); platform paper: "
            "https://doi.org/10.1038/s41592-018-0255-0")
    graph=ProtocolGraph("BD Rhapsody WTA",authority=source)
    cells=graph.input("cells","single-cell suspension",MaterialKind.BIOLOGICAL,primary=True)
    beads=graph.input("capture_beads","proprietary cell-labelled capture beads",
                      MaterialKind.ADAPTER,state=tuple(bead_capture_segments()),
                      inferred=True,source=source+" (exact bead sequence proprietary)")
    microwells=graph.input("microwells","BD microwell array",MaterialKind.CONDITION)
    captured=graph.product("captured_rna","microwell-confined bead-captured poly(A) RNA",
                           MaterialKind.MOLECULAR,state=capture_scene(),inferred=True,
                           source=source+" (bead bases proprietary)")
    empty=graph.product("empty_wells","unproductive cell/bead loading outcomes",
                        MaterialKind.WASTE)
    graph.react("partition_capture","Load cells and beads into microwells",StepKind.PARTITION,
                ((cells,InputRole.TEMPLATE),(beads,InputRole.STRUCTURAL),
                 (microwells,InputRole.CONDITION)),(captured,empty))
    graph.discard(empty,"wells without a productive cell-bead pairing do not yield library")
    rt=graph.input("reverse_transcriptase","reverse transcriptase",MaterialKind.CATALYST)
    first=graph.product("first_strand",bead_cdna().name,MaterialKind.MOLECULAR,
                        state=bead_cdna(),inferred=True,
                        source=source+" (bead-proximal bases proprietary)")
    graph.react("reverse_transcribe","Synthesize bead-barcoded first-strand cDNA",
                StepKind.REVERSE_TRANSCRIPTION,
                ((captured,InputRole.TEMPLATE),(beads,InputRole.STRUCTURAL),
                 (rt,InputRole.CATALYST)),(first,))
    exonuclease=graph.input("exonuclease","Exonuclease I",MaterialKind.CATALYST)
    alkali=graph.input("rna_removal","RNA-removal condition",MaterialKind.CONDITION)
    cleaned=graph.product("cleaned_first_strand","RNA-free bead-bound first-strand cDNA",
                          MaterialKind.MOLECULAR,state=bead_cdna(),inferred=True,
                          source=source+" (bead-proximal bases proprietary)")
    cleanup_waste=graph.product("cleanup_waste","unused capture oligos and removed RNA",
                                MaterialKind.WASTE)
    graph.react("cleanup","Remove unused capture oligos and RNA",StepKind.CLEANUP,
                ((first,InputRole.TEMPLATE),(exonuclease,InputRole.CATALYST),
                 (alkali,InputRole.CONDITION)),(cleaned,cleanup_waste))
    graph.discard(cleanup_waste,"unused oligos and RNA are excluded from RPE products")
    randomer=graph.input("randomer","published BD Randomer",MaterialKind.OLIGO,
                         state=tuple(randomer_segments()))
    klenow=graph.input("klenow","Klenow exo-minus polymerase",MaterialKind.CATALYST)
    extended=graph.product("rpe_product",rpe_product().name,MaterialKind.MOLECULAR,
                           state=rpe_product(),inferred=True,
                           source=source+" (universal bead-arm complement proprietary)")
    graph.react("random_extension","Extend the random primer toward the bead",
                StepKind.EXTENSION,
                ((cleaned,InputRole.TEMPLATE),(randomer,InputRole.STRUCTURAL),
                 (klenow,InputRole.CATALYST)),(extended,))
    universal=graph.input("universal_oligo","proprietary Universal Oligo",MaterialKind.OLIGO,
                          inferred=True,source=source+" (sequence proprietary)")
    rpe_reverse=graph.input("rpe_reverse","RPE reverse primer",MaterialKind.OLIGO,
                            inferred=True,source=source+" (sequence proprietary)")
    rpe_pol=graph.input("rpe_polymerase","RPE PCR polymerase",MaterialKind.CATALYST)
    amplified=graph.product("amplified_rpe","amplified cell-labelled RPE products",
                            MaterialKind.MOLECULAR,state=rpe_product(),inferred=True,
                            source=source+" (primer boundaries proprietary)")
    graph.react("rpe_pcr","Amplify random-primer-extension products",StepKind.PCR,
                ((extended,InputRole.TEMPLATE),(universal,InputRole.STRUCTURAL),
                 (rpe_reverse,InputRole.STRUCTURAL),(rpe_pol,InputRole.CATALYST)),
                (amplified,),geometry=PCRGeometry.unavailable(
                    rpe_product(),reason="BD does not publish the Universal Oligo or both RPE primer-binding sequences."))
    library_forward=graph.input("library_forward","proprietary library forward primer",
                                MaterialKind.OLIGO,inferred=True,
                                source=source+" (sequence proprietary)")
    library_reverse=graph.input("library_reverse","proprietary indexed library reverse primer",
                                MaterialKind.OLIGO,inferred=True,
                                source=source+" (sequence proprietary)")
    library_pol=graph.input("library_polymerase","library PCR polymerase",MaterialKind.CATALYST)
    final=final_library()
    library=graph.product("library",final.name,MaterialKind.LIBRARY,state=final,
                          inferred=True,source=source+" (outer arms proprietary)")
    graph.react("library_pcr","Add the indexed sequencing-library arms",StepKind.PCR,
                ((amplified,InputRole.TEMPLATE),(library_forward,InputRole.STRUCTURAL),
                 (library_reverse,InputRole.STRUCTURAL),(library_pol,InputRole.CATALYST)),
                (library,),geometry=PCRGeometry.unavailable(
                    final,reason="BD does not publish the proprietary library-primer sequences."))
    graph.sequence(library,SequencingPlan(
        acquired_roles=("Read 1","Index 1 (i7)","Read 2"),
        read_lengths={"Read 1":101,"Read 2":101},
        read_uses={"Read 1":"extract","Index 1 (i7)":"demultiplex","Read 2":"include"},
        limitation="BD does not publish the proprietary sequencing-primer sites."))
    graph.declare_reviewed_scope((
        "cells","capture_beads","microwells","partition_capture","reverse_transcriptase",
        "reverse_transcribe","exonuclease","rna_removal","cleanup","randomer","klenow",
        "random_extension","universal_oligo","rpe_reverse","rpe_polymerase","rpe_pcr",
        "library_forward","library_reverse","library_polymerase","library_pcr","library"),source)
    return graph.require_publishable()
