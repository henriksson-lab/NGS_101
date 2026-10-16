"""Construct model for the documented BD Rhapsody WTA alpha workflow.

Provenance.  The platform paper (Yoon et al., Nat Methods 2019) is a kit user and prints
no Rhapsody oligo.  The one real sequence drawn here, `RANDOMER_HANDLE`, is printed by BD
itself in the Whole Transcriptome Analysis Alpha Protocol 23-21179-00 (12/2018), a vendor
document mirrored in the scg_lib_structs repository under `data/BD/` (the manifest records
the mirrored copy and the hand-fetched one as md5-identical, so it is an ordinary vendor
original whichever name it sits under).

The remaining bead details are unresolved in the checked sources.  The BD instrument user guide (Doc ID
214062 Rev. 2.0) gives only the cartoon order `Univ | CL | UMI | poly(T)` with the legend
"Univ: universal oligo; CL: cell label; UMI: Unique Molecular Identifier" -- an order, not
bases and not lengths.  The checked guides give no bead bases, cell-label or UMI length, no
Universal Oligo sequence and no library-primer sequence, so those are symbolic placeholders created
`inferred=True`; the drawing therefore claims no width for them.

The scg_lib_structs drawn page is **not** a source here: it is inspiration only, nothing
is taken from it, and agreement with it is not verification.  Its own author describes its
bead and primer sequences as an educational guess reconstructed from read data, which is
not a vendor or author statement; an earlier revision of the research note reproduced
them, and they are no longer asserted anywhere.
"""
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
    """BD's published cartoon order only; every unresolved region is symbolic.

    BD gives neither bases nor lengths for the universal region, cell label or UMI, so
    each is a symbolic placeholder: the drawing must not imply a cycle count BD has not
    published.
    """
    return [_s("bead universal region","[UNIVERSAL OLIGO]",inferred=True,placeholder=True),
            _s("cell label","[CELL LABEL]","cbc",inferred=True,placeholder=True),
            _s("UMI","[UMI]","umi",inferred=True,placeholder=True),
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
                      _s("UMI'","[UMI]","umi",inferred=True,placeholder=True),
                      _s("cell label'","[CELL LABEL]","cbc",inferred=True,placeholder=True),
                      _s("bead universal complement","[UNIVERSAL OLIGO]",inferred=True,
                         placeholder=True)],name="released RPE product")

def final_library():
    return Construct([
        _s("unresolved P5/Read-1 arm","[P5 / READ 1 ARM]",inferred=True,placeholder=True),
        _s("cell label","[CELL LABEL]","cbc",inferred=True,placeholder=True),
        _s("UMI","[UMI]","umi",inferred=True,placeholder=True),
        _s("poly(T)/cDNA","TTTXXXXXXXX...XXXXXXXX",placeholder=True),
        _s("published Randomer handle",RANDOMER_HANDLE,"r2"),
        _s("unresolved index/P7 arm","[i7 / P7 ARM]",inferred=True,placeholder=True),
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
    beads=graph.input("capture_beads","cell-labelled capture beads (bases unverified)",
                      MaterialKind.ADAPTER,state=tuple(bead_capture_segments()),
                      inferred=True,source=source+" (exact bead sequence unresolved in checked guides)")
    microwells=graph.input("microwells","BD microwell array",MaterialKind.CONDITION)
    captured=graph.product("captured_rna","microwell-confined bead-captured poly(A) RNA",
                           MaterialKind.MOLECULAR,state=capture_scene(),inferred=True,
                           source=source+" (bead bases unresolved in checked guides)")
    empty=graph.product("empty_wells","unproductive cell/bead loading outcomes",
                        MaterialKind.WASTE)
    graph.react("partition_capture","INFERRED — Load cells and beads into microwells (the checked guide gives bead region order, not bases or lengths)",StepKind.PARTITION,
                ((cells,InputRole.TEMPLATE),(beads,InputRole.STRUCTURAL),
                 (microwells,InputRole.CONDITION)),(captured,empty))
    graph.discard(empty,"wells without a productive cell-bead pairing do not yield library")
    rt=graph.input("reverse_transcriptase","reverse transcriptase",MaterialKind.CATALYST)
    first=graph.product("first_strand",bead_cdna().name,MaterialKind.MOLECULAR,
                        state=bead_cdna(),inferred=True,
                        source=source+" (bead-proximal bases unresolved)")
    graph.react("reverse_transcribe","INFERRED — Synthesize bead-barcoded first-strand cDNA",
                StepKind.REVERSE_TRANSCRIPTION,
                ((captured,InputRole.TEMPLATE),(beads,InputRole.STRUCTURAL),
                 (rt,InputRole.CATALYST)),(first,))
    exonuclease=graph.input("exonuclease","Exonuclease I",MaterialKind.CATALYST)
    alkali=graph.input("rna_removal","RNA-removal condition",MaterialKind.CONDITION)
    cleaned=graph.product("cleaned_first_strand","RNA-free bead-bound first-strand cDNA",
                          MaterialKind.MOLECULAR,state=bead_cdna(),inferred=True,
                          source=source+" (bead-proximal bases unresolved)")
    cleanup_waste=graph.product("cleanup_waste","unused capture oligos and removed RNA",
                                MaterialKind.WASTE)
    graph.react("cleanup","INFERRED — Remove unused capture oligos and RNA",StepKind.CLEANUP,
                ((first,InputRole.TEMPLATE),(exonuclease,InputRole.CATALYST),
                 (alkali,InputRole.CONDITION)),(cleaned,cleanup_waste))
    graph.discard(cleanup_waste,"unused oligos and RNA are excluded from RPE products")
    randomer=graph.input("randomer","published BD Randomer",MaterialKind.OLIGO,
                         state=tuple(randomer_segments()))
    klenow=graph.input("klenow","Klenow exo-minus polymerase",MaterialKind.CATALYST)
    extended=graph.product("rpe_product",rpe_product().name,MaterialKind.MOLECULAR,
                           state=rpe_product(),inferred=True,
                           source=source+" (universal bead-arm complement unresolved)")
    graph.react("random_extension","INFERRED — Extend the random primer toward the bead (the Randomer is printed; copied bead bases remain unverified)",
                StepKind.EXTENSION,
                ((cleaned,InputRole.TEMPLATE),(randomer,InputRole.STRUCTURAL),
                 (klenow,InputRole.CATALYST)),(extended,))
    universal=graph.input("universal_oligo","Universal Oligo (bases unverified)",MaterialKind.OLIGO,
                          inferred=True,source=source+" (sequence unresolved)")
    rpe_reverse=graph.input("rpe_reverse","RPE reverse primer",MaterialKind.OLIGO,
                            inferred=True,source=source+" (sequence unresolved)")
    rpe_pol=graph.input("rpe_polymerase","RPE PCR polymerase",MaterialKind.CATALYST)
    amplified=graph.product("amplified_rpe","amplified cell-labelled RPE products",
                            MaterialKind.MOLECULAR,state=rpe_product(),inferred=True,
                            source=source+" (primer boundaries unresolved)")
    graph.react("rpe_pcr","INFERRED — Amplify random-primer-extension products",StepKind.PCR,
                ((extended,InputRole.TEMPLATE),(universal,InputRole.STRUCTURAL),
                 (rpe_reverse,InputRole.STRUCTURAL),(rpe_pol,InputRole.CATALYST)),
                (amplified,),geometry=PCRGeometry.unavailable(
                    rpe_product(),reason="The checked BD alpha guide does not provide the Universal Oligo or both RPE primer-binding sequences."))
    library_forward=graph.input("library_forward","library forward primer (bases unverified)",
                                MaterialKind.OLIGO,inferred=True,
                                source=source+" (sequence unresolved)")
    library_reverse=graph.input("library_reverse","indexed library reverse primer (bases unverified)",
                                MaterialKind.OLIGO,inferred=True,
                                source=source+" (sequence unresolved)")
    library_pol=graph.input("library_polymerase","library PCR polymerase",MaterialKind.CATALYST)
    final=final_library()
    library=graph.product("library",final.name,MaterialKind.LIBRARY,state=final,
                          inferred=True,source=source+" (outer arms unresolved)")
    graph.react("library_pcr","INFERRED — Add the indexed sequencing-library arms",StepKind.PCR,
                ((amplified,InputRole.TEMPLATE),(library_forward,InputRole.STRUCTURAL),
                 (library_reverse,InputRole.STRUCTURAL),(library_pol,InputRole.CATALYST)),
                (library,),geometry=PCRGeometry.unavailable(
                    final,reason="The checked BD alpha guide does not provide the library-primer sequences."))
    graph.sequence(library,SequencingPlan(
        acquired_roles=("Read 1","Index 1 (i7)","Read 2"),
        read_lengths={"Read 1":101,"Read 2":101},
        read_uses={"Read 1":"extract","Index 1 (i7)":"demultiplex","Read 2":"include"},
        limitation="The checked BD sources do not establish the sequencing-primer sites."))
    graph.declare_reviewed_scope((
        "cells","capture_beads","microwells","partition_capture","reverse_transcriptase",
        "reverse_transcribe","exonuclease","rna_removal","cleanup","randomer","klenow",
        "random_extension","universal_oligo","rpe_reverse","rpe_polymerase","rpe_pcr",
        "library_forward","library_reverse","library_polymerase","library_pcr","library"),source)
    return graph.require_publishable()
