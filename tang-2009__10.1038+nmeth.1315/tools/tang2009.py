"""Sequence model for Tang 2009 single-cell mRNA sequencing.

The defining 2009 method used an unbarcoded SOLiD fragment library.  The 16-plex
P2 adapter belongs to the 2010 Nature Protocols extension and is deliberately not
part of :func:`final_library`.
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))

from chemdraw import Construct, Scene, Segment, revcomp, complement_segments
from adapter_ligation import AdapterLigationGeometry
from pcr import PCRGeometry, PCRPrimer
from protocol_graph import (Disposition, Terminal, InputRole,MaterialKind,ProtocolGraph,SequencingPlan,StepKind)
UP1_HANDLE = "ATATGGATCCGGCGCGCCGTCGAC"
UP2_HANDLE = "ATATCTCGAGGGCGCGCCGGATCC"
DT_LEN = 24
UP1 = UP1_HANDLE + "T" * DT_LEN
UP2 = UP2_HANDLE + "T" * DT_LEN
AUP1, AUP2 = UP1, UP2


def amplified_cdna():
    return Construct([Segment("UP2 handle", UP2_HANDLE), Segment("poly(T)", "T"*24),
                      role("cDNA", "transcript-derived DNA"), Segment("poly(A)", "A"*24),
                      Segment("UP1 handle'", revcomp(UP1_HANDLE))], name="Tang native cDNA PCR daughter")


def role(name, text):
    return Segment(name, "["+text+"]", placeholder=True)


def native_state(stage):
    sc=Scene()
    if stage == "sample":
        sc.strand("rna", [role("RNA", "cellular RNA; poly(A) transcripts are RT substrates")])
        sc.strand("other", [role("other", "other lysate polymers; recovery not class-resolved")])
        return sc
    if stage in ("rt", "exo", "tail"):
        parts=[Segment("UP1",UP1),role("copy","first-strand DNA copy of RNA")]
        if stage == "tail": parts += [role("TdT tail","variable 3′ poly(dA) tail; length unassigned")]
        sc.strand("first",parts,label="native first-strand cDNA")
    elif stage == "second":
        sc.strand("second",list(amplified_cdna()),label="actual UP2-extended second-strand DNA; not yet PCR-amplified")
        sc.strand("original",[Segment("UP1",UP1),role("first copy","original first-strand DNA"),role("tail","variable TdT poly(dA)")],label="original first strand; annealed extent unassigned")
    elif stage in ("pcr", "blocked"):
        parts=list(amplified_cdna())
        sc.strand("daughter",parts,label="PCR daughter",mod5="NH2-C6" if stage=="blocked" else "")
        sc.anneal("daughter_complement",complement_segments(parts),to="daughter",pair=("UP2 handle'","UP2 handle"),unpaired=("cDNA'",),mod5="NH2-C6" if stage=="blocked" else "")
        sc.strand("originals",[role("originals","possible co-carried original cDNA / prior PCR templates; amine status is lineage-specific")])
    elif stage in ("fragment", "repaired"):
        sc.strand("fragments",[role("fragments","80–130 bp sheared DNA; internal fragments plus possible NH2-C6 terminal fragments")])
        sc.strand("opposite",[role("opposite","opposite fragment strands; exact repair geometry unassigned" if stage=="repaired" else "opposite fragment strands / shearing ends unassigned")])
        sc.strand("blocked_end",[role("blocked WTA end","separate NH2-C6-primer-derived terminal fragment; blocking survives shearing and is not automatically removed by end repair")],mod5="NH2-C6")
    elif stage in ("ligated", "library"):
        if stage=="ligated":
            sc.strand("top",[Segment("P1",P1_TOP),role("insert","sheared DNA fragment"),Segment("native P2 upper including distal TT",P2_TOP)],label="native ligated strand; nick/fill extent unassigned")
            sc.strand("opposite",[Segment("native P2 lower",P2_BOTTOM),role("opposite insert","opposite fragment strand"),Segment("P1 lower",P1_BOTTOM)],label="opposite native strand; not completed PCR duplex")
        else:
            parts=list(final_library());sc.strand("daughter",parts,label="SOLiD PCR daughter")
            sc.anneal("opposite",complement_segments(parts),to="daughter",pair=("SOLiD P1'","SOLiD P1"))
        sc.strand("residuals",[role("residuals","possible unligated / blocked terminal fragments and original ligated templates; recovery unassigned")])
    sc.strand("rna",[role("RNA residual","original RNA or RNA fragments; this capped source does not assign its chemical/physical fate")])
    sc.strand("other",[role("lysate residual","possible co-carried lysate polymers; recoveries not class-resolved")])
    if stage not in ("rt",):
        sc.strand("primer_fragments",[role("primer fragments","ExoI products / surviving oligo; digestion is not physical removal")])
    return sc



# Unbarcoded SOLiD adapters, verbatim from the 2009 Nature Methods Supplementary Table
# (Supplementary Data, `41592_2009_BFnmeth1315_MOESM306_ESM.xls`), which prints
# "P1 Adaptor (Upper/Lower strand)", "P2 Adaptor (Upper/Lower strand)" and
# "Library PCR Primer 1/2" in spaced triplets; the lower strands are written 3'->5'
# there and are stored here 5'->3'.
P1_TOP = "CCACTACGCCTCCGCTTTCCTCTCTATGGGCAGTCGGTGAT"
P1_BOTTOM = "ATCACCGACTGCCCATAGAGAGGAAAGCGGAGGCGTAGTGGTT"
P2_TOP = "AGAGAATGAGGAACCCGGGGCAGTT"
P2_BOTTOM = "CTGCCCCGGGTTCCTCATTCTCT"
LIBRARY_PCR_1 = "CCACTACGCCTCCGCTTTCCTCTCTATG"
LIBRARY_PCR_2 = "CTGCCCCGGGTTCCTCATTCT"

# The later 2010 protocol adds a 19-nt internal adapter and one of sixteen 6-nt
# barcodes to P2.  These constants identify the extension without presenting it
# as part of the defining 2009 construct.  Source: Tang et al. 2010, Nat Protoc
# 5:516-535 (PMC3847604), which prints all sixteen "CGCCTTGGCCGTACAGCAG<6-nt
# barcode>AGAGAATGAGGAACCCGGGGCAGTT" P2 adapters and the two 16Barcode library
# PCR primers.
BARCODE2010_INTERNAL = "CGCCTTGGCCGTACAGCAG"
BARCODE2010_LEN = 6
BARCODE2010_COUNT = 16


def final_library() -> Construct:
    """2009 unbarcoded SOLiD library, shown in the orientation P1 -> insert -> P2."""
    return Construct([
        Segment("SOLiD P1", P1_TOP, "p5"),
        Segment("cDNA fragment", "XXXXXXXX...XXXXXXXX", placeholder=True),
        # PCR primer 2 leaves this 23-nt P2-derived end in the amplified library.
        Segment("SOLiD P2", revcomp(P2_BOTTOM), "p7"),
    ], name="Tang 2009 unbarcoded SOLiD library")


def solid_adapter_geometry():
    scenes=[]
    for name,top,bottom in (("P1",P1_TOP,P1_BOTTOM),("P2",P2_TOP,P2_BOTTOM)):
        sc=Scene();sc.strand("top",[Segment("top",top)])
        if name=="P1":
            lower=[Segment("paired",bottom[:-2]),Segment("distal TT",bottom[-2:])]
            sc.anneal("lower",lower,to="top",pair=("paired","top"))
        else:
            sc=Scene();sc.strand("top",[Segment("paired",top[:-2]),Segment("distal TT",top[-2:])])
            sc.anneal("lower",[Segment("paired lower",bottom)],to="top",pair=("paired lower","paired"))
        sc.ligation_end("top","top" if name=="P1" else "paired",kind="blunt",owner="adapter",text="fragment-facing blunt end; distal 3′ TT retained")
        scenes.append(sc)
    return AdapterLigationGeometry.known("blunt",*scenes,statement="Published Table9 gives both unbarcoded adapter strands; fragment-facing blunt end and distal 3′ TT are distinct ends. No unreported adapter phosphorylation is assigned.")


def protocol_graph():
    source="https://doi.org/10.1038/nmeth.1315; Online Methods Single-cell cDNA preparation and mRNA-Seq library preparation; published Table9"
    g=ProtocolGraph("Tang 2009",authority=source)
    M,R,K=MaterialKind,InputRole,StepKind
    def inp(i,label,kind=M.CATALYST,state=None):return g.input(i,label,kind,state=state)
    def prod(i,label,stage,kind=M.MOLECULAR):return g.product(i,label,kind,state=native_state(stage))
    def cleanup(i,label,x):
        reagent=inp(i+"_selection",label+" recovery; exact chemistry unassigned",M.SELECTION)
        y=g.product(i+"_retained",label+" retained mixture",x.kind,state=x.state)
        waste=g.product(i+"_waste","unretained fraction; class composition unassigned",M.WASTE)
        g.react(i,label,K.CLEANUP,((x,R.TEMPLATE),(reagent,R.SELECTION)),(y,waste));g.discard(waste,"physically unretained fraction; no complete RNA/oligo clearance asserted")
        return y
    cell=g.input("cell","single-cell lysate",M.BIOLOGICAL,state=native_state("sample"),primary=True)
    up1=inp("up1","UP1",M.OLIGO,Construct([Segment("UP1",UP1)]));up2=inp("up2","UP2",M.OLIGO,Construct([Segment("UP2",UP2)]))
    rt=inp("rt","reverse transcriptase; enzyme identity deferred to cited predecessor")
    first=prod("first","UP1-primed ss cDNA plus original RNA","rt")
    g.react("reverse_transcribe","UP1 reverse transcription 50°C30min; 70°C15min inactivation",K.REVERSE_TRANSCRIPTION,((cell,R.TEMPLATE),(up1,R.STRUCTURAL),(rt,R.CATALYST)),(first,))
    exo=inp("exoi","Exonuclease I");clean=prod("primer_digested","cDNA mixture plus primer hydrolysates","exo")
    g.react("remove_up1","Digest nonreactive UP1 with Exonuclease I",K.DIGESTION,((first,R.TEMPLATE),(exo,R.CATALYST)),(clean,))
    tdt=inp("tdt","TdT");datp=inp("datp","dATP",M.CONDITION);tail=prod("tailed","poly(dA)-tailed native first strand","tail")
    g.react("tail","TdT poly(dA) tailing; native tail length unassigned",K.TAILING,((clean,R.TEMPLATE),(tdt,R.CATALYST),(datp,R.CONDITION)),(tail,))
    pol=inp("cdna_polymerase","cDNA DNA polymerase; cited predecessor conditions not reread")
    second=prod("second","actual UP2-extended strands plus original first strands","second")
    g.react("second_strand","Extend UP2 on poly(dA)-tailed cDNA",K.EXTENSION,((tail,R.TEMPLATE),(up2,R.STRUCTURAL),(pol,R.CATALYST)),(second,))
    amp=prod("cdna_amplified","20-cycle cDNA daughters plus residual originals","pcr")
    g.react("cdna_pcr","UP1/UP2 PCR20cycles (95°C30s,67°C1min,72°C6min+6s/cycle)",K.PCR,((second,R.TEMPLATE),(up1,R.STRUCTURAL),(up2,R.STRUCTURAL),(pol,R.CATALYST)),(amp,),geometry=PCRGeometry.unavailable(second.state,reason="Actual UP2-extended native strands, first-strand originals and residual RNA enter PCR; variable tail/landing positions are not assigned exact geometry."))
    recovered=cleanup("cdna_cleanup","Purify first PCR product",amp)
    astate=Scene();astate.strand("aup1",[Segment("AUP1",AUP1)],mod5="NH2-C6")
    bstate=Scene();bstate.strand("aup2",[Segment("AUP2",AUP2)],mod5="NH2-C6")
    a1=inp("aup1","NH2-C6-UP1",M.OLIGO,astate);a2=inp("aup2","NH2-C6-UP2",M.OLIGO,bstate)
    blocked=prod("blocked","amine-primer-derived daughters; prior templates distinct","blocked")
    g.react("blocked_pcr","NH2-C6-UP1/UP2 PCR9cycles",K.PCR,((recovered,R.TEMPLATE),(a1,R.STRUCTURAL),(a2,R.STRUCTURAL),(pol,R.CATALYST)),(blocked,),geometry=PCRGeometry.unavailable(recovered.state,reason="The amine resides on newly incorporated primer ends; original template strands are not retroactively amine-modified."))
    recovered=cleanup("blocked_cleanup","Purify amine-primer PCR products",blocked)
    shear=inp("covaris","Covaris S2; 200–500ng cDNA (0.5–3kb) to80–130bp",M.CONDITION)
    fragments=prod("fragments","sheared internal and end-blocked fragments","fragment")
    g.react("fragment","Shear cDNA to80–130bp",K.FRAGMENTATION,((recovered,R.TEMPLATE),(shear,R.CONDITION)),(fragments,))
    repair=inp("end_repair_mix","SOLiD low-input end repair; vendor implementation unread")
    repaired=prod("repaired","repaired native fragments; detailed end chemistry bounded","repaired")
    g.react("repair","Repair fragment ends; low-input kit details deferred",K.END_REPAIR,((fragments,R.TEMPLATE),(repair,R.CATALYST)),(repaired,))
    adapters=inp("adapters","published P1/P2 strands; no index",M.ADAPTER,(P1_TOP,P1_BOTTOM,P2_TOP,P2_BOTTOM))
    ligase=inp("ligase","DNA ligase");ligated=prod("ligated","native P1/P2 ligated strands and residual classes","ligated")
    g.react("ligate","Ligate unbarcoded P1/P2 adapters",K.LIGATION,((repaired,R.TEMPLATE),(adapters,R.STRUCTURAL),(ligase,R.CATALYST)),(ligated,),geometry=solid_adapter_geometry())
    gel=inp("gel","6%PAGE;150–200bp fraction",M.SELECTION);selected=prod("selected","selected native ligated population","ligated");w=g.product("offsize","unselected gel fractions",M.WASTE)
    g.react("size_select","Excise150–200bp gel fraction",K.SELECTION,((ligated,R.TEMPLATE),(gel,R.SELECTION)),(selected,w));g.discard(w,"outside excised gel fraction")
    f=inp("lib_primer_1","Library PCR Primer1",M.OLIGO,Construct([Segment("primer1",LIBRARY_PCR_1)]));r=inp("lib_primer_2","Library PCR Primer2",M.OLIGO,Construct([Segment("primer2",LIBRARY_PCR_2)]))
    lib=prod("amplified_library","SOLiD daughters plus possibly co-carried native substrates","library",M.LIBRARY)
    g.react("library_pcr","Library PCR8–10cycles; cycle selection by FlashGel aliquot",K.PCR,((selected,R.TEMPLATE),(f,R.STRUCTURAL),(r,R.STRUCTURAL),(pol,R.CATALYST)),(lib,),geometry=PCRGeometry.unavailable(selected.state,reason="Actual ligated strands enter PCR; native nick/fill state and exact heterogeneous fragment boundaries remain unassigned."))
    final=cleanup("library_cleanup","Purify SOLiD fragment library",lib)
    g.sequence(final,SequencingPlan(acquired_roles=("Read 1",),read_uses={"Read 1":"include"},read_lengths={"Read 1":50},limitation="Online Methods declares50-base SOLiD; supplementary runs also35mer. Emulsion PCR uses500pg library and1.6billion P1-primed beads; bead-primer/linkage/full emulsion chemistry and loaded sequencing oligo remain outside this library-preparation graph."))
    g.declare_reviewed_scope(("blocked_cleanup_selection","cdna_cleanup_selection","cdna_polymerase","covaris","datp","end_repair_mix","exoi","gel","library_cleanup_retained","library_cleanup_selection","ligase","rt","tdt","cell","up1","up2","reverse_transcribe","remove_up1","tail","second_strand","cdna_pcr","cdna_cleanup","aup1","aup2","blocked_pcr","blocked_cleanup","fragment","repair","adapters","ligate","size_select","lib_primer_1","lib_primer_2","library_pcr","library_cleanup"),source)
    return g.require_publishable()
