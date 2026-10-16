"""Sequence model for Tang 2009 single-cell mRNA sequencing.

The defining 2009 method used an unbarcoded SOLiD fragment library.  The 16-plex
P2 adapter belongs to the 2010 Nature Protocols extension and is deliberately not
part of :func:`final_library`.
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))

from chemdraw import Construct, Scene, Segment, revcomp
from adapter_ligation import AdapterLigationGeometry
from pcr import PCRGeometry, PCRPrimer
from protocol_graph import (Disposition, Terminal, InputRole,MaterialKind,ProtocolGraph,SequencingPlan,StepKind)
from tang_cdna import (AUP1, AUP2, DT_LEN, UP1, UP1_HANDLE, UP2, UP2_HANDLE,
                       amplified_cdna)


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


def solid_adapter_geometry() -> AdapterLigationGeometry:
    scenes=[]
    for name,top,bottom in (("P1",P1_TOP,P1_BOTTOM),("P2",P2_TOP,P2_BOTTOM)):
        sc=Scene();sc.strand(f"{name} top",[Segment(f"{name} top strand",top)],label=name)
        sc.anneal(f"{name} bottom",[Segment(f"{name} bottom strand",bottom)],
                  to=f"{name} top",pair=(f"{name} bottom strand",f"{name} top strand"),
                  label=name,unpaired=(f"{name} bottom strand",))
        sc.ligation_end(f"{name} top",f"{name} top strand",kind="blunt",owner="adapter",
                        text=f"{name} fragment-ligating blunt end")
        scenes.append(sc)
    return AdapterLigationGeometry.known(
        "blunt",*scenes,
        statement="The 2009 supplementary table prints both strands of the unbarcoded SOLiD P1 and P2 adapters used after end repair.")


def protocol_graph() -> ProtocolGraph:
    source="https://doi.org/10.1038/nmeth.1315; Online Methods and Supplementary Table 9"
    graph=ProtocolGraph("Tang 2009 single-cell mRNA sequencing",authority=source)
    cell=graph.input("single_cell","single cell lysate",MaterialKind.BIOLOGICAL,primary=True)
    up1=graph.input("up1","UP1 oligo-dT reverse-transcription primer",MaterialKind.OLIGO,
                    state=(Segment("UP1",UP1),))
    rt=graph.input("reverse_transcriptase","reverse transcriptase",MaterialKind.CATALYST)
    first_state=Construct([Segment("UP1",UP1),Segment("first-strand cDNA","X"*24,
                                                      placeholder=True)],
                          name="UP1-primed first-strand cDNA")
    first=graph.product("first_strand",first_state.name,MaterialKind.MOLECULAR,
                        state=first_state)
    graph.react("reverse_transcribe","Reverse-transcribe poly(A) RNA with UP1 (50 C, 30 min)",
                StepKind.REVERSE_TRANSCRIPTION,
                ((cell,InputRole.TEMPLATE),(up1,InputRole.STRUCTURAL),(rt,InputRole.CATALYST)),
                (first,))
    exo=graph.input("exonuclease","Exonuclease I",MaterialKind.CATALYST)
    cleaned=graph.product("cleaned_first_strand","first-strand cDNA after free-primer removal",
                          MaterialKind.MOLECULAR,state=first_state)
    primer_waste=graph.product("primer_waste","digested unused UP1",MaterialKind.MOLECULAR)
    graph.react("remove_up1","Digest unused UP1",StepKind.DIGESTION,
                ((first,InputRole.TEMPLATE),(exo,InputRole.CATALYST)),(cleaned,primer_waste))
    graph.mark_terminal(Terminal(primer_waste.id,Disposition.UNRESOLVED,"UP1 digestion does not establish physical removal; the fragment fate is unassigned."))
    tdt=graph.input("tdt","terminal deoxynucleotidyl transferase",MaterialKind.CATALYST)
    datp=graph.input("datp","dATP tailing substrate",MaterialKind.CONDITION)
    tailed_state=Construct([*first_state.segments,Segment("TdT poly(A)","A"*12)],
                           name="poly(A)-tailed first-strand cDNA")
    tailed=graph.product("tailed_first_strand",tailed_state.name,MaterialKind.MOLECULAR,
                         state=tailed_state)
    graph.react("tail","Add the second poly(A) tail with TdT",StepKind.TAILING,
                ((cleaned,InputRole.TEMPLATE),(datp,InputRole.CONDITION),
                 (tdt,InputRole.CATALYST)),(tailed,))
    up2=graph.input("up2","UP2 oligo-dT second-strand primer",MaterialKind.OLIGO,
                    state=(Segment("UP2",UP2),))
    polymerase=graph.input("second_strand_polymerase","second-strand polymerase",
                           MaterialKind.CATALYST)
    duplex=graph.product("double_stranded_cdna",amplified_cdna().name,MaterialKind.MOLECULAR,
                         state=amplified_cdna())
    graph.react("second_strand","Synthesize the second strand from UP2",StepKind.EXTENSION,
                ((tailed,InputRole.TEMPLATE),(up2,InputRole.STRUCTURAL),
                 (polymerase,InputRole.CATALYST)),(duplex,))
    pcr_pol=graph.input("pcr_polymerase","cDNA PCR polymerase",MaterialKind.CATALYST)
    amplified=graph.product("amplified_cdna","UP1/UP2-amplified cDNA",MaterialKind.MOLECULAR,
                            state=amplified_cdna())
    graph.react("cdna_pcr","Amplify cDNA with UP1 and UP2 (20 cycles)",StepKind.PCR,
                ((duplex,InputRole.TEMPLATE),(up1,InputRole.STRUCTURAL),
                 (up2,InputRole.STRUCTURAL),(pcr_pol,InputRole.CATALYST)),(amplified,),
                geometry=PCRGeometry.unavailable(
                    amplified_cdna(),reason="The source prints UP1/UP2, but the representative insert does not assert transcript-specific primer landing coordinates."))
    aup1=graph.input("aup1","5-prime-amine-blocked AUP1",MaterialKind.OLIGO,
                     state=(Segment("AUP1",AUP1),))
    aup2=graph.input("aup2","5-prime-amine-blocked AUP2",MaterialKind.OLIGO,
                     state=(Segment("AUP2",AUP2),))
    blocked=graph.product("amine_blocked_cdna","5-prime-amine-blocked amplified cDNA",
                          MaterialKind.MOLECULAR,state=amplified_cdna())
    graph.react("blocked_pcr","Re-amplify with amine-blocked AUP1/AUP2 (9 cycles)",StepKind.PCR,
                ((amplified,InputRole.TEMPLATE),(aup1,InputRole.STRUCTURAL),
                 (aup2,InputRole.STRUCTURAL),(pcr_pol,InputRole.CATALYST)),(blocked,),
                geometry=PCRGeometry.unavailable(
                    amplified_cdna(),reason="The source prints AUP1/AUP2, but the representative insert does not assert transcript-specific primer landing coordinates."))
    shear=graph.input("shearing","mechanical shearing condition",MaterialKind.CONDITION)
    fragments=graph.product("fragments","80–130 bp cDNA fragments",MaterialKind.MOLECULAR)
    graph.react("fragment","Shear amplified cDNA",StepKind.FRAGMENTATION,
                ((blocked,InputRole.TEMPLATE),(shear,InputRole.CONDITION)),(fragments,))
    repair=graph.input("end_repair_mix","end-repair enzyme mix",MaterialKind.CATALYST)
    repaired=graph.product("repaired_fragments","end-repaired cDNA fragments",
                           MaterialKind.MOLECULAR)
    graph.react("end_repair","End-repair the sheared cDNA",StepKind.END_REPAIR,
                ((fragments,InputRole.TEMPLATE),(repair,InputRole.CATALYST)),(repaired,))
    adapters=graph.input("solid_adapters","unbarcoded SOLiD P1/P2 adapter duplexes",
                         MaterialKind.ADAPTER,
                         state=(P1_TOP,P1_BOTTOM,P2_TOP,P2_BOTTOM),
                         source=source+"; 2009 Supplementary Data adapter table")
    ligase=graph.input("ligase","DNA ligase",MaterialKind.CATALYST)
    ligated=graph.product("ligated_library","P1/P2-ligated cDNA fragments",
                          MaterialKind.MOLECULAR,state=final_library())
    graph.react("adapter_ligation","Ligate unbarcoded SOLiD P1/P2 adapters",
                StepKind.LIGATION,
                ((repaired,InputRole.TEMPLATE),(adapters,InputRole.STRUCTURAL),
                 (ligase,InputRole.CATALYST)),(ligated,),geometry=solid_adapter_geometry())
    gel=graph.input("gel_selection","6% polyacrylamide gel, 150–200 bp fraction",MaterialKind.SELECTION)
    selected=graph.product("selected_ligated_library","150–200 bp adapter-ligated fraction",
                           MaterialKind.MOLECULAR,state=final_library())
    rejected=graph.product("offsize_library","excluded gel fractions",MaterialKind.WASTE)
    graph.react("size_select","Select 150–200 bp adapter-ligated DNA",StepKind.SELECTION,
                ((ligated,InputRole.TEMPLATE),(gel,InputRole.SELECTION)),(selected,rejected))
    graph.discard(rejected,"outside the source-selected gel fraction")
    lib1=graph.input("library_primer_1","SOLiD library PCR primer 1",MaterialKind.OLIGO,
                     state=(Segment("Library PCR Primer 1",LIBRARY_PCR_1),))
    lib2=graph.input("library_primer_2","SOLiD library PCR primer 2",MaterialKind.OLIGO,
                     state=(Segment("Library PCR Primer 2",LIBRARY_PCR_2),))
    final=final_library()
    library=graph.product("library",final.name,MaterialKind.LIBRARY,state=final)
    graph.react("library_pcr","Amplify the adapter-ligated SOLiD library (8–10 cycles)",StepKind.PCR,
                ((selected,InputRole.TEMPLATE),(lib1,InputRole.STRUCTURAL),
                 (lib2,InputRole.STRUCTURAL),(pcr_pol,InputRole.CATALYST)),(library,),
                geometry=PCRGeometry.exact(
                    final,PCRPrimer("library PCR 1",LIBRARY_PCR_1),
                    PCRPrimer("library PCR 2",LIBRARY_PCR_2),
                    statement="The 2009 supplementary table prints both SOLiD library PCR primers and adapter sites."))
    graph.sequence(library,SequencingPlan(
        acquired_roles=("Read 1",),read_uses={"Read 1":"include"},
        limitation="Online Methods specifies 50-base SOLiD output; supplements also label35mer runs. Loaded sequencing-primer identity is outside this model evidence scope."))
    graph.declare_reviewed_scope((
        "single_cell","up1","reverse_transcriptase","reverse_transcribe","exonuclease",
        "remove_up1","tdt","datp","tail","up2","second_strand_polymerase",
        "second_strand","pcr_polymerase","cdna_pcr","aup1","aup2","blocked_pcr",
        "shearing","fragment","end_repair_mix","end_repair","solid_adapters","ligase",
        "adapter_ligation","gel_selection","selected_ligated_library","size_select","library_primer_1","library_primer_2","library_pcr","library"),source)
    return graph.require_publishable()

