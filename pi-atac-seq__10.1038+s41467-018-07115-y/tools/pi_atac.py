"""Molecular model for Pi-ATAC-seq; unmapped primer-pair choices infer by default."""
from __future__ import annotations
import illumina as il
import nextera as nx
import seqprimers as sp
from chemdraw import Construct,Scene,Segment,complement_segments,feature
from pcr import PCRGeometry
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind)

def seg(name,top,tag=None,**kw): return Segment(name=name,top=top,tag=tag,**kw)
def inferred(name,top,tag=None,**kw): return Segment(name=name,top=top,tag=tag,inferred=True,**kw)

def tagmented_scene():
 gap=nx.TAGMENTATION_GAP; core=seg("genomic core","X"*24,placeholder=True)
 top=[inferred("s5",nx.S5,"s5"),inferred("left ME",nx.ME,"me"),seg("left gap","X"*gap,placeholder=True),core]
 bot=[inferred("s7",nx.S7,"s7"),inferred("right ME",nx.ME,"me"),seg("right gap","X"*gap,placeholder=True),*complement_segments([core])]
 sc=Scene(); sc.strand("top",top,label="fixed nuclear DNA"); sc.anneal("bottom",bot,to="top",pair=("genomic core'","genomic core"),label="fixed nuclear DNA")
 sc.anneal("left ME bottom",[inferred("left ME'",nx.ME_RC,"me")],to="top",pair=("left ME'","left ME"),label="",mod5="p")
 sc.anneal("right ME bottom",[inferred("right ME'",nx.ME_RC,"me")],to="bottom",pair=("right ME'","right ME"),label="",mod5="p",above=True); return sc

SEQ_PRIMERS=tuple(sp.NEXTERA[k] for k in ("R1","I1","I2","R2"))
def final_library(insert_nt=34):
 lib=Construct([inferred("P5",il.P5,"p5"),inferred("i5","I"*8,"cbc",placeholder=True,feature=feature("sample_i5","sample_index","unknown")),inferred("s5",nx.S5,"s5"),inferred("ME",nx.ME,"me"),
 seg("accessible genomic DNA","X"*insert_nt,placeholder=True),inferred("ME reverse complement",nx.ME_RC,"me"),inferred("s7 reverse complement",nx.S7_RC,"s7"),
 inferred("i7 reverse complement","J"*8,"cbc",placeholder=True,feature=feature("sample_i7","sample_index","unknown")),inferred("P7 reverse complement",il.P7_RC,"p7")],name="inferred Pi-ATAC-seq library")
 problems=sp.verify(lib,SEQ_PRIMERS)
 if problems: raise ValueError("invalid inferred Pi-ATAC-seq library: "+"; ".join(problems))
 return lib
def primer_landings():
 lib=final_library(); return [(p,sp.locate(lib,p)) for p in SEQ_PRIMERS]

def protocol_graph():
 source="Defining source: https://doi.org/10.1038/s41467-018-07115-y"
 inferred_source=source+" (Pi-ATAC well-to-primer assignment unavailable; linked Buenrostro Ad1/Ad2 family does not identify the chosen pairs)"
 graph=ProtocolGraph("Pi-ATAC-seq",authority=source)
 cells=graph.input("cells","fixed cells",MaterialKind.BIOLOGICAL,primary=True)
 antibodies=graph.input("antibodies","fluorescent protein-detection antibodies",MaterialKind.MOLECULAR)
 stained=graph.product("stained_cells","antibody-stained fixed cells",MaterialKind.BIOLOGICAL)
 unbound=graph.product("unbound_antibody","unbound antibody wash fraction",MaterialKind.WASTE)
 graph.react("antibody_stain","Bind fluorescent antibodies to cellular proteins",StepKind.BINDING,
             ((cells,InputRole.TEMPLATE),(antibodies,InputRole.ANALYTE)),(stained,unbound))
 graph.discard(unbound,"unbound antibodies are removed before cytometric sorting")
 transposome=graph.input("transposome","Nextera-compatible Tn5 transposome",MaterialKind.ADAPTER,
                         inferred=True,source=inferred_source)
 tn5=graph.input("tn5","Tn5 transposase",MaterialKind.CATALYST)
 tagged=graph.product("tagmented_cells","fixed cells containing tagmented chromatin",
                      MaterialKind.MOLECULAR,state=tagmented_scene(),inferred=True,source=inferred_source)
 twaste=graph.product("tagmentation_waste","non-library tagmentation products",MaterialKind.WASTE)
 graph.react("tagment","Tagment antibody-stained cells in bulk",StepKind.TAGMENTATION,
             ((stained,InputRole.TEMPLATE),(transposome,InputRole.STRUCTURAL),
              (tn5,InputRole.CATALYST)),(tagged,twaste))
 graph.discard(twaste,"fragments lacking usable paired adapter ends are not amplified")
 sorter=graph.input("index_sort","one-cell-per-indexed-well FACS rule",MaterialKind.SELECTION)
 wells=graph.product("single_cell_wells","one protein-profiled cell per indexed well",MaterialKind.MOLECULAR)
 excluded=graph.product("excluded_events","empty wells and multiplet events",MaterialKind.WASTE)
 graph.react("sort","Record protein fluorescence and index-sort one cell per well",StepKind.SELECTION,
             ((tagged,InputRole.TEMPLATE),(sorter,InputRole.SELECTION)),(wells,excluded))
 graph.discard(excluded,"empty and multiplet sorting events are excluded")
 proteinase=graph.input("proteinase_k","proteinase K reverse-crosslink reagent",MaterialKind.CATALYST)
 released=graph.product("released_fragments","reverse-crosslinked tagmented fragments",MaterialKind.MOLECULAR)
 graph.react("reverse_crosslinks","Reverse crosslinks and release Tn5",StepKind.DEPROTEINIZATION,
             ((wells,InputRole.TEMPLATE),(proteinase,InputRole.CATALYST)),(released,))
 gap_pol=graph.input("gap_fill_polymerase","gap-fill and PCR polymerase",MaterialKind.CATALYST)
 filled=graph.product("gap_filled","gap-filled tagmented fragments",MaterialKind.MOLECULAR)
 graph.react("gap_fill","Fill transposition gaps",StepKind.EXTENSION,
             ((released,InputRole.TEMPLATE),(gap_pol,InputRole.CATALYST)),(filled,))
 i5=graph.input("i5_primer","well-specific i5 barcoding primer",MaterialKind.OLIGO,
                inferred=True,source=inferred_source)
 i7=graph.input("i7_primer","well-specific i7 barcoding primer",MaterialKind.OLIGO,
                inferred=True,source=inferred_source)
 final=final_library()
 library=graph.product("library",final.name,MaterialKind.LIBRARY,state=final,
                       inferred=True,source=inferred_source)
 graph.react("indexed_pcr","Amplify with the well-specific barcoding-primer pair",StepKind.PCR,
             ((filled,InputRole.TEMPLATE),(i5,InputRole.STRUCTURAL),
              (i7,InputRole.STRUCTURAL),(gap_pol,InputRole.CATALYST)),(library,),
             geometry=PCRGeometry.unavailable(final,
                 reason="The paper reports 96 x 90 primer combinations but not their well-to-primer assignments; the cited predecessor prints related Ad1/Ad2 primer families."))
 graph.sequence(library,SequencingPlan(
     library=final,primers=SEQ_PRIMERS,acquired_roles=tuple(p.role for p in SEQ_PRIMERS),
     read_lengths={"Read 1":75,"Index 1 (i7)":8,"Index 2 (i5)":8,"Read 2":75},
     read_uses={"Read 1":"include","Index 1 (i7)":"demultiplex",
                "Index 2 (i5)":"demultiplex","Read 2":"include"},
     limitation="Pi-ATAC well-to-primer assignments are not specified; the cited predecessor prints related Ad1/Ad2 primer families."))
 graph.declare_reviewed_scope((
     "cells","antibodies","antibody_stain","transposome","tn5","tagment",
     "index_sort","sort","proteinase_k","reverse_crosslinks","gap_fill_polymerase",
     "gap_fill","i5_primer","i7_primer","indexed_pcr","library"),source)
 return graph.require_publishable()
