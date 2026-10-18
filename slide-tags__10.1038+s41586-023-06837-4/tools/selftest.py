#!/usr/bin/env python3
from pathlib import Path
import sys
H=Path(__file__).resolve().parent;sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import slide_tags as M
from checks import Check


def ancestors(graph, material_id):
    """Step ids that can contribute to one terminal library."""
    seen_materials=set()
    seen_steps=set()
    pending=[material_id]
    while pending:
        material=pending.pop()
        if material in seen_materials:
            continue
        seen_materials.add(material)
        producer=graph.producer_of(material)
        if producer is not None:
            seen_steps.add(producer.id)
            pending.extend(port.material for port in producer.inputs)
    return seen_steps


g=M.protocol_graph()
c=Check()
c("TAGS3′,SLACbrain3′ andMultiome retain source-distinct outputs",len(g.sequencing),7)
c("3′ RNA bypasses Multiome tagmentation","multi_transpose" in ancestors(g,"tags_rna_library"),False)
c("3′ tag bypasses Multiome tagmentation","multi_transpose" in ancestors(g,"tags_tag_library"),False)
c("Multiome ATAC follows tagmentation","multi_transpose" in ancestors(g,"multi_atac_library"),True)
c("3′spatialtagPCR follows actual cDNA/SPRI stages",{"tags_cdna_pcr","tags_spri_split"}.issubset(ancestors(g,"tags_tag_library")),True)
c("Multiome spatialtagPCR follows mixed preamp/SPRI without GEXindexPCR",("multi_preamp_spri" in ancestors(g,"multi_tag_library"),"multi_gex_index_pcr" in ancestors(g,"multi_tag_library")),(True,False))
c.report()
