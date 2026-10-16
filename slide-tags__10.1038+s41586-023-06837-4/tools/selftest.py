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
c("both assays have separate RNA and spatial-tag readouts",len(g.sequencing),5)
c("3′ RNA bypasses Multiome tagmentation","tagment" in ancestors(g,"rna_3p_library"),False)
c("3′ tag bypasses Multiome tagmentation","tagment" in ancestors(g,"tag_3p_library"),False)
c("Multiome ATAC follows tagmentation","tagment" in ancestors(g,"atac_library"),True)
c.report()
