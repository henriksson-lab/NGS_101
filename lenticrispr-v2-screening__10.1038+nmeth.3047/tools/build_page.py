#!/usr/bin/env python3
"""Build the lentiCRISPR v2 screening page."""
from pathlib import Path
import sys
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "lib"))
import pooled_crispr_page as P
import pooled_crispr as M
from chemdraw import Construct, Segment, oligo
from pcr import PCRGeometry
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, StepKind
from protocol_model import molecule, role

OUT = HERE.parent / "lenticrispr-v2-screening.html"
SOURCE = 'Defining source: <a href="https://doi.org/10.1038/nmeth.3047">Sanjana, Shalem &amp; Zhang, <i>Nature Methods</i> (2014)</a>. The final-library example uses the separately cited Joung one-PCR readout.'


# Sanjana 2014 Supplementary Methods, “GeCKO library pooled synthesis and cloning”.
ARRAY_F = "TAACTTGAAAGTATTTCGATTTCTTGGCTTTATATATCTTGTGGAAAGGACGAAACACCG"
ARRAY_R = "ACTTTTTCAAGTTGATAACGGACTAGCCTTATTTTAACTTGCTATTTCTAGCTCTAAAAC"


def protocol_graph() -> ProtocolGraph:
    source = SOURCE
    graph = ProtocolGraph("lentiCRISPR v2 screening", authority=source)
    vector = graph.input("vector", "lentiCRISPR v2 transfer vector", MaterialKind.MOLECULAR,
                         state=molecule("lentiCRISPR v2 plasmid", (role("BsmBI cloning cassette"),)),
                         primary=True)
    pool_state = molecule("74-nt synthesized GeCKO v2 oligo pool",
                          (role("array oligo containing variable guide and shared flanks"),))
    pool = graph.input("array_pool", "74-nt synthesized guide-oligo pool",
                       MaterialKind.OLIGO, state=pool_state, primary=True)
    forward = graph.input("array_f", "ArrayF pool-construction primer", MaterialKind.OLIGO,
                          state=Construct([Segment("ArrayF", ARRAY_F)], name="ArrayF"))
    reverse = graph.input("array_r", "ArrayR pool-construction primer", MaterialKind.OLIGO,
                          state=Construct([Segment("ArrayR", ARRAY_R)], name="ArrayR"))
    polymerase = graph.input("array_polymerase", "pool-amplification polymerase",
                             MaterialKind.CATALYST)
    insert = graph.product("array_insert", "140-bp amplified pool inserts", MaterialKind.MOLECULAR,
                           state=molecule("guide-pool insert", (role("vector-homology arm"),
                                          role("variable guide"), role("opposite vector-homology arm"))))
    graph.react("amplify_pool", "Amplify the synthesized guide pool with ArrayF and ArrayR",
                StepKind.PCR,
                ((pool, InputRole.TEMPLATE), (forward, InputRole.STRUCTURAL),
                 (reverse, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)),
                (insert,), geometry=PCRGeometry.unavailable(Construct([role("74-nt synthesis-pool member")], name="GeCKO synthesis pool"),
                    reason="Both primers are printed in the original supplement; this model has not transcribed an individual 74-nt synthesis-pool member and therefore does not claim exact template coordinates."))
    nuclease = graph.input("bsmbi", "BsmBI restriction enzyme", MaterialKind.CATALYST)
    opened = graph.product("opened_vector", "BsmBI-opened lentiCRISPR v2 backbone",
                           MaterialKind.MOLECULAR, state=vector.state)
    stuffer = graph.product("excised_stuffer", "excised cloning cassette", MaterialKind.WASTE)
    graph.react("open_vector", "Digest the transfer vector with BsmBI", StepKind.DIGESTION,
                ((vector, InputRole.TEMPLATE), (nuclease, InputRole.CATALYST)), (opened, stuffer))
    graph.discard(stuffer, "excluded when the vector backbone is gel-purified")
    phosphatase = graph.input("fastap", "FastAP phosphatase", MaterialKind.CATALYST)
    dephosphorylated = graph.product("dephosphorylated_vector", "dephosphorylated vector backbone",
                                     MaterialKind.MOLECULAR, state=opened.state)
    graph.react("dephosphorylate", "Dephosphorylate the opened backbone with FastAP",
                StepKind.DEPHOSPHORYLATION,
                ((opened, InputRole.TEMPLATE), (phosphatase, InputRole.CATALYST)), (dephosphorylated,))
    assembly = graph.input("gibson_mix", "Gibson assembly enzyme mixture", MaterialKind.CATALYST)
    cloned = graph.product("cloned_vector", "guide-bearing pooled transfer vectors", MaterialKind.MOLECULAR,
                           state=M.cloned_guide_cassette())
    graph.react("clone_guide", "Gibson-assemble the amplified pool into the opened vector",
                StepKind.MAGIC,
                ((dephosphorylated, InputRole.TEMPLATE), (insert, InputRole.STRUCTURAL),
                 (assembly, InputRole.CATALYST)), (cloned,),
                reason="The defining pooled-library construction uses Gibson assembly; the exact synthesized-member junctions are not yet reconstructed here.")
    packaging = graph.input("packaging_system", "lentiviral packaging system", MaterialKind.CONDITION)
    virus = graph.product("virus", "guide-library lentiviral particles", MaterialKind.MOLECULAR,
                          state=cloned.state)
    graph.react("package", "Package the pooled transfer-vector library", StepKind.GENERIC,
                ((cloned, InputRole.TEMPLATE), (packaging, InputRole.CONDITION)),
                (virus,), changes_molecule=False)
    cells = graph.input("cells", "screen target-cell population", MaterialKind.BIOLOGICAL,
                        state=molecule("screen cells", (role("cellular genomes"),)), primary=True)
    integration = graph.input("integration_catalyst", "lentiviral integration machinery",
                              MaterialKind.CATALYST)
    screened = graph.product("screened_cells", "selected cells carrying integrated guides",
                             MaterialKind.MOLECULAR, state=cells.state)
    graph.react("integrate", "Transduce cells and integrate one guide cassette", StepKind.INTEGRATION,
                ((cells, InputRole.TEMPLATE), (virus, InputRole.STRUCTURAL),
                 (integration, InputRole.CATALYST)), (screened,))
    extraction = graph.input("extraction_condition", "genomic-DNA extraction condition",
                             MaterialKind.CONDITION)
    gdna = graph.product("genomic_dna", "screen genomic DNA containing guide cassettes",
                         MaterialKind.MOLECULAR, state=screened.state)
    graph.react("extract", "Extract genomic DNA after pooled selection", StepKind.EXTRACTION,
                ((screened, InputRole.TEMPLATE), (extraction, InputRole.CONDITION)), (gdna,))
    reviewed = ["vector", "array_pool", "array_f", "array_r", "array_polymerase",
                "amplify_pool", "bsmbi", "open_vector", "fastap", "dephosphorylate",
                "gibson_mix", "clone_guide", "packaging_system",
                "package", "cells", "integration_catalyst", "integrate", "extraction_condition", "extract"]
    # The Joung library is an explicitly cited compatible example, not Sanjana et al.'s original readout.
    M._pcr(graph, "readout", gdna, M.one_pcr_library(), M.one_pcr_sequencing_primers(),
           "Compatible one-PCR guide-cassette readout", reviewed)
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()


def render() -> str:
    source = SOURCE
    body = f'''<h2>Transfer vector</h2>
{P.map_panel("v2")}
<p>v2 moves cPPT/CTS upstream of the hU6 guide cassette. A cPPT reverse primer that gives a short v1 product therefore must not be assumed to work on v2.</p>
<h2>Pooled guide-library construction</h2>
<p>The defining GeCKO v2 library starts with a synthesized 74-nt oligo pool. ArrayF/ArrayR PCR makes 140-bp inserts; BsmBI digestion and FastAP treatment prepare the vector backbone for Gibson assembly (Sanjana 2014 Supplementary Methods, “GeCKO library pooled synthesis and cloning”). The graph shows this pooled construction.</p>
<seq>{oligo("ArrayF", [Segment("ArrayF", ARRAY_F)])}
{oligo("ArrayR", [Segment("ArrayR", ARRAY_R)])}</seq>
<p>The original supplement describes a separate Illumina-adapter PCR for validation of pool representation. Its historical primer and barcode choices remain under source review. The final library below is the explicitly cited later Joung genomic readout.</p>
{P.one_pcr_readout("Compatible one-PCR genomic readout")}
<p>The Moffat LCV2::TKOv3 two-step indexed readout is shown on its own page; it is not an extra step in the Joung one-PCR protocol above.</p>'''
    return P.page("lentiCRISPR v2 screening", "01_lenticrispr-v2.html", source, body,
                  graph=protocol_graph())


def main() -> None:
    OUT.write_text(render(), encoding="utf-8")
    print(f"wrote {OUT}  ({OUT.stat().st_size:,} bytes)")


if __name__ == "__main__": main()
