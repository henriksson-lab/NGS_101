#!/usr/bin/env python3
"""Build the lentiCRISPR v2 screening page."""
from pathlib import Path
import sys
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "lib"))
import pooled_crispr_page as P
import pooled_crispr as M
from adapter_ligation import AdapterLigationGeometry
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, StepKind
from protocol_model import molecule, role

OUT = HERE.parent / "lenticrispr-v2-screening.html"
SOURCE = 'Defining source: <a href="https://doi.org/10.1038/nmeth.3047">Sanjana, Shalem &amp; Zhang, <i>Nature Methods</i> (2014)</a>. The final-library example uses the separately cited Joung one-PCR readout.'


def protocol_graph() -> ProtocolGraph:
    source = SOURCE
    graph = ProtocolGraph("lentiCRISPR v2 screening", authority=source)
    vector = graph.input("vector", "lentiCRISPR v2 transfer vector", MaterialKind.MOLECULAR,
                         state=molecule("lentiCRISPR v2 plasmid", (role("BsmBI cloning cassette"),)),
                         primary=True)
    guide = graph.input("guide_duplex", "annealed BsmBI-compatible guide oligos", MaterialKind.OLIGO)
    ligase = graph.input("cloning_ligase", "guide-cloning DNA ligase", MaterialKind.CATALYST)
    cloned = graph.product("cloned_vector", "guide-bearing transfer vector", MaterialKind.MOLECULAR,
                           state=M.cloned_guide_cassette())
    graph.react("clone_guide", "Ligate the guide duplex into the BsmBI-opened cassette",
                StepKind.LIGATION,
                ((vector, InputRole.TEMPLATE), (guide, InputRole.STRUCTURAL),
                 (ligase, InputRole.CATALYST)), (cloned,),
                geometry=AdapterLigationGeometry.unavailable(
                    "The two BsmBI-directed seams are checked in cloning_scene; a single adapter-end object does not represent a two-ended plasmid insertion."))
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
    reviewed = ["vector", "guide_duplex", "cloning_ligase", "clone_guide", "packaging_system",
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
{P.cloning()}
{P.one_pcr_readout("Compatible one-PCR genomic readout")}
<p>The Moffat LCV2::TKOv3 two-step indexed readout is shown on its own page; it is not an extra step in the Joung one-PCR protocol above.</p>'''
    return P.page("lentiCRISPR v2 screening", "01_lenticrispr-v2.html", source, body,
                  graph=protocol_graph())


def main() -> None:
    OUT.write_text(render(), encoding="utf-8")
    print(f"wrote {OUT}  ({OUT.stat().st_size:,} bytes)")


if __name__ == "__main__": main()
