#!/usr/bin/env python3
"""Build the lentiCRISPR v1 screening page."""
from pathlib import Path
import sys
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "lib"))
import pooled_crispr_page as P
import pooled_crispr as M
from adapter_ligation import AdapterLigationGeometry
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, StepKind
from protocol_model import molecule, role

OUT = HERE.parent / "lenticrispr-v1-screening.html"
SOURCE = 'Defining source: <a href="https://doi.org/10.1126/science.1247005">Shalem et al., <i>Science</i> (2014)</a>. The final-library example uses the separately cited Joung one-PCR readout.'


def protocol_graph() -> ProtocolGraph:
    source = SOURCE
    graph = ProtocolGraph("lentiCRISPR v1 screening", authority=source)
    vector = graph.input("vector", "lentiCRISPR v1 transfer vector", MaterialKind.MOLECULAR,
                         state=molecule("lentiCRISPR v1 plasmid", (role("BsmBI cloning cassette"),)),
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
    # The Joung library is an explicitly cited compatible example, not Shalem's nested readout.
    M._pcr(graph, "readout", gdna, M.one_pcr_library(), M.one_pcr_sequencing_primers(),
           "Compatible one-PCR guide-cassette readout", reviewed)
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()


def render() -> str:
    source = SOURCE
    body = f'''<h2>The screen the vector was built for</h2>
<p>Shalem et al. deliver the GeCKO library — 64,751 unique guide sequences against 18,080
genes, 3–4 guides per gene — as lentiCRISPR particles at a multiplicity of infection of
0.3, so that most transduced cells carry a single integrated cassette, and read the guide
distribution out of genomic DNA 3 and 14 days after transduction (A375 melanoma and
HUES62 cells; the vemurafenib arm adds 14 days of drug). Those are the only construction
parameters the accessible article text states.</p>
<h2>Transfer vector</h2>
{P.map_panel("v1")}
<p>The guide cassette is upstream of cPPT/CTS in v1. That order matters when a reverse readout primer is placed in cPPT.</p>
{P.cloning()}
{P.one_pcr_readout("Compatible one-PCR genomic readout")}
<p>The readout the original screen ran on genomic DNA is described only in the paper's
Supplementary Methods, which are not reachable here; its details are kept in the research
notes and are not blended into the compatible Joung one-PCR example above.</p>
<p class="caveat"><b>Source availability.</b> Everything this paper publishes about how the
screen libraries were amplified and sequenced sits in its Supplementary Information. The
only full text reachable here is the PubMed Central author manuscript (PMC4089965), which
carries the Results but no Methods and no supplement; the PMC supplementary PDF is behind
a download gate, Europe PMC rejects the record as not open access, and the live
<code>genome-engineering.org/gecko/</code> host still does not respond. The Shalem
Methods' own readout primer pair, its barcodes and the acquired read lengths therefore
remain unverified here, rather than established as unpublished.</p>
<p class="caveat"><b>The author-hosted GeCKO primer workbook (recovered 2026-10-17).</b>
The Zhang-lab workbook <code>GeCKO-plasmid-readout-primers-July2014.xlsx</code> was
recovered from the Internet Archive copy of that page and read. It states that
lentiCRISPR&nbsp;v1 (and lentiGuide-Puro) are amplified in a <b>single PCR</b> with its
forward and reverse Illumina primers — a mix of forward primers supplying the stagger —
while lentiCRISPRv2 needs an extra adaptor PCR first. Its twelve forward primers carry
P5, TruSeq Read&nbsp;1, a 1–9&nbsp;nt stagger, an <b>8-nt inline barcode</b> and the hU6
3′ priming site used by the drawn forward primer; its reverse primers carry P7, a second
8-nt inline barcode, TruSeq Read&nbsp;2, a stagger and a 25-nt site downstream of the
guide scaffold (in cPPT, which sits after the cassette in v1 only). Sample identity is
therefore read inside Read&nbsp;1/Read&nbsp;2 rather than as the i7 index of the library
drawn above. Two limits: the workbook is headed <i>plasmid</i> readout and refers
genomic-DNA readout back to the gated Supplementary Methods, and it prints no read
lengths or cycling conditions — so it is a v1-compatible primer set, not a transcript of
the 2014 screen's own protocol.</p>'''
    return P.page("lentiCRISPR v1 screening", "01_lenticrispr-v1.html", source, body,
                  graph=protocol_graph())


def main() -> None:
    OUT.write_text(render(), encoding="utf-8")
    print(f"wrote {OUT}  ({OUT.stat().st_size:,} bytes)")


if __name__ == "__main__": main()
