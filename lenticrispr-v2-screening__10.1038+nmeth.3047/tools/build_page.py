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
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind)
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
    polymerase = graph.input("array_polymerase", "Phusion HS Flex polymerase",
                             MaterialKind.CATALYST)
    cycling = graph.input("array_cycling", "63 °C annealing, 15 s extension, 20 cycles",
                          MaterialKind.CONDITION)
    amplicon = graph.product("array_amplicon", "pool PCR product", MaterialKind.MOLECULAR,
                             state=molecule("guide-pool insert", (role("vector-homology arm"),
                                            role("variable guide"), role("opposite vector-homology arm"))))
    graph.react("amplify_pool", "Amplify the synthesized guide pool with ArrayF and ArrayR",
                StepKind.PCR,
                ((pool, InputRole.TEMPLATE), (forward, InputRole.STRUCTURAL),
                 (reverse, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST),
                 (cycling, InputRole.CONDITION)),
                (amplicon,), geometry=PCRGeometry.unavailable(Construct([role("74-nt synthesis-pool member")], name="GeCKO synthesis pool"),
                    reason="Both primers are printed in the original supplement; this model has not transcribed an individual 74-nt synthesis-pool member and therefore does not claim exact template coordinates."))
    gel = graph.input("insert_gel", "2% agarose E-Gel EX size selection", MaterialKind.CONDITION)
    insert = graph.product("array_insert", "140-bp size-selected pool inserts", MaterialKind.MOLECULAR,
                           state=amplicon.state)
    offsize = graph.product("array_offsize", "off-size pool PCR products", MaterialKind.WASTE)
    graph.react("size_select_insert", "Size-select the 140-bp amplicon on a 2% E-Gel EX",
                StepKind.CLEANUP,
                ((amplicon, InputRole.TEMPLATE), (gel, InputRole.CONDITION)),
                (insert, offsize), changes_molecule=False)
    graph.discard(offsize, "cut away from the excised 140-bp band")
    nuclease = graph.input("bsmbi", "BsmBI restriction enzyme (Fermentas)", MaterialKind.CATALYST)
    opened = graph.product("opened_vector", "BsmBI-opened lentiCRISPR v2 backbone",
                           MaterialKind.MOLECULAR, state=vector.state)
    stuffer = graph.product("excised_stuffer", "excised 1.8-kb filler", MaterialKind.WASTE)
    graph.react("open_vector", "Digest the transfer vector with BsmBI", StepKind.DIGESTION,
                ((vector, InputRole.TEMPLATE), (nuclease, InputRole.CATALYST)), (opened, stuffer))
    graph.discard(stuffer, "removed when the backbone is gel-purified on a 1% E-Gel EX")
    phosphatase = graph.input("fastap", "FastAP phosphatase, 37 °C for 2 h", MaterialKind.CATALYST)
    dephosphorylated = graph.product("dephosphorylated_vector", "dephosphorylated vector backbone",
                                     MaterialKind.MOLECULAR, state=opened.state)
    graph.react("dephosphorylate", "Dephosphorylate the opened backbone with FastAP",
                StepKind.DEPHOSPHORYLATION,
                ((opened, InputRole.TEMPLATE), (phosphatase, InputRole.CATALYST)), (dephosphorylated,))
    assembly = graph.input("gibson_mix", "Gibson assembly enzyme mixture", MaterialKind.CATALYST)
    cloned = graph.product("cloned_vector", "guide-bearing pooled transfer vectors", MaterialKind.MOLECULAR,
                           state=M.cloned_guide_cassette())
    graph.react("clone_guide", "Gibson-assemble 10 ng of insert into 25 ng of opened vector (20 µl, NEB)",
                StepKind.MAGIC,
                ((dephosphorylated, InputRole.TEMPLATE), (insert, InputRole.STRUCTURAL),
                 (assembly, InputRole.CATALYST)), (cloned,),
                reason="The defining pooled-library construction uses Gibson assembly; the exact synthesized-member junctions are not yet reconstructed here.")
    bacteria = graph.input("electrocompetent_cells",
                           "electrocompetent cells (Lucigen) with 50 µg/ml carbenicillin plates",
                           MaterialKind.SELECTION)
    colonies = graph.product("transformed_colonies",
                             "carbenicillin-resistant colonies at 30–200× library coverage",
                             MaterialKind.MOLECULAR, state=cloned.state)
    untransformed = graph.product("untransformed_cells", "cells without an assembled plasmid",
                                  MaterialKind.WASTE)
    graph.react("transform", "Electroporate parallel transformations and select on carbenicillin",
                StepKind.SELECTION,
                ((cloned, InputRole.TEMPLATE), (bacteria, InputRole.SELECTION)),
                (colonies, untransformed), changes_molecule=False)
    graph.discard(untransformed, "not recovered from the carbenicillin plates")
    maxiprep = graph.input("maxiprep", "Endotoxin-Free Plasmid Maxiprep (Qiagen)",
                           MaterialKind.CONDITION)
    plasmid_library = graph.product("plasmid_library", "maxiprepped pooled plasmid library",
                                    MaterialKind.MOLECULAR, state=cloned.state)
    graph.react("maxiprep_library", "Scrape and combine the colonies, then extract plasmid DNA",
                StepKind.EXTRACTION,
                ((colonies, InputRole.TEMPLATE), (maxiprep, InputRole.CONDITION)),
                (plasmid_library,), changes_molecule=False)
    qc_primers = graph.input("qc_adapter_primers",
                             "barcoded Illumina-adapter primers (sequences not printed)",
                             MaterialKind.OLIGO)
    qc_polymerase = graph.input("qc_polymerase", "Phusion Flash High-Fidelity polymerase, 16 cycles",
                                MaterialKind.CATALYST)
    qc_amplicon = graph.product("qc_amplicon",
                                "duplicate barcoded library-representation amplicons",
                                MaterialKind.LIBRARY)
    graph.react("validate_library",
                "Amplify the plasmid library in duplicate for representation sequencing",
                StepKind.PCR,
                ((plasmid_library, InputRole.TEMPLATE), (qc_primers, InputRole.STRUCTURAL),
                 (qc_polymerase, InputRole.CATALYST)),
                (qc_amplicon,),
                geometry=PCRGeometry.unavailable(
                    Construct([role("cloned guide cassette in the pooled plasmid")],
                              name="pooled plasmid library"),
                    reason="The supplement says 16-cycle Phusion Flash PCR adds Illumina adaptors with duplicate barcoding, but prints neither primer strand nor the barcodes.",
                    symmetric_species=True))
    graph.sequence(qc_amplicon, SequencingPlan(
        limitation="MiSeq v3 150-cycle kit with 10% PhiX (supplement, 'Library sequencing and validation'). "
                   "The adapter-adding primers, their barcodes and the read split are not printed, so no "
                   "library construct or primer site is drawn for this QC endpoint."))
    packaging = graph.input("packaging_system", "lentiviral packaging system", MaterialKind.CONDITION)
    virus = graph.product("virus", "guide-library lentiviral particles", MaterialKind.MOLECULAR,
                          state=plasmid_library.state)
    graph.react("package", "Package the pooled transfer-vector library", StepKind.GENERIC,
                ((plasmid_library, InputRole.TEMPLATE), (packaging, InputRole.CONDITION)),
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
                "array_cycling", "amplify_pool", "insert_gel", "size_select_insert",
                "bsmbi", "open_vector", "fastap", "dephosphorylate",
                "gibson_mix", "clone_guide", "electrocompetent_cells", "transform",
                "maxiprep", "maxiprep_library", "qc_adapter_primers", "qc_polymerase",
                "validate_library", "qc_amplicon", "packaging_system",
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
<p>The defining GeCKO v2 library starts with a synthesized 74-nt oligo pool. ArrayF/ArrayR PCR with Phusion HS Flex (63&nbsp;°C annealing, 15&nbsp;s extension, 20 cycles) makes the 140-bp insert, which is size-selected on a 2% E-Gel EX; the backbone is cut with BsmBI, dephosphorylated with FastAP for 2&nbsp;h at 37&nbsp;°C and gel-purified, and 10&nbsp;ng of insert is joined to 25&nbsp;ng of vector in a 20&nbsp;µl Gibson reaction (Sanjana 2014 Supplementary Methods, “GeCKO library pooled synthesis and cloning”). Parallel electroporations plated under carbenicillin give 30–200× library coverage; the colonies are scraped together and maxiprepped before packaging. The graph shows this pooled construction.</p>
<seq>{oligo("ArrayF", [Segment("ArrayF", ARRAY_F)])}
{oligo("ArrayR", [Segment("ArrayR", ARRAY_R)])}</seq>
<p>The paper's own sequencing endpoint is library QC, not a screen readout: the plasmid library is amplified for 16 cycles with Phusion Flash and primers that add Illumina adaptors, in duplicate and barcoded to measure PCR bias, and run on a MiSeq v3 150-cycle kit with 10% PhiX. That endpoint is now in the graph. The supplement prints neither those adapter primers nor the barcodes, so no library construct or primer site is drawn for it; that is an unresolved source question, not a statement that the authors published nothing. The final library below is the explicitly cited later Joung genomic readout, shown because this vector is used for screens, and it is not part of Sanjana et al.'s own experiment.</p>
{P.one_pcr_readout("Compatible one-PCR genomic readout")}
<p>The Moffat LCV2::TKOv3 two-step indexed readout is shown on its own page; it is not an extra step in the Joung one-PCR protocol above.</p>'''
    return P.page("lentiCRISPR v2 screening", "01_lenticrispr-v2.html", source, body,
                  graph=protocol_graph())


def main() -> None:
    OUT.write_text(render(), encoding="utf-8")
    print(f"wrote {OUT}  ({OUT.stat().st_size:,} bytes)")


if __name__ == "__main__": main()
