#!/usr/bin/env python3
from pathlib import Path
import sys
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "lib"))
import pooled_crispr as M
import pooled_crispr_page as P
import seqprimers as sp
from pcr import PCRGeometry
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind
from protocol_model import molecule, role

OUT = HERE.parent / "moffat-two-step-crispr-screen-readout.html"
SOURCE = 'Source: Moffat Lab <a href="https://media.addgene.org/cms/filer_public/44/9f/449f22c6-f818-4af6-8af4-2bab1d20f00a/moffat_lab_crispr_screen_protocol_v2.pdf">Pooled CRISPR Screen Protocol v2</a> (2019), LCV2::TKOv3 branch.'


def protocol_graph() -> ProtocolGraph:
    """PCR1 enrichment and PCR2 indexing as separate local reactions."""
    source = SOURCE
    graph = ProtocolGraph("Moffat two-step CRISPR screen readout", authority=source)
    gdna = graph.input("genomic_dna", "genomic DNA from the LCV2::TKOv3 screen",
                       MaterialKind.BIOLOGICAL,
                       state=molecule("screen genomic DNA", (role("integrated guide cassette"),)),
                       primary=True)
    first_forward = graph.input("pcr1_forward", "PCR1 locus forward primer", MaterialKind.OLIGO)
    first_reverse = graph.input("pcr1_reverse", "PCR1 locus reverse primer", MaterialKind.OLIGO)
    first_polymerase = graph.input("pcr1_polymerase", "PCR1 polymerase", MaterialKind.CATALYST)
    intermediate = M.two_step_intermediate()
    first = graph.product("pcr1_product", "approximately 600-bp guide-locus enrichment product",
                          MaterialKind.MOLECULAR, state=intermediate)
    graph.react("pcr1", "PCR1 enriches the integrated guide locus", StepKind.PCR,
                ((gdna, InputRole.TEMPLATE), (first_forward, InputRole.STRUCTURAL),
                 (first_reverse, InputRole.STRUCTURAL), (first_polymerase, InputRole.CATALYST)),
                (first,), geometry=PCRGeometry.unavailable(
                    intermediate, reason="The surrounding integrated-vector/genomic template is schematic."))
    second_forward = graph.input("pcr2_forward", "PCR2 adds the dual-index sequencing arms forward primer",
                                 MaterialKind.OLIGO)
    second_reverse = graph.input("pcr2_reverse", "PCR2 adds the dual-index sequencing arms reverse primer",
                                 MaterialKind.OLIGO)
    second_polymerase = graph.input("pcr2_polymerase", "PCR2 adds the dual-index sequencing arms polymerase",
                                   MaterialKind.CATALYST)
    final = M.two_step_library()
    library = graph.product("pcr2_library", "PCR2 adds the dual-index sequencing arms sequencing library",
                            MaterialKind.LIBRARY, state=final)
    graph.react("pcr2_pcr", "PCR2 adds the dual-index sequencing arms", StepKind.PCR,
                ((first, InputRole.TEMPLATE), (second_forward, InputRole.STRUCTURAL),
                 (second_reverse, InputRole.STRUCTURAL), (second_polymerase, InputRole.CATALYST)),
                (library,), geometry=PCRGeometry.unavailable(
                    final, reason="The genomic/vector template is larger than the exact final amplicon construct."))
    graph.sequence(library, SequencingPlan(library=final,
                                           primers=tuple(M.two_step_sequencing_primers())))
    graph.declare_reviewed_scope(
        ("genomic_dna", "pcr1_forward", "pcr1_reverse", "pcr1_polymerase", "pcr1",
         "pcr2_forward", "pcr2_reverse", "pcr2_polymerase", "pcr2_pcr", "pcr2_library"), source)
    return graph.require_publishable()


def render() -> str:
    source = SOURCE
    intermediate, final = M.two_step_intermediate(), M.two_step_library()
    body = f'''<h2>PCR1 — enrich the integrated guide locus</h2>
{P.construct_panel(intermediate, "PCR1 uses ordinary locus primers and produces the ~600-bp enrichment product. No flow-cell adapters are present.")}
<h2>PCR2 — add dual-index TruSeq arms</h2>
{P.construct_panel(final, "PCR2 uses the nested hU6 and scaffold sites and yields the ~200-bp sequencing library.")}
{sp.section(final, M.two_step_sequencing_primers(), heading="Sequencing primers", intro="The protocol uses single-end sequencing with 21 dark cycles, then 26 imaged guide cycles, plus 8-cycle i7 and i5 reads.", required_roles=("Read 1", "Index 1 (i7)", "Index 2 (i5)"))}'''
    return P.page("Moffat two-step CRISPR screen readout", "01_moffat-two-step.html",
                  source, body, graph=protocol_graph())


def main() -> None:
    OUT.write_text(render(), encoding="utf-8")
    print(f"wrote {OUT}  ({OUT.stat().st_size:,} bytes)")


if __name__ == "__main__": main()
