#!/usr/bin/env python3
from pathlib import Path
import sys
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "lib"))
import pooled_crispr_page as P
import pooled_crispr as M
from pcr import PCRGeometry
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind
from protocol_model import molecule, role

OUT = HERE.parent / "joung-one-pcr-crispr-screen-readout.html"
SOURCE = 'Source: <a href="https://doi.org/10.1038/nprot.2017.016">Joung et al., <i>Nature Protocols</i> (2017)</a>. This is one PCR, not the PCR2 half of a nested workflow.'


def protocol_graph() -> ProtocolGraph:
    """The single amplification belongs to this published readout."""
    source = SOURCE
    graph = ProtocolGraph("Joung one-PCR CRISPR screen readout", authority=source)
    gdna = graph.input("genomic_dna", "genomic DNA from the pooled screen",
                       MaterialKind.BIOLOGICAL,
                       state=molecule("screen genomic DNA", (role("integrated guide cassette"),)),
                       primary=True)
    forward = graph.input("readout_forward", "One-PCR guide-cassette readout forward primer",
                          MaterialKind.OLIGO)
    reverse = graph.input("readout_reverse", "One-PCR guide-cassette readout reverse primer",
                          MaterialKind.OLIGO)
    polymerase = graph.input("readout_polymerase", "One-PCR guide-cassette readout polymerase",
                             MaterialKind.CATALYST)
    final = M.one_pcr_library()
    library = graph.product("readout_library", "One-PCR guide-cassette readout sequencing library",
                            MaterialKind.LIBRARY, state=final)
    graph.react("readout_pcr", "One-PCR guide-cassette readout", StepKind.PCR,
                ((gdna, InputRole.TEMPLATE), (forward, InputRole.STRUCTURAL),
                 (reverse, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)),
                (library,), geometry=PCRGeometry.unavailable(
                    final, reason="The genomic/vector template is larger than the exact final amplicon construct."))
    graph.sequence(library, SequencingPlan(library=final,
                                           primers=tuple(M.one_pcr_sequencing_primers())))
    graph.declare_reviewed_scope(
        ("genomic_dna", "readout_forward", "readout_reverse", "readout_polymerase",
         "readout_pcr", "readout_library"), source)
    return graph.require_publishable()


def render() -> str:
    source = SOURCE
    body = P.one_pcr_readout("Genomic DNA to sequencing library — one PCR")
    return P.page("Joung one-PCR CRISPR screen readout", "01_joung-one-pcr.html",
                  source, body, graph=protocol_graph())


def main() -> None:
    OUT.write_text(render(), encoding="utf-8")
    print(f"wrote {OUT}  ({OUT.stat().st_size:,} bytes)")


if __name__ == "__main__": main()
