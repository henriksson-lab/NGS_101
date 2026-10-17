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
from chemdraw import revcomp

OUT = HERE.parent / "joung-one-pcr-crispr-screen-readout.html"
SOURCE = 'Source: <a href="https://doi.org/10.1038/nprot.2017.016">Joung et al., <i>Nature Protocols</i> (2017)</a>. This is one PCR, not the PCR2 half of a nested workflow.'


# Table 3 (PMC5526071 JATS lines 2164/2184): representative members as printed.
NGS_LIB_FWD_1 = ("AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT"
                 "TAAGTAGAGGCTTTATATATCTTGTGGAAAGGACGAAACACC")
NGS_LIB_KO_REV_1 = ("CAAGCAGAAGACGGCATACGAGATTCGCCTTGGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT"
                    "CCGACTCGGTGCCACTTTTTCAA")
READ_LENGTHS = {"Read 1": 80, "Index 1 (i7)": 8}  # Step 37


def _checked_library():
    final = M.one_pcr_library()
    top = final.top()
    if not (top.startswith(NGS_LIB_FWD_1) and top.endswith(revcomp(NGS_LIB_KO_REV_1))):
        raise ValueError("one-PCR library ends do not match the Table 3 primer pair")
    return final


def protocol_graph() -> ProtocolGraph:
    """The single amplification belongs to this published readout."""
    source = SOURCE
    graph = ProtocolGraph("Joung one-PCR CRISPR screen readout", authority=source)
    gdna = graph.input("genomic_dna", "genomic DNA from the pooled screen (up to 2.5 µg per 50 µl reaction, Step 58)",
                       MaterialKind.BIOLOGICAL,
                       state=molecule("screen genomic DNA", (role("integrated guide cassette"),)),
                       primary=True)
    forward = graph.input("readout_forward",
                          "NGS-Lib-Fwd-1 (Table 3; one of ten staggered forward primers, one per reaction, all ten required)",
                          MaterialKind.OLIGO, state=NGS_LIB_FWD_1)
    reverse = graph.input("readout_reverse",
                          "NGS-Lib-KO-Rev-1 (Table 3; one of eight i7-barcoded GeCKO reverse primers, one per screen sample)",
                          MaterialKind.OLIGO, state=NGS_LIB_KO_REV_1)
    polymerase = graph.input("readout_polymerase", "NEBNext High-Fidelity 2× PCR Master Mix (NEB M0541)",
                             MaterialKind.CATALYST)
    final = _checked_library()
    library = graph.product("readout_library", "One-PCR guide-cassette readout sequencing library",
                            MaterialKind.LIBRARY, state=final)
    graph.react("readout_pcr",
                "One PCR per forward primer (Steps 32–33, 58): 95 °C 5 min; 12 cycles 98 °C 20 s, 60 °C 15 s, 72 °C 15 s; 72 °C 1 min",
                StepKind.PCR,
                ((gdna, InputRole.TEMPLATE), (forward, InputRole.STRUCTURAL),
                 (reverse, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)),
                (library,), geometry=PCRGeometry.unavailable(
                    final, reason="The genomic/vector template is larger than the exact final amplicon construct."))
    column = graph.input("zymo_spin_v", "Zymo-Spin V column with reservoir (Steps 59–63)",
                         MaterialKind.SELECTION)
    clean = graph.product("purified_library", "Pooled, column-purified readout library (Step 64)",
                          MaterialKind.LIBRARY, state=final)
    flow = graph.product("column_flowthrough", "Zymo-Spin V flow-through and washes",
                         MaterialKind.WASTE)
    graph.react("pool_purify", "Pool PCRs and purify on Zymo-Spin V columns", StepKind.CLEANUP,
                ((library, InputRole.TEMPLATE), (column, InputRole.SELECTION)), (clean, flow))
    graph.discard(flow, "unbound PCR-reaction material removed by the column")
    primers = tuple(M.one_pcr_sequencing_primers())
    graph.sequence(clean, SequencingPlan(library=final, primers=primers,
                                         acquired_roles=tuple(p.role for p in primers),
                                         read_lengths=READ_LENGTHS))
    graph.declare_reviewed_scope(
        ("genomic_dna", "readout_forward", "readout_reverse", "readout_polymerase",
         "readout_pcr", "readout_library", "zymo_spin_v", "pool_purify",
         "purified_library"), source)
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
