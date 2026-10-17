"""Architecture model for 10x Chromium Single Cell 3' Gene Expression v1."""
from __future__ import annotations

import illumina as il
import seqprimers as sp
from chemdraw import Construct, Scene, Segment, feature, revcomp
import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[2]))
from adapter_ligation import AdapterLigationGeometry
from pcr import PCRGeometry
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind
from protocol_model import molecule, role

CELL_BARCODE = feature("cell_barcode", "cell_barcode", "whitelist",
                       whitelist="10x-chromium-3prime-v1")
UMI = feature("umi", "umi", "random")
I5 = feature("sample_index_i5", "sample_index", "unknown")
FEATURES = {"cell barcode": CELL_BARCODE, "cell barcode reverse complement": CELL_BARCODE,
            "UMI": UMI, "i5 sample index": I5}


def seg(name: str, top: str, tag: str | None = None, **kw) -> Segment:
    kw.setdefault("feature", FEATURES.get(name))
    return Segment(name=name, top=top, tag=tag, **kw)


def inferred(name: str, sequence: str, tag: str | None = None, **kw) -> Segment:
    """Mark reconstructed adapter constants as inferred from the v1 architecture."""
    return seg(name, sequence, tag, inferred=True, **kw)


def bead_oligo() -> Construct:
    return Construct([
        inferred("P7", il.P7, "p7"),
        seg("cell barcode", "B" * 14, "cbc", placeholder=True),
        inferred("TruSeq Read 2", il.TRUSEQ_READ2, "r2"),
        seg("UMI", "N" * 10, "umi", placeholder=True), seg("poly(T)", "T" * 30),
    ], name="10x 3-prime v1 bead oligo architecture")


def rt_scene(insert_nt: int = 28) -> Scene:
    mrna = [seg("transcript", "X" * insert_nt, placeholder=True), seg("poly(A)", "A" * 30)]
    sc = Scene()
    sc.strand("mRNA", mrna, label="mRNA", mod3="poly(A)")
    sc.anneal("bead", list(bead_oligo()), to="mRNA", pair=("poly(T)", "poly(A)"),
              label="released gel-bead primer")
    sc.arrow("bead", "reverse transcription and template switching")
    return sc


def first_strand(insert_nt: int = 28) -> Construct:
    return Construct([*bead_oligo(), seg("antisense cDNA", "X" * insert_nt,
                                        placeholder=True),
                      inferred("template-switch end", "N" * 18, "tso",
                               placeholder=True)],
                     name="10x 3-prime v1 first-strand architecture")


def amplified_cdna(insert_nt: int = 28) -> Construct:
    return Construct(list(first_strand(insert_nt)),
                     name="10x 3-prime v1 amplified-cDNA architecture")


SEQ_PRIMERS = tuple(sp.TRUSEQ[k] for k in ("R1", "I1", "I2", "R2"))


def final_library(insert_nt: int = 28) -> Construct:
    lib = Construct([
        inferred("P5", il.P5, "p5"), seg("i5 sample index", "I" * 8, "cbc", placeholder=True),
        inferred("TruSeq Read 1", il.TRUSEQ_READ1, "r1"),
        seg("3-prime cDNA", "X" * insert_nt, placeholder=True),
        seg("poly(A)", "A" * 30), seg("UMI", "N" * 10, "umi", placeholder=True),
        inferred("TruSeq Read-2 reverse complement", revcomp(il.TRUSEQ_READ2), "r2"),
        seg("cell barcode reverse complement", "B" * 14, "cbc", placeholder=True),
        inferred("P7 reverse complement", il.P7_RC, "p7"),
    ], name="10x Chromium 3-prime v1 library")
    problems = sp.verify(lib, SEQ_PRIMERS)
    if problems:
        raise ValueError("invalid 10x 3-prime v1 library: " + "; ".join(problems))
    return lib


def final_scene() -> Scene:
    return Scene.duplex(list(final_library()), label="library")


def _cleanup(graph, substrate, step_id, label, kept_id, kept_label, waste_label, *,
             kind=MaterialKind.MOLECULAR, state=None):
    """One bead cleanup: retained eluate plus the discarded supernatant/wash fraction."""
    beads = graph.input(step_id + "_beads", label, MaterialKind.SELECTION)
    kept = graph.product(kept_id, kept_label, kind,
                         state=substrate.state if state is None else state,
                         inferred=substrate.inferred, source=substrate.source)
    waste = graph.product(step_id + "_waste", waste_label, MaterialKind.WASTE)
    graph.react(step_id, label, StepKind.CLEANUP,
                ((substrate, InputRole.TEMPLATE), (beads, InputRole.SELECTION)), (kept, waste))
    graph.discard(waste, "removed with the supernatant and ethanol washes")
    return kept


def protocol_graph():
    """Local v1 route (paper Methods; CG00026 Rev B Steps 1-3).

    Physical Covaris shearing, not enzymatic fragmentation, and ligation of the
    R1 Adaptor Mix (Read 1 side), not a Read 2 adapter: the bead supplies P7/Read 2.
    """
    authority = ('Defining source: <a href="https://doi.org/10.1038/ncomms14049">'
                 'Zheng et al., <i>Nature Communications</i> (2017)</a>; protocol '
                 'authority: 10x Genomics CG00026 Rev B.')
    graph = ProtocolGraph("10x Chromium Single Cell 3-prime Gene Expression v1", authority=authority)
    inferred_src = authority + " (intermediate structure inferred)"
    first_state = first_strand()
    final = final_library()
    cells = graph.input("cells", "cells containing poly(A) RNA", MaterialKind.BIOLOGICAL,
                        state=molecule("cells", (role("poly(A) RNA"),)), primary=True)
    gems = graph.input("gems", "Chromium GEM partitioning condition", MaterialKind.CONDITION)
    occupied = graph.product("occupied", "occupied GEMs containing cells", MaterialKind.MOLECULAR,
                             state=cells.state)
    failed = graph.product("failed", "empty, multiplet or unresolved GEM products",
                           MaterialKind.WASTE)
    graph.react("partition", "Partition cells with gel beads", StepKind.PARTITION,
                ((cells, InputRole.TEMPLATE), (gems, InputRole.CONDITION)), (occupied, failed),
                changes_molecule=False)
    graph.discard(failed, "does not retain a resolvable single-cell barcode")
    bead = graph.input("rt_primer", "cell-barcode/UMI oligo(dT) gel-bead primer",
                       MaterialKind.OLIGO, inferred=True,
                       source=authority + " (architecture only; complete bead oligo unresolved)")
    rt_primer_mix = graph.input("rt_primer_mix", "kit RT Primer (lyophilized; identity unstated)",
                                MaterialKind.OLIGO, inferred=True,
                                source="CG00026 Rev B GEM master mix; role and sequence not given")
    rt = graph.input("rt", "RT Enzyme (identity not named)", MaterialKind.CATALYST)
    first = graph.product("first_strand", "first-strand transcript cDNA", MaterialKind.MOLECULAR,
                          state=first_state, inferred=True, source=inferred_src)
    graph.react("reverse_transcribe", "GEM-RT 55 °C 2 h, 85 °C 5 min", StepKind.REVERSE_TRANSCRIPTION,
                ((occupied, InputRole.TEMPLATE), (bead, InputRole.STRUCTURAL),
                 (rt_primer_mix, InputRole.STRUCTURAL), (rt, InputRole.CATALYST)), (first,))
    tso = graph.input("tso", "template-switch oligo (bases unresolved)", MaterialKind.OLIGO,
                      inferred=True, source=authority + " (TSO bases not asserted)")
    switched = graph.product("switched", "template-switched cell-barcoded cDNA",
                             MaterialKind.MOLECULAR, state=first_state, inferred=True,
                             source=inferred_src)
    graph.react("template_switch", "INFERRED — install the template-switch handle",
                StepKind.TEMPLATE_SWITCHING,
                ((first, InputRole.TEMPLATE), (tso, InputRole.STRUCTURAL),
                 (rt, InputRole.CATALYST)), (switched,))
    recovery = graph.input("gem_recovery", "Recovery Agent (GEM break)", MaterialKind.CONDITION)
    broken = graph.product("broken_gems", "aqueous post-GEM-RT phase", MaterialKind.MOLECULAR,
                           state=first_state, inferred=True, source=inferred_src)
    oil = graph.product("recovery_oil", "Recovery Agent/Partitioning Oil phase", MaterialKind.WASTE)
    graph.react("recover", "Break GEMs with Recovery Agent; remove the oil phase",
                StepKind.EXTRACTION, ((switched, InputRole.TEMPLATE), (recovery, InputRole.CONDITION)),
                (broken, oil), changes_molecule=False)
    graph.discard(oil, "Recovery Agent/Partitioning Oil removed from the tube bottom")
    silane = _cleanup(graph, broken, "silane_cleanup", "DynaBeads MyOne Silane cleanup",
                      "silane_cdna", "silane-purified first-strand cDNA",
                      "leftover GEM-RT reagents and primers")
    recovered = _cleanup(graph, silane, "post_rt_spri", "SPRIselect 0.6x cleanup",
                         "recovered", "recovered cell-barcoded cDNA",
                         "unbound short species in the 0.6x supernatant")
    paper_primers = "paper: primers complementary to the switch oligos and sequencing adapters"
    pre_tso = graph.input("preamp_forward", "cDNA Primer Mix: switch-oligo-side primer (not printed)",
                          MaterialKind.OLIGO, inferred=True, source=paper_primers)
    pre_adapter = graph.input("preamp_reverse",
                              "cDNA Primer Mix: sequencing-adapter-side primer (not printed)",
                              MaterialKind.OLIGO, inferred=True, source=paper_primers)
    pre_pol = graph.input("preamp_polymerase", "Amplification Master Mix", MaterialKind.CATALYST)
    amplified = graph.product("amplified", "amplified cell-barcoded cDNA", MaterialKind.MOLECULAR,
                              state=amplified_cdna(), inferred=True, source=inferred_src)
    graph.react("preamp", "Amplify cDNA: 98 °C 3 min; 10–14 × (98/67/72 °C) by targeted recovery",
                StepKind.PCR, ((recovered, InputRole.TEMPLATE), (pre_tso, InputRole.STRUCTURAL),
                               (pre_adapter, InputRole.STRUCTURAL),
                               (pre_pol, InputRole.CATALYST)), (amplified,),
                geometry=PCRGeometry.unavailable(first_state,
                    reason="cDNA Primer Mix sequences are not printed in the paper or CG00026 Rev B."))
    cleaned = _cleanup(graph, amplified, "cdna_spri", "SPRIselect 0.6x cDNA cleanup",
                       "cleaned_cdna", "purified amplified cDNA",
                       "short products and PCR reagents in the 0.6x supernatant")
    shearing = graph.input("covaris", "Covaris M220 acoustic shearing (~200 bp target)",
                           MaterialKind.CONDITION)
    fragments = graph.product("fragments", "sheared cDNA fragments", MaterialKind.MOLECULAR)
    graph.react("fragment", "Shear amplified cDNA to ~200 bp (Covaris)", StepKind.FRAGMENTATION,
                ((cleaned, InputRole.TEMPLATE), (shearing, InputRole.CONDITION)), (fragments,))
    size_beads = graph.input("shear_size_selection_beads",
                             "Double-sided SPRIselect (30 µl, keep supernatant; +10 µl, keep beads)",
                             MaterialKind.SELECTION)
    sized = graph.product("sized_fragments", "size-selected sheared fragments", MaterialKind.MOLECULAR)
    large = graph.product("shear_size_selection_large", "large fragments bound by the first bead addition",
                          MaterialKind.WASTE)
    small = graph.product("shear_size_selection_small", "small fragments left in the second supernatant",
                          MaterialKind.WASTE)
    graph.react("shear_size_selection", "Post-shearing double-sided SPRIselect",
                StepKind.SELECTION, ((fragments, InputRole.TEMPLATE), (size_beads, InputRole.SELECTION)),
                (sized, large, small))
    graph.discard(large, "first-bead pellet is discarded; the supernatant is carried forward")
    graph.discard(small, "second-bead supernatant is removed; bead-bound fraction is eluted")
    end_mix = graph.input("end_mix", "End Repair & A-tailing enzymes", MaterialKind.CATALYST)
    tailed = graph.product("tailed", "end-repaired dA-tailed cDNA fragments", MaterialKind.MOLECULAR)
    graph.react("end_prep", "End-repair and A-tail: 20 °C 30 min, 65 °C 30 min",
                StepKind.END_REPAIR, ((sized, InputRole.TEMPLATE), (end_mix, InputRole.CATALYST)),
                (tailed,))
    adapter = graph.input("r1_adaptor", "R1 Adaptor Mix (Read 1 side; sequence not printed)",
                          MaterialKind.ADAPTER, inferred=True,
                          source="CG00026 Rev B Adaptor Ligation Mix; P5/R1/sample index added in library construction")
    ligase = graph.input("ligase", "DNA Ligase + Ligation Buffer", MaterialKind.CATALYST)
    adapted = graph.product("adapted", "R1-adaptor-ligated cDNA fragments", MaterialKind.MOLECULAR)
    graph.react("ligate", "Ligate the R1 Adaptor Mix at 20 °C, 15 min", StepKind.LIGATION,
                ((tailed, InputRole.TEMPLATE), (adapter, InputRole.STRUCTURAL),
                 (ligase, InputRole.CATALYST)), (adapted,),
                geometry=AdapterLigationGeometry.unavailable(
                    "CG00026 Rev B names the R1 Adaptor Mix but prints neither its oligos nor "
                    "its ligating end geometry (overhang, polarity, fork)."))
    lig1 = _cleanup(graph, adapted, "ligation_cleanup_1", "Post-ligation SPRIselect cleanup 1 (88 µl)",
                    "ligated_1", "purified ligated fragments",
                    "unligated adaptor and ligation reagents in the supernatant")
    ligated = _cleanup(graph, lig1, "ligation_cleanup_2", "Post-ligation SPRIselect cleanup 2 (50 µl)",
                       "ligated", "twice-purified ligated fragments",
                       "residual adaptor in the supernatant")
    si_primer = graph.input("index_forward", "SI-PCR Primer (common; sequence not printed)",
                            MaterialKind.OLIGO, inferred=True,
                            source="CG00026 Rev B Sample Index PCR Mix")
    si_plate = graph.input("index_reverse",
                           "Single Cell 3' Sample Index Plate well (4 sequences per set; not printed)",
                           MaterialKind.OLIGO, inferred=True,
                           source="CG00026 Rev B Sample Index PCR")
    index_pol = graph.input("index_polymerase", "Amplification Master Mix", MaterialKind.CATALYST)
    indexed = graph.product("indexed", "sample-index PCR product", MaterialKind.LIBRARY, state=final)
    graph.react("index_pcr", "Sample index PCR: 98 °C 45 s; 10 × (98/60/72 °C)", StepKind.PCR,
                ((ligated, InputRole.TEMPLATE), (si_primer, InputRole.STRUCTURAL),
                 (si_plate, InputRole.STRUCTURAL), (index_pol, InputRole.CATALYST)), (indexed,),
                geometry=PCRGeometry.unavailable(final,
                    reason="SI-PCR Primer and Sample Index Plate oligo sequences are not printed."))
    si1 = _cleanup(graph, indexed, "si_cleanup_1", "Post-SI-PCR SPRIselect cleanup 1 (100 µl)",
                   "si_clean_1", "purified indexed library", "PCR reagents and primers in the supernatant",
                   kind=MaterialKind.LIBRARY, state=final)
    library = _cleanup(graph, si1, "si_cleanup_2", "Post-SI-PCR SPRIselect cleanup 2 (50 µl)",
                       "library", "indexed Chromium 3' v1 library", "residual primers in the supernatant",
                       kind=MaterialKind.LIBRARY, state=final)
    graph.sequence(library, SequencingPlan(library=final, primers=SEQ_PRIMERS,
                   acquired_roles=tuple(p.role for p in SEQ_PRIMERS),
                   read_lengths={"Read 1": 98, "Index 1 (i7)": 14, "Index 2 (i5)": 8, "Read 2": 10}))
    # Checked 2026-10-17 against paper Methods (JATS line 252) and CG00026 Rev B Steps 1-3.
    graph.declare_reviewed_scope([
        "cells", "gems", "partition", "rt_primer", "rt_primer_mix", "rt", "reverse_transcribe",
        "tso", "template_switch", "gem_recovery", "recover", "silane_cleanup",
        "silane_cleanup_beads", "post_rt_spri", "post_rt_spri_beads", "preamp_forward", "preamp_reverse",
        "preamp_polymerase", "preamp", "cdna_spri", "cdna_spri_beads", "covaris", "fragment",
        "shear_size_selection_beads", "shear_size_selection", "end_mix", "end_prep",
        "r1_adaptor", "ligase", "ligate", "ligation_cleanup_1", "ligation_cleanup_1_beads",
        "ligation_cleanup_2", "ligation_cleanup_2_beads", "index_forward", "index_reverse",
        "index_polymerase", "index_pcr", "si_cleanup_1", "si_cleanup_1_beads", "si_cleanup_2",
        "si_cleanup_2_beads", "library"], authority)
    return graph.require_publishable()
