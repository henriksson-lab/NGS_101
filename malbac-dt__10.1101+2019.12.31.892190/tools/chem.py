"""Source-owned protocol definition and reaction graph.

Authority (published version first): Chapman, Lee et al., "Correlated gene modules
uncovered by high-precision single-cell transcriptomics", PNAS 119, e2206938119 (2022),
doi:10.1073/pnas.2206938119, PMC9907105, Materials and Methods "MALBAC-DT Protocol" and
"Automated MALBAC-DT Protocol" and "Sequence Processing".  The bioRxiv preprint
10.1101/2019.12.31.892190 (the directory's defining DOI) is superseded and not used for
active claims.  Neither the article nor its SI Appendix prints a sequence for RT-An,
RT-Bn, GAT-7N, GAT-COM, Tru2-Gn-RT or Ix-Tru2, so their regions are role tokens and
their order inside each oligo is inferred from the oligo names only.
"""
from __future__ import annotations

import illumina as il
import nextera as nx
import seqprimers as sp
from batch_ngs import seg
from chemdraw import Construct, feature
from pcr import PCRGeometry
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind)
from protocol_model import molecule, role
from protocol_specs import Definition, get

PNAS = "Chapman, Lee et al., PNAS 119, e2206938119 (2022), PMC9907105"
SEQUENCING = (
    "The original plate (U2OS, HEK293T) was tagmented with the Nextera kit and amplified "
    "with a Nextera S5XX i5 primer and an Ix-Tru2 primer substituted for the S7XX primer, "
    "selecting the barcoded 3′ ends. The automated plate (K562) used the Vazyme TruePrep "
    "kit and NovaSeq 2 × 150 bp with a custom Read 2 sequencing primer, so the barcode/UMI "
    "end is read from the Tru2 side; transcript sequence comes from the Tn5 (Read 1) end.")

DEFINITION = Definition(
    "MALBAC-DT — high-detectability single-cell transcriptomics", "10.1073/pnas.2206938119",
    PNAS + "; preprint bioRxiv 10.1101/2019.12.31.892190",
    "A barcoded poly(T)/UMI primer reverse-transcribes each sorted cell; excess primer is "
    "removed by Exo I; random-primer MALBAC looping cycles and PCR amplify the first strand, "
    "a Tru2-Gn-RT primer finishes amplification, and pooled products are tagmented and "
    "amplified with an Ix-Tru2 primer that selects the barcoded 3′ ends.",
    "nextera", ("transcript cDNA",), SEQUENCING, "RNA")

CAVEAT = (
    "<b>Oligo sequences and barcode order are not published in the inspected sources.</b> "
    "PNAS Methods name RT-An, RT-Bn, GAT-7N, GAT-COM, Tru2-Gn-RT and Ix-Tru2 but print no "
    "sequence; SI Appendix contains only performance and module tables. The UMI pattern "
    "(HBDV)<sub>5</sub> (20 nt) is from Sequence Processing. The order of poly(T), UMI, "
    "inner barcode and RT handle inside RT-An, the Tru2/Gn/RT order (from the primer name) "
    "and the P7 end added by Ix-Tru2 are inferred. The role of RT-Bn, added with the Exo I "
    "mix and used as the (VDBH)<sub>5</sub> pattern that consensus UMIs must differ from, is "
    "not stated. The automated-plate oligos (96 well-specific RT3-An, one Tru2-G-RT) are "
    "drawn with the same roles; the cited automated-protocol reference (33) is not "
    "resolvable in the reference list.")


def _final() -> tuple[Construct, tuple]:
    well = "malbac_dt_well"
    lib = Construct([
        seg("P5", il.P5, "p5"),
        seg("i5", "J" * 8, "cbc", placeholder=True,
            feature=feature("sample_index_i5", "sample_index", "unknown")),
        seg("S5", nx.S5, "s5"),
        seg("left mosaic end", nx.ME, "me"),
        role("transcript cDNA"),
        role("poly(A)/poly(T)", inferred=True),
        role("UMI (HBDV)5", tag="umi", length_bp=20, inferred=True,
             molecular_feature=feature("umi", "umi", "degenerate",
                                       note="(HBDV)5 RT-An pattern; consensus must differ "
                                            "from the (VDBH)5 RT-Bn pattern at >=3 bases")),
        role("inner RT-An barcode", tag="cbc", inferred=True,
             molecular_feature=feature("inner_barcode", "cell_barcode", "unknown",
                                       group=well, part="inner (RT-An row or well)")),
        role("RT handle", inferred=True),
        role("outer Gn barcode", tag="cbc", inferred=True,
             molecular_feature=feature("outer_barcode", "cell_barcode", "unknown",
                                       group=well, part="outer (Tru2-Gn-RT column)")),
        role("Tru2 handle", inferred=True),
        role("Ix index", tag="cbc", inferred=True,
             molecular_feature=feature("sample_index_ix", "sample_index", "unknown")),
        seg("P7 reverse complement", il.P7_RC, "p7", inferred=True),
    ], name=DEFINITION.title + " final library")
    primers = (sp.NEXTERA["R1"],)
    problems = sp.verify(lib, primers, required_roles=("Read 1",))
    if problems:
        raise ValueError("; ".join(problems))
    return lib, primers


FINAL_LIBRARY, SEQ_PRIMERS = _final()
SPEC = get("malbac-dt", definition=DEFINITION, final_override=(FINAL_LIBRARY, SEQ_PRIMERS),
           sequencing_unavailable=(
               "Read 2 uses a custom primer on the Tru2 side (automated plate) whose sequence, "
               "like the Tru2 handle and Ix index, is not printed; index-read primers and the "
               "original plate's run configuration are not stated."))
SPEC.CAVEAT = CAVEAT


def _parts(names, inferred):
    return tuple(role(n, inferred=n in inferred) for n in names)


def _state(label, *names, inferred=()):
    return molecule(label, _parts(names, inferred))


def _local_graph(d, final, primers):
    source = f"{PNAS}; doi:{d.doi}"
    graph = ProtocolGraph(d.title, authority=source)
    rna = graph.input("rna", "single sorted cell in 3 µL SSIV/ICA-630/betaine lysis buffer",
                      MaterialKind.BIOLOGICAL,
                      state=_state("single-cell poly(A) RNA", "poly(A) mRNA"), primary=True)
    rt_primer = graph.input("rt_primer", "RT-An primer (poly-T, cell barcode, (HBDV)5 UMI; "
                            "sequence not printed), 0.2 µL 50 µM", MaterialKind.OLIGO)
    rt = graph.input("reverse_transcriptase", "SuperScript IV; 72 °C 3 min, then 55 °C 10 min",
                     MaterialKind.CATALYST)
    cdna = graph.product("cdna", "barcoded first-strand cDNA", MaterialKind.MOLECULAR,
                         state=_state("barcoded first-strand cDNA (antisense)", "RT handle",
                                      "inner RT-An barcode", "UMI", "poly(T)",
                                      "antisense transcript cDNA",
                                      inferred=("poly(T)", "UMI", "inner RT-An barcode",
                                                "RT handle")))
    graph.react("reverse_transcribe", "Reverse-transcribe with the barcoded poly-T RT-An primer",
                StepKind.REVERSE_TRANSCRIPTION,
                ((rna, InputRole.TEMPLATE), (rt_primer, InputRole.STRUCTURAL),
                 (rt, InputRole.CATALYST)), (cdna,))
    exo = graph.input("exo_i", "Exo I (0.6 µL in 1 µL mix); 37 °C 30 min, 80 °C 20 min",
                      MaterialKind.CATALYST)
    rt_b = graph.input("rt_bn", "RT-Bn primer added with the Exo I mix (role not stated; "
                       "sequence not printed)", MaterialKind.OLIGO)
    exo_cdna = graph.product("exo_cdna", "first-strand cDNA after excess RT-An digestion",
                             MaterialKind.MOLECULAR, state=cdna.state)
    graph.react("exo_digest", "Degrade excess RT primers with exonuclease I",
                StepKind.DIGESTION,
                ((cdna, InputRole.TEMPLATE), (exo, InputRole.CATALYST),
                 (rt_b, InputRole.STRUCTURAL)), (exo_cdna,), changes_molecule=False)
    random_primer = graph.input("malbac_primer", "GAT-7N and GAT-COM primers, 0.4 µL 50 µM "
                                "each (sequences not printed)", MaterialKind.OLIGO)
    polymerase = graph.input("malbac_polymerase", "Deep Vent (exo-) in ThermoPol buffer",
                             MaterialKind.CATALYST)
    amplified = graph.product("amplified", "MALBAC-amplified barcoded cDNA",
                              MaterialKind.MOLECULAR,
                              state=_state("MALBAC-DT amplicon", "GAT common sequence",
                                           "transcript cDNA", "poly(T)", "UMI",
                                           "inner RT-An barcode", "RT handle",
                                           inferred=("GAT common sequence", "poly(T)", "UMI",
                                                     "inner RT-An barcode", "RT handle")))
    graph.react("malbac", "MALBAC: 10 looping cycles (4–65 °C ramp, 95 °C, 58 °C), then "
                "PCR steps 12–16 (written 'Goto 13 17×')", StepKind.EXTENSION,
                ((exo_cdna, InputRole.TEMPLATE), (random_primer, InputRole.STRUCTURAL),
                 (polymerase, InputRole.CATALYST)), (amplified,))
    tru2 = graph.input("tru2_gn_rt", "Tru2-Gn-RT primer, 0.4 µL 50 µM (12 column primers in "
                       "the original plate; one in the automated plate)", MaterialKind.OLIGO)
    finished_names = ("GAT common sequence", "transcript cDNA", "poly(T)", "UMI",
                      "inner RT-An barcode", "RT handle", "outer Gn barcode", "Tru2 handle")
    finished_inferred = tuple(n for n in finished_names if n != "transcript cDNA")
    finished_state = _state("Tru2-tagged amplicon", *finished_names,
                            inferred=finished_inferred)
    finished = graph.product("finished", "Tru2-Gn-RT-extended amplicons", MaterialKind.MOLECULAR,
                             state=finished_state)
    graph.react("tru2_pcr", "Add Tru2-Gn-RT and run five more PCR cycles (steps 12–15)",
                StepKind.PCR,
                ((amplified, InputRole.TEMPLATE), (tru2, InputRole.STRUCTURAL),
                 (random_primer, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)),
                (finished,),
                geometry=PCRGeometry.unavailable(
                    Construct(list(_parts(finished_names, finished_inferred)),
                              name="Tru2-tagged amplicon"),
                    reason="Tru2-Gn-RT, GAT-COM and the RT-handle target are not printed; the GAT primers remaining in the same well are assumed to prime the opposite end."))
    beads = graph.input("ampure", "0.8× AMPure beads", MaterialKind.SELECTION)
    pooled = graph.product("pooled", "pooled, bead-purified plate amplicons (1 µL per well)",
                           MaterialKind.MOLECULAR, state=finished_state)
    bead_waste = graph.product("ampure_supernatant", "bead supernatant with short products, "
                               "primers and nucleotides", MaterialKind.WASTE)
    graph.react("pool_cleanup", "Pool 1 µL per well and purify with 0.8× AMPure",
                StepKind.CLEANUP,
                ((finished, InputRole.TEMPLATE), (beads, InputRole.SELECTION)),
                (pooled, bead_waste))
    graph.discard(bead_waste, "bead supernatant discarded")
    transposome = graph.input("transposome", "Nextera transposome (original plate) or Vazyme "
                              "TruePrep (automated plate)", MaterialKind.ADAPTER)
    tn5 = graph.input("tn5", "kit Tn5 transposase", MaterialKind.CATALYST)
    tagged = graph.product("tagged", "tagmented amplicons", MaterialKind.MOLECULAR,
                           state=_state("tagmented barcoded 3′ fragment", "Nextera S5/ME end",
                                        "transcript cDNA", "poly(T)", "UMI",
                                        "inner RT-An barcode", "RT handle", "outer Gn barcode",
                                        "Tru2 handle",
                                        inferred=("poly(T)", "UMI", "inner RT-An barcode",
                                                  "RT handle", "outer Gn barcode",
                                                  "Tru2 handle")))
    waste = graph.product("tn5_waste", "internal fragments without the Tru2 end",
                          MaterialKind.WASTE)
    graph.react("tagment", "Tagment the pooled amplicons", StepKind.TAGMENTATION,
                ((pooled, InputRole.TEMPLATE), (transposome, InputRole.STRUCTURAL),
                 (tn5, InputRole.CATALYST)), (tagged, waste))
    graph.discard(waste, "lacks the Tru2 site, so is not amplified by the Ix-Tru2 primer and "
                  "is not physically separated in the inspected Methods")
    forward = graph.input("library_forward", "Nextera S5XX i5 index primer", MaterialKind.OLIGO)
    reverse = graph.input("library_reverse", "Ix-Tru2 primer, substituted for S7XX (sequence "
                          "not printed)", MaterialKind.OLIGO)
    pcr_pol = graph.input("library_polymerase", "kit PCR polymerase", MaterialKind.CATALYST)
    library = graph.product("library", "MALBAC-DT 3′-end sequencing library",
                            MaterialKind.LIBRARY, state=final)
    graph.react("library_pcr", "Amplify barcoded 3′ ends with S5XX and Ix-Tru2", StepKind.PCR,
                ((tagged, InputRole.TEMPLATE), (forward, InputRole.STRUCTURAL),
                 (reverse, InputRole.STRUCTURAL), (pcr_pol, InputRole.CATALYST)), (library,),
                geometry=PCRGeometry.unavailable(
                    final, reason="Ix-Tru2 sequence and its Tru2 target are not printed."))
    graph.sequence(library, SequencingPlan(library=final, primers=tuple(primers)))
    reviewed = ["rna", "rt_primer", "reverse_transcriptase", "reverse_transcribe", "exo_i",
                "rt_bn", "exo_digest", "malbac_primer", "malbac_polymerase", "malbac",
                "tru2_gn_rt", "tru2_pcr", "ampure", "pool_cleanup", "transposome", "tn5",
                "tagment", "library_forward", "library_reverse", "library_polymerase",
                "library_pcr", "library"]
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()


def protocol_graph():
    return _local_graph(DEFINITION, SPEC.FINAL_LIBRARY, SPEC.SEQ_PRIMERS)


SPEC.protocol_graph = protocol_graph
