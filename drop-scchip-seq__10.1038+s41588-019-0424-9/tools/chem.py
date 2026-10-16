"""Source-backed Drop-scChIP-seq bead assembly and sequencing library."""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))

import illumina as il
import seqprimers as sp
from adapter_ligation import CohesiveDuplex, cohesive_duplex
from adapter_ligation import AdapterLigationGeometry
from chemdraw import (Construct, MolecularState, Scene, Segment, Workflow,
                      duplex_rows, feature, revcomp)
from pcr import PCRGeometry
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph,
                            SequencingPlan, StepKind)


TITLE = "Drop-scChIP-seq — droplet single-cell ChIP-seq"
NOTES = "01_drop-scchip-seq.html"
SOURCE = ('Defining source: <a href="https://doi.org/10.1038/s41588-019-0424-9">'
          'Grosselin et al., <i>Nature Genetics</i> (2019)</a>.')
SUMMARY = ("MNase-released nucleosomes receive a bead identity assembled by three "
           "96-way cohesive-end ligations. After immunoprecipitation, T7 transcription, "
           "random-primed reverse transcription and PCR produce a library in which "
           "Read 1 is genomic DNA and Read 2 traverses the bead barcode.")
CAVEAT = ("The drawings use the published A01 member from each 96-oligo barcode plate. "
          "The barcode is three plate identities, not a UMI. The final 12-bp linker and "
          "76-bp barcode dimensions follow Supplementary Figure 4; the source does not "
          "resolve individual bases retained at the PacI processing junction beyond "
          "those dimensions. The article Methods cite Supplementary Table 4 for SEQ1, "
          "SEQ2, SEQ4 and SEQ5 (local article lines 1880–1890, 1943–1947), but those "
          "oligos are in Supplementary Table 5 (workbook lines 4750–4762); Table 4 "
          "contains the three barcode plates. The article's Table 5 citation for "
          "barcode structure at line 1969 is likewise reversed. The model uses the "
          "workbook locations. "
          "The first four Read 2 cycles were dark (Methods lines 1948–1951), so "
          "diagram positions describe molecule order, not an uninterrupted base-called read.")

BC1 = feature("drop_scchip_bc1", "cell_barcode", "whitelist",
              group="drop_scchip_cell", part="round 1",
              whitelist="Grosselin et al. Supplementary Table 4, Index 1")
BC2 = feature("drop_scchip_bc2", "cell_barcode", "whitelist",
              group="drop_scchip_cell", part="round 2",
              whitelist="Grosselin et al. Supplementary Table 4, Index 2")
BC3 = feature("drop_scchip_bc3", "cell_barcode", "whitelist",
              group="drop_scchip_cell", part="round 3",
              whitelist="Grosselin et al. Supplementary Table 4, Index 3")
I7 = feature("drop_scchip_i7", "sample_index", "unknown")

# Supplementary Table 5 common oligos, written 5′→3′ without ordering syntax.
SEQ1_TOP = "TTAAGAATTTAATACGACTCACTATAGGGAGAGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT"
SEQ1_BOTTOM = "CGAAAGATCGGAAGAGCACACGTCTGAACTCCAGTCACTCTCCCTATAGTGAGTCGTATTAAATTCTTAA"
SEQ2_TOP = "CAACGTGATTGCTTGTGACTTAA"
SEQ2_BOTTOM = "TTAAGTCACAAGCAATCAC"
SEQ3_TOP = "GATACCGTCGAC"       # preceded by one of TTCG / TGAC / TCCC linkers
SEQ3_BOTTOM = "GTCGACGGTATC"
RT_PRIMER = "TACACGACGCTCTTCCGATCT" + "N" * 6
FORWARD_PCR = il.TRUSEQ_P5_FULL

# Supplementary Table 4, A01 representative of each 96-member plate.
ROUND_WRITTEN = (
    ("round 1", "TTCGGCTTGGTTGCGCTTGGTTGC", "GTCAGCAACCAAGCGCAACCAAGC", BC1),
    ("round 2", "TGACACTTGGCGCTACTTGGCGCT", "GGGAAGCGCCAAGTAGCGCCAAGT", BC2),
    ("round 3", "TCCCGCCATTGGACGCCATTGGAC", "GTTGGTCCAATGGCGTCCAATGGC", BC3),
)


def _plate_duplex(name: str, top: str, bottom: str, identity) -> CohesiveDuplex:
    return cohesive_duplex(name, top, bottom, core_name=f"{name} A01 identity",
                           core_feature=identity, top_mod5="/5Phos/",
                           bottom_mod5="/5Phos/")


SEQ1 = cohesive_duplex("SEQ1 bead flank", SEQ1_TOP, SEQ1_BOTTOM,
                       core_name="PacI half-site + T7 promoter + Read-2 flank",
                       top_mod5="/5PCBio/", bottom_mod3="/3C3/")
ROUNDS = tuple(_plate_duplex(*row) for row in ROUND_WRITTEN)
SEQ2 = cohesive_duplex("SEQ2 terminal flank", SEQ2_TOP, SEQ2_BOTTOM,
                       core_name="terminal common flank", top_mod5="/5Phos/",
                       bottom_mod5="/5Phos/")

# Compatibility is an assembly invariant, not a prose assertion.
_END_CHAIN = (SEQ1.bottom_overhang, ROUNDS[0].top_overhang,
              ROUNDS[0].bottom_overhang, ROUNDS[1].top_overhang,
              ROUNDS[1].bottom_overhang, ROUNDS[2].top_overhang,
              ROUNDS[2].bottom_overhang, SEQ2.top_overhang)
for upstream, downstream in zip(_END_CHAIN[::2], _END_CHAIN[1::2]):
    if revcomp(upstream) != downstream:
        raise ValueError(f"incompatible Drop-scChIP ligation: {upstream}/{downstream}")


def oligos():
    rows = [f"<b>SEQ1 top</b> — 5′-/5PCBio/{SEQ1_TOP}-3′",
            f"<b>SEQ1 bottom</b> — 5′-{SEQ1_BOTTOM}/3C3/-3′"]
    for name, top, bottom, _identity in ROUND_WRITTEN:
        rows.extend((f"<b>{name} A01 top</b> — 5′-/5Phos/{top}-3′",
                     f"<b>{name} A01 bottom</b> — 5′-/5Phos/{bottom}-3′"))
    rows.extend((f"<b>SEQ2 top</b> — 5′-/5Phos/{SEQ2_TOP}-3′",
                 f"<b>SEQ2 bottom</b> — 5′-/5Phos/{SEQ2_BOTTOM}-3′",
                 f"<b>Random RT primer</b> — 5′-{RT_PRIMER}-3′",
                 f"<b>Forward PCR primer</b> — 5′-{FORWARD_PCR}-3′",
                 "<b>Reverse PCR primer</b> — 5′-CAAGCAGAAGACGGCATACGAGAT[index]GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT-3′"))
    return tuple(rows)


def _barcode_segments(through: int = 3) -> list[Segment]:
    out = []
    for round_ in ROUNDS[:through]:
        if round_.top_overhang:
            out.append(Segment(f"{round_.name} incoming cohesive end",
                               round_.top_overhang))
        feat = next(s.feature for s in round_.scene.strands["top"].segs if s.feature)
        out.append(Segment(f"{round_.name} A01 identity", round_.paired_core,
                           "cbc", feature=feat))
    return out


def _assembled_barcode_scene(through: int) -> Scene:
    parts = [Segment("bead flank", SEQ1_TOP), *_barcode_segments(through)]
    sc = Scene(); sc.strand("barcode strand", parts,
                            label=f"bead barcode through round {through}",
                            mod5="/5PCBio/")
    previous = "bead flank"
    for round_ in ROUNDS[:through]:
        incoming = f"{round_.name} incoming cohesive end"
        sc.junction("barcode strand", previous, incoming, f"{round_.name} ligation")
        previous = f"{round_.name} A01 identity"
    sc.labels("barcode strand")
    return sc


def _completed_barcode_scene() -> Scene:
    parts = [Segment("bead flank", SEQ1_TOP), *_barcode_segments(),
             Segment("terminal CAAC cohesive end", SEQ2.top_overhang),
             Segment("terminal common flank", SEQ2.paired_core)]
    sc = Scene(); sc.strand("barcode strand", parts, label="completed bead barcode",
                            mod5="/5PCBio/")
    previous = "bead flank"
    for round_ in ROUNDS:
        incoming = f"{round_.name} incoming cohesive end"
        sc.junction("barcode strand", previous, incoming, f"{round_.name} ligation")
        previous = f"{round_.name} A01 identity"
    sc.junction("barcode strand", previous, "terminal CAAC cohesive end",
                "terminal-flank ligation")
    sc.note("barcode strand", "76-bp three-round barcode followed by the terminal common flank")
    sc.labels("barcode strand")
    return sc


def _barcoded_nucleosome_scene() -> Scene:
    parts = [Segment("T7 promoter", "TAATACGACTCACTATAGGGAGA"),
             Segment("Read-2 primer site", il.TRUSEQ_READ2, "r2"),
             *_barcode_segments(),
             Segment("terminal CAAC barcode boundary", SEQ2.top_overhang),
             Segment("published 12-bp linker", "GTGATTGCTTGT"),
             Segment("nucleosomal DNA", "X" * 40, placeholder=True)]
    sc = Scene(); sc.strand("template", parts, label="barcoded nucleosome template")
    sc.junction("template", "published 12-bp linker", "nucleosomal DNA",
                "barcode-to-nucleosome junction")
    sc.labels("template")
    return sc


def _rt_scene() -> Scene:
    target = [Segment("barcode-derived RNA", "N" * 28, placeholder=True),
              Segment("RNA landing site", "N" * 6, placeholder=True),
              Segment("RNA beyond landing site", "N" * 18, placeholder=True)]
    primer = [Segment("Read-1 partial arm", RT_PRIMER[:-6], "r1"),
              Segment("random hexamer", RT_PRIMER[-6:], placeholder=True)]
    sc = Scene(); sc.strand("IVT RNA", target, label="linearly amplified RNA")
    sc.anneal("RT primer", primer, to="IVT RNA",
              pair=("random hexamer", "RNA landing site"),
              unpaired=("Read-1 partial arm",), label="SEQ4 random RT primer")
    sc.arrow("RT primer", "reverse transcription")
    sc.labels("RT primer")
    return sc


# Supplementary Figure 4 draws the final amplified strand from P7 to P5.
P5_UNIQUE = il.P5[:-4]
FINAL_LIBRARY = Construct([
    Segment("P7", il.P7, "p7"),
    Segment("i7 sample index", "I" * 8, "cbc", placeholder=True, feature=I7),
    Segment("Read-2 primer site", il.TRUSEQ_READ2, "r2"),
    *_barcode_segments(),
    Segment("terminal CAAC barcode boundary", SEQ2.top_overhang),
    Segment("published 12-bp linker", "GTGATTGCTTGT"),
    Segment("nucleosomal DNA", "X" * 150, placeholder=True),
    Segment("Read-1 primer site reverse complement", revcomp(il.TRUSEQ_READ1), "r1"),
    Segment("P5 unique portion reverse complement", revcomp(P5_UNIQUE), "p5"),
], name="Drop-scChIP-seq library")

ACQUIRED_ROLES = ("Read 1", "Index 1 (i7)", "Read 2")
SEQ_PRIMERS = sp.standard_primers(FINAL_LIBRARY, ACQUIRED_ROLES)
READ_LENGTHS = {"Read 1": 50, "Index 1 (i7)": 8, "Read 2": 100}
READ_USES = {"Read 1": "include", "Index 1 (i7)": "demultiplex",
             "Read 2": "extract_then_include"}
FINAL_CAPTION = ("Source-oriented P7→P5 strand. Read 1 enters nucleosomal DNA; Read 2 "
                 "reads the 76-bp bead-barcode region, the 12-bp linker and then genomic DNA.")
SEQUENCING_INTRO = ("Standard TruSeq primer sites are inferred from the completed library. "
                    "The defining study acquired 50 cycles for genomic Read 1 and 100 cycles "
                    "for barcode-side Read 2, with its first four cycles dark; the i7 primer "
                    "is present in the exact PCR oligo.")


def workflow() -> Workflow:
    start_rows = [*SEQ1.scene.rows(),
                  *Scene.duplex([Segment("MNase-released nucleosomal DNA", "X" * 34,
                                         placeholder=True)],
                                label="MNase-released chromatin").rows()]
    wf = Workflow(MolecularState("Bead flank and MNase-released chromatin", tuple(start_rows)))
    for number, round_ in enumerate(ROUNDS, 1):
        wf.react(f"Split, ligate {round_.name} A01 representative, pool and resplit",
                 [*round_.scene.rows(), *_assembled_barcode_scene(number).rows()],
                 name=f"Bead barcode after {round_.name}",
                 note="Source-written strands determine both 4-nt 5′ cohesive ends; compatibility with adjacent rounds is checked during construction.")
    wf.react("Ligate the terminal SEQ2 common flank and cap incomplete products",
             [*SEQ2.scene.rows(), *_completed_barcode_scene().rows()],
             name="Completed bead barcode",
             note="SEQ3 blocking oligos carry a 5′ C3 spacer and 3′ inverted dT; they cap failed chains rather than contribute to complete products.")
    wf.react("Photorelease bead barcodes and barcode nucleosomes in fused droplets",
             _barcoded_nucleosome_scene().rows(), name="Barcoded nucleosomes",
             note="The final source schematic resolves a 76-bp barcode and 12-bp linker next to nucleosomal DNA. PacI processing closes and reopens the flanking half-sites.")
    wf.react("Immunoprecipitate the chosen histone mark and transcribe from T7",
             _barcoded_nucleosome_scene().rows(), name="Linearly amplified barcode–DNA RNA",
             note="Antibody enrichment selects the chromatin target; it does not add an inline target barcode. T7 IVT amplifies the recovered templates.")
    wf.react("Random-prime reverse transcription with SEQ4", _rt_scene().rows(),
             name="Random-primed cDNA",
             note="The six random bases are a priming tract, not a UMI. The 5′ arm supplies the partial Read-1 site.")
    wf.react("PCR with the exact P5 and indexed P7 primers",
             duplex_rows(FINAL_LIBRARY, label="sequencing library"),
             name="Final Illumina library",
             note="PCR completes the Read-1 site and adds P5, i7 and P7 while retaining the bead barcode beside the Read-2 site.")
    return wf


def protocol_graph() -> ProtocolGraph:
    source = "Defining source: https://doi.org/10.1038/s41588-019-0424-9"
    graph = ProtocolGraph("Drop-scChIP-seq", authority=source)
    chromatin = graph.input("chromatin", "MNase-released nucleosomal chromatin",
                            MaterialKind.BIOLOGICAL, primary=True)
    bead_flank = graph.input("bead_flank", "SEQ1 bead flank", MaterialKind.ADAPTER,
                             state=SEQ1.scene)
    current = bead_flank
    for number, round_ in enumerate(ROUNDS, 1):
        adapter = graph.input(f"round_{number}_adapter", f"barcode round {number} duplex",
                              MaterialKind.ADAPTER, state=round_.scene)
        assembled = graph.product(f"barcode_round_{number}",
                                  f"bead barcode after ligation round {number}",
                                  MaterialKind.MOLECULAR,
                                  state=_assembled_barcode_scene(number))
        graph.react(f"barcode_ligation_{number}",
                    f"Ligate barcode round {number}, pool and resplit", StepKind.LIGATION,
                    ((current, InputRole.TEMPLATE), (adapter, InputRole.STRUCTURAL)),
                    (assembled,),
                    geometry=AdapterLigationGeometry.known(
                        "cohesive", round_.scene,
                        statement="Supplementary Table 4 supplies the source-written cohesive duplex."))
        current = assembled
    terminal_flank = graph.input("terminal_flank", "SEQ2 terminal common flank",
                                 MaterialKind.ADAPTER, state=SEQ2.scene)
    completed = graph.product("completed_barcode", "completed three-round bead barcode",
                              MaterialKind.MOLECULAR, state=_completed_barcode_scene())
    graph.react("terminal_ligation", "Ligate the terminal SEQ2 common flank",
                StepKind.LIGATION,
                ((current, InputRole.TEMPLATE), (terminal_flank, InputRole.STRUCTURAL)),
                (completed,), geometry=AdapterLigationGeometry.known(
                    "cohesive", SEQ2.scene,
                    statement="Supplementary Table 5 supplies the terminal cohesive duplex."))
    droplet_condition = graph.input("droplet_condition", "bead–nucleosome droplet fusion",
                                    MaterialKind.CONDITION)
    barcoded = graph.product("barcoded_nucleosomes", "barcoded nucleosome templates",
                             MaterialKind.MOLECULAR, state=_barcoded_nucleosome_scene())
    graph.react("droplet_barcoding", "Photorelease bead barcodes and barcode nucleosomes",
                StepKind.ADAPTER_ATTACHMENT,
                ((chromatin, InputRole.TEMPLATE), (completed, InputRole.STRUCTURAL),
                 (droplet_condition, InputRole.CONDITION)), (barcoded,))
    antibody = graph.input("antibody", "histone-mark antibody", MaterialKind.SELECTION)
    selected = graph.product("selected_templates", "immunoprecipitated barcoded nucleosomes",
                             MaterialKind.MOLECULAR, state=_barcoded_nucleosome_scene())
    rejected = graph.product("unselected_templates", "unselected chromatin", MaterialKind.WASTE)
    graph.react("immunoprecipitation", "Immunoprecipitate the chosen histone mark",
                StepKind.AFFINITY_CAPTURE,
                ((barcoded, InputRole.TEMPLATE), (antibody, InputRole.SELECTION)),
                (selected, rejected))
    graph.discard(rejected, "chromatin not retained by the antibody is washed away")
    t7 = graph.input("t7_polymerase", "T7 RNA polymerase", MaterialKind.CATALYST)
    rna = graph.product("ivt_rna", "linearly amplified barcode–DNA RNA",
                        MaterialKind.MOLECULAR, state=_barcoded_nucleosome_scene())
    graph.react("ivt", "Transcribe selected templates from the T7 promoter",
                StepKind.TRANSCRIPTION,
                ((selected, InputRole.TEMPLATE), (t7, InputRole.CATALYST)), (rna,))
    rt_primer = graph.input("random_rt_primer", "SEQ4 random RT primer", MaterialKind.OLIGO,
                            state=_rt_scene())
    reverse_transcriptase = graph.input("reverse_transcriptase", "reverse transcriptase",
                                       MaterialKind.CATALYST)
    cdna = graph.product("cdna", "random-primed cDNA", MaterialKind.MOLECULAR,
                         state=_rt_scene())
    graph.react("reverse_transcription", "Random-prime reverse transcription with SEQ4",
                StepKind.REVERSE_TRANSCRIPTION,
                ((rna, InputRole.TEMPLATE), (rt_primer, InputRole.STRUCTURAL),
                 (reverse_transcriptase, InputRole.CATALYST)), (cdna,))
    rna_digestion = graph.input("rna_digestion", "RNA digestion after reverse transcription", MaterialKind.CONDITION)
    dna_only = graph.product("rna_depleted_cdna", "cDNA after source-described RNA digestion", MaterialKind.MOLECULAR)
    graph.react("digest_rna", "Digest RNA before PCR", StepKind.MAGIC,
                ((cdna, InputRole.TEMPLATE), (rna_digestion, InputRole.CONDITION)), (dna_only,),
                reason="The inspected Methods specify RNA digestion after RT but its exact reagent and products remain unassigned in this pass.")
    p5 = graph.input("p5_primer", "exact P5 PCR primer", MaterialKind.OLIGO)
    p7 = graph.input("p7_primer", "indexed P7 PCR primer", MaterialKind.OLIGO)
    polymerase = graph.input("pcr_polymerase", "PCR polymerase", MaterialKind.CATALYST)
    library = graph.product("library", FINAL_LIBRARY.name, MaterialKind.LIBRARY,
                            state=FINAL_LIBRARY)
    graph.react("library_pcr", "PCR with the exact P5 and indexed P7 primers", StepKind.PCR,
                ((dna_only, InputRole.TEMPLATE), (p5, InputRole.STRUCTURAL),
                 (p7, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)), (library,),
                geometry=PCRGeometry.unavailable(
                    FINAL_LIBRARY, reason="The final construct is exact, while this page does not expose both PCR-primer landing objects."))
    graph.sequence(library, SequencingPlan(
        library=FINAL_LIBRARY, primers=SEQ_PRIMERS, acquired_roles=ACQUIRED_ROLES,
        read_lengths=READ_LENGTHS, read_uses=READ_USES))
    graph.declare_reviewed_scope((
        "chromatin", "bead_flank", "round_1_adapter", "barcode_ligation_1",
        "round_2_adapter", "barcode_ligation_2", "round_3_adapter",
        "barcode_ligation_3", "terminal_flank", "terminal_ligation",
        "droplet_condition", "droplet_barcoding", "antibody", "immunoprecipitation",
        "t7_polymerase", "ivt", "random_rt_primer", "reverse_transcriptase",
        "reverse_transcription", "rna_digestion", "digest_rna", "p5_primer", "p7_primer", "pcr_polymerase",
        "library_pcr", "library"), source)
    return graph.require_publishable()
