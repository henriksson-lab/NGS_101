from chemdraw import Workflow
from types import SimpleNamespace
from rna_cell_assays import AdapterLigationGeometry, Construct, Definition, InputRole, MaterialKind, PCRGeometry, ProtocolGraph, SequencingPlan, StepKind, custom_library, single, token
from rna_cell_assays import feature

from rna_cell_assays import branches, chain
from chemdraw import Segment
from illumina import P5, P7
from protocol_graph import SourceClaim, SourceDiscrepancy

def lineage_workflow() -> Workflow:
    return chain("integrated scGESTALT transgene", single("lineage cassette", (
        token("DsRed transcript"), token("CRISPR target array"), token("poly(A)"))), (
        ("Induce Cas9 and guide pools at successive developmental times",
         single("edited lineage cassette", (token("DsRed transcript"),
                                              token("cumulative target-array indels"), token("poly(A)"))),
         "Repair-derived insertions and deletions form a permanent, cumulative lineage record."),
        ("Capture cassette RNA and transcriptome by inDrops",
         single("cell-indexed cDNA", (token("transcript or lineage cDNA"),
                                      token("cell-barcode part 1", tag="cbc", length=8),
                                      token("cell-barcode part 2", tag="cbc", length=8),
                                      token("UMI", tag="umi", length=6))),
         "The inDrops cell identifier is the join key between expression and lineage libraries."),
        ("IVT-amplify, then split fragmented and unfragmented RNA",
         branches("same indexed amplified RNA", (
             ("transcriptome", (token("fragmented amplified RNA"),)),
             ("lineage", (token("unfragmented lineage RNA"),)))),
         "Fragmented RNA makes the expression library; unfragmented RNA preserves the full target array."),
        ("Random-hexamer RT, two nested PCRs, then final R1-PCRix/R2-PCR indexing PCR",
         single("scGESTALT lineage library", (token("lineage barcode array"),
                                               token("cell-barcode part 1", tag="cbc", length=8),
                                               token("cell-barcode part 2", tag="cbc", length=8),
                                               token("UMI", tag="umi", length=6), token("library index", tag="cbc", length=8))),
         "The lineage and transcriptome libraries use the same multiplexing index so cells can be matched."),
    ))

DEFINITION = Definition("scGESTALT — single-cell lineage recording", "10.1038/s41596-018-0058-x",
        "Raj &amp; Gagnon, <i>Nature Protocols</i> (2018)", "01_scgestalt.html",
        "Cumulative CRISPR indels in an expressed target array record ancestry; inDrops captures that lineage transcript alongside the cell transcriptome and supplies their shared cell identifier.",
        "The lineage run (MiSeq 300-cycle kit, 20% PhiX, standard Illumina primers) uses Read 1 = 250 cycles for the edited array, Index 1 = 8 cycles for cell-barcode part 1, Index 2 = 8 cycles for library index, and Read 2 = 14 cycles for cell-barcode part 2 plus 6-nt UMI.")

d = DEFINITION
final = custom_library('scGESTALT lineage library', (token('lineage barcode array', length=44), token('cell-barcode part 1', tag='cbc', length=8, molecular_feature=feature('cell_barcode_1', 'cell_barcode', 'combinatorial', group='cell_barcode', part='part 1')), token('cell-barcode part 2', tag='cbc', length=8, molecular_feature=feature('cell_barcode_2', 'cell_barcode', 'combinatorial', group='cell_barcode', part='part 2')), token('UMI', tag='umi', length=6, molecular_feature=feature('umi', 'umi', 'random')), token('library index', tag='cbc', length=8, molecular_feature=feature('sample_index', 'sample_index', 'fixed'))))
primers = ()
reads = None
workflow_fn = lineage_workflow
SPEC = SimpleNamespace(
    TITLE=d.title, NOTES=d.notes,
    SOURCE=f'Defining source: <a href="https://doi.org/{d.doi}">{d.citation}</a>.',
    SUMMARY=d.summary, CAVEAT=d.caveat, FINAL_LIBRARY=final,
    SEQ_PRIMERS=primers, FINAL_CAPTION=d.endpoint, SEQUENCING_INTRO=d.endpoint,
    READ_LENGTHS=reads, SEQUENCING_ENDING=d.endpoint,
    SEQUENCING_UNAVAILABLE=("" if primers else "Exact sequencing-primer sequences are not fixed by the defining protocol or were not exposed in its accessible primary text; read entry and cycle roles are stated explicitly."),
    workflow=workflow_fn, oligos=lambda: oligo_rows(),
)


_CLEANUP_WASTE = "AMPure XP bead supernatant (primers, reagents and fragments below the cut-off)"


def _ampure(graph, ident, label, material, ratio, *, state=None, kind=MaterialKind.MOLECULAR,
            kept_label=None, source=None):
    """An AMPure XP bead cleanup with the retained eluate and the removed supernatant."""
    beads = graph.input(ident + "_beads", f"AMPure XP beads, {ratio}", MaterialKind.SELECTION,
                        source=source)
    kept = graph.product(ident + "_eluate", kept_label or label + " (eluate)", kind,
                         state=material.state if state is None else state, source=source)
    waste = graph.product(ident + "_supernatant", _CLEANUP_WASTE, MaterialKind.WASTE,
                          source=source)
    graph.react(ident, label, StepKind.CLEANUP,
                ((material, InputRole.TEMPLATE), (beads, InputRole.SELECTION)),
                (kept, waste), changes_molecule=False, source=source)
    graph.discard(waste, "supernatant removed from the magnet-held AMPure XP beads")
    return kept


def _oligo(name, *parts):
    return Construct([Segment(label, bases, tag, feature=feat) for label, bases, tag, feat in parts],
                     name=name)


# Raj & Gagnon 2018 Nat. Protoc. PMC6279253 Reagents primer table (JATS .txt lines 263-273)
# and Supplementary Table 2 (MOESM1 .txt lines 88-112).  GP12 prints nine random N.
GP6 = _oligo("inDrops_GP6", ("DsRed-transgene forward", "GAGGACTACACCATCGTGGAG", None, None))
PE1SA = _oligo("inDrops_PE1Sa", ("PE1Sa (inDrops adaptor overlap)",
                                 "CTCTTTCCCTACACGACGCTGGGTGTCGGGTGCAG", None, None))
GP12 = _oligo("inDrops_GP12", ("GP12 Nextera-type 5' arm", "TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG", None, None),
              ("9 random bases (diversity)", "NNNNNNNNN", None, None),
              ("barcode-array forward", "TCGAGCTCAAGCTTCGGAC", None, None))
PE1SB = _oligo("inDrops_PE1Sb", ("PE1Sb (inDrops adaptor overlap)", "CTCTTTCCCTACACGACGCT", None, None))
R1_PCR_IX1 = _oligo("R1-PCR-ix1", ("P5", P5, "p5", None),
                    ("library index ix1", "CTCTCTAT", "cbc",
                     feature("sample_index", "sample_index", "fixed")),
                    ("R1-PCR 3' end", "TCGTCGGCAGCGTC", None, None))
R2_PCR = _oligo("R2-PCR", ("P7", P7, "p7", None), ("PE1 overlap", "GGGTGTCGGGTGCAG", None, None))
_TABLE = "Raj &amp; Gagnon 2018, Reagents primer table"
OLIGO_SOURCES = ((GP6, _TABLE), (PE1SA, _TABLE), (GP12, _TABLE), (PE1SB, _TABLE),
                 (R1_PCR_IX1, "Supplementary Table 2 (first of 24 R1-PCR-ix members)"),
                 (R2_PCR, _TABLE))


def scgestalt_graph(definition, final: Construct) -> ProtocolGraph:
    """scGESTALT graph with matched transcriptome and lineage-library outcomes."""
    source = f"{definition.citation}; doi:{definition.doi}"
    graph = ProtocolGraph(definition.title, authority=source)
    cell = graph.input(
        "transgenic_cell", "cell with integrated expressed lineage cassette",
        MaterialKind.BIOLOGICAL,
        state=single("lineage cassette and transcriptome", (
            token("DsRed transcript"), token("CRISPR target array"),
            token("poly(A)"), token("cellular transcriptome"))), primary=True)
    guides = graph.input("guide_pools", "sgRNAs 1–4 injected at one cell; U6 sgRNAs 5–9 transgenic",
                         MaterialKind.OLIGO)
    cas9 = graph.input("cas9", "injected Cas9 protein, then heat-shock-induced hsp:Cas9",
                       MaterialKind.CATALYST)
    edited_cell = graph.product(
        "edited_cell", "cell with cumulatively edited lineage cassette",
        MaterialKind.MOLECULAR,
        state=single("edited lineage cassette and transcriptome", (
            token("DsRed transcript"), token("cumulative target-array indels"),
            token("poly(A)"), token("cellular transcriptome"))))
    graph.react("edit_lineage", "Induce cumulative lineage-cassette editing",
                StepKind.GENOME_EDITING,
                ((cell, InputRole.TEMPLATE), (guides, InputRole.STRUCTURAL),
                 (cas9, InputRole.CATALYST)), (edited_cell,))
    host_rnap = graph.input("host_rna_polymerase", "cellular RNA polymerase (heat-shock-induced cassette)",
                            MaterialKind.CATALYST)
    cell_rna = graph.product(
        "cellular_rna", "polyadenylated edited-cassette RNA plus cellular mRNA",
        MaterialKind.MOLECULAR,
        state=single("polyadenylated cellular RNA", (
            token("edited lineage transcript"), token("cellular transcripts"))))
    graph.react("transcribe_cassette", "Transcribe and polyadenylate the edited cassette",
                StepKind.TRANSCRIPTION,
                ((edited_cell, InputRole.TEMPLATE),
                 (host_rnap, InputRole.CATALYST)), (cell_rna,))

    indrops_primer = graph.input("indrops_capture_primer",
                                 "inDrops V3 hydrogel-bead primer (T7–PE1–cell ID–UMI–oligo(dT))",
                                 MaterialKind.OLIGO,
                                 source=source + "; Fig. 3b legend; bases in inDrops ref. 73, not printed")
    rtase = graph.input("reverse_transcriptase", "inDrops in-droplet reverse transcriptase (ref. 73)",
                        MaterialKind.CATALYST)
    first_strand = graph.product(
        "barcoded_first_strand", "cell-indexed lineage and transcriptome cDNA",
        MaterialKind.MOLECULAR,
        state=single("cell-indexed cDNA", (
            token("transcript or lineage cDNA"),
            token("cell-barcode part 1", tag="cbc", length=8),
            token("cell-barcode part 2", tag="cbc", length=8),
            token("UMI", tag="umi", length=6))))
    graph.react("reverse_transcription", "Capture lineage RNA and transcriptome by inDrops",
                StepKind.REVERSE_TRANSCRIPTION,
                ((cell_rna, InputRole.TEMPLATE),
                 (indrops_primer, InputRole.STRUCTURAL),
                 (rtase, InputRole.CATALYST)), (first_strand,))
    second_primer = graph.input("second_strand_primer", "inDrops second-strand synthesis reagents (ref. 73)",
                                MaterialKind.OLIGO, inferred=True)
    second_polymerase = graph.input("second_strand_polymerase", "second-strand polymerase",
                                    MaterialKind.CATALYST)
    duplex = graph.product("barcoded_cdna_duplex", "barcoded second-strand cDNA",
                           MaterialKind.MOLECULAR, state=first_strand.state)
    graph.react("second_strand", "Synthesize the second cDNA strand", StepKind.EXTENSION,
                ((first_strand, InputRole.TEMPLATE),
                 (second_primer, InputRole.STRUCTURAL),
                 (second_polymerase, InputRole.CATALYST)), (duplex,))
    t7 = graph.input("t7_rna_polymerase", "T7 RNA polymerase", MaterialKind.CATALYST)
    amplified_rna = graph.product("amplified_rna", "IVT-amplified barcoded RNA",
                                  MaterialKind.MOLECULAR, state=duplex.state)
    graph.react("ivt", "Linearly amplify barcoded cDNA by T7 IVT", StepKind.TRANSCRIPTION,
                ((duplex, InputRole.TEMPLATE), (t7, InputRole.CATALYST)),
                (amplified_rna,))
    purified_rna = _ampure(graph, "post_ivt_cleanup",
                           "Purify IVT RNA (1.3× AMPure XP; 25 µl RNA elution buffer)",
                           amplified_rna, "1.3×")
    transcriptome_rna = graph.product("transcriptome_rna", "9 µl RNA for transcriptome library",
                                      MaterialKind.MOLECULAR,
                                      state=single("transcriptome branch", (
                                          token("amplified transcriptome RNA"),)))
    lineage_rna = graph.product("lineage_rna", "remaining unfragmented lineage-cassette RNA",
                                MaterialKind.MOLECULAR,
                                state=single("lineage branch", (
                                    token("unfragmented lineage RNA"),)))
    graph.react("split_ivt", "Split purified IVT RNA: 9 µl to transcriptome, rest kept unfragmented",
                StepKind.SPLIT, ((purified_rna, InputRole.TEMPLATE),),
                (transcriptome_rna, lineage_rna))

    # Transcriptome branch (Steps 83-84 defer to inDrops ref. 73; Fig. 3b legend summarizes it).
    fragmentation = graph.input("rna_fragmentation_reagent", "chemical RNA fragmentation",
                                MaterialKind.CONDITION)
    fragmented = graph.product("fragmented_transcriptome_rna", "fragmented amplified RNA",
                               MaterialKind.MOLECULAR, state=transcriptome_rna.state)
    graph.react("fragment_transcriptome", "Chemically fragment the transcriptome branch",
                StepKind.FRAGMENTATION,
                ((transcriptome_rna, InputRole.TEMPLATE),
                 (fragmentation, InputRole.CONDITION)), (fragmented,))
    transcriptome_rt_primer = graph.input("transcriptome_rt_primer",
                                          "random-hexamer RT primer with PE2 adaptor (ref. 73)",
                                          MaterialKind.OLIGO, inferred=True,
                                          source=source + "; Fig. 3b legend; bases not printed")
    expression_cdna = graph.product("expression_cdna", "PE2-tagged transcriptome ssDNA",
                                    MaterialKind.MOLECULAR, state=fragmented.state,
                                    inferred=True)
    graph.react("transcriptome_rt", "Random-hexamer–PE2 RT of the fragmented RNA",
                StepKind.REVERSE_TRANSCRIPTION,
                ((fragmented, InputRole.TEMPLATE),
                 (transcriptome_rt_primer, InputRole.STRUCTURAL),
                 (rtase, InputRole.CATALYST)), (expression_cdna,))
    kapa = graph.input("kapa_hifi", "2× KAPA HiFi HotStart PCR mix", MaterialKind.CATALYST)
    r1 = graph.input("r1_pcrix", "R1-PCRix (24 indexed members; ix1 shown)", MaterialKind.OLIGO,
                     state=R1_PCR_IX1, source=source + "; Supplementary Table 2")
    r2 = graph.input("r2_pcr", "R2-PCR", MaterialKind.OLIGO, state=R2_PCR,
                     source=source + "; Reagents primer table")
    expression_final = custom_library("scGESTALT matched transcriptome library", (
        token("cell barcode", tag="cbc", length=16), token("UMI", tag="umi", length=6),
        token("transcript cDNA", length=44)))
    expression_library = graph.product(
        "expression_library", "matched scGESTALT transcriptome library",
        MaterialKind.LIBRARY, state=expression_final, inferred=True)
    graph.react("expression_pcr", "PCR the transcriptome branch with the same R1-PCRix/R2-PCR pair",
                StepKind.PCR,
                ((expression_cdna, InputRole.TEMPLATE),
                 (r1, InputRole.STRUCTURAL), (r2, InputRole.STRUCTURAL),
                 (kapa, InputRole.CATALYST)), (expression_library,),
                geometry=PCRGeometry.unavailable(
                    expression_final,
                    reason="The PE2 RT-primer and bead PE1 bases are in inDrops ref. 73, not this protocol, so primer landing is not assembled."))
    graph.sequence(expression_library, SequencingPlan(
        limitation="NextSeq 75-cycle high output, standard Illumina primers.",
        read_lengths={"Read 1": 61, "Index 1": 8, "Index 2": 8, "Read 2": 14}))

    # Lineage branch, Steps 85-99.
    lineage_rt_primer = graph.input("lineage_random_hexamer", "random hexamers (3.75 µM)",
                                    MaterialKind.OLIGO)
    primescript = graph.input("primescript", "PrimeScript reverse transcriptase + RNaseOut",
                              MaterialKind.CATALYST)
    lineage_cdna = graph.product("lineage_cdna", "full target-array lineage cDNA",
                                 MaterialKind.MOLECULAR, state=lineage_rna.state)
    graph.react("lineage_rt", "Random-hexamer RT of unfragmented RNA (70 °C 3 min; 30 °C 10 min, 42 °C 1 h, 70 °C 15 min)",
                StepKind.REVERSE_TRANSCRIPTION,
                ((lineage_rna, InputRole.TEMPLATE),
                 (lineage_rt_primer, InputRole.STRUCTURAL),
                 (primescript, InputRole.CATALYST)), (lineage_cdna,))
    lineage_cdna = _ampure(graph, "lineage_rt_cleanup", "Purify RT product (1.2× AMPure XP)",
                           lineage_cdna, "1.2×")
    q5 = graph.input("q5", "Q5 DNA polymerase with Q5 enhancer", MaterialKind.CATALYST)
    nested1_forward = graph.input("lineage_nested1_forward", "inDrops_GP6 (DsRed)",
                                  MaterialKind.OLIGO, state=GP6)
    nested1_reverse = graph.input("lineage_nested1_reverse", "inDrops_PE1Sa",
                                  MaterialKind.OLIGO, state=PE1SA)
    nested = graph.product("lineage_nested_amplicon", "long DsRed–array–cell-ID amplicon",
                           MaterialKind.MOLECULAR, state=lineage_cdna.state)
    graph.react("lineage_nested_pcr1", "PCR round 1: GP6 + PE1Sa, Q5, 15 cycles (61 °C anneal)",
                StepKind.PCR,
                ((lineage_cdna, InputRole.TEMPLATE),
                 (nested1_forward, InputRole.STRUCTURAL),
                 (nested1_reverse, InputRole.STRUCTURAL),
                 (q5, InputRole.CATALYST)), (nested,),
                geometry=PCRGeometry.unavailable(
                    final, reason="The DsRed/array cassette and inDrops PE1 bead bases are not printed, so primer landing is not assembled."))
    nested = _ampure(graph, "pcr1_cleanup", "Purify round-1 PCR (0.6× AMPure XP)", nested, "0.6×")
    phusion = graph.input("phusion", "Phusion High-Fidelity DNA polymerase (HF buffer)",
                          MaterialKind.CATALYST)
    nested2_forward = graph.input("lineage_nested2_forward", "inDrops_GP12 (N9 diversity)",
                                  MaterialKind.OLIGO, state=GP12)
    nested2_reverse = graph.input("lineage_nested2_reverse", "inDrops_PE1Sb",
                                  MaterialKind.OLIGO, state=PE1SB)
    short = graph.product("lineage_short_amplicon", "short array–cell-ID amplicon with GP12 arm",
                          MaterialKind.MOLECULAR, state=nested.state)
    graph.react("lineage_nested_pcr2", "PCR round 2: GP12 + PE1Sb, Phusion (60 °C anneal; cycle count ambiguous)",
                StepKind.PCR,
                ((nested, InputRole.TEMPLATE),
                 (nested2_forward, InputRole.STRUCTURAL),
                 (nested2_reverse, InputRole.STRUCTURAL),
                 (phusion, InputRole.CATALYST)), (short,),
                geometry=PCRGeometry.unavailable(
                    final, reason="The barcode-array and PE1 bead bases are not printed, so primer landing is not assembled."))
    short = _ampure(graph, "pcr2_cleanup", "Purify round-2 PCR (0.6× AMPure XP)", short, "0.6×")
    lineage_library = graph.product("lineage_library", "scGESTALT lineage library",
                                    MaterialKind.LIBRARY, state=final, inferred=True)
    graph.react("lineage_final_pcr", "Final PCR on 1:10 diluted eluate: R1-PCRix + R2-PCR, KAPA HiFi (2 cycles 55 °C, 6–8 cycles 65 °C)",
                StepKind.PCR,
                ((short, InputRole.TEMPLATE),
                 (r1, InputRole.STRUCTURAL), (r2, InputRole.STRUCTURAL),
                 (kapa, InputRole.CATALYST)), (lineage_library,),
                geometry=PCRGeometry.unavailable(
                    final, reason="R1-PCRix overlaps the GP12 arm and R2-PCR the PE1 overlap, but the bead PE1/cell-ID bases are in inDrops ref. 73."))
    lineage_library = _ampure(graph, "final_cleanup", "Purify final library (0.6× AMPure XP)",
                              lineage_library, "0.6×", kind=MaterialKind.LIBRARY,
                              kept_label="sequencing-ready scGESTALT lineage library")
    graph.sequence(lineage_library, SequencingPlan(
        limitation=definition.endpoint,
        read_lengths={"Read 1": 250, "Index 1": 8, "Index 2": 8, "Read 2": 14}))
    graph.annotate(SourceDiscrepancy(
        "r1_pcrix", "Which supplementary table lists the R1-PCRix primers",
        (SourceClaim("Reagents primer table, R1-PCRix row (PMC6279253)",
                     "Index sequences: see Supplementary Table 2."),
         SourceClaim("Reagent setup, 'R1-PCRix/R2-PCR primer mix'",
                     "See Supplementary Table 1 for a list of 24 R1-PCRix primer sequences."),
         SourceClaim("Supplementary Information (MOESM1)",
                     "Table 1 lists sgRNA target sites; Table 2 lists R1-PCR-ix1 to ix24.")),
        "The model takes R1-PCRix members from Supplementary Table 2 and shows ix1 as one example.",
        "The supplement itself puts the 24 R1-PCR-ix sequences in Table 2.",
        "Only the table label is in question; no primer sequence differs between the passages."))
    graph.annotate(SourceDiscrepancy(
        "lineage_nested_pcr2", "PCR round-2 cycle numbering",
        (SourceClaim("Procedure, PCR round 2 thermocycler table",
                     "Cycle rows read '1', '2–9', then '16' (final extension) and '17' (hold)."),
         SourceClaim("Procedure, PCR round 1 thermocycler table",
                     "The analogous table reads '1', '2–16', '17', '18' (15 amplification cycles).")),
        "The model gives no round-2 cycle count; it states only the 60 °C anneal.",
        "Rows 2–9 imply 8 cycles, but '16' for the next row implies 2–15 (14 cycles).",
        "Which count the authors intended is not stated; Step 96 says to lower input or cycles if yield is high."))
    graph.declare_reviewed_scope(tuple(graph.materials) + tuple(graph.steps), source)
    return graph.require_publishable()


def oligo_rows():
    """Key-oligo rows rendered from the graph's source-transcribed oligo states."""
    return tuple(f"<b>{c.name}</b> ({where}): 5′-{c.top()}-3′" for c, where in OLIGO_SOURCES)


def protocol_graph():
    d = DEFINITION
    final = SPEC.FINAL_LIBRARY
    primers = SPEC.SEQ_PRIMERS
    key = 'scgestalt'
    return scgestalt_graph(d, final)
SPEC.protocol_graph = protocol_graph
