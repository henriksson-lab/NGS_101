"""FLASH-seq standard branch, declared from its own article and author protocol."""
from __future__ import annotations

import os
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "lib"))

import illumina as il
import nextera as nx
import seqprimers as sp
from chemdraw import Construct, Scene, Segment, feature, oligo, panel, revcomp
from checks import Check, Source, have
from page import head, info
from pcr import PCRGeometry, PCRPrimer
from protocol_graph import (AnnotationSlot, Highlight, InputRole, MaterialKind,
                            ProtocolGraph, SequencingPlan, SourceClaim,
                            SourceDiscrepancy, StepKind)
from protocol_views import assemble_page, WorkflowSection

PAGE_FILE = "flash-seq.html"
SOURCE = "Hahaut et al. 2022 doi:10.1038/s41587-022-01312-3, Methods and Supplementary Fig.1/Note2"
HANDLE = "AAGCAGTGGTATCAACGCAGAGTAC"
DT_LENGTH = 30
RT_SEQUENCE = HANDLE + "T" * DT_LENGTH + "VN"
TSO_AS_WRITTEN = HANDLE + "rGrGrG"
TSO_SEQUENCE = HANDLE + "GGG"
# Own author table, article-linked standard v3 pp5-6. An example, not a well map.
I5_S502 = "CTCTCTAT"
I7_N714 = "TCATGAGC"
INDEX_FORWARD = il.P5 + I5_S502 + nx.S5
INDEX_REVERSE = il.P7 + I7_N714 + nx.S7
SEQ_PRIMERS = tuple(sp.NEXTERA[k] for k in ("R1", "I1", "I2"))
READ_LENGTHS = {"Read 1": 75, "Index 1 (i7)": 8, "Index 2 (i5)": 8}


def rt_primer():
    return Construct([Segment("FS handle", HANDLE), Segment("dT30", "T" * DT_LENGTH),
                      Segment("VN anchor", "VN", placeholder=True)], "FS dT30VN")


def tso():
    return Construct([Segment("FS handle", HANDLE), Segment("RNA G tail", "GGG")], "FS TSO")


def amplified_cdna():
    """Primer-derived duplex reference; RNA primer residues use the DNA alphabet."""
    return Construct([*tso().segments, Segment("cDNA insert", "X" * 32, placeholder=True),
                      Segment("anchor complement", revcomp("VN"), placeholder=True),
                      Segment("polyA", "A" * DT_LENGTH),
                      Segment("FS handle reverse complement", revcomp(HANDLE))],
                     "FS semi-suppressive amplification product reference")


def _commercial_me(name, sequence):
    # The author prints these arms for its separate in-house option. Their complete
    # commercial transposome/junction assignment has not been source-checked here.
    return Segment(name, sequence, "me", inferred=True)


def final_library():
    library = Construct([
        Segment("P5", il.P5, "p5"),
        Segment("illustrative S502 i5", I5_S502, "cbc",
                feature=feature("sample_index_i5", "sample_index", "fixed")),
        Segment("s5", nx.S5, "s5"), _commercial_me("canonical ME", nx.ME),
        Segment("cDNA insert", "X" * 32, placeholder=True),
        _commercial_me("canonical ME reverse complement", nx.ME_RC),
        Segment("s7 reverse complement", nx.S7_RC, "s7"),
        Segment("illustrative N714 i7 reverse complement", revcomp(I7_N714), "cbc",
                feature=feature("sample_index_i7", "sample_index", "fixed")),
        Segment("P7 reverse complement", il.P7_RC, "p7")],
        "FS commercial-XT endpoint; illustrative S502/N714, inferred ME assignment")
    index_geometry(library)
    problems = sp.verify(library, SEQ_PRIMERS, required_roles=tuple(p.role for p in SEQ_PRIMERS))
    if problems:
        raise ValueError("; ".join(problems))
    return library


def index_geometry(library):
    return PCRGeometry.exact(library,
        PCRPrimer("S502", INDEX_FORWARD, anneal_nt=len(INDEX_FORWARD)),
        PCRPrimer("N714", INDEX_REVERSE, anneal_nt=len(INDEX_REVERSE)),
        statement="Complete source-table index primers face inward on this endpoint reference; exact input-to-product transposome geometry remains unassigned.")


def protocol_graph():
    graph = ProtocolGraph("FLASH-seq — standard commercial Nextera XT branch", authority=SOURCE)
    rna = graph.input("rna", "Single-cell poly(A) RNA", MaterialKind.MOLECULAR,
        state=Construct([Segment("transcript", "X" * 32, placeholder=True),
                         Segment("polyA", "A" * DT_LENGTH)]), primary=True)
    dt = graph.input("dt", "Biotinylated FS dT30VN", MaterialKind.OLIGO, state=rt_primer())
    rt = graph.input("rt", "SuperScript IV reverse transcriptase", MaterialKind.CATALYST)
    first = graph.product("first", "First-strand cDNA; terminal cytidine addition not resolved base by base", MaterialKind.MOLECULAR)
    graph.react("reverse_transcription", "Reverse-transcribe in the combined RT-PCR mix", StepKind.REVERSE_TRANSCRIPTION,
                ((rna, InputRole.TEMPLATE), (dt, InputRole.STRUCTURAL), (rt, InputRole.CATALYST)), (first,))
    switch = graph.input("tso", "Biotinylated FS TSO ending in rGrGrG", MaterialKind.OLIGO, state=tso())
    switched = graph.product("switched", "Template-switched cDNA; terminal pairing geometry unassigned", MaterialKind.MOLECULAR)
    graph.react("template_switch", "Template-switch with the FS RNA-G tail", StepKind.TEMPLATE_SWITCHING,
                ((first, InputRole.TEMPLATE), (switch, InputRole.STRUCTURAL), (rt, InputRole.CATALYST)), (switched,))
    kapa = graph.input("kapa", "KAPA HiFi HotStart in the same RT-PCR mix", MaterialKind.CATALYST)
    preamp = amplified_cdna()
    amplified = graph.product("amplified", "Semi-suppressive amplified FS cDNA", MaterialKind.MOLECULAR, state=preamp)
    geometry = PCRGeometry.exact(preamp,
        PCRPrimer("FS TSO, DNA-alphabet reference", TSO_SEQUENCE, anneal_nt=len(TSO_SEQUENCE)),
        PCRPrimer("FS dT30VN", RT_SEQUENCE, anneal_nt=len(RT_SEQUENCE)),
        statement="Source Fig.1 uses FS TSO and FS-dT for amplification and crosses out ISPCR. Sites are checked on a primer-derived product reference; this does not prove RNA-primer extension or the complete first-cycle mechanism.")
    graph.react("preamp", "Semi-suppressive PCR in the same vessel; no separate ISPCR addition", StepKind.PCR,
                ((switched, InputRole.TEMPLATE), (switch, InputRole.STRUCTURAL),
                 (dt, InputRole.STRUCTURAL), (kapa, InputRole.CATALYST)), (amplified,), geometry=geometry)
    beads = graph.input("beads", "SeraMag cleanup beads", MaterialKind.SELECTION)
    clean = graph.product("clean", "Purified cDNA for QC and normalization", MaterialKind.MOLECULAR, state=preamp)
    rejected = graph.product("cleanup_waste", "Material outside the retained cDNA fraction", MaterialKind.WASTE)
    graph.react("cleanup", "Clean up, quantify and normalize standard FS cDNA", StepKind.CLEANUP,
                ((amplified, InputRole.TEMPLATE), (beads, InputRole.SELECTION)), (clean, rejected))
    graph.discard(rejected, "Removed in the standard-FS cDNA cleanup; no quantitative recovery is asserted.")
    adapter = graph.input("transposome", "Commercial Nextera XT ATM transposome", MaterialKind.ADAPTER)
    tn5 = graph.input("tn5", "Tn5 activity in ATM", MaterialKind.CATALYST)
    tagged = graph.product("tagged", "Tagmented cDNA; exact commercial end state unassigned", MaterialKind.MOLECULAR)
    graph.react("tagment", "Tagment using the commercial Nextera XT option", StepKind.TAGMENTATION,
                ((clean, InputRole.TEMPLATE), (adapter, InputRole.STRUCTURAL), (tn5, InputRole.CATALYST)), (tagged,))
    nt = graph.input("nt", "Nextera NT neutralization buffer", MaterialKind.CONDITION)
    neutralized = graph.product("neutralized", "Neutralized tagmentation mixture", MaterialKind.MOLECULAR)
    graph.react("neutralize", "Neutralize with NT", StepKind.QUENCHING,
                ((tagged, InputRole.TEMPLATE), (nt, InputRole.CONDITION)), (neutralized,))
    npm = graph.input("npm", "Nextera PCR Master mix (NPM)", MaterialKind.CATALYST)
    forward = graph.input("s502", "Biotinylated S502 primer with terminal phosphorothioate", MaterialKind.OLIGO,
                          state=Construct([Segment("S502", INDEX_FORWARD)]))
    reverse = graph.input("n714", "Biotinylated N714 primer with terminal phosphorothioate", MaterialKind.OLIGO,
                          state=Construct([Segment("N714", INDEX_REVERSE)]))
    filled = graph.product("filled", "Initial extension intermediate; exact ends unassigned", MaterialKind.MOLECULAR)
    graph.react("gap_fill", "Initial extension in the NPM/index-primer mixture", StepKind.EXTENSION,
                ((neutralized, InputRole.TEMPLATE), (forward, InputRole.STRUCTURAL),
                 (reverse, InputRole.STRUCTURAL), (npm, InputRole.CATALYST)), (filled,))
    final = final_library()
    library = graph.product("library", "Illustrative indexed FS library", MaterialKind.LIBRARY, state=final)
    graph.react("index_pcr", "INFERRED — primer-derived endpoint with canonical commercial ME assignment", StepKind.PCR,
                ((filled, InputRole.TEMPLATE), (forward, InputRole.STRUCTURAL),
                 (reverse, InputRole.STRUCTURAL), (npm, InputRole.CATALYST)), (library,), geometry=index_geometry(final))
    pooled = graph.product("pooled_library", "Pooled and bead-purified FS sequencing library", MaterialKind.LIBRARY, state=final)
    final_waste = graph.product("final_waste", "Material removed in pooled-library cleanup", MaterialKind.WASTE)
    graph.react("pool_cleanup", "Pool cells and clean the indexed library", StepKind.CLEANUP,
                ((library, InputRole.TEMPLATE), (beads, InputRole.SELECTION)), (pooled, final_waste))
    graph.discard(final_waste, "Removed during the source pooled-library bead cleanup.")
    graph.sequence(pooled, SequencingPlan(library=final, primers=SEQ_PRIMERS,
        acquired_roles=tuple(p.role for p in SEQ_PRIMERS), read_lengths=READ_LENGTHS))
    graph.declare_reviewed_scope(tuple(graph.steps) + ("rna", "dt", "rt", "tso", "kapa", "beads", "transposome", "tn5", "nt", "npm", "s502", "n714", "pooled_library"), SOURCE)
    graph.annotate(Highlight("preamp", "FLASH-seq removes the separate ISPCR reagent", "Original Supplementary Fig.1 visibly crosses out ISPCR and labels the amplifying oligos TSO and FS-dT. The article and detailed standard recipe list the same two oligos in one combined SSIV/KAPA reaction."))
    graph.annotate(SourceDiscrepancy(
        "tagment", "Commercial tagmentation recipe has inconsistent volumes",
        (SourceClaim("Author standard protocol v3, p21, step10.1 reagent table",
                     "ATM 0.250 µl + TD 2.000 µl + water 0.500 µl sums to 2.750 µl."),
         SourceClaim("Same v3 protocol, p22, continuation of step10.1",
                     "The table declares 3.000 µl total and instructs dispensing 3 µl mix before adding 1 µl cDNA."),
         SourceClaim("Defining article, standard tagmentation Methods; later author v4 pp21–22",
                     "The article specifies 0.5 µl ATM + 2 µl TD + 1 µl cDNA, with water to 4 µl. The later v4 table retains the 0.250/2.000/0.500 entries and 3.000 µl total.")),
        "The graph follows the commercial XT/NT/NPM reagent route without assigning reagent volumes or silently correcting a row.",
        "Both readers inspected the historical v3 original, including its printed table and continuation; the article independently identifies this route.",
        "The protocol table does not establish which component supplies its missing 0.250 µl. A volume or concentration cannot be selected from that inconsistent table alone."))
    return graph.require_publishable()


def render_page():
    graph = protocol_graph()
    ordered = [oligo("FS dT30VN", rt_primer().segments, mods="/5Biosg/"),
        oligo("FS TSO", [Segment("FS handle", HANDLE), Segment("RNA G tail", "rGrGrG", placeholder=True)], mods="/5Biosg/")]
    for name, sequence in (("S502, illustrative", INDEX_FORWARD), ("N714, illustrative", INDEX_REVERSE)):
        ordered.append(oligo(name, [Segment("primer", sequence[:-1]), Segment("terminal bond", "*" + sequence[-1], placeholder=True)], mods="/5Biosg/"))
    return assemble_page(graph, [head("FLASH-seq library chemistry"),
        '<div class="wrap" data-protocol-complete="true"><h1>FLASH-seq — standard full-length single-cell RNA</h1>',
        info('Defining source: <a href="https://doi.org/10.1038/s41587-022-01312-3">Hahaut et al. (2022)</a>; index-order example from its linked <a href="https://doi.org/10.17504/protocols.io.kxygxzkrwv8j/v3">author standard protocol v3 (2022)</a>.'),
        info("This page selects standard FS with cDNA cleanup and the commercial Nextera XT/NT/NPM route. FS-LA, FS-UMI, in-house Tn5/KAPA enrichment and plexWell are separate published options. The recovered original v3 supplies this ordering example; it is not an experimental well assignment."),
        info("INFERRED — canonical Nextera ME regions in the commercial endpoint are model assignments, marked in the sequence. Their complete commercial transposome and intermediate end geometry has not been independently source-checked here. Exact printed FS and index-primer inputs remain distinct from that inference."),
        info("The tagmentation model follows the represented library route; other commercial product structures and fates remain unassigned. A single represented class does not assert complete recovery or physical discard of all other products."),
        AnnotationSlot.INTRO, "<h2>Source-defined oligos</h2><seq>" + "".join(ordered) + "</seq>",
        WorkflowSection("Library generation"), AnnotationSlot.AFTER_WORKFLOW,
        "<h2>Final library — illustrative S502/N714</h2>", panel(Scene.duplex(final_library().segments).rows(), cls="long"),
        sp.section(final_library(), SEQ_PRIMERS, heading="Source standard run: 75 / 8 / 8 / 0", read_lengths=READ_LENGTHS,
                   required_roles=tuple(p.role for p in SEQ_PRIMERS), show_read_content=True),
        info("Read1 is transcript sequence; both eight-cycle index reads serve sample demultiplexing. The defining Methods report NextSeq550 standard 75-8-8-0. No UMI is assigned to this branch. Full physical primer sites on the inferred ME regions are conditional on that endpoint assignment."), "</div>"])


def selftest():
    c = Check()
    bundle = Path(os.environ.get("CHEM_DATA", HERE.parents[1] / "_data")) / "sources" / HERE.parent.name
    original = Source(bundle / "s41587-022-01312-3_PMC9546769.xml.txt", "https://doi.org/10.1038/s41587-022-01312-3")
    if have(c, original, label="FLASH defining Methods transcription"):
        text = original.path.read_text()
        dt = text.split("FS-dT30VN (5′Bio-", 1)[1].split("-3′", 1)[0]
        dt = re.sub(r"T(\d+)", lambda m: "T" * int(m[1]), dt)
        c("FS dT30VN literal source expansion", RT_SEQUENCE, dt)
        source_tso = text.split("FS-TSO (5′Bio-", 1)[1].split("-3′", 1)[0]
        c("FS TSO modifications and RNA G tail", TSO_AS_WRITTEN, source_tso)
        c("standard source acquisition", "75-8-8-0" in text.split("Sample pooling and sequencing (all protocols)", 1)[1])
    author = Source(bundle / "protocolsio_flash-seq-protocol-kxygxzkrwv8j_v3.pdf.txt",
                    "https://doi.org/10.17504/protocols.io.kxygxzkrwv8j/v3")
    if have(c, author, label="article-linked v3 modified index oligos"):
        compact = re.sub(r"\s+", "", author.path.read_text())
        for name, sequence in (("S502", INDEX_FORWARD), ("N714", INDEX_REVERSE)):
            c(name + " full modified primer in original v3 table",
              "/5Biosg/" + sequence[:-1] + "*" + sequence[-1] in compact)
    protocol_graph()
    raise SystemExit(c.report())
