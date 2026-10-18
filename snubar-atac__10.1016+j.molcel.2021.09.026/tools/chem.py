"""Source-owned protocol definition and reaction graph."""
from __future__ import annotations

from types import SimpleNamespace
from protocol_model import molecule, role
from pcr import PCRGeometry, PCRPrimer
from protocol_graph import Disposition, Terminal, SequencingPlan
from chemdraw import Segment, Scene, feature, oligo
import illumina as il
import nextera as nx
import seqprimers as sp
from chemdraw import Construct, revcomp
from protocol_graph import Highlight, InputRole, SourceClaim, SourceDiscrepancy
from protocol_graph import MaterialKind
from protocol_graph import ProtocolGraph
from protocol_graph import StepKind

DEFINITION = SimpleNamespace(
    caveat="The separate ATAC and SBO libraries retain their own PCR templates and read plans; the reaction graph also builds the separate SBO barcode library from the printed SBO, TruSeq RPIX and P5 oligos and the spiked custom Read 1 and I5 index primers. The 16-bp cell barcode junctions (P5 directly to N16 directly to the 14-nt Nextera S5 tail), the 8-bp ATAC i7 index and the 50/8/16/50 read lengths follow the authors' Figure S1A and Methods S1 step 10; the 10x guides CG000168 Rev A / CG000209 Rev D themselves were not inspected. Methods S1 names no Index 1 or Read 2 primer for the SBO library, whose RPIX i7 index is 6 nt while Index 1 runs 8 cycles, so those acquired mixed-run reads retain an unresolved SBO priming assignment. Individual 18-base sample barcodes are listed in Table S1; the model keeps the generic N18 window.",
    title='SNuBar-ATAC — oligo sample multiplexing for scATAC',
    doi='10.1016/j.molcel.2021.09.026',
    source_label='Wang et al., Molecular Cell (2021)',
    summary='An Illumina TDE1 transposome tags each sample, then a sample-barcode oligo anneals to its universal tail before nuclei are pooled for 10x scATAC.',
    sequencing='After 10x GEM cleanup, separate ATAC and SBO barcode templates are PCR-amplified and the two libraries are mixed for sequencing. ATAC reads report genomic accessibility and a cell barcode; the SBO branch reports sample barcode and UMI.',
)


# Sequences printed in the original STAR Methods (PMC8611914 JATS lines 187 and 193;
# key-resources table lines 769-777; published STAR Methods, Document S2 lines
# 1672-1697) and in Methods S1 step 7(2) / Table S1 (Document S1).
SBO_HANDLE = "CCTTGGCACCCGAGAATTCCA"
# Table S1 caption: 19-bp universal sequence then the 14-bp complementary sequence of
# the transposome universal tail (CST).  The spiked I5 index primer prints both.
SBO_UNIVERSAL = "CTACTCTGCGTTGATACCA"
SBO_CST = "GACGCTGCCGACGA"
SBO_TAIL = SBO_UNIVERSAL + SBO_CST
RPIX_PRIMER = "CAAGCAGAAGACGGCATACGAGATNNNNNNGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"
RPIX_P7_ARM = il.P7  # Printed RPIX primer P7 arm, STAR Methods.
# RPIX bases between the index and SBO_HANDLE, derived from the printed primer.
RPIX_INNER_ARM = RPIX_PRIMER[len(RPIX_P7_ARM) + 6:-len(SBO_HANDLE)]
if RPIX_PRIMER != RPIX_P7_ARM + "N" * 6 + RPIX_INNER_ARM + SBO_HANDLE:
    raise AssertionError("printed RPIX primer no longer splits into P7 arm, N6, arm and handle")
P5_ADAPTOR = il.P5  # Printed P5 Adaptor, Methods S1 step 7(2).
SBO_READ1_PRIMER = "TCGTCGGCAGCGTCTGGTATCAACGCAGAGTAG"
SBO_I5_INDEX_PRIMER = "CTACTCTGCGTTGATACCAGACGCTGCCGACGA"
_SRC = ("Wang et al. 2021 STAR Methods (PMC8611914 JATS line 193; key resources line 775-777; "
        "published Document S2 lines 1697-1698); Methods S1 step 10")
if revcomp(SBO_TAIL) != SBO_READ1_PRIMER or SBO_I5_INDEX_PRIMER != SBO_TAIL:
    raise AssertionError("printed SBO tail and spiked sequencing primers no longer interlock")
if revcomp(SBO_CST) != nx.S5:
    raise AssertionError("Table S1 CST no longer complements the Nextera S5 tail")


def _cell_barcode():
    return Segment("10x cell barcode", "N" * 16, "cbc", placeholder=True, length_bp=16,
                   feature=feature("snubar_cell", "cell_barcode", "unknown"),
                   note="16 bp, read as Index 2 (STAR Methods line 205). Junctions as in the "
                        "authors' Figure S1A cell-barcode primer "
                        "AATGATACGGCGACCACCGAGATCTACAC-N16-TCGTCGGCAGCGTC; 10x guide "
                        "CG000168 Rev A / CG000209 Rev D not inspected.")


def atac_library():
    """ATAC branch: 10x cell barcode in the i5 position, single Chromium i7 sample index."""
    lib = Construct([
        Segment("P5", il.P5, "p5"),
        _cell_barcode(),
        Segment("S5", nx.S5, "s5"),
        Segment("left mosaic end", nx.ME, "me"),
        Segment("accessible DNA", "X" * 28, placeholder=True),
        Segment("right mosaic end reverse complement", nx.ME_RC, "me"),
        Segment("S7 reverse complement", nx.S7_RC, "s7"),
        Segment("Chromium i7 sample index (Set A)", "I" * 8, "cbc", placeholder=True,
                length_bp=8,
                feature=feature("snubar_atac_i7", "sample_index", "unknown"),
                note="8 bp: Figure S1A ATAC library primer and final library (-I8-), "
                     "Methods S1 step 10 (Index 1: 8bp)."),
        Segment("P7 reverse complement", il.P7_RC, "p7"),
    ], name="SNuBar-ATAC accessible-DNA library")
    primers = (sp.NEXTERA["R1"], sp.NEXTERA["I1"], sp.NEXTERA["I2"], sp.NEXTERA["R2"])
    problems = sp.verify(lib, primers, required_roles=tuple(p.role for p in primers))
    if problems:
        raise ValueError("; ".join(problems))
    return lib, primers


def sbo_library():
    """SBO barcode branch, assembled from the printed SBO, RPIX and P5 oligos."""
    lib = Construct([
        Segment("P5 adaptor", P5_ADAPTOR, "p5"),
        _cell_barcode(),
        Segment("CST complement (custom Read 1 site, part 1)", revcomp(SBO_CST), "r1"),
        Segment("universal sequence complement (custom Read 1 site, part 2)",
                revcomp(SBO_UNIVERSAL), "r1"),
        Segment("sample barcode", "N" * 18, "cbc", placeholder=True, length_bp=18,
                feature=feature("snubar_sample", "sample_index", "unknown")),
        Segment("UMI", "N" * 8, "umi", placeholder=True, length_bp=8,
                feature=feature("snubar_umi", "umi", "random")),
        Segment("PCR handle complement", revcomp(SBO_HANDLE), "r2"),
        Segment("RPIX inner arm complement", revcomp(RPIX_INNER_ARM), "r2"),
        Segment("RPIX i7 index", "I" * 6, "cbc", placeholder=True, length_bp=6,
                feature=feature("snubar_pcr_index", "sample_index", "unknown")),
        Segment("P7 reverse complement", revcomp(RPIX_P7_ARM), "p7"),
    ], name="SNuBar-ATAC SBO barcode library")
    primers = (
        sp.custom("Read 1", "SNuBar custom Read 1", SBO_READ1_PRIMER, _SRC,
                  "Spiked into the NextSeq run; reads the 18-bp sample barcode then 8-bp UMI."),
        sp.custom("Index 2 (i5)", "SNuBar I5 index primer", SBO_I5_INDEX_PRIMER, _SRC,
                  "Spiked into the NextSeq run; reads the 16-bp cell barcode."),
    )
    problems = sp.verify(lib, primers, required_roles=tuple(p.role for p in primers))
    if problems:
        raise ValueError("; ".join(problems))
    return lib, primers


def snubar_oligos():
    """Published ATAC-only SBO family and barcode-library primers (STAR Methods)."""
    return (
        oligo("SNuBar-ATAC SBO family", [
            Segment("PCR handle", SBO_HANDLE),
            Segment("UMI", "N" * 8, "umi", placeholder=True, length_bp=8,
                    feature=feature("snubar_umi", "umi", "random")),
            Segment("sample barcode", "N" * 18, "cbc", placeholder=True, length_bp=18,
                    feature=feature("snubar_sample", "sample_index", "unknown")),
            Segment("universal sequence", SBO_UNIVERSAL),
            Segment("CST", SBO_CST),
        ]),
        oligo("TruSeq RPIX barcode-library primer", [
            Segment("P7 and index-primer arm", RPIX_P7_ARM),
            Segment("sample index", "N" * 6, "cbc", placeholder=True, length_bp=6,
                    feature=feature("snubar_pcr_index", "sample_index", "unknown")),
            Segment("barcode-library arm", RPIX_INNER_ARM + SBO_HANDLE),
        ]),
        oligo("P5 barcode-library primer", [
            Segment("P5 primer", P5_ADAPTOR),
        ]),
    )


# Each branch owns its source-stated endpoint, primers and acquired cycle counts.
READ_LENGTHS = {"Read 1": 50, "Index 1 (i7)": 8, "Index 2 (i5)": 16, "Read 2": 50}
ATAC_FORWARD = il.P5[:22]  # Figure S1A, ATAC library primer (AATG...GAGA).
ATAC_REVERSE = il.P7 + "N" * 8 + nx.S7  # Figure S1A, Set A i7 primer family.
SOURCE_OLIGO_EXCEPTIONS = {
    "atac_reverse": "Figure S1A prints P7-I8-S7 rather than eight literal N bases. The local source test checks the original I8 family and both flanking arms; individual Set A index choices remain unassigned.",
}
ATAC_FINAL, ATAC_SEQ_PRIMERS = atac_library()
SPEC = SimpleNamespace(
    TITLE=DEFINITION.title, NOTES="01_snubar_atac.html",
    SOURCE='<a href="https://doi.org/10.1016/j.molcel.2021.09.026">Wang et al. 2021</a>; Figure S1A and Methods S1.',
    SUMMARY=DEFINITION.summary, CAVEAT=DEFINITION.caveat,
    FINAL_LIBRARY=ATAC_FINAL, SEQ_PRIMERS=ATAC_SEQ_PRIMERS,
    FINAL_CAPTION=DEFINITION.sequencing, SEQUENCING_INTRO=DEFINITION.sequencing,
    SEQUENCING_ENDING=DEFINITION.sequencing, SEQUENCING_UNAVAILABLE="",
    READ_LENGTHS=READ_LENGTHS, oligos=snubar_oligos,
)

def _atac_template():
    # Figure S1A step 3: P5/CB16 is already copied during single-primer GEM
    # amplification. The indexed P7 tail is added only in step 4.
    return Construct(list(ATAC_FINAL)[:-2], name="cell-barcoded ATAC pre-index-PCR template")

def _sbo_template():
    # Figure S1A step 3, before RPIX adds its P7/index/inner arm.
    return Construct(list(sbo_library()[0])[:7], name="cell-barcoded SBO pre-PCR template")


def _sbo_oligo():
    return Construct([
        Segment("PCR handle", SBO_HANDLE),
        Segment("UMI", "N" * 8, "umi", placeholder=True, length_bp=8,
                feature=feature("snubar_umi", "umi", "random")),
        Segment("sample barcode", "N" * 18, "cbc", placeholder=True, length_bp=18,
                feature=feature("snubar_sample", "sample_index", "unknown")),
        Segment("universal sequence", SBO_UNIVERSAL), Segment("CST", SBO_CST),
    ], name="SNuBar-ATAC SBO family")

def _tagged_chromatin(with_sbo=False):
    # Figure S1A step 1 shows SBO hybridization to the S5 tail, rather
    # than covalent insertion of the sample tag into the genomic strand.
    state = Scene()
    state.strand("genomic", list(_atac_template())[2:-1],
                 label="tagmented genomic strand (opposite strand not reconstructed)")
    if with_sbo:
        state.anneal("SBO", list(_sbo_oligo()), to="genomic", pair=("CST", "S5"), label="annealed SBO",
                     unpaired=("PCR handle", "UMI", "sample barcode", "universal sequence"))
        state.note("SBO", "Only the 14-nt CST pairs with S5; the upstream SBO tail is unpaired.")
    return state

def _gap_filled_products():
    # Figure S1A step 2: a completed ATAC template and a distinct SBO.
    state = Scene()
    state.strand("ATAC", list(_atac_template())[2:], label="gap-filled ATAC (one strand shown)")
    state.strand("SBO", list(_sbo_oligo()), label="separate SBO template")
    return state

def _gem_products():
    state = Scene()
    state.strand("ATAC", list(_atac_template()), label="ATAC branch (one strand shown)")
    state.strand("SBO", list(_sbo_template()), label="SBO branch (one strand shown)")
    return state


def _local_graph(d, final, primers):
    source = f"{d.source_label}; doi:{d.doi}; published Figure S1A steps 1–4 and Methods S1, ATAC steps 3–10"
    graph = ProtocolGraph(d.title, authority=source)
    nuclei = graph.input("nuclei", "sample-separated nuclei containing accessible chromatin",
                         MaterialKind.BIOLOGICAL,
                         state=molecule("sample-specific nuclei", (role("accessible chromatin"),)),
                         primary=True)
    transposome = graph.input("transposome", "Illumina TDE1 loaded transposome",
                              MaterialKind.ADAPTER)
    tn5 = graph.input("tn5", "Tn5 transposase", MaterialKind.CATALYST)
    tagged = graph.product("tagged", "tagmented chromatin carrying the universal tail",
                           MaterialKind.MOLECULAR,
                           state=_tagged_chromatin())
    waste = graph.product("tn5_waste", "other tagmentation products (not individually resolved)", MaterialKind.MOLECULAR)
    graph.react("tagment", "Tagment each sample with Illumina TDE1",
                StepKind.TAGMENTATION,
                ((nuclei, InputRole.TEMPLATE), (transposome, InputRole.STRUCTURAL),
                 (tn5, InputRole.CATALYST)), (tagged, waste))
    graph.mark_terminal(Terminal(waste.id, Disposition.UNRESOLVED,
        "The inspected original does not establish class-specific physical removal of other tagmentation outcomes; amplification and recovery are unresolved."))
    snubar = graph.input("snubar", "sample-barcode/UMI SNuBar oligo",
                         MaterialKind.OLIGO, state=_sbo_oligo(), source=_SRC)
    sample_tagged = graph.product("sample_tagged", "ATAC fragments with hybridized, noncovalent SBO",
                                  MaterialKind.MOLECULAR, state=_tagged_chromatin(with_sbo=True))
    graph.react("anneal_snubar", "Anneal the SNuBar oligo to the transposome tail",
                StepKind.HYBRIDIZATION,
                ((tagged, InputRole.TEMPLATE), (snubar, InputRole.STRUCTURAL)),
                (sample_tagged,))
    droplets = graph.input("droplets", "10x single-cell ATAC droplet condition",
                           MaterialKind.CONDITION)
    barcoded = graph.product("barcoded", "partitioned sample-tagged nuclei",
                             MaterialKind.MOLECULAR,
                             state=sample_tagged.state)
    collisions = graph.product("collisions", "empty or ambiguous droplets", MaterialKind.MOLECULAR)
    graph.react("partition", "Pool samples and partition nuclei into 10x droplets",
                StepKind.PARTITION,
                ((sample_tagged, InputRole.TEMPLATE), (droplets, InputRole.CONDITION)),
                (barcoded, collisions))
    graph.mark_terminal(Terminal(collisions.id, Disposition.UNRESOLVED,
        "Barcode ambiguity is a computational selection criterion; physical material is not shown to be removed here."))
    gem_enzyme = graph.input("gem_enzyme", "10x Barcoding Enzyme PN-2000125/2000139", MaterialKind.CATALYST, source="Methods S1 step 5")
    gap_filled = graph.product("gap_filled", "gap-filled ATAC molecules and separate SBO templates", MaterialKind.MOLECULAR,
        state=_gap_filled_products())
    graph.react("gem_gap_fill", "Gap-fill in GEMs (72 °C, 5 min)", StepKind.EXTENSION,
        ((barcoded, InputRole.TEMPLATE), (gem_enzyme, InputRole.CATALYST)), (gap_filled,))
    cell_primer = graph.input("cell_barcode_primer", "10x cell-barcode primer family, Figure S1A", MaterialKind.OLIGO,
        state=Construct([Segment("P5", P5_ADAPTOR), _cell_barcode(), Segment("S5", nx.S5)], name="GEM cell-barcode primer"),
        source="Figure S1A step 3; Methods S1 step 5")
    amplified = graph.product("gem_amplified", "cell-barcoded ATAC and SBO products of GEM linear amplification", MaterialKind.MOLECULAR,
        state=_gem_products())
    graph.react("gem_linear_amplification", "Single-primer linear amplification in GEMs (12 cycles, 98/59/72 °C)", StepKind.EXTENSION,
        ((gap_filled, InputRole.TEMPLATE), (cell_primer, InputRole.STRUCTURAL), (gem_enzyme, InputRole.CATALYST)), (amplified,))
    atac_template = graph.product("atac_template", "cell-barcoded ATAC template",
                                  MaterialKind.MOLECULAR,
                                  state=_atac_template())
    barcode_template = graph.product("barcode_template", "cell-barcoded SBO sample-tag template",
                                     MaterialKind.MOLECULAR,
                                     state=_sbo_template())
    graph.react("separate_templates", "Separate templates by 1.2X AMPure after Dynabead cleanup: ATAC on beads, SBO barcodes in supernatant (then 1.8X)",
                StepKind.SPLIT, ((amplified, InputRole.TEMPLATE),),
                (atac_template, barcode_template), changes_molecule=False)
    reviewed = ["nuclei", "transposome", "tn5", "tagment", "snubar", "anneal_snubar",
                "droplets", "partition", "separate_templates"]
    # Figure S1A gives these primers and the actual pre-PCR templates.
    # Exact geometry checks only inward primer landing, not reaction efficiency.
    atac_forward = graph.input("atac_forward", "SI-PCR Primer B (Figure S1A ATAC forward primer)", MaterialKind.OLIGO,
        state=Construct([Segment("SI-PCR Primer B", ATAC_FORWARD)]), source="Figure S1A step 4; Methods S1 step 7(1)")
    atac_reverse = graph.input("atac_reverse", "Chromium i7 Sample Index N, Set A primer family", MaterialKind.OLIGO,
        state=Construct([Segment("P7", il.P7), Segment("i7 sample index", "N" * 8, "cbc", placeholder=True, length_bp=8,
            feature=feature("snubar_atac_i7", "sample_index", "unknown")), Segment("S7", nx.S7)]), source="Figure S1A step 4")
    atac_mix = graph.input("atac_polymerase", "10x Amp Mix PN-2000047/2000103", MaterialKind.CATALYST, source="Methods S1 step 7(1)")
    atac_pcr = graph.product("atac_pcr_product", "ATAC sample-index PCR product", MaterialKind.MOLECULAR, state=final)
    graph.react("atac_pcr", "ATAC sample-index PCR (9–12 cycles)", StepKind.PCR,
        ((atac_template, InputRole.TEMPLATE), (atac_forward, InputRole.STRUCTURAL), (atac_reverse, InputRole.STRUCTURAL),
         (atac_mix, InputRole.CATALYST)), (atac_pcr,), geometry=PCRGeometry.exact(atac_template.state,
            PCRPrimer("SI-PCR Primer B", ATAC_FORWARD, len(ATAC_FORWARD)), PCRPrimer("Chromium i7 Set A family", ATAC_REVERSE, len(nx.S7)),
            statement="Figure S1A steps 3–4 gives the pre-PCR template and primer families; Methods S1 step 7(1) specifies Amp Mix and 9–12 cycles."))
    atac_library_product = graph.product("atac_library", "SNuBar-ATAC accessible-DNA sequencing library", MaterialKind.LIBRARY, state=final)
    atac_selection = graph.input("atac_selection", "SPRIselect 0.4X/1.2X double-sided selection; optional extra 1.2X", MaterialKind.SELECTION)
    atac_waste = graph.product("atac_cleanup_waste", "ATAC size-selection excluded fractions and washes", MaterialKind.WASTE)
    graph.react("atac_cleanup", "Size-select the ATAC library", StepKind.SELECTION,
        ((atac_pcr, InputRole.TEMPLATE), (atac_selection, InputRole.SELECTION)), (atac_library_product, atac_waste), changes_molecule=False)
    graph.discard(atac_waste, "excluded fractions and ethanol washes are removed during SPRIselect cleanup (Methods S1 step 8)")
    graph.sequence(atac_library_product, SequencingPlan(library=final, primers=tuple(primers), read_lengths=READ_LENGTHS))
    barcode_final, barcode_primers = sbo_library()
    sbo_forward = graph.input("sbo_forward", "P5 Adaptor", MaterialKind.OLIGO,
        state=Construct([Segment("P5 Adaptor", P5_ADAPTOR)]), source=_SRC)
    sbo_reverse = graph.input("sbo_reverse", "TruSeq RPIX primer family", MaterialKind.OLIGO,
        state=Construct([Segment("P7 arm", RPIX_P7_ARM), Segment("RPIX i7", "N" * 6, "cbc", placeholder=True, length_bp=6,
            feature=feature("snubar_pcr_index", "sample_index", "unknown")), Segment("inner arm and handle", RPIX_INNER_ARM + SBO_HANDLE)]), source=_SRC)
    sbo_mix = graph.input("sbo_polymerase", "2X KAPA HiFi HotStart ReadyMix", MaterialKind.CATALYST, source="Methods S1 step 7(2)")
    sbo_pcr = graph.product("sbo_pcr_product", "SBO barcode PCR product", MaterialKind.MOLECULAR, state=barcode_final)
    graph.react("sbo_pcr", "SBO barcode-library PCR (8–10 cycles)", StepKind.PCR,
        ((barcode_template, InputRole.TEMPLATE), (sbo_forward, InputRole.STRUCTURAL), (sbo_reverse, InputRole.STRUCTURAL),
         (sbo_mix, InputRole.CATALYST)), (sbo_pcr,), geometry=PCRGeometry.exact(barcode_template.state,
            PCRPrimer("P5 Adaptor", P5_ADAPTOR, len(P5_ADAPTOR)), PCRPrimer("TruSeq RPIX", RPIX_PRIMER, len(SBO_HANDLE)),
            statement="Figure S1A step 3 supplies the pre-PCR SBO template; the printed RPIX and P5 primers (Methods S1 step 7(2)) add the index/P7 arm."))
    barcode_library = graph.product("sbo_library", "SNuBar-ATAC SBO barcode sequencing library", MaterialKind.LIBRARY, state=barcode_final)
    sbo_selection = graph.input("sbo_selection", "SPRIselect 1.5X", MaterialKind.SELECTION)
    sbo_waste = graph.product("sbo_cleanup_waste", "SBO cleanup supernatant and ethanol washes", MaterialKind.WASTE)
    graph.react("sbo_cleanup", "Purify the SBO barcode library", StepKind.SELECTION,
        ((sbo_pcr, InputRole.TEMPLATE), (sbo_selection, InputRole.SELECTION)), (barcode_library, sbo_waste), changes_molecule=False)
    graph.discard(sbo_waste, "supernatant and ethanol washes are removed (Methods S1 step 9)")
    graph.sequence(barcode_library, SequencingPlan(library=barcode_final, primers=tuple(barcode_primers),
        read_lengths={"Read 1": 50, "Index 2 (i5)": 16},
        unplaced_reads=(
            sp.read_spec("Index 1 (i7)", 8, "RPIX i7 has 6 bases; SBO primer/site unassigned", "undetermined"),
            sp.read_spec("Read 2", 50, "SBO primer/site unassigned in this mixed run", "undetermined"),
        ),
        limitation="Methods S1 step 10 acquires R1/I1/I2/R2 at 50/8/16/50 cycles for the mixed libraries. SBO R1 and I2 sites are source-defined; SBO Index 1 and Read 2 primer assignments remain unresolved. Its RPIX index is only 6 nt."))
    reviewed.extend(("gem_enzyme", "gem_gap_fill", "cell_barcode_primer", "gem_linear_amplification",
        "atac_template", "barcode_template", "atac_forward", "atac_reverse", "atac_polymerase", "atac_pcr",
        "atac_selection", "atac_cleanup", "atac_library", "sbo_forward", "sbo_reverse", "sbo_polymerase",
        "sbo_pcr", "sbo_selection", "sbo_cleanup", "sbo_library"))
    graph.declare_reviewed_scope(reviewed, source)
    graph.annotate(Highlight(
        "atac_pcr", "ATAC sample-index PCR",
        "Methods S1 step 7(1): 10x Amp Mix, SI-PCR Primer B and one Chromium i7 Sample "
        "Index N, Set A, 9&ndash;12 cycles, then 0.4X/1.2X double-sided SPRIselect "
        "selection plus an optional extra 1.2X cleanup. Figure S1A draws the i7 primer as "
        "CAAGCAGAAGACGGCATACGAGAT-I8-GTCTCGTGGGCTCGG."))
    graph.annotate(Highlight(
        "sbo_pcr", "SBO barcode-library PCR and run",
        "Methods S1 step 7(2): 2X KAPA HiFi HotStart ReadyMix, 10 &micro;l barcode "
        "template, 1.5 &micro;l each of 10 &micro;M TruSeq RPIX and P5 Adaptor, 8&ndash;10 "
        "cycles, then 1.5X SPRIselect. ATAC and barcode libraries are mixed 8:2 (STAR "
        "Methods) or 4:1&ndash;20:1 (Methods S1 step 10) for NextSeq with the custom Read 1 "
        "and I5 index primers spiked in; run lengths Read 1 50, Index 1 8, Index 2 16, "
        "Read 2 50. No Index 1 or Read 2 primer is named for this library and its RPIX "
        "index is 6 nt, so those reads are not drawn on it."))
    graph.annotate(SourceDiscrepancy(
        "sbo_pcr", "RPIX junction drawn reversed in Figure S1A",
        (SourceClaim("Wang et al. 2021 STAR Methods and Methods S1 step 7(2); Figure S1A "
                     "barcode-library primer",
                     "RPIX primer CAAGCAGAAGACGGCATACGAGAT-N6-GTGACTGGAGTTCCTTGGCACCC"
                     "GAGAATTCCA, so the library top strand reads "
                     "...CCAAGGAACTCCAGTCAC-I6-ATCTCG..."),
         SourceClaim("Figure S1A, final barcode library",
                     "Top strand ...TGGAATTCTCGGGTGCCAAGGTTGAGGTCAGTG-I6-ATCTCG...: the 12 "
                     "bases after the handle are GTGACTGGAGTT written in reverse.")),
        "The model takes the junction from the printed RPIX primer.",
        "The primer is printed identically in the STAR Methods, Methods S1 and Figure S1A's "
        "own primer list; the final-library drawing is the only source with the reversed "
        "12 bases.",
        "None for the sequence; the figure line is treated as a drawing inconsistency."))
    return graph.require_publishable()


def protocol_graph():
    return _local_graph(DEFINITION, SPEC.FINAL_LIBRARY, SPEC.SEQ_PRIMERS)

SPEC.protocol_graph = protocol_graph
# Methods S1 step 10 (Document S1): "Read 1: 50bp, Read 2: 50bp, Index 1: 8bp, Index 2: 16bp".
SPEC.READ_LENGTHS = {"Read 1": 50, "Index 1 (i7)": 8, "Index 2 (i5)": 16, "Read 2": 50}
