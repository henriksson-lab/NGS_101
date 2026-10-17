"""Source-owned protocol definition and reaction graph."""
from __future__ import annotations

from protocol_specs import Definition, get
from chemdraw import Segment, feature, oligo
import illumina as il
import nextera as nx
import seqprimers as sp
from chemdraw import Construct, revcomp
from protocol_graph import Highlight, InputRole, SourceClaim, SourceDiscrepancy
from protocol_graph import MaterialKind
from protocol_graph import ProtocolGraph
from protocol_graph import StepKind
from rna_library_operations import finish_pcr
from rna_library_operations import molecule

DEFINITION = Definition("SNuBar-ATAC — oligo sample multiplexing for scATAC", "10.1016/j.molcel.2021.09.026", "Wang et al., Molecular Cell (2021)", "An Illumina TDE1 transposome tags each sample, then a sample-barcode oligo anneals to its universal tail before nuclei are pooled for 10x scATAC.", "nextera", ("cell barcode", "accessible DNA"), "After 10x GEM cleanup, separate ATAC and SBO barcode templates are PCR-amplified and the two libraries are mixed for sequencing. ATAC reads report genomic accessibility and a cell barcode; the SBO branch reports sample barcode and UMI.", "DNA", "The displayed final molecule is the ATAC branch; the reaction graph also builds the separate SBO barcode library from the printed SBO, TruSeq RPIX and P5 oligos and the spiked custom Read 1 and I5 index primers. The 16-bp cell barcode junctions (P5 directly to N16 directly to the 14-nt Nextera S5 tail), the 8-bp ATAC i7 index and the 50/8/16/50 read lengths follow the authors' Figure S1A and Methods S1 step 10; the 10x guides CG000168 Rev A / CG000209 Rev D themselves were not inspected. Methods S1 names no Index 1 or Read 2 primer for the SBO library, whose RPIX i7 index is 6 nt while Index 1 runs 8 cycles, so no SBO Index 1 or Read 2 read is drawn. Individual 18-base sample barcodes are listed in Table S1; the model keeps the generic N18 window.")


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
RPIX_P7_ARM = "CAAGCAGAAGACGGCATACGAGAT"
# RPIX bases between the index and SBO_HANDLE, derived from the printed primer.
RPIX_INNER_ARM = RPIX_PRIMER[len(RPIX_P7_ARM) + 6:-len(SBO_HANDLE)]
if RPIX_PRIMER != RPIX_P7_ARM + "N" * 6 + RPIX_INNER_ARM + SBO_HANDLE:
    raise AssertionError("printed RPIX primer no longer splits into P7 arm, N6, arm and handle")
P5_ADAPTOR = "AATGATACGGCGACCACCGAGATCTACAC"
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


from protocol_specs import workflow, two_branch


def _source_workflow():
    return workflow(molecule("sample-specific nuclei", ("chromatin",)), (
        ("Tagment each sample with Illumina TDE1",
         molecule("tagmented chromatin", ("genomic fragment", "transposome universal tail")),
         "The sample barcode is not carried by Tn5 itself."),
        ("Anneal a SNuBar oligo to the transposome tail",
         molecule("sample-tagged fragment", ("genomic fragment", "sample barcode", "UMI", "PCR handle"), joined=True),
         "The SBO transfers sample identity before samples are pooled."),
        ("Pool samples, run 10x GEM barcoding, then separate ATAC and SBO barcode templates",
         two_branch("same cell", "ATAC library", "SBO barcode library"),
         "The two products are separately amplified, then mixed for sequencing."),
    ))

SPEC = get('snubar-atac', definition=DEFINITION,
           workflow_fn=_source_workflow, oligos_fn=snubar_oligos,
           final_override=atac_library())


def _local_graph(d, final, primers):
    source = f"{d.source_label}; doi:{d.doi}"
    graph = ProtocolGraph(d.title, authority=source)
    nuclei = graph.input("nuclei", "sample-separated nuclei containing accessible chromatin",
                         MaterialKind.BIOLOGICAL,
                         state=molecule("sample-specific nuclei", ("accessible chromatin",)),
                         primary=True)
    transposome = graph.input("transposome", "Illumina TDE1 loaded transposome",
                              MaterialKind.ADAPTER)
    tn5 = graph.input("tn5", "Tn5 transposase", MaterialKind.CATALYST)
    tagged = graph.product("tagged", "tagmented chromatin carrying the universal tail",
                           MaterialKind.MOLECULAR,
                           state=molecule("tagmented chromatin", (
                               "genomic fragment", "transposome universal tail")))
    waste = graph.product("tn5_waste", "non-library tagmentation products", MaterialKind.WASTE)
    graph.react("tagment", "Tagment each sample with Illumina TDE1",
                StepKind.TAGMENTATION,
                ((nuclei, InputRole.TEMPLATE), (transposome, InputRole.STRUCTURAL),
                 (tn5, InputRole.CATALYST)), (tagged, waste))
    graph.discard(waste, "does not retain an accessible-DNA insert with the universal tail")
    snubar = graph.input("snubar", "sample-barcode/UMI SNuBar oligo",
                         MaterialKind.OLIGO)
    sample_tagged = graph.product("sample_tagged", "SNuBar-tagged ATAC fragments",
                                  MaterialKind.MOLECULAR,
                                  state=molecule("sample-tagged fragment", (
                                      "genomic fragment", "sample barcode", "UMI", "PCR handle")))
    graph.react("anneal_snubar", "Anneal the SNuBar oligo to the transposome tail",
                StepKind.HYBRIDIZATION,
                ((tagged, InputRole.TEMPLATE), (snubar, InputRole.STRUCTURAL)),
                (sample_tagged,))
    droplets = graph.input("droplets", "10x single-cell ATAC droplet condition",
                           MaterialKind.CONDITION)
    barcoded = graph.product("barcoded", "cell- and sample-barcoded ATAC molecules",
                             MaterialKind.MOLECULAR,
                             state=molecule("SNuBar-ATAC products", (
                                 "10x cell barcode", "sample barcode", "UMI", "ATAC fragment")))
    collisions = graph.product("collisions", "empty or ambiguous droplets", MaterialKind.WASTE)
    graph.react("partition", "Pool samples and partition nuclei into 10x droplets",
                StepKind.PARTITION,
                ((sample_tagged, InputRole.TEMPLATE), (droplets, InputRole.CONDITION)),
                (barcoded, collisions))
    graph.discard(collisions, "does not retain an interpretable cell barcode")
    atac_template = graph.product("atac_template", "cell-barcoded ATAC template",
                                  MaterialKind.MOLECULAR,
                                  state=molecule("ATAC template", ("10x cell barcode", "ATAC fragment")))
    barcode_template = graph.product("barcode_template", "cell-barcoded SBO sample-tag template",
                                     MaterialKind.MOLECULAR,
                                     state=molecule("SBO barcode template", (
                                         "10x cell barcode", "sample barcode", "UMI")))
    graph.react("separate_templates", "Separate templates by 1.2X AMPure after Dynabead cleanup: ATAC on beads, SBO barcodes in supernatant (then 1.8X)",
                StepKind.SPLIT, ((barcoded, InputRole.TEMPLATE),),
                (atac_template, barcode_template), changes_molecule=False)
    reviewed = ["nuclei", "transposome", "tn5", "tagment", "snubar", "anneal_snubar",
                "droplets", "partition", "separate_templates"]
    finish_pcr(graph, "atac", atac_template, final, primers, d.sequencing, reviewed,
               label="SNuBar-ATAC accessible-DNA")
    barcode_final, barcode_primers = sbo_library()
    finish_pcr(graph, "sbo", barcode_template, barcode_final, barcode_primers,
               "separate SBO sample-barcode read", reviewed, label="SNuBar-ATAC SBO barcode")
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
