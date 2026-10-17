"""Robertson et al. 2007 ChIP-seq and original single-read Solexa library."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))

import illumina as il
import seqprimers as sp
from chemdraw import Construct, Row, Scene, Segment, revcomp
from chromatin_epigenetics import antibody_enzyme_rows
from adapter_ligation import AdapterLigationGeometry
from pcr import PCRGeometry
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            SourceClaim, SourceDiscrepancy, StepKind)
from protocol_model import molecule, role


def seg(name, top, tag=None, **kw):
    return Segment(name, top, tag, **kw)


TITLE = "ChIP-seq — conventional immunoprecipitation library"
NOTES = "01_chip-seq.html"
SOURCE = 'Defining protocol: <a href="https://doi.org/10.1038/nmeth1068">Robertson et al., <i>Nature Methods</i> 4:651–657 (2007)</a>.'
SUPPLEMENT = "Robertson et al. 2007 Supplementary Methods (Springer MOESM305), ChIP DNA preparation"
# Robertson et al. 2007 publisher PDF, Methods "Illumina library construction and
# sequencing" (nmeth1068.pdf.txt lines 751-786), read 2026-10-17.
ARTICLE = ("Robertson et al. 2007, Methods \u201cIllumina library construction and sequencing\u201d")
GUIDE = ("Illumina, Preparing Samples for ChIP Sequencing of DNA, part 11257047 Rev. A (2007); "
         "Robertson et al. follow Illumina's \u201cmanufacturer's instructions\u201d for end repair, "
         "A-tailing, ligation and PCR without naming the document, so the guide's reagent details "
         "are assumed to be those instructions")
ADAPTER_INFERENCE = ("INFERRED — no inspected original prints the Adapter Oligo Mix (part 1000521) "
                     "or PCR primer 1.1/2.1 (parts 1000537/1000538) bases; they are reconstructed "
                     "from the single-read reamplification primers 1.2/2.2 that Schmidt et al. 2009 "
                     "(a later, different-laboratory protocol) print; equality with 1.1/2.1 is assumed")
SUMMARY = "Crosslink and sonicate chromatin, immunoprecipitate target-associated DNA, reverse crosslinks, and prepare the recovered fragments for the original single-read Solexa Genome Analyzer."
CAVEAT = ("This 2007 protocol predates TruSeq and sample indexing. Chromatin preparation follows Robertson et al.’s Supplementary Methods; "
          "the library order, gel, cleanups, PCR cycles and 27 sequencing cycles follow the article’s Methods. For enzyme mixes the article "
          "defers to Illumina’s instructions, so reagent details shown come from Illumina’s 2007 ChIP-Seq Sample Prep guide "
          "(part 11257047), whose identity with those instructions is assumed. Neither source prints the Adapter Oligo Mix or PCR "
          "primer 1.1/2.1 bases, so the adapter, primer and sequencing-primer bases (<inf>marked inferred</inf>) are reconstructed. "
          "No sample index is described; an unindexed read is inferred.")

GENOMIC_ADAPTER_BOTTOM = il.TRUSEQ_READ1
# INFERRED: Schmidt et al. 2009 (PMC4052679, reamplification step) prints these bases as
# "Primer 2.2"; no inspected original prints PCR primer 2.1 (part 1000538) itself.
GENOMIC_PCR2 = "CAAGCAGAAGACGGCATACGAGCTCTTCCGATCT"
# The primer's 3'-terminal T pairs with the insert's dA junction; the rest determines
# the inferred phosphorylated adapter oligo rather than duplicating it independently.
GENOMIC_ADAPTER_TOP = revcomp(GENOMIC_PCR2)[1:]


def adapter_scene() -> Scene:
    sc = Scene()
    sc.strand("bottom", [
        seg("Read 1 fork arm", GENOMIC_ADAPTER_BOTTOM[:-13], "r1", inferred=True),
        seg("genomic stem complement", GENOMIC_ADAPTER_BOTTOM[-13:-1], "r1", inferred=True),
        seg("3′ T", GENOMIC_ADAPTER_BOTTOM[-1], "r1", inferred=True),
    ], label="Genomic Adapter Oligo")
    sc.anneal("top", [
        seg("genomic stem", GENOMIC_ADAPTER_TOP[:12], "p7", inferred=True),
        seg("genomic fork arm", GENOMIC_ADAPTER_TOP[12:], "p7", inferred=True),
    ], to="bottom", pair=("genomic stem", "genomic stem complement"),
              label="Genomic Adapter Oligo", mod5="p",
              unpaired=("genomic fork arm",))
    sc.ligation_end("bottom", "3′ T", kind="cohesive", polarity="3-prime",
                    owner="adapter", text="3′-T ligation overhang")
    return sc


SEQ_PRIMERS = (sp.custom(
    "Read 1", "Genomic DNA sequencing primer", il.TRUSEQ_READ1,
    "INFERRED: 3′ part of Schmidt et al. 2009 reamplification primer 1.2; no inspected 2007 original prints a sequencing primer"),)


def final_library() -> Construct:
    lib = Construct([
        seg("P5 / PCR primer 1.1", il.P5, "p5", inferred=True),
        seg("Read 1 remainder", il.TRUSEQ_READ1[4:], "r1", inferred=True),
        seg("ChIP-enriched DNA", "X" * 36, placeholder=True),
        seg("dA junction", "A"),
        seg("Genomic adapter / PCR primer 2.1 arm", GENOMIC_ADAPTER_TOP, "p7", inferred=True),
    ], name="2007 ChIP-seq library")
    problems = sp.verify(lib, SEQ_PRIMERS, required_roles=("Read 1",))
    if problems:
        raise ValueError("invalid original ChIP-seq library: " + "; ".join(problems))
    return lib


FINAL_LIBRARY = final_library()
# Article Methods line 785 ("27 cycles of sequencing on the Illumina cluster station
# and 1G analyzer"); also Supplementary Data 1 README (MOESM306) lines 18-21.
FINAL_READ_LENGTH = 27
FINAL_CAPTION = ("Unindexed single-read Solexa library from the Adapter Oligo Mix and PCR primers 1.1/2.1. "
                 "INFERRED — adapter and primer bases are reconstructed; no inspected 2007 original prints them.")
SEQUENCING_INTRO = ("One single-end genomic read of 27 cycles on the Illumina cluster station and 1G analyzer (article "
                    "Methods: \u201c27 cycles of sequencing\u201d; each 27-bp read \u201cdirectionally represents one end of a DNA "
                    "fragment\u201d). No sample index is described; the unindexed layout is inferred.")


def sections():
    return [
        ("Crosslink and fragment chromatin",
         [Row(chunks=[("protein—DNA chromatin → formaldehyde crosslink → sonicated fragments", None, False)])],
         "Mechanical shearing defines fragment ends independently of the target protein."),
        ("Immunoprecipitate the target",
         antibody_enzyme_rows("bead capture", "retain target-associated fragments"),
         "A target-specific antibody enriches crosslinked protein–DNA complexes; crosslinks are then reversed."),
        ("Repair and add one 3-prime dA to each end",
         Scene.duplex([seg("3-prime dA", "A"), seg("ChIP-enriched DNA", "X" * 36, placeholder=True), seg("3-prime dA", "A")], label="dA-tailed ChIP DNA").rows(),
         "The 3-prime dA overhangs are shown before the forked adapter is introduced."),
        ("Ligate the Genomic Adapter", adapter_scene().rows(),
         "The original Illumina single-read fork has a 3′-T ligation overhang. PCR primers 1.1 and 2.1 complete the flow-cell arms."),
    ]


def protocol_graph():
    graph = ProtocolGraph(TITLE, authority=SOURCE)
    cells = graph.input("cells", "input cells", MaterialKind.BIOLOGICAL,
                        state=molecule("cells", (role("native protein-DNA chromatin"),)),
                        primary=True)
    formaldehyde = graph.input("formaldehyde", "1% formaldehyde, 10 min room temperature; 125 mM glycine stop",
                               MaterialKind.CONDITION, source=SUPPLEMENT)
    fixed = graph.product("fixed_chromatin", "formaldehyde-cross-linked chromatin",
                          MaterialKind.MOLECULAR, state=cells.state)
    graph.react("crosslink", "Cross-link protein-DNA chromatin",
                StepKind.CROSS_LINKING,
                ((cells, InputRole.TEMPLATE), (formaldehyde, InputRole.CONDITION)),
                (fixed,))
    sonication = graph.input("sonication", "nuclear lysis in RIPA; Branson 250 sonication, seven 30-s pulses, ~700 bp average",
                             MaterialKind.CONDITION, source=SUPPLEMENT)
    fragments = graph.product("chromatin_fragments", "sonicated cross-linked chromatin",
                              MaterialKind.MOLECULAR, state=fixed.state)
    graph.react("fragment", "Sonicate cross-linked chromatin", StepKind.FRAGMENTATION,
                ((fixed, InputRole.TEMPLATE), (sonication, InputRole.CONDITION)),
                (fragments,))
    antibody = graph.input("target_antibody", "anti-STAT1 rabbit polyclonal (Santa Cruz sc-345) overnight 4 °C, Protein A agarose beads, RIPA/PBS washes",
                           MaterialKind.SELECTION, source=SUPPLEMENT)
    selected = graph.product("selected_chromatin", "target-associated chromatin fragments",
                             MaterialKind.MOLECULAR, state=fragments.state)
    flowthrough = graph.product("ip_flowthrough", "unbound chromatin fragments",
                                MaterialKind.WASTE)
    graph.react("immunoprecipitate", "Immunoprecipitate target-associated chromatin",
                StepKind.AFFINITY_CAPTURE,
                ((fragments, InputRole.TEMPLATE), (antibody, InputRole.SELECTION)),
                (selected, flowthrough))
    graph.discard(flowthrough, "is washed away from the Protein A agarose immunoprecipitate")
    reversal = graph.input("crosslink_reversal", "elution in TE/1% SDS; cross-link reversal overnight at 65 °C",
                           MaterialKind.CONDITION, source=SUPPLEMENT)
    purification = graph.input(
        "dna_purification", "RNase A 1 h 37 °C, proteinase K 2 h 45 °C, phenol:chloroform:isoamyl alcohol, ethanol precipitation",
        MaterialKind.CATALYST, source=SUPPLEMENT)
    removed = graph.product("removed_protein_rna", "digested RNA, protein and organic/aqueous extraction waste",
                            MaterialKind.WASTE, source=SUPPLEMENT)
    dna = graph.product("recovered_dna", "de-cross-linked ChIP-enriched DNA",
                        MaterialKind.MOLECULAR, state=selected.state)
    graph.react("reverse_crosslinks", "Reverse cross-links and recover ChIP DNA",
                StepKind.DEPROTEINIZATION,
                ((selected, InputRole.TEMPLATE), (reversal, InputRole.CONDITION),
                 (purification, InputRole.CATALYST)),
                (dna, removed), source=SUPPLEMENT)
    graph.discard(removed, "is removed by phenol:chloroform extraction and ethanol precipitation")
    page_gel = graph.input("page_size_selection",
                           "12% PAGE of 5–50 ng immunoprecipitate; excise 100–300 bp; overnight 4 °C "
                           "elution in LoTE/7.5 M ammonium acetate (5:1); QIAquick",
                           MaterialKind.SELECTION, source=ARTICLE)
    sized = graph.product("size_selected", "100–300 bp ChIP DNA fragments",
                          MaterialKind.MOLECULAR, state=dna.state)
    rejected = graph.product("rejected_sizes", "ChIP DNA outside the excised 100–300 bp slice",
                             MaterialKind.WASTE, source=ARTICLE)
    graph.react("size_select", "Size-fractionate the immunoprecipitated DNA before library preparation",
                StepKind.SELECTION,
                ((dna, InputRole.TEMPLATE), (page_gel, InputRole.SELECTION)),
                (sized, rejected), source=ARTICLE)
    graph.discard(rejected, "remains in the polyacrylamide gel outside the excised 100–300 bp slice")
    repair = graph.input("end_repair_mix",
                         "1:5 T4:Klenow DNA polymerase mix (Illumina), 30 min 20 °C; phenol–chloroform–isoamyl "
                         "alcohol (phase-lock tube), ethanol/glycogen/NH4OAc precipitation. Guide mix also "
                         "contains T4 PNK, dNTPs and ATP ligase buffer",
                         MaterialKind.CATALYST, source=ARTICLE + "; " + GUIDE)
    repaired = graph.product("end_repaired", "end-repaired phosphorylated ChIP DNA",
                             MaterialKind.MOLECULAR, state=dna.state)
    graph.react("end_repair", "Repair and phosphorylate recovered DNA ends",
                StepKind.END_REPAIR,
                ((sized, InputRole.TEMPLATE), (repair, InputRole.CATALYST)), (repaired,),
                source=ARTICLE)
    tailing = graph.input("tailing_polymerase",
                          "Klenow exo− (3′ and 5′ exo minus; Illumina), 30 min 37 °C; phenol–chloroform "
                          "extraction, precipitation, two 70% ethanol washes, 10 µL LoTE",
                          MaterialKind.CATALYST, source=ARTICLE)
    datp = graph.input("datp", "1 mM dATP tailing nucleotide", MaterialKind.MOLECULAR, source=GUIDE)
    tailed = graph.product(
        "da_tailed", "3′-dA-tailed ChIP DNA", MaterialKind.MOLECULAR,
        state=Scene.duplex([seg("3-prime dA", "A"),
                            seg("ChIP-enriched DNA", "X" * 36, placeholder=True),
                            seg("opposite 3-prime dA", "A")], label="dA-tailed ChIP DNA"))
    graph.react("da_tail", "Add one 3′ dA to each repaired fragment end",
                StepKind.TAILING,
                ((repaired, InputRole.TEMPLATE), (datp, InputRole.STRUCTURAL),
                 (tailing, InputRole.CATALYST)), (tailed,), source=ARTICLE)
    adapter = graph.input("genomic_adapter", "Illumina adaptor oligonucleotide mix (Adapter Oligo Mix, part 1000521), diluted 1/10",
                          MaterialKind.ADAPTER, state=adapter_scene(),
                          source=ADAPTER_INFERENCE, inferred=True)
    ligase = graph.input("adapter_ligase",
                         "DNA ligase per Illumina instructions (guide: 15 min room temperature); "
                         "QIAquick, 30 µL elution buffer at 50 °C",
                         MaterialKind.CATALYST, source=ARTICLE + "; " + GUIDE)
    ligated = graph.product("adapter_ligated", "Genomic-Adapter-ligated ChIP DNA",
                            MaterialKind.MOLECULAR, state=FINAL_LIBRARY)
    graph.react("adapter_ligation", "Ligate the original Illumina Genomic Adapter",
                StepKind.LIGATION,
                ((tailed, InputRole.TEMPLATE), (adapter, InputRole.STRUCTURAL),
                 (ligase, InputRole.CATALYST)), (ligated,), source=ARTICLE,
                geometry=AdapterLigationGeometry.known(
                    "cohesive", adapter_scene(),
                    statement=("The article ligates Illumina adapters to single adenine-tailed DNA; the 2007 Illumina "
                               "guide states that the adapters carry a single 3′ T overhang that pairs with the "
                               "fragment's 3′ dA. " + ADAPTER_INFERENCE + ".")))
    primer1 = graph.input("pcr_primer_1_1", "PCR primer 1.1 (Illumina, part 1000537)",
                          MaterialKind.OLIGO, source=ARTICLE)
    primer2 = graph.input("pcr_primer_2_1", "PCR primer 2.1 (Illumina, part 1000538)",
                          MaterialKind.OLIGO, source=ARTICLE)
    polymerase = graph.input("pcr_polymerase",
                             "Phusion polymerase (Finnzymes); separate 10- and 15-cycle reactions; "
                             "QIAquick MinElute, 15 µL elution buffer at 50 °C; Agilent DNA 1000 QC",
                             MaterialKind.CATALYST, source=ARTICLE)
    library = graph.product("sequencing_library", FINAL_LIBRARY.name,
                            MaterialKind.LIBRARY, state=FINAL_LIBRARY)
    graph.react("library_pcr", "Amplify with Genomic PCR primers 1.1 and 2.1",
                StepKind.PCR,
                ((ligated, InputRole.TEMPLATE), (primer1, InputRole.STRUCTURAL),
                 (primer2, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)),
                (library,), source=ARTICLE, geometry=PCRGeometry.unavailable(
                    FINAL_LIBRARY,
                    reason="The adapter-ligated pre-PCR fork is represented through its inferred adapter end geometry rather than as one exact PCR template."))
    graph.annotate(SourceDiscrepancy(
        "size_select", "Size selection and PCR: article versus Illumina 2007 guide",
        (SourceClaim("Robertson et al. 2007, Methods \u201cIllumina library construction and sequencing\u201d",
                     "5–50 ng immunoprecipitate size-fractionated by 12% PAGE (100–300 bp) before end repair; "
                     "no post-ligation gel described; separate 10- and 15-cycle Phusion PCRs."),
         SourceClaim("Illumina part 11257047 Rev. A (2007), gel purification and PCR sections",
                     "2% agarose gel after adapter ligation, excising 200 ± 25 bp; 18 PCR cycles.")),
        "The page follows the article: PAGE size selection of the ChIP DNA first, QIAquick after "
        "ligation, then 10- or 15-cycle PCR. Guide details are used only for reagent compositions "
        "the article leaves to the manufacturer's instructions.",
        "The article describes its own library steps and cites Illumina instructions only for the "
        "enzymatic reactions; it does not name the guide.",
        "Whether the authors also performed the guide's post-ligation gel is not stated; the "
        "article does not say which of the 10- and 15-cycle products was sequenced."))
    graph.sequence(library, SequencingPlan(
        library=FINAL_LIBRARY, primers=tuple(SEQ_PRIMERS),
        acquired_roles=("Read 1",), read_lengths={"Read 1": FINAL_READ_LENGTH}))
    # Reviewed against Robertson et al. Supplementary Methods (ChIP DNA preparation).
    graph.declare_reviewed_scope((
        "cells", "formaldehyde", "crosslink", "sonication", "fragment",
        "target_antibody", "immunoprecipitate", "crosslink_reversal",
        "dna_purification", "reverse_crosslinks"), SUPPLEMENT)
    # Reviewed against the article Methods (lines 751-786); the adapter and primer
    # bases and the sequencing primer are not in the reviewed scope (no inspected
    # original prints them).
    graph.declare_reviewed_scope({
        "article_page_size_selection": "page_size_selection", "article_size_select": "size_select",
        "article_end_repair_mix": "end_repair_mix", "article_end_repair": "end_repair",
        "article_tailing_polymerase": "tailing_polymerase", "article_da_tail": "da_tail",
        "article_adapter_ligase": "adapter_ligase", "article_adapter_ligation": "adapter_ligation",
        "article_pcr_primer_1_1": "pcr_primer_1_1", "article_pcr_primer_2_1": "pcr_primer_2_1",
        "article_pcr_polymerase": "pcr_polymerase", "article_library_pcr": "library_pcr",
        # The article names the adaptor oligonucleotide mix (1/10) and the resulting
        # library, not their bases (the material is inferred).
        "article_adapter_oligo_mix": "genomic_adapter",
        "article_chip_seq_library": "sequencing_library"}, ARTICLE)
    # Guide-only reagent detail: the 1 mM dATP tailing nucleotide.
    graph.declare_reviewed_scope({"guide_datp": "datp"}, GUIDE)
    return graph.require_publishable()
