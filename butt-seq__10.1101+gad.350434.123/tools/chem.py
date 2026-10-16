"""Butt-seq source oligos, blocked TGIRT starters and reaction graph."""
import re

from chemdraw import Construct, Scene, Segment, feature, oligo, revcomp
import illumina as il
import nextera as nx
import seqprimers as sp
from pcr import PCRGeometry, PCRPrimer
from protocol_graph import (AnnotationSlot, DetailView, InputRole, MaterialKind,
                            ProtocolGraph, SequencingPlan, SourceClaim,
                            SourceDiscrepancy, StepKind)
from protocol_model import molecule, role
from extended_protocols import Definition, UMI

# Main article, Materials and Methods, and Supplemental_Primers.xlsx, TGIRT Primers.
SCR2_RNA = "rArArGrArUrCrGrGrArArGrArGrCrArCrArCrGrUrCrUrGrArArCrUrCrCrArGrUrCrArC/3SpC3/"
R1ME_RNA = "rCrUrGrUrCrUrCrUrUrArUrArCrArCrArUrCrUrGrArCrGrCrUrGrCrCrGrArCrGrA/3SpC3/"
SCR2_DNA_CORE = "CAAGCAGAAGACGGCATACGAGATNNNNNNNNGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTT"
R1ME_DNA_CORE = "TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG"
SCR2_DNA_POOL = tuple(SCR2_DNA_CORE + base for base in "ATCG")
R1ME_DNA_POOL = tuple(R1ME_DNA_CORE + base for base in "ATCG")
SCP7 = "CAAGCAGAAGACGGCATACGAG"
AD1_60 = "AATGATACGGCGACCACCGAGATCTACACTACTGACCTCGTCGGCAGCGTCAGATGTGTAT"


def _rna_bases(written):
    """The RNA base sequence used for pairing; keep modifications in ordering text."""
    return "".join(re.findall(r"r([ACGU])", written)).replace("U", "T")


def starter_scene(rna_written, dna_core, name):
    """Check the published 3′ blocked RNA against a representative A-overhang DNA."""
    rna_bases = _rna_bases(rna_written)
    if not rna_written.endswith("/3SpC3/"):
        raise ValueError("TGIRT starter RNA lacks its source-published 3′ spacer")
    if not dna_core.endswith(revcomp(rna_bases)):
        raise ValueError("starter DNA does not pair with its RNA donor")
    prefix = dna_core[:-len(rna_bases)]
    if prefix not in ("", il.P7 + "N" * 8):
        raise ValueError("unexpected unpaired starter DNA handle")
    tag = "r2" if prefix else "r1"
    dna = ([Segment("P7", il.P7, "p7"),
            Segment("8-base UMI", "N" * 8, "umi", placeholder=True,
                    feature=feature("butt_seq.umi", "umi", "random"))]
           if prefix else [])
    dna.extend((Segment("paired DNA core", revcomp(rna_bases), tag),
                Segment("one-base 3′ overhang", "A")))
    sc = Scene()
    sc.strand("RNA donor", [Segment("blocked RNA donor", rna_bases, tag)],
              label=name + " RNA, 5′→3′", mod3="3SpC3")
    sc.anneal("DNA primer", dna,
              to="RNA donor", pair=("paired DNA core", "blocked RNA donor"),
              label=name + " DNA, A variant")
    sc.labels("RNA donor")
    sc.labels("DNA primer")
    return sc


def first_starter():
    return starter_scene(SCR2_RNA, SCR2_DNA_CORE, "SCR2R")


def second_starter():
    return starter_scene(R1ME_RNA, R1ME_DNA_CORE, "R1MER")

DEFINITION = Definition('Butt-seq', '10.1101/gad.350434.123', 'Yu and Rosbash, Genes & Development (2023)', 'TGIRT template switching captures native nascent-RNA 3′ termini without RNA fragmentation, followed by second-strand synthesis and indexed PCR.', 'The terminal insert base reports the nascent RNA 3′ end; a separate UMI read supports molecule deduplication.', '01_butt_seq.html')
TITLE = DEFINITION.title
NOTES = DEFINITION.notes
SOURCE = f'Defining source: <a href="https://doi.org/{DEFINITION.doi}">{DEFINITION.citation}</a>.'
SUMMARY = DEFINITION.summary
CAVEAT = ("The starter drawings show the A-ended member of each source-published "
          "A/T/C/G equimolar DNA pool. Both RNA donors have a 3′ /3SpC3/ blocker; "
          "the 3′ DNA base remains a one-nucleotide overhang. The final duplex "
          "places source-published handles and an example Ad1.60 i5 index; insert "
          "sequence and the two overhang bases vary by molecule. Primer identities "
          "are sequence-compatible standard Illumina entries, not instrument oligos "
          "individually named by the paper.")


def _known_library_parts():
    """P5-facing strand implied by the two source-printed starters and PCR primers."""
    donor = SCR2_DNA_CORE[len(il.P7) + 8:]
    if (not SCR2_DNA_CORE.startswith(il.P7 + "N" * 8)
            or donor != revcomp(_rna_bases(SCR2_RNA))
            or R1ME_DNA_CORE != nx.S5 + nx.ME
            or AD1_60[:len(il.P5)] != il.P5
            or AD1_60[len(il.P5) + 8:] != (nx.S5 + nx.ME)[:24]
            or SCP7 != il.P7[:len(SCP7)]):
        raise ValueError("Butt-seq source oligos no longer support the inferred library ends")
    return [Segment("Nextera S5", nx.S5, "s5"),
            Segment("mosaic end", nx.ME, "me"),
            Segment("second-starter 3′ overhang", "N", placeholder=True, length_bp=1),
            role("nascent-RNA insert"),
            Segment("first-starter 3′ overhang", "N", placeholder=True, length_bp=1),
            Segment("SCR2 read-2-side donor complement", revcomp(donor), "r2"),
            Segment("8-base UMI complement", "N" * 8, "umi", placeholder=True,
                    feature=feature("butt_seq.umi", "umi", "random")),
            Segment("P7 complement", revcomp(il.P7), "p7")]


def pre_pcr_library():
    return Construct(_known_library_parts(), "Butt-seq second-strand PCR template")


def final_library():
    index = AD1_60[len(il.P5):len(il.P5) + 8]
    return Construct([Segment("P5", il.P5, "p5"),
                      Segment("Ad1.60 i5 index", index, "i5",
                              feature=feature("butt_seq.i5", "sample_index", "whitelist",
                                              whitelist="Butt Supplemental_Primers.xlsx v2_Ad1 panel"))]
                     + _known_library_parts(), TITLE + " final library (Ad1.60 example)")


FINAL_LIBRARY = final_library()
SEQ_PRIMERS = (sp.NEXTERA["R1"], sp.TRUSEQ["I1"],
               sp.NEXTERA["I2"], sp.TRUSEQ["R2"])
FINAL_CAPTION = DEFINITION.endpoint
SEQUENCING_INTRO = ("The first index acquisition is an eight-cycle UMI read, output "
                   "as FASTQ R2 after the author changes RunInfo.xml to non-indexed. "
                   "The second index acquisition reads the Ad1 i5 sample index. "
                   "Biological Read 1 and Read 2 are FASTQ R1 and R3, respectively. "
                   "The table places standard primer sequences on the constructed "
                   "duplex; the paper does not name the physical instrument oligos.")
READ_LENGTHS = {"Index 1 (i7)": 8, "Index 2 (i5)": 8}
def oligos():
    rows = []
    for name, sequence in (("SCR2R_RNA, 3′ blocked", SCR2_RNA),
                           ("R1MER_RNA, 3′ blocked; supplemental workbook version", R1ME_RNA)):
        rows.append(oligo(name, [Segment(name, sequence, placeholder=True)], three=""))
    for family, pool in (("SCR2_DNA", SCR2_DNA_POOL), ("R1ME", R1ME_DNA_POOL)):
        for base, sequence in zip("ATCG", pool):
            rows.append(oligo(f"{family}_{base}", [Segment(family, sequence,
                                                            placeholder="N" in sequence)]))
    rows.append(oligo("SCP7 PCR primer", [Segment("SCP7", SCP7)]))
    rows.append(oligo("v2_Ad1.60 indexed PCR primer, example", [Segment("Ad1.60", AD1_60)]))
    return tuple(rows)


def protocol_graph() -> ProtocolGraph:
    d = DEFINITION
    final = FINAL_LIBRARY
    source = f"{d.citation}; doi:{d.doi}"
    graph = ProtocolGraph(d.title, authority=source)
    rna = graph.input(
        "nascent_rna", "purified chromatin-associated nascent RNA",
        MaterialKind.BIOLOGICAL,
        state=molecule("purified chromatin-associated nascent RNA", (
            role("nascent RNA 5-prime body"), role("native RNA 3-prime end"))),
        primary=True)
    starter = graph.input("starter_duplex", "SCR2R blocked RNA/DNA starter; A/T/C/G overhang pool",
                          MaterialKind.OLIGO, state=first_starter())
    graph.annotate(DetailView(
        "starter_duplex", "First TGIRT starter: published SCR2R duplex",
        "The RNA has a 3′ C3 spacer. A is one representative DNA overhang; "
        "the experiment mixes A, T, C and G equally."))
    tgirt = graph.input("tgirt", "TGIRT reverse transcriptase", MaterialKind.CATALYST)
    first = graph.product(
        "first_strand", "first-strand reaction mixture with nascent-RNA cDNA",
        MaterialKind.MOLECULAR,
        state=molecule("first-strand cDNA", (
            UMI, role("nascent RNA 3-prime-end cDNA"))))
    graph.react("template_switch", "Template-switch TGIRT onto the native RNA 3-prime end",
                StepKind.TEMPLATE_SWITCHING,
                ((rna, InputRole.TEMPLATE), (starter, InputRole.STRUCTURAL),
                 (tgirt, InputRole.CATALYST)), (first,))
    exonuclease = graph.input("exonuclease_iii", "exonuclease III", MaterialKind.CATALYST)
    exo_cleaned = graph.product("exo_cleaned", "first-strand reaction after unused-starter digestion",
                                MaterialKind.MOLECULAR, state=first.state)
    primer_fragments = graph.product("unused_starter_fragments", "digested unused starter primers",
                                     MaterialKind.WASTE)
    graph.react("digest_unused_starters", "Partially digest unused starter primers with exonuclease III",
                StepKind.DIGESTION,
                ((first, InputRole.TEMPLATE), (exonuclease, InputRole.CATALYST)),
                (exo_cleaned, primer_fragments))
    graph.discard(primer_fragments, "unused starters are partially digested before second-strand synthesis")
    alkali = graph.input("alkali", "1 N NaOH followed by 1 N HCl neutralization", MaterialKind.CONDITION)
    graph.annotate(SourceDiscrepancy(
        "alkali", "Preprint and published alkaline-hydrolysis concentrations differ",
        (SourceClaim("bioRxiv v1 10.1101/2023.03.04.531084, library-prep Methods lines 522–523",
                     "3 µL of 1.3 M NaOH followed by 3 µL of 1.3 M HCl."),
         SourceClaim("Published article PMC10270195 Materials and Methods, RNA-removal paragraph; "
                     "Supplemental_Detailed_Protocol.docx lines 289–301",
                     "3 µL of 1 M (1 N) NaOH followed by 3 µL of 1 M (1 N) HCl.")),
        "Use the published article and detailed-protocol 1 N reagents.",
        "The published article and its practical protocol agree on 1 N for both reagents.",
        "The sources do not say whether the earlier 1.3 M condition was tested or changed during revision."))
    dna_only = graph.product("rna_removed", "first-strand cDNA after alkaline RNA removal",
                             MaterialKind.MOLECULAR, state=first.state)
    removed_rna = graph.product("removed_rna", "hydrolyzed RNA and dissociated TGIRT",
                                MaterialKind.WASTE)
    graph.react("remove_rna", "Hydrolyze RNA and dissociate TGIRT with alkali",
                StepKind.CHEMICAL_CONVERSION,
                ((exo_cleaned, InputRole.TEMPLATE), (alkali, InputRole.CONDITION)),
                (dna_only, removed_rna))
    graph.discard(removed_rna, "RNA is removed before the second TGIRT reaction")
    second_primer = graph.input("second_strand_primer", "R1MER blocked RNA/DNA starter; A/T/C/G overhang pool",
                                MaterialKind.OLIGO, state=second_starter())
    graph.annotate(DetailView(
        "second_strand_primer", "Second TGIRT starter: supplemental R1MER duplex",
        "The supplemental primer workbook gives this full-length, 3′ C3-blocked RNA. "
        "A is one representative DNA overhang of the equimolar A/T/C/G pool."))
    graph.annotate(SourceDiscrepancy(
        "second_strand_primer", "Published second-starter sequences disagree",
        (SourceClaim("Article Materials and Methods, PMC10270195 paragraph beginning 'To prepare MER1R primer'",
                     "R1ME RNA is 27 nt before /3SpC3/ and the paired DNA is GCAGCGTCAGATGTGTATAAGAGACAGN (28 nt)."),
         SourceClaim("Supplemental_Primers.xlsx, TGIRT Primers, R1ME_A/T/C/G and R1MER_RNA rows",
                     "R1MER RNA is 33 nt before /3SpC3/ and the DNA family is 34 nt, with a TCGTCG 5′ extension.")),
        "Draw the full-length supplemental workbook variant.",
        "The detailed supplemental protocol explicitly directs using the primer-table "
        "sequences and hand-mixing the four terminal DNA variants.",
        "The article does not explain whether its shorter second starter was an earlier "
        "working variant or a truncation in the Methods text."))
    double = graph.product(
        "double_stranded", "double-stranded Butt-seq product",
        MaterialKind.MOLECULAR,
        state=pre_pcr_library())
    graph.react("second_strand", "Second TGIRT template switch and strand synthesis",
                StepKind.TEMPLATE_SWITCHING,
                ((dna_only, InputRole.TEMPLATE), (second_primer, InputRole.STRUCTURAL),
                 (tgirt, InputRole.CATALYST)), (double,))
    reviewed = ["nascent_rna", "starter_duplex", "tgirt", "template_switch",
                "exonuclease_iii", "digest_unused_starters", "alkali", "remove_rna",
                "second_strand_primer", "second_strand"]
    forward = graph.input("pcr_forward", "SCP7 universal PCR primer", MaterialKind.OLIGO,
                          state=SCP7)
    reverse = graph.input("pcr_reverse", "v2_Ad1 indexed PCR primer family (Ad1.60 shown)",
                          MaterialKind.OLIGO, state=AD1_60)
    polymerase = graph.input("pcr_polymerase", "library PCR polymerase", MaterialKind.CATALYST)
    amplified = graph.product("amplified_library", "amplified Butt-seq library",
                              MaterialKind.MOLECULAR, state=final)
    graph.react("library_pcr", "PCR-amplify and index", StepKind.PCR,
                ((double, InputRole.TEMPLATE), (forward, InputRole.STRUCTURAL),
                 (reverse, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)),
                (amplified,), geometry=PCRGeometry.exact(
                    double.state, PCRPrimer("v2_Ad1.60", AD1_60, anneal_nt=24),
                    PCRPrimer("SCP7", SCP7),
                    statement="Primer binding is derived from the published workbook; "
                              "Ad1.60 is one member of the 25-member index panel, and "
                              "both one-base TGIRT overhangs vary by molecule."))
    size_rule = graph.input("size_selection", "Butt-seq library-size selection rule",
                            MaterialKind.SELECTION)
    library = graph.product("sequencing_library", "size-selected Butt-seq library",
                            MaterialKind.LIBRARY, state=final)
    rejected = graph.product("rejected_sizes", "off-size Butt-seq products",
                             MaterialKind.WASTE)
    graph.react("select_size", "Size-select the indexed library", StepKind.SELECTION,
                ((amplified, InputRole.TEMPLATE), (size_rule, InputRole.SELECTION)),
                (library, rejected))
    graph.discard(rejected, "falls outside the selected nascent-RNA library size")
    graph.sequence(library, SequencingPlan(limitation=d.endpoint))
    graph.annotate(SourceDiscrepancy(
        "sequencing_library", "First Read 2 base depends on an unreported run primer",
        (SourceClaim("Yu and Rosbash 2023, Results paragraph beginning 'To analyze transcriptional dynamics'",
                     "The first sequenced Read 2 base should reflect the last nascent-RNA base."),
         SourceClaim("Supplemental_Primers.xlsx SCR2_DNA family and Illumina Adapter Sequences "
                     "#1000000002694 TruSeq Read 2 sequencing primer",
                     "The SCR2 donor is 35 bases; the standard TruSeq Read 2 primer pairs "
                     "with 34 bases and would read a fixed T before the variable overhang.")),
        "Show the standard-primer-compatible landing and its fixed-T first base.",
        "The starter sequence is source-printed and the 34-base canonical primer is "
        "sequence-compatible with the constructed library.",
        "The paper does not identify the physical Read 2 oligo or provide Primers.xml. "
        "A 35-base custom primer could skip the fixed T; the reported first-base "
        "interpretation cannot be confirmed from the available run records.",
        slot=AnnotationSlot.SEQUENCING))
    reviewed.extend(("pcr_forward", "pcr_reverse", "pcr_polymerase", "library_pcr",
                     "size_selection", "select_size", "sequencing_library"))
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()
