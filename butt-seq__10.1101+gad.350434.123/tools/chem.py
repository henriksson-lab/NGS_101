"""Butt-seq source oligos, blocked TGIRT starters and reaction graph."""
import re

from chemdraw import Scene, Segment, feature, oligo, revcomp
import illumina as il
from pcr import PCRGeometry
from protocol_graph import (AnnotationSlot, DetailView, InputRole, MaterialKind,
                            ProtocolGraph, SequencingPlan, SourceClaim,
                            SourceDiscrepancy, StepKind)
from protocol_model import molecule, role, role_library
from extended_protocols import Definition, UMI, SAMPLE

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
PAYLOAD = (UMI, role("nascent-RNA terminal insert"))

TITLE = DEFINITION.title
NOTES = DEFINITION.notes
SOURCE = f'Defining source: <a href="https://doi.org/{DEFINITION.doi}">{DEFINITION.citation}</a>.'
SUMMARY = DEFINITION.summary
CAVEAT = ("The starter drawings show the A-ended member of each source-published "
          "A/T/C/G equimolar DNA pool. Both RNA donors have a 3′ /3SpC3/ blocker; "
          "the 3′ DNA base remains a one-nucleotide overhang.")
FINAL_LIBRARY = role_library(TITLE + " final library roles", PAYLOAD, index=SAMPLE)
SEQ_PRIMERS = ()
FINAL_CAPTION = DEFINITION.endpoint
SEQUENCING_INTRO = DEFINITION.endpoint
READ_LENGTHS = None
SEQUENCING_ENDING = (
    "The author analysis pipeline expects three FASTQ files: R1 is Read 1, "
    "R2 carries the eight-cycle UMI, and R3 is the biological Read 2. Its README "
    "instructs changing the second RunInfo.xml read from indexed to non-indexed "
    "before demultiplexing. <a href=\"https://github.com/albertdyu/BuTTSeq\">"
    "Author analysis pipeline</a>. Exact sequencing-primer landing sites remain "
    "unplaced in this role-level construct.")
SEQUENCING_UNAVAILABLE = ("The source publishes the starter and PCR oligos, but this "
                          "role-level final state does not yet establish the complete "
                          "post-PCR sequence or exact sequencing-primer landing sites.")

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
    alkali = graph.input("alkali", "NaOH followed by HCl neutralization", MaterialKind.CONDITION)
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
        state=molecule("double-stranded Butt-seq product", (
            UMI, role("nascent-RNA insert"))))
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
                (amplified,), geometry=PCRGeometry.unavailable(
                    final,
                    reason="The primer oligos are published, but the complete pre-PCR duplex "
                           "and terminal overhang choice are not yet represented at base level."))
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
    reviewed.extend(("pcr_forward", "pcr_reverse", "pcr_polymerase", "library_pcr",
                     "size_selection", "select_size", "sequencing_library"))
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()
