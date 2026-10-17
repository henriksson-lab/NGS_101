"""Safe-SeqS exogenous-UID library model (Kinde et al. 2011, CTNNB1 assay).

Oligo bases are transcribed from SI Table S4 (page 10 of 10, "CTNNB1 mutations in
DNA from normal human cells") as reproduced in the original-facsimile PTAB exhibit
(local download-documents.pdf, exhibit pages 0015-0016).  The table is an image, so
these strings were read visually from the rendered page (2026-10-17 bounded pass)
and are not reachable by the text-based source-oligo guard; the graph inputs are
therefore left without exact states.
"""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))

from batch_ngs import seg
from chemdraw import Construct, Scene, feature, inline_sequence, revcomp
import illumina as il
import seqprimers as sp

TITLE = "Safe-SeqS — single-strand UID-family consensus"
NOTES = "01_safe_seqs.html"
SOURCE = ('Defining source: <a href="https://doi.org/10.1073/pnas.1105422108">'
          'Kinde et al., <i>PNAS</i> 108:9530 (2011)</a>, Methods and SI Materials and '
          'Methods; oligos from SI Table S4 (CTNNB1 assay).')
SUMMARY = ("A random UID is installed by two cycles of target-specific PCR before bulk "
           "amplification. Redundant reads with one UID form a molecular family, so late "
           "PCR and sequencing errors can be rejected before rare variants are counted.")
CAVEAT = ("Oligo bases are a visual transcription of the image-only SI Table S4 (as "
          "reproduced in a patent-office exhibit of the original PNAS SI); they have one "
          "reading and are not reachable by the repository's text-based transcription "
          "guard. The page draws the CTNNB1 exogenous-UID assay (14-nt UID, index on the "
          "reverse library primer). Other Table S4 assays differ: the polymerase-fidelity "
          "and mitochondrial assays place the index on the forward (P5-side) primer, and "
          "the endogenous-UID capture and inverse-PCR libraries use ligated adapters. "
          "The CTNNB1 amplicon interior between the two target arms is not given in the "
          "inspected sources and is drawn as a placeholder.")

# --- SI Table S4, CTNNB1 block, 5'->3' as printed ----------------------------------
FWD_TAG = "CGACGTAAAACGACGGCCAGT"          # black "universal sequence", UID primer 5' end
REV_TAG = "CACACAGGAAACAGCTATGACCATG"      # black universal sequence on the reverse primer
CTNNB1_FWD_ARM = "GCAGCAACAGTCTTACCTGGACT"  # blue: complementary to template
CTNNB1_REV_ARM = "TCCACATCCTCTTCCTCAGGATT"  # blue: complementary to template
UID_LENGTH = 14                             # green N run (12 or 14 N across the study)
# Library amplification rev-1 ... rev-10: P7 graft + purple experiment index + tag.
REV_INDICES = ("ATCACG", "CGATGT", "TGACCA", "GCCAAT", "CAGATC",
               "ACTTGA", "GATCAG", "TAGCTT", "GGCTAC", "CTTGTA")
SEQ_PRIMER = "CGACGTAAAACGACGGCCAGT"       # "to read UID and internal sequences"
INDEX_PRIMER = "CATGGTCATAGCTGTTTCCTGTGTG"  # "to read experiment indexes"
# Written forms with the two 3'-terminal phosphorothioate linkages (*).
LIB_FOR_WRITTEN = il.P5 + FWD_TAG[:-2] + "*" + FWD_TAG[-2] + "*" + FWD_TAG[-1]
LIB_REV_WRITTEN = (il.P7 + "{index}" + REV_TAG[:-2] + "*" + REV_TAG[-2]
                   + "*" + REV_TAG[-1])
if SEQ_PRIMER != FWD_TAG:
    raise ValueError("Table S4 sequencing primer should equal the forward universal tag")
if INDEX_PRIMER != revcomp(REV_TAG):
    raise ValueError("Table S4 index primer should be the reverse complement of the reverse tag")

UID = feature("safe_seqs_uid", "umi", "random",
              note="14 random nucleotides 5' of the CTNNB1 forward arm (Table S4)")
SAMPLE_INDEX = feature("safe_seqs_experiment_index", "sample_index", "fixed",
                       note="one of ten 6-nt indices on library amplification rev-1..rev-10")
_SRC = "Kinde et al. 2011 SI Table S4 (CTNNB1 block)"


def _uid():
    return seg(f"{UID_LENGTH}-nt UID", "N" * UID_LENGTH, "umi", placeholder=True, feature=UID)


def _interior():
    return seg("CTNNB1 amplicon interior (not given)", "X" * 24, placeholder=True)


def final_library() -> Construct:
    return Construct([
        seg("P5 grafting sequence", il.P5, "p5"),
        seg("forward universal tag", FWD_TAG, "r1"),
        _uid(),
        seg("CTNNB1 forward arm", CTNNB1_FWD_ARM, "me"),
        _interior(),
        seg("CTNNB1 reverse arm (rc)", revcomp(CTNNB1_REV_ARM), "me"),
        seg("reverse universal tag (rc)", revcomp(REV_TAG), "r2"),
        seg("experiment index (rc)", "I" * 6, "cbc", placeholder=True, feature=SAMPLE_INDEX),
        seg("P7 grafting sequence (rc)", il.P7_RC, "p7"),
    ], name="Safe-SeqS CTNNB1 sequencing library")


FINAL_LIBRARY = final_library()
SEQ_PRIMERS = (
    sp.custom("Read 1", "Safe-SeqS sequencing primer", SEQ_PRIMER, _SRC,
              "Reads the UID and then the CTNNB1 target."),
    sp.custom("Index 1 (i7)", "Safe-SeqS index primer", INDEX_PRIMER, _SRC,
              "Reads the experiment index adjacent to the P7 graft."),
)
_problems = sp.verify(FINAL_LIBRARY, SEQ_PRIMERS,
                      required_roles=tuple(p.role for p in SEQ_PRIMERS))
if _problems:
    raise ValueError("; ".join(_problems))
FINAL_CAPTION = ("Library as installed by the Table S4 primers: P5 graft, forward universal "
                 "tag, 14-nt UID, CTNNB1 amplicon, reverse universal tag, 6-nt experiment "
                 "index and P7 graft. The amplicon interior is a placeholder.")
SEQUENCING_INTRO = ("Illumina GA IIx; the SI reports read lengths of 36–73 bases across "
                    "experiments without assigning one to this assay, so no cycle count is "
                    "drawn. Table S4 lists a custom sequencing primer (equal to the forward "
                    "tag) and an index primer (reverse complement of the reverse tag); no "
                    "second genomic read is described.")


def oligos():
    rev1 = LIB_REV_WRITTEN.format(index=REV_INDICES[0])
    return (
        "<p>UID assignment amplification, for: "
        + inline_sequence(FWD_TAG + "N" * UID_LENGTH + CTNNB1_FWD_ARM) + "</p>",
        "<p>UID assignment amplification, rev: "
        + inline_sequence(REV_TAG + CTNNB1_REV_ARM) + "</p>",
        "<p>Library amplification, for: " + inline_sequence(LIB_FOR_WRITTEN) + "</p>",
        "<p>Library amplification, rev-1: " + inline_sequence(rev1)
        + " (rev-2…rev-10 replace the index " + REV_INDICES[0] + " with "
        + ", ".join(REV_INDICES[1:]) + ")</p>",
        "<p>Sequencing primer: " + inline_sequence(SEQ_PRIMER)
        + " · Index primer: " + inline_sequence(INDEX_PRIMER) + "</p>",
        "<p>* phosphorothioate linkage (Table S4 legend).</p>",
    )


def assignment_scene() -> Scene:
    template = [seg("CTNNB1 target body", "X" * 16, placeholder=True),
                seg("forward-arm site", revcomp(CTNNB1_FWD_ARM), "me")]
    primer = [seg("forward universal tag", FWD_TAG, "r1"), _uid(),
              seg("CTNNB1 forward arm", CTNNB1_FWD_ARM, "me")]
    sc = Scene()
    sc.strand("template", template, label="one CTNNB1 template strand")
    sc.anneal("UID primer", primer, to="template",
              pair=("CTNNB1 forward arm", "forward-arm site"),
              label="UID assignment primer (for)",
              unpaired=("forward universal tag", f"{UID_LENGTH}-nt UID"), above=True)
    sc.arrow("UID primer", "2 cycles, Phusion Hot Start II")
    sc.mark("UID primer", f"{UID_LENGTH}-nt UID", "one UID per template strand")
    sc.labels("UID primer")
    return sc


def family_scene() -> Scene:
    parts = [seg("forward universal tag", FWD_TAG, "r1"), _uid(),
             seg("CTNNB1 forward arm", CTNNB1_FWD_ARM, "me"), _interior(),
             seg("CTNNB1 reverse arm (rc)", revcomp(CTNNB1_REV_ARM), "me"),
             seg("reverse universal tag (rc)", revcomp(REV_TAG), "r2")]
    sc = Scene.duplex(parts, label="UID-tagged CTNNB1 fragment")
    sc.mark("top", f"{UID_LENGTH}-nt UID", "copied into every family member")
    sc.labels("top")
    return sc


from pcr import PCRGeometry
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind
from protocol_model import molecule, role


def safe_seqs_graph(authority: str, final: Construct, primers) -> ProtocolGraph:
    graph = ProtocolGraph(TITLE, authority=authority)
    target = graph.input("target_dna", "genomic DNA (normal human cells), quantified by digital PCR",
                         MaterialKind.BIOLOGICAL,
                         state=molecule("target DNA", (role("CTNNB1 target strands"),)), primary=True)
    uid_primer = graph.input("uid_primer", "UID assignment amplification, for (tag-N14-CTNNB1 arm)",
                             MaterialKind.OLIGO)
    opposing = graph.input("opposing_primer", "UID assignment amplification, rev (tag-CTNNB1 arm)",
                           MaterialKind.OLIGO)
    assignment_pol = graph.input("assignment_polymerase", "Phusion Hot Start II (2 U, 1x HF buffer)",
                                 MaterialKind.CATALYST)
    assigned = graph.product("uid_assigned", "UID-tagged CTNNB1 fragments",
                             MaterialKind.MOLECULAR, state=assignment_scene())
    graph.react("assign_uid",
                "Assign UIDs: 98 °C 30 s, then 2 cycles of 98 °C 10 s / 61 °C 120 s / 72 °C 10 s",
                StepKind.PCR,
                ((target, InputRole.TEMPLATE), (uid_primer, InputRole.STRUCTURAL),
                 (opposing, InputRole.STRUCTURAL),
                 (assignment_pol, InputRole.CATALYST)), (assigned,),
                geometry=PCRGeometry.unavailable(
                    final, reason="The CTNNB1 amplicon interior is not given in the inspected sources."))
    exoi = graph.input("exonuclease_i", "Exonuclease I (60 U, Enzymatics)",
                       MaterialKind.CATALYST)
    cleared = graph.product("primer_cleared", "primer-cleared UID-tagged fragments",
                            MaterialKind.MOLECULAR, state=family_scene())
    excess = graph.product("digested_primers", "digested excess UID-assignment primers",
                           MaterialKind.WASTE)
    graph.react("remove_excess_primers",
                "Digest single-stranded primers with Exonuclease I, 37 °C 1 h; 98 °C 5 min",
                StepKind.DIGESTION,
                ((assigned, InputRole.TEMPLATE), (exoi, InputRole.CATALYST)),
                (cleared, excess))
    graph.discard(excess, "single-stranded primers are degraded in the well; no purification follows")
    universal_f = graph.input("universal_forward",
                              "Library amplification, for (P5-tag, 3′ phosphorothioates)",
                              MaterialKind.OLIGO)
    universal_r = graph.input("universal_reverse",
                              "Library amplification, rev-1..10 (P7-index-tag, 3′ phosphorothioates)",
                              MaterialKind.OLIGO)
    family_pol = graph.input("family_polymerase", "Phusion polymerase carried over from UID assignment",
                             MaterialKind.CATALYST)
    amplified = graph.product("amplified", "redundantly amplified UID families",
                              MaterialKind.MOLECULAR, state=final)
    graph.react("amplify_families",
                "Add 0.5 µM library primers; 25 cycles of 98 °C 10 s / 72 °C 15 s",
                StepKind.PCR,
                ((cleared, InputRole.TEMPLATE), (universal_f, InputRole.STRUCTURAL),
                 (universal_r, InputRole.STRUCTURAL), (family_pol, InputRole.CATALYST)),
                (amplified,), geometry=PCRGeometry.unavailable(
                    final, reason="The CTNNB1 amplicon interior is not given in the inspected sources."))
    gel = graph.input("page_gel", "QIAquick column, then 10% TBE polyacrylamide gel",
                      MaterialKind.SELECTION)
    library = graph.product("library", "size-selected Safe-SeqS CTNNB1 library",
                            MaterialKind.LIBRARY, state=final)
    rejected = graph.product("gel_rejected", "fragments outside the excised band",
                             MaterialKind.WASTE)
    graph.react("purify_library",
                "Pool wells, QIAquick purify, excise expected-size band from 10% TBE PAGE",
                StepKind.CLEANUP,
                ((amplified, InputRole.TEMPLATE), (gel, InputRole.SELECTION)),
                (library, rejected))
    graph.discard(rejected, "left in the gel; size composition not represented")
    graph.sequence(library, SequencingPlan(
        library=final, primers=tuple(primers),
        acquired_roles=tuple(p.role for p in primers)))
    graph.declare_reviewed_scope([
        "target_dna", "uid_primer", "opposing_primer", "assignment_polymerase",
        "assign_uid", "exonuclease_i", "primer_cleared", "digested_primers",
        "remove_excess_primers", "universal_forward", "universal_reverse",
        "family_polymerase", "amplify_families", "amplified", "page_gel",
        "purify_library", "gel_rejected", "library"], authority)
    return graph.require_publishable()


def protocol_graph():
    return safe_seqs_graph(SOURCE, FINAL_LIBRARY, SEQ_PRIMERS)
