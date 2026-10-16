"""Safe-SeqS targeted exogenous-UID library model."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))

from batch_ngs import duplex_fragment, seg, truseq_amplicon_library
from chemdraw import MolecularState, Workflow, feature
from molecular_consensus import SingleStrandUID

TITLE = "Safe-SeqS — single-strand UID-family consensus"
NOTES = "01_safe_seqs.html"
SOURCE = ('Defining source: <a href="https://doi.org/10.1073/pnas.1105422108">'
          'Kinde et al., <i>PNAS</i> (2011)</a>.')
SUMMARY = ("A random UID is installed before bulk target amplification. Redundant reads "
           "with one UID form a molecular family, so late PCR and sequencing errors can be "
           "rejected before rare variants are counted.")
CAVEAT = ("The original supplement lists the study's ordered oligos in Table S4, but its "
          "sequence cells have not yet been independently transcribed. This representative "
          "14-nt-UID amplicon uses an INFERRED — TruSeq-compatible sequencing shell; its "
          "adapter bases and read-primer landing sites are not verified for Safe-SeqS. "
          "The assay-specific target arms are shown by molecular role.")

UID = feature("safe_seqs_uid", "umi", "random")
FINAL_LIBRARY, SEQ_PRIMERS = truseq_amplicon_library([
    seg("14-nt UID", "N" * 14, "umi", placeholder=True, feature=UID),
    seg("target amplicon", "X" * 42, placeholder=True)],
    "Safe-SeqS representative sequencing library", dual_index=False,
    inferred_adapters=True)
FINAL_CAPTION = ("INFERRED — in this representative orientation, the first 14 insert "
                 "bases identify one amplified template strand; the remaining insert "
                 "is the target amplicon used for family consensus.")
SEQUENCING_INTRO = ("In this inferred representative orientation, Read 1 encounters the "
                    "14-nt UID before the target. Primer sites are computed from the "
                    "modeled shell, not yet verified against source Table S4.")


def workflow():
    chemistry = SingleStrandUID(14)
    w = Workflow(MolecularState("Target DNA", tuple(duplex_fragment("target DNA").rows())))
    w.react("Assign a random UID with two target-specific PCR cycles",
            chemistry.assignment().rows(), name="UID-assigned target",
            note="Each original strand receives a primer-borne random identifier before bulk amplification.")
    w.react("Digest excess UID-assignment primers with Exonuclease I",
            chemistry.assignment().rows(), name="primer-cleared UID-assigned target",
            note="The original supplemental method uses Exonuclease I before universal amplification; the exact oligo geometry still needs Table S4 transcription.")
    w.react("Amplify through the universal handles",
            chemistry.family().rows(), name="UID family",
            note="All redundant descendants of one early tagged strand carry the same UID.")
    w.react("Group reads and form a family consensus",
            chemistry.family().rows(), name="single-strand consensus family",
            note="The defining analysis calls a supermutant only when at least 95% of the family agrees.")
    return w


def oligos():
    return ()



from pcr import PCRGeometry
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind
from protocol_model import molecule, role

def safe_seqs_graph(authority: str, final: Construct, primers) -> ProtocolGraph:
    graph = ProtocolGraph("Safe-SeqS — single-strand UID-family consensus",
                          authority=authority)
    target = graph.input("target_dna", "target DNA containing assay primer sites",
                         MaterialKind.BIOLOGICAL,
                         state=molecule("target DNA", (role("target strands"),)), primary=True)
    uid_primer = graph.input("uid_primer", "target-specific primer carrying a random UID",
                             MaterialKind.OLIGO)
    opposing = graph.input("opposing_primer", "opposing target-specific primer",
                           MaterialKind.OLIGO)
    assignment_pol = graph.input("assignment_polymerase", "UID-assignment PCR polymerase",
                                 MaterialKind.CATALYST)
    assigned = graph.product("uid_assigned", "UID-assigned target strands",
                             MaterialKind.MOLECULAR, state=SingleStrandUID(14).assignment())
    graph.react("assign_uid", "Assign a random UID in two target-specific PCR cycles",
                StepKind.PCR,
                ((target, InputRole.TEMPLATE), (uid_primer, InputRole.STRUCTURAL),
                 (opposing, InputRole.STRUCTURAL),
                 (assignment_pol, InputRole.CATALYST)), (assigned,),
                geometry=PCRGeometry.unavailable(
                    final, reason="The assay-specific target template and opposing site vary by panel."))
    exoi = graph.input("exonuclease_i", "Exonuclease I for excess primer removal",
                       MaterialKind.CATALYST)
    cleared = graph.product("primer_cleared", "primer-cleared UID-assigned strands",
                            MaterialKind.MOLECULAR, state=SingleStrandUID(14).assignment())
    excess = graph.product("digested_primers", "digested excess UID-assignment primers",
                           MaterialKind.WASTE)
    graph.react("remove_excess_primers", "Digest excess UID primers with Exonuclease I",
                StepKind.DIGESTION,
                ((assigned, InputRole.TEMPLATE), (exoi, InputRole.CATALYST)),
                (cleared, excess))
    graph.discard(excess, "excess primers are degraded before universal amplification")
    universal_f = graph.input("universal_forward", "universal-handle forward primer",
                              MaterialKind.OLIGO)
    universal_r = graph.input("universal_reverse", "universal-handle reverse primer",
                              MaterialKind.OLIGO)
    family_pol = graph.input("family_polymerase", "UID-family amplification polymerase",
                             MaterialKind.CATALYST)
    library = graph.product("library", "redundantly amplified Safe-SeqS UID families",
                            MaterialKind.LIBRARY, state=final)
    graph.react("amplify_families", "Amplify through the universal handles",
                StepKind.PCR,
                ((cleared, InputRole.TEMPLATE), (universal_f, InputRole.STRUCTURAL),
                 (universal_r, InputRole.STRUCTURAL), (family_pol, InputRole.CATALYST)),
                (library,), geometry=PCRGeometry.unavailable(
                    final, reason="The source Table S4 shell and target-specific insert coordinates remain unverified."))
    graph.sequence(library, SequencingPlan(library=final, primers=tuple(primers)))
    graph.declare_reviewed_scope([
        "target_dna", "uid_primer", "opposing_primer", "assignment_polymerase",
        "assign_uid", "exonuclease_i", "primer_cleared", "digested_primers",
        "remove_excess_primers", "universal_forward", "universal_reverse",
        "family_polymerase", "amplify_families", "library"], authority)
    return graph.require_publishable()


def protocol_graph():
    return safe_seqs_graph(SOURCE, FINAL_LIBRARY, SEQ_PRIMERS)
