"""Safe-SeqS targeted exogenous-UID library model."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))

from batch_ngs import duplex_fragment, seg, truseq_library
from chemdraw import MolecularState, Workflow, feature
from molecular_consensus import SingleStrandUID

TITLE = "Safe-SeqS — single-strand UID-family consensus"
NOTES = "01_safe_seqs.html"
SOURCE = ('Defining source: <a href="https://doi.org/10.1073/pnas.1105422108">'
          'Kinde et al., <i>PNAS</i> (2011)</a>.')
SUMMARY = ("A random UID is installed before bulk target amplification. Redundant reads "
           "with one UID form a molecular family, so late PCR and sequencing errors can be "
           "rejected before rare variants are counted.")
CAVEAT = ("The defining paper allows assay-specific target primers and illustrates a 14-nt "
          "exogenous UID rather than defining one universal target oligo set. Variable arms "
          "are therefore shown by molecular role.")

UID = feature("safe_seqs_uid", "umi", "random")
FINAL_LIBRARY, SEQ_PRIMERS = truseq_library([
    seg("14-nt UID", "N" * 14, "umi", placeholder=True, feature=UID),
    seg("target amplicon", "X" * 42, placeholder=True)],
    "Safe-SeqS sequencing library", dual_index=False)
FINAL_CAPTION = ("The first 14 insert bases identify one amplified template strand; the "
                 "remaining insert is the target amplicon used for family consensus.")
SEQUENCING_INTRO = ("Read 1 encounters the 14-nt UID before the target. Standard Illumina "
                    "primer sites are derived from the final modeled library.")


def workflow():
    chemistry = SingleStrandUID(14)
    w = Workflow(MolecularState("Target DNA", tuple(duplex_fragment("target DNA").rows())))
    w.react("Assign a random UID with two target-specific PCR cycles",
            chemistry.assignment().rows(), name="UID-assigned target",
            note="Each original strand receives a primer-borne random identifier before bulk amplification.")
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
                ((assigned, InputRole.TEMPLATE), (universal_f, InputRole.STRUCTURAL),
                 (universal_r, InputRole.STRUCTURAL), (family_pol, InputRole.CATALYST)),
                (library,), geometry=PCRGeometry.unavailable(
                    final, reason="The final TruSeq shell is exact, while target-specific insert coordinates vary."))
    graph.sequence(library, SequencingPlan(library=final, primers=tuple(primers)))
    graph.declare_reviewed_scope([
        "target_dna", "uid_primer", "opposing_primer", "assignment_polymerase",
        "assign_uid", "universal_forward", "universal_reverse", "family_polymerase",
        "amplify_families", "library"], authority)
    return graph.require_publishable()


def protocol_graph():
    return safe_seqs_graph(SOURCE, FINAL_LIBRARY, SEQ_PRIMERS)
