"""Protocol-owned source, oligos, molecular workflow and reaction graph."""
from adapter_ligation import AdapterLigationGeometry
from chemdraw import Construct, Segment, Workflow, feature, revcomp
from circular import CircularizationGeometry
import illumina as il
from pcr import PCRGeometry, PCRPrimer
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind
from protocol_model import chain, molecule, role
import seqprimers as sp
from rna_end_readout import (StopReadout, BaseConversionReadout, SelectiveCleavageReadout,
                             _finish_pcr, _apply_linear_step)


TITLE = 'PIP-seq — protein interaction profile sequencing'
NOTES = '01_protein_pip_seq.html'
SOURCE = 'Silverman et al., Genome Biology (2014), <a href="https://doi.org/10.1186/gb-2014-15-1-r3">doi:10.1186/gb-2014-15-1-r3</a>.'
SUMMARY = 'Four matched libraries combine single- and double-strand-specific RNases with protein-present and protein-removed branches, revealing both RNA accessibility and protein-protected sites.'
CAVEAT = 'This is protein-interaction-profile PIP-seq, not the unrelated particle-templated single-cell RNA-seq protocol with the same acronym.'
FINAL_LIBRARY = Construct((role('P5'), role('RNase-protected strand-specific cDNA'), role('sample index', tag='cbc', molecular_feature=feature('sample_index', 'sample_index', 'unknown')), role('P7')), name='protein-interaction-profile PIP-seq library')
FINAL_CAPTION = 'One strand-specific sequencing library is made for each protein-state × nuclease-specificity branch.'
SEQ_PRIMERS = ()
READ_LENGTHS = None
SEQUENCING_INTRO = ''
SEQUENCING_ENDING = 'The defining study sequenced strand-specific libraries from all four branches; branch identity must be retained through sample indexing.'
SEQUENCING_UNAVAILABLE = 'Exact commercial library-kit primer sequences are not represented in the defining schematic.'

def oligos():
    return ()

def workflow() -> Workflow:
    return chain(molecule("cross-linked cellular ribonucleoproteins", (role("RNA"),
                  role("bound proteins"))), (
        ("Split material into footprint and structure-control branches",
         molecule("four PIP-seq branches", (role("protein-bound ssRNase branch"),
                  role("protein-bound dsRNase branch"), role("protein-removed ssRNase branch"),
                  role("protein-removed dsRNase branch"))),
         "Proteinase K precedes nuclease only in controls; footprint branches retain protein during digestion."),
        ("Digest exposed single- or double-stranded RNA",
         molecule("RNase-protected fragments", (role("ssRNase-resistant RNA"),
                  role("dsRNase-resistant RNA"))),
         "The paired nuclease specificities distinguish RNA structure while the branch order reveals protein protection."),
        ("Reverse cross-links, purify RNA and prepare strand-specific libraries",
         molecule("strand-specific PIP-seq libraries", (role("5-prime sequencing adapter"),
                  role("RNase-protected RNA-derived cDNA"), role("3-prime sequencing adapter"))),
         "Rehybridization and duplex-specific nuclease normalization reduce abundant-transcript domination."),
        ("Sequence and compare footprint with protein-depleted controls",
         molecule("PIP-seq library quartet", (role("P5"),
                  role("branch-specific protected-fragment cDNA"), role("sample index"), role("P7"))),
         "Enrichment in footprint libraries marks protein-protected sites; control branches report RNA accessibility."),
    ))


def protocol_graph():
    source = SOURCE
    graph = ProtocolGraph(TITLE, authority=source)
    rnps = graph.input("input", "crosslinked cellular ribonucleoproteins",
                       MaterialKind.BIOLOGICAL,
                       state=molecule("crosslinked RNPs", (role("RNA"), role("bound proteins"))),
                       primary=True)
    branches = []
    for identifier, label in (
            ("bound_ss", "protein-bound single-strand-RNase branch"),
            ("bound_ds", "protein-bound double-strand-RNase branch"),
            ("removed_ss", "protein-removed single-strand-RNase control"),
            ("removed_ds", "protein-removed double-strand-RNase control")):
        branches.append(graph.product(identifier, label, MaterialKind.MOLECULAR,
                                      state=rnps.state))
    graph.react("split", "Split into the four protein-state by nuclease branches",
                StepKind.SPLIT, ((rnps, InputRole.TEMPLATE),), tuple(branches),
                changes_molecule=False)
    reviewed = ["input", "split"]
    for branch in branches:
        prefix = branch.id
        current = branch
        if prefix.startswith("removed"):
            current = _apply_linear_step(graph, current, prefix + "_deproteinize",
                                         "Remove protein before nuclease digestion",
                                         StepKind.DEPROTEINIZATION, reviewed)
        current = _apply_linear_step(graph, current, prefix + "_digest",
                                     "Digest with the branch-specific RNA nuclease",
                                     StepKind.DIGESTION, reviewed)
        if prefix.startswith("bound"):
            current = _apply_linear_step(graph, current, prefix + "_deproteinize",
                                         "Reverse crosslinks after footprint digestion",
                                         StepKind.DEPROTEINIZATION, reviewed)
        current = _apply_linear_step(graph, current, prefix + "_rt",
                                     "Reverse-transcribe protected RNA fragments",
                                     StepKind.REVERSE_TRANSCRIPTION, reviewed)
        _finish_pcr(graph, prefix, current, FINAL_LIBRARY,
                    branch.label + ": " + FINAL_CAPTION, reviewed,
                    label="PIP-seq " + branch.label)
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()
