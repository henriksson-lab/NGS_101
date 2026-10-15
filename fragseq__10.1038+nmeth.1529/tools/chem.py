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


TITLE = 'FragSeq — nuclease accessibility profiling of RNA'
NOTES = '01_fragseq.html'
SOURCE = 'Underwood et al., Nature Methods (2010), <a href="https://doi.org/10.1038/nmeth.1529">doi:10.1038/nmeth.1529</a>.'
SUMMARY = 'Partial nuclease P1 digestion converts exposed single-stranded RNA into mapped fragment termini, directly profiling RNA accessibility across many transcripts.'
CAVEAT = 'FragSeq probes purified, refolded RNA. It reports nuclease accessibility, not native intracellular accessibility or abundance alone.'
FINAL_LIBRARY = Construct((role('SOLiD P1 adapter'), role('P1-cleavage-defined cDNA'), role('SOLiD P2 adapter')), name='FragSeq library')
FINAL_CAPTION = 'The informative coordinate is a nuclease-P1-defined RNA fragment end.'
SEQ_PRIMERS = ()
READ_LENGTHS = None
SEQUENCING_INTRO = ''
SEQUENCING_ENDING = 'The defining experiment converted treated and control RNA fragments with the SOLiD small-RNA workflow and sequenced their ends.'
SEQUENCING_UNAVAILABLE = 'Exact proprietary SOLiD adapter and sequencing-primer sequences are not represented.'

def oligos():
    return ()

def workflow() -> Workflow:
    return chain(molecule("purified nuclear RNA", (role("native folded RNA"),)), (
        ("Refold RNA and split treated and undigested-control aliquots",
         molecule("refolded RNA", (role("RNA secondary structure"),)),
         "The control measures fragmentation and library-preparation background."),
        ("Partially digest exposed single-stranded RNA with nuclease P1",
         molecule("P1 accessibility fragments", (role("5-prime phosphate at accessible site"),
                  role("downstream RNA fragment"))),
         "Nuclease P1 preferentially cleaves single-stranded nucleic acid; fragment starts encode accessibility."),
        ("Size-select and ligate small-RNA sequencing adapters",
         molecule("adapter-ligated FragSeq RNA", (role("5-prime sequencing adapter"),
                  role("P1-cleavage-defined RNA"), role("3-prime sequencing adapter"))),
         "Treated and control branches are converted independently so their endpoint enrichments can be compared."),
        ("Reverse-transcribe, amplify and sequence fragment ends",
         molecule("FragSeq library", (role("SOLiD P1 end"),
                  role("P1-cleavage-defined cDNA"), role("SOLiD P2 end"))),
         "Mapped fragment termini provide nucleotide-resolution nuclease-accessibility scores."),
    ))


def protocol_graph():
    source = SOURCE
    graph = ProtocolGraph(TITLE, authority=source)
    rna = graph.input("input", "purified refolded nuclear RNA", MaterialKind.BIOLOGICAL,
                      state=molecule("refolded nuclear RNA", (role("native folded RNA"),)),
                      primary=True)
    treated = graph.product("treated_branch", "nuclease-P1-treated aliquot",
                            MaterialKind.MOLECULAR, state=rna.state)
    control = graph.product("control_branch", "matched undigested control aliquot",
                            MaterialKind.MOLECULAR, state=rna.state)
    graph.react("split", "Split treated and undigested-control aliquots",
                StepKind.SPLIT, ((rna, InputRole.TEMPLATE),), (treated, control),
                changes_molecule=False)
    reviewed = ["input", "split"]
    treated = _apply_linear_step(graph, treated, "p1_digest",
                                 "Partially digest accessible RNA with nuclease P1",
                                 StepKind.DIGESTION, reviewed)
    for prefix, template, label in (("treated", treated, "P1-treated"),
                                    ("control", control, "undigested control")):
        ligated = _apply_linear_step(graph, template, prefix + "_ligate",
                                     "Ligate small-RNA sequencing adapters to the " + label + " branch",
                                     StepKind.LIGATION, reviewed)
        cdna = _apply_linear_step(graph, ligated, prefix + "_rt",
                                  "Reverse-transcribe the " + label + " branch",
                                  StepKind.REVERSE_TRANSCRIPTION, reviewed)
        _finish_pcr(graph, prefix, cdna, FINAL_LIBRARY,
                    label + ": " + FINAL_CAPTION, reviewed,
                    label="FragSeq " + label)
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()
