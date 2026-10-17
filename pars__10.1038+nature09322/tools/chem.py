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


TITLE = 'PARS — parallel analysis of RNA structure'
NOTES = '01_pars.html'
SOURCE = 'Kertesz et al., Nature (2010), <a href="https://doi.org/10.1038/nature09322">doi:10.1038/nature09322</a>; detailed protocol: Wan et al., Nature Protocols (2013), <a href="https://doi.org/10.1038/nprot.2013.045">doi:10.1038/nprot.2013.045</a>.'
SUMMARY = 'Matched RNase V1 and S1 digests expose double- and single-stranded positions; selective capture of nuclease-created 5′ phosphates turns branch-specific read starts into a per-base structure score.'
CAVEAT = 'The defining experiment used SOLiD. Detailed adapter/recovery steps still require verification against its supplementary methods; the later 2013 protocol is not a substitute for that check.'
FINAL_LIBRARY = Construct((role('SOLiD P1 adapter'), role('V1- or S1-cleavage-defined cDNA'), role('SOLiD P2 adapter')), name='PARS library')
FINAL_CAPTION = 'One library is made for each nuclease branch; the defining experiment used the SOLiD platform.'
SEQ_PRIMERS = ()
READ_LENGTHS = None
SEQUENCING_INTRO = ''
SEQUENCING_ENDING = 'The defining libraries used Applied Biosystems SOLiD; reads were truncated to 35 bases for mapping (not an asserted acquisition length); the vendor guide publishes library PCR primers but not the complete adaptor mix or sequencing primer.'
SEQUENCING_UNAVAILABLE = 'Adaptor Mix A remains a degenerate kit oligo set without complete published bases; the exact sequencing primer is also unresolved.'

def oligos():
    return ('SOLiD Small RNA kit common P1 / 5′ PCR primer — 5′-CCACTACGCCTCCGCTTTCCTCTCTATGGGCAGTCGGTGAT-3′', "SOLiD Small RNA kit indexed 3′ PCR primer 1 (kit example; experiment's primer set unstated) — 5′-CTGCCCCGGGTTCCTCATTCTCTAAGCCCCTGCTGTACGGCCAAGGCG-3′")

def workflow() -> Workflow:
    return chain(molecule("poly(A)-selected RNA", (role("purified transcriptome"),)), (
        ("Refold RNA in vitro and split it into matched aliquots",
         molecule("refolded transcriptome", (role("RNA secondary structures"),)),
         "The two branches start from the same refolded RNA pool."),
        ("Digest double-stranded or single-stranded regions",
         molecule("paired V1/S1 digests", (role("RNase V1: double-strand cleavage"), role("RNase S1: single-strand cleavage"))),
         "RNase V1 and S1 generate separate structure-selective 5′-phosphate endpoints."),
        ("Randomly fragment, then ligate the 5′ adapter selectively",
         molecule("cleavage-selected RNA", (role("5-prime adapter"), role("RNA downstream of nuclease cut"))),
         "Random fragmentation yields 5′-OH ends; only the nuclease-created 5′-P ends accept the first adapter."),
        ("Dephosphorylate 3′ ends and ligate the second adapter",
         molecule("dual-adapter PARS RNA", (role("5-prime adapter"), role("structure-defined RNA fragment"), role("3-prime adapter"))),
         "Alkaline phosphatase makes the second ligation compatible with nuclease and fragmentation end chemistries."),
        ("Reverse-transcribe, PCR-amplify and sequence both branches",
         molecule("PARS library pair", (role("SOLiD P1 end"), role("V1 or S1 cleavage-defined cDNA"), role("SOLiD P2 end"))),
         "The log2 normalized V1/S1 read-start ratio is the PARS score for the upstream nucleotide."),
    ))


def protocol_graph():
    source = SOURCE
    graph = ProtocolGraph(TITLE, authority=source)
    rna = graph.input("input", "refolded poly(A)-selected transcriptome",
                      MaterialKind.BIOLOGICAL,
                      state=molecule("refolded transcriptome", (role("RNA secondary structures"),)),
                      primary=True)
    v1 = graph.product("v1_branch", "RNase-V1 double-strand branch",
                       MaterialKind.MOLECULAR, state=rna.state)
    s1 = graph.product("s1_branch", "RNase-S1 single-strand branch",
                       MaterialKind.MOLECULAR, state=rna.state)
    graph.react("split", "Split refolded RNA into matched V1 and S1 aliquots",
                StepKind.SPLIT, ((rna, InputRole.TEMPLATE),), (v1, s1),
                changes_molecule=False)
    reviewed = ["input", "split"]
    for prefix, template, label in (("v1", v1, "RNase V1"), ("s1", s1, "RNase S1")):
        digested = _apply_linear_step(graph, template, prefix + "_digest",
                                      "Digest with " + label, StepKind.DIGESTION, reviewed)
        first = _apply_linear_step(graph, digested, prefix + "_ligate_5",
                                   "Select nuclease-created 5-prime phosphates by adapter ligation",
                                   StepKind.LIGATION, reviewed)
        dephos = _apply_linear_step(graph, first, prefix + "_dephosphorylate",
                                    "Dephosphorylate remaining RNA ends",
                                    StepKind.DEPHOSPHORYLATION, reviewed)
        second = _apply_linear_step(graph, dephos, prefix + "_ligate_3",
                                    "Ligate the second adapter", StepKind.LIGATION, reviewed)
        cdna = _apply_linear_step(graph, second, prefix + "_rt",
                                  "Reverse-transcribe the " + label + " branch",
                                  StepKind.REVERSE_TRANSCRIPTION, reviewed)
        _finish_pcr(graph, prefix, cdna, FINAL_LIBRARY,
                    label + " branch: " + FINAL_CAPTION, reviewed,
                    label="PARS " + label)
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()
