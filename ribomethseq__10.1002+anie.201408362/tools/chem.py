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


TITLE = 'RiboMethSeq — 2′-O-methylation by protected cleavage'
NOTES = '01_ribomethseq.html'
SOURCE = 'Birkedal et al., Angewandte Chemie (2015), <a href="https://doi.org/10.1002/anie.201408362">doi:10.1002/anie.201408362</a>; Illumina implementation: Marchand et al., Nucleic Acids Research (2016), <a href="https://doi.org/10.1093/nar/gkw547">doi:10.1093/nar/gkw547</a>.'
SUMMARY = 'Partial alkaline hydrolysis avoids phosphodiester bonds immediately 3′ of 2′-O-methylated riboses; depletion of mapped fragment ends identifies and quantifies Nm sites.'
CAVEAT = 'The displayed NEBNext/Illumina workflow is the 2016 Marchand implementation. The catalogue’s defining 2015 Birkedal article has not been checked in this pass; the two implementations must not be assumed to use identical library chemistry.'
FINAL_LIBRARY = Construct((role('P5'), role('hydrolysis-end cDNA'), role('sample index', tag='cbc', molecular_feature=feature('sample_index', 'sample_index', 'unknown')), role('P7')), name='RiboMethSeq library')
FINAL_CAPTION = 'Standard small-RNA library retaining the exact alkaline-cleavage endpoints.'
SEQ_PRIMERS = ()
READ_LENGTHS = None
SEQUENCING_INTRO = ''
SEQUENCING_ENDING = 'The Illumina implementation used paired-end MiSeq or single-read 50-cycle HiSeq 1000 sequencing; standard kit primers enter from the ligated adapters.'
SEQUENCING_UNAVAILABLE = 'The method paper names the commercial NEBNext Small RNA kit but does not reproduce a distinct sequencing-primer oligo.'

def oligos():
    return ('NEBNext 3′ small-RNA adapter — 5′-/rApp/AGATCGGAAGAGCACACGTCT/3AmMO/-3′',)

def workflow() -> Workflow:
    return chain(molecule("total or enriched RNA", (role("RNA containing 2-prime-O-methylated riboses"),)), (
        ("Partially hydrolyse RNA under mild alkaline conditions",
         molecule("alkaline RNA fragments", (role("random fragments with gaps beside protected Nm sites"),)),
         "A 2′-O-methyl ribose protects its 3′-adjacent phosphodiester bond from cleavage."),
        ("Convert cyclic-phosphate/5′-OH ends to ligatable termini",
         molecule("end-repaired RNA fragments", (role("5-prime phosphate"), role("RNA fragment"), role("3-prime hydroxyl"))),
         "Antarctic phosphatase followed by T4 PNK creates the ends required by the standard small-RNA kit."),
        ("Ligate 3′ and 5′ small-RNA adapters",
         molecule("dual-adapter RNA fragments", (role("NEBNext 5-prime small-RNA adapter"), role("hydrolysis fragment"), role("NEBNext 3-prime small-RNA adapter"))),
         "Direct ligation preserves the exact hydrolysis endpoints used for methylation scoring."),
        ("Reverse-transcribe, index by PCR and sequence fragment ends",
         molecule("RiboMethSeq library", (role("P5"), role("hydrolysis-end cDNA"), role("sample index"), role("P7"))),
         "Paired-end MiSeq can recover both ends; SR50 HiSeq analysis can use the less biased 5′-end profile."),
    ))


PCR_GEOMETRY = None
ACQUIRED_ROLES = ()
READ_USES = None

RECIPE = (('hydrolyse', 'Partially hydrolyse RNA under mild alkaline conditions', StepKind.CHEMICAL_CONVERSION), ('dephosphorylate', 'Dephosphorylate RNA 3-prime ends with Antarctic phosphatase', StepKind.DEPHOSPHORYLATION), ('phosphorylate', 'Phosphorylate RNA 5-prime ends with T4 PNK and ATP', StepKind.PHOSPHORYLATION), ('ligate_3', 'Ligate the 3-prime small-RNA adapter', StepKind.LIGATION), ('ligate_5', 'Ligate the 5-prime small-RNA adapter', StepKind.LIGATION), ('rt', 'Reverse-transcribe hydrolysis fragments', StepKind.REVERSE_TRANSCRIPTION))

def protocol_graph() -> ProtocolGraph:
    graph = ProtocolGraph(TITLE, authority=SOURCE)
    start = graph.input("input", TITLE + " source RNA", MaterialKind.BIOLOGICAL,
                        state=molecule("source RNA", (role("input RNA"),)), primary=True)
    reviewed = ["input"]
    current = start
    for identifier, label, kind in RECIPE:
        current = _apply_linear_step(graph, current, identifier, label, kind, reviewed)
    _finish_pcr(
        graph, "library", current, FINAL_LIBRARY,
        SEQUENCING_ENDING or FINAL_CAPTION, reviewed, label=TITLE,
        geometry=PCR_GEOMETRY,
        primers=SEQ_PRIMERS,
        acquired_roles=ACQUIRED_ROLES,
        read_lengths=READ_LENGTHS,
        read_uses=READ_USES)
    graph.declare_reviewed_scope(reviewed, SOURCE)
    return graph.require_publishable()
