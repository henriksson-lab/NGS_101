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


TITLE = 'AlkAniline-Seq — positive selection of modification-derived RNA cleavage'
NOTES = '01_alkaniline_seq.html'
SOURCE = 'Marchand et al., Angewandte Chemie (2018), <a href="https://doi.org/10.1002/anie.201810946">doi:10.1002/anie.201810946</a>.'
SUMMARY = 'Alkali-labile m7G, m3C and related residues form abasic sites; after background dephosphorylation, aniline cleavage creates the only 5′ phosphates accepted by the first sequencing adapter.'
CAVEAT = 'AlkAniline-Seq also responds to dihydrouridine and some other abasic-site-forming lesions. A signal identifies the chemistry, not uniquely one modification class.'
FINAL_LIBRARY = Construct((role('P5'), role('cDNA beginning N+1 after aniline-sensitive modification'), role('sample index', tag='cbc', molecular_feature=feature('sample_index', 'sample_index', 'unknown')), role('P7')), name='AlkAniline-Seq library')
FINAL_CAPTION = 'The informative coordinate is the cDNA/read start one nucleotide downstream of the modified base.'
SEQ_PRIMERS = ()
READ_LENGTHS = None
SEQUENCING_INTRO = ''
SEQUENCING_ENDING = 'The defining workflow used single-read 50-nt Illumina sequencing; a standard small-RNA read primer enters from the selectively ligated 5′ adapter.'
SEQUENCING_UNAVAILABLE = 'NEB confirms standard Illumina sequencing primers for this kit family; the role-level final molecule has not yet reconstructed their exact landing site.'

def oligos():
    return ('NEBNext E7300/E7330 3′ SR adapter — 5′-rAppAGATCGGAAGAGCACACGTCT-NH2-3′', 'NEBNext E7300/E7330 5′ SR adapter — 5′-rGrUrUrCrArGrArGrUrUrCrUrArCrArGrUrCrCrGrArCrGrArUrC-3′', 'NEBNext E7300/E7330 SR RT primer — 5′-AGACGTGTGCTCTTCCGATCT-3′', 'NEBNext E7300/E7330 SR PCR primer — 5′-AATGATACGGCGACCACCGAGATCTACACGTTCAGAGTTCTACAGTCCG*A-3′', "NEBNext E7300/E7330 index 1 PCR primer (kit example; experiment's index choice unstated) — 5′-CAAGCAGAAGACGGCATACGAGATCGTGATGTGACTGGAGTTCAGACGTGTGCTCTTCCGATC*T-3′")

def workflow() -> Workflow:
    readout = SelectiveCleavageReadout(("m7G", "m3C", "D"), "aniline", "5-prime phosphate", 1)
    targets = ", ".join(readout.target_modifications)
    return chain(molecule("total or enriched RNA", (role(f"RNA containing {targets}"),)), (
        ("Partially hydrolyse RNA in bicarbonate at 96 °C",
         molecule("alkali-treated RNA", (role("size-compatible RNA fragments"), role("abasic sites at alkali-labile modifications"))),
         "Ring opening and depurination create abasic intermediates at susceptible modified bases."),
        ("Remove every pre-existing terminal phosphate",
         molecule("dephosphorylated RNA", (role("5-prime hydroxyl"), role("RNA with abasic sites"), role("3-prime hydroxyl"))),
         "Antarctic phosphatase erases background 5′-P and fragmentation-derived 3′ phosphates."),
        ("Cleave abasic sites with 1 M aniline",
         molecule("aniline-cleaved RNA", (role("upstream cleavage fragment"), role("downstream fragment with new 5-prime phosphate"))),
         "Fifteen minutes at 60 °C and pH 4.5 creates the sole ligation-competent 5′ ends."),
        ("Ligate NEBNext small-RNA adapters",
         molecule("selected AlkAniline RNA", (
             Segment("NEBNext 5-prime SR adapter", "GTTCAGAGTTCTACAGTCCGACGATC", "r1"),
             role("RNA beginning N+1 after modification"),
             Segment("NEBNext 3-prime SR adapter", "AGATCGGAAGAGCACACGTCT", "r2"))),
         "The 5′ adapter positively selects the aniline-generated phosphate; the read start is N+1. The current NEB E7300/E7330 manual identifies both adapter sequences."),
        ("Reverse-transcribe, barcode by PCR and single-end sequence",
         molecule("AlkAniline-Seq library", (role("P5"), role("cleavage-defined cDNA"), role("sample index"), role("P7"))),
         "Only read 5′ ends are required for normalized-cleavage scoring."),
    ))


PCR_GEOMETRY = None
ACQUIRED_ROLES = ()
READ_USES = None

RECIPE = (('hydrolyse', 'Partially hydrolyse RNA under alkaline heat', StepKind.CHEMICAL_CONVERSION), ('dephosphorylate', 'Remove pre-existing terminal phosphates', StepKind.DEPHOSPHORYLATION), ('aniline', 'Cleave abasic sites with aniline', StepKind.CHEMICAL_CONVERSION), ('ligate', 'Ligate small-RNA adapters to newly phosphorylated fragments', StepKind.LIGATION), ('rt', 'Reverse-transcribe cleavage-selected RNA', StepKind.REVERSE_TRANSCRIPTION))

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
