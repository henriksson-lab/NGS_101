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

SS2_HAIRPIN = Segment("Structure-seq2 hairpin donor", "TGAAGAGCCTAGTCGCTGTTCANNNNNNCTGCCCATAGAG", "r1", placeholder=True)

SS2_RT_SEQUENCE = "CAGACGTGTGCTCTTCCGATC" + "N" * 6

SS2_FORWARD_SEQUENCE = "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCTTGAACAGCGACTAGGCTCTTCA"

SS2_REVERSE_SEQUENCE = "CAAGCAGAAGACGGCATACGAGAT" + "N" * 6 + "GTGACTGGAGTTCAGACGTGTGCTCTTCCGATC"

# Reagent states for the graph inputs.  The RT primer is drawn only in Supplementary
# Figure S1 (step 1A, 3'-NNNNNNCTAGCCTTCTCGTGTGCAGAC-5', 27 nt as stated in Figure S4);
# the hairpin donor and both PCR primers are printed in the Methods.
SS2_RT_PRIMER = Construct((
    Segment("RT-primer fixed handle", SS2_RT_SEQUENCE[:-6], "r2"),
    Segment("random hexamer", SS2_RT_SEQUENCE[-6:],
            note="degenerate priming bases; they anneal at random RNA positions"),
), name="Structure-seq2 RT primer")

SS2_HAIRPIN_DONOR = Construct((
    Segment("hairpin 5' arm and loop", "TGAAGAGCCTAGTCGCTGTTCA", "r1",
            note="5' phosphate (Methods: '5'-p')"),
    Segment("degenerate six-base region", "NNNNNN",
            note="pairs with diverse cDNA 3' ends"),
    Segment("hairpin 3' arm", "CTGCCCATAGAG", "r1",
            note="3' C3 spacer (Methods: '3'-Spacer')"),
), name="Structure-seq2 hairpin donor")

SS2_REVERSE_PRIMER = Construct((
    Segment("P7", il.P7, "p7"),
    Segment("sample barcode", "NNNNNN", "cbc", placeholder=True,
            feature=feature("sample_index", "sample_index", "unknown",
                            note="Methods state a 6-8-nt barcode; the six-base member is shown")),
    Segment("TruSeq Read 2 / Index 1 arm", SS2_REVERSE_SEQUENCE[-33:], "r2"),
), name="Structure-seq2 indexed reverse primer (six-base member)")

if SS2_REVERSE_PRIMER.top() != SS2_REVERSE_SEQUENCE:
    raise AssertionError("indexed reverse primer state differs from its printed sequence")

# Read by tools/selftest_source_oligos.py.
SOURCE_OLIGO_EXCEPTIONS = {
    "rt_primer": ("Printed only as artwork in Supplementary Figure S1 (File006 p. S2, "
                  "step 1A/1B: 3'-NNNNNNCTAGCCTTCTCGTGTGCAGAC-5'); the extracted text "
                  "has no bases, so it was transcribed from the rendered figure."),
}

SS2_FORWARD = Segment("Structure-seq2 TruSeq forward", SS2_FORWARD_SEQUENCE, "p5")

SS2_REVERSE = Segment("Structure-seq2 indexed reverse", SS2_REVERSE_SEQUENCE, "p7", placeholder=True,
                      feature=feature("sample_index", "sample_index", "unknown"))

SS2_CUSTOM = Segment("Structure-seq2 custom sequencing primer", "TCTTCCGATCTTGAACAGCGACTAGGCTCTTCA", "r1")

SS2_PCR_TEMPLATE = Construct((
    Segment("hairpin forward-primer site", SS2_FORWARD_SEQUENCE[-22:], "r1"),
    Segment("stop-defined cDNA", "X" * 42, placeholder=True),
    Segment("RT-primer-derived reverse site", revcomp(SS2_REVERSE_SEQUENCE[-21:]), "r2"),
), name="Structure-seq2 PCR-ready template")

SS2_PCR_GEOMETRY = PCRGeometry.exact(
    SS2_PCR_TEMPLATE,
    PCRPrimer("Structure-seq2 TruSeq forward", SS2_FORWARD_SEQUENCE),
    PCRPrimer("Structure-seq2 indexed reverse", SS2_REVERSE_SEQUENCE, anneal_nt=21),
    statement=("The main Methods print both PCR primers, and Supplementary Figure S1 "
               "draws their hairpin- and RT-primer-derived landing regions."))

SS2_FINAL_LIBRARY = Construct((
    Segment("forward PCR primer", SS2_FORWARD_SEQUENCE, "p5"),
    Segment("stop-defined cDNA", "X" * 42, placeholder=True),
    Segment("Index 1 primer site", revcomp(SS2_REVERSE_SEQUENCE[-33:]), "r2"),
    Segment("sample index complement", "N" * 6, "cbc", placeholder=True,
            feature=feature("sample_index", "sample_index", "unknown")),
    Segment("P7 reverse complement", il.P7_RC, "p7"),
), name="Structure-seq2 six-base-index library")

SS2_SEQ_PRIMERS = (
    sp.custom("Read 1", "Structure-seq2 custom Read 1", SS2_CUSTOM.top,
              "Ritchey et al. Supplementary Figure S1"),
    sp.TRUSEQ["I1"],
)


TITLE = 'Structure-seq2 — low-bias in-vivo RNA structure probing'
NOTES = '01_structure_seq2.html'
SOURCE = 'Ritchey et al., Nucleic Acids Research (2017); <a href="https://doi.org/10.1093/nar/gkx533">doi:10.1093/nar/gkx533</a>.'
SUMMARY = 'Random-primed DMS-stop cDNA receives a degenerate hairpin adapter by T4 DNA ligase; PAGE or biotin–streptavidin purification removes the primer–adapter by-product.'
CAVEAT = 'This graph follows the PAGE purification branch. The biotin/streptavidin variation and the original Structure-seq comparison libraries are separate experimental variants.'
FINAL_LIBRARY = SS2_FINAL_LIBRARY
PCR_GEOMETRY = SS2_PCR_GEOMETRY
FINAL_CAPTION = 'Representative six-base-index PCR product. The custom primer abuts the unique insert, bypassing the 22-nt low-diversity prefix.'
SEQ_PRIMERS = SS2_SEQ_PRIMERS
ACQUIRED_ROLES = ('Read 1', 'Index 1 (i7)')
READ_LENGTHS = {'Read 1': 150, 'Index 1 (i7)': 6}
READ_USES = {'Read 1': 'include', 'Index 1 (i7)': 'demultiplex'}
SEQUENCING_INTRO = 'The custom primer and standard multiplexing-index primer both land on source-exact PCR-derived sites.'
SEQUENCING_ENDING = ''
SEQUENCING_UNAVAILABLE = ''

def oligos():
    return ('RT primer — 5′-CAGACGTGTGCTCTTCCGATCNNNNNN-3′', 'hairpin donor — 5′-/5Phos/TGAAGAGCCTAGTCGCTGTTCANNNNNNCTGCCCATAGAG/3SpC3/-3′', 'TruSeq forward — 5′-AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCTTGAACAGCGACTAGGCTCTTCA-3′', 'indexed reverse — 5′-CAAGCAGAAGACGGCATACGAGAT[6–8-nt barcode]GTGACTGGAGTTCAGACGTGTGCTCTTCCGATC-3′', 'custom sequencing primer — 5′-TCTTCCGATCTTGAACAGCGACTAGGCTCTTCA-3′')

def workflow() -> Workflow:
    return chain(molecule("poly(A)-selected RNA", (role("RNA folded in vivo"),)), (
        ("Probe exposed A and C in vivo with DMS",
         molecule("DMS-probed RNA", (role("RNA with structure-dependent DMS adducts"),)),
         "A matched −DMS sample measures background reverse-transcription stops."),
        ("Random-prime reverse transcription to each DMS stop",
         molecule("stop-defined cDNA", (role("TruSeq handle plus random hexamer"), role("cDNA ending before modified base"))),
         "The optional branch incorporates biotin-dCTP for streptavidin purification instead of two PAGE steps."),
        ("Ligate the degenerate hairpin donor with T4 DNA ligase",
         molecule("hairpin-adapted cDNA", (SS2_HAIRPIN, role("stop-defined cDNA")),
                  junctions=(("Structure-seq2 hairpin donor", "stop-defined cDNA", "T4 DNA ligase"),)),
         "Six random bases let the donor base-pair at diverse cDNA ends and reduce ligation bias."),
        ("PCR with TruSeq forward and indexed reverse primers",
         molecule("Structure-seq2 library", (SS2_FORWARD, role("stop-defined cDNA"), SS2_REVERSE)),
         "A final 220–600-nt PAGE selection removes primer dimers and residual by-product."),
    ))



REAGENTS = {
    "rt": dict(reagent_label="RT primer: TruSeq handle plus random hexamer",
               reagent_state=SS2_RT_PRIMER,
               reagent_source="Ritchey et al. (2017) Supplementary Figure S1, step 1A"),
    "ligate": dict(reagent_label="Structure-seq2 hairpin donor (5'-phosphate, 3'-C3 spacer)",
                   reagent_state=SS2_HAIRPIN_DONOR,
                   reagent_source="Ritchey et al. (2017), Methods, hairpin ligation (PMC5737731)"),
}

RECIPE = (('probe', 'Probe exposed A and C in vivo with DMS', StepKind.CHEMICAL_CONVERSION), ('rt', 'Random-prime DMS-stop cDNA with SuperScript III (55°C, 50 min)', StepKind.REVERSE_TRANSCRIPTION), ('digest_rna', 'Digest residual RNA with RNase H (37°C, 20 min)', StepKind.DIGESTION), ('ligate', 'Ligate the degenerate hairpin donor', StepKind.LIGATION))

def protocol_graph() -> ProtocolGraph:
    graph = ProtocolGraph(TITLE, authority=SOURCE)
    start = graph.input("input", TITLE + " source RNA", MaterialKind.BIOLOGICAL,
                        state=molecule("source RNA", (role("input RNA"),)), primary=True)
    reviewed = ["input"]
    current = start
    for identifier, label, kind in RECIPE:
        current = _apply_linear_step(graph, current, identifier, label, kind, reviewed,
                                     **REAGENTS.get(identifier, {}))
    _finish_pcr(
        graph, "library", current, FINAL_LIBRARY,
        SEQUENCING_ENDING or FINAL_CAPTION, reviewed, label=TITLE,
        step_label="Amplify and index with Illumina TruSeq forward and barcoded reverse primers (Q5)",
        forward_label="Illumina TruSeq forward primer",
        reverse_label="Illumina TruSeq reverse primer with 6–8-nt sample barcode",
        polymerase_label="Q5 High-Fidelity DNA polymerase (NEB)",
        forward_state=SS2_FORWARD,
        reverse_state=SS2_REVERSE_PRIMER,
        primer_source=("Ritchey et al. (2017), Methods, 'Library amplification by PCR' "
                       "(PMC5737731)."),
        geometry=PCR_GEOMETRY,
        primers=SEQ_PRIMERS,
        acquired_roles=ACQUIRED_ROLES,
        read_lengths=READ_LENGTHS,
        read_uses=READ_USES)
    graph.declare_reviewed_scope(reviewed, SOURCE)
    return graph.require_publishable()
