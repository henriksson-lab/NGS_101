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


TITLE = 'RiboMethSeq (Illumina) — 2′-O-methylation by protected cleavage'
NOTES = '01_ribomethseq.html'
SOURCE = 'Marchand et al., Nucleic Acids Research (2016), <a href="https://doi.org/10.1093/nar/gkw547">doi:10.1093/nar/gkw547</a>; original Ion Torrent method (not drawn): Birkedal et al., Angewandte Chemie (2015), <a href="https://doi.org/10.1002/anie.201408362">doi:10.1002/anie.201408362</a>.'
SUMMARY = 'Partial alkaline hydrolysis avoids phosphodiester bonds immediately 3′ of 2′-O-methylated riboses; depletion of mapped fragment ends identifies and quantifies Nm sites.'
CAVEAT = ('This page documents the 2016 Marchand Illumina implementation (NEBNext Small RNA library), which is its defining source. The original 2015 Birkedal RiboMethSeq protocol is a distinct earlier Ion Torrent implementation and is not drawn; nothing below is transferred onto the drawn construct. '
          'The Birkedal Angewandte article itself has still not been obtained; its Wiley Supporting Information (anie_201408362_sm_miscellaneous_information.pdf, obtained 2026-10-17) was read. '
          'According to that SI (Materials and Methods, pp. S3–S5; Figure S4), RNA at 1 µg/µl was mixed with an equal volume of 100 mM Na₂CO₃ pH 9.0 and heated at 90 °C for 1–10 min; 20–40-nt fragments were excised from a 10% denaturing gel and eluted near neutral pH to keep their 2′,3′-cyclic phosphates. '
          'Fragments (5′-OH, 2′,3′-cyclic phosphate) were cloned with a kinase-dead Arabidopsis tRNA ligase (AtRNL S701A, “AtRNLΔPNK”): first a T4 PNK-phosphorylated 3′ RNA adaptor corresponding to the Ion Torrent A sequence with a 4-nt barcode (CR30 5′-CGCACUGAGUCGGAG or CR35 5′-CGAGCUGAGUCGGAG; 30 °C 1 h), then one T4 PNK reaction to phosphorylate the ligation product and remove cyclic phosphates of unligated fragments, gel purification, ligation of an in vitro transcribed, HDV-ribozyme-cleaved 5′ adaptor (P1 IVT, 5′-triphosphate, 43 nt, 2′,3′-cyclic phosphate end), gel purification, FastAP removal of the junction 2′-phosphates, and reverse transcription with a short Ion A primer carrying barcode #1 (CGCA); libraries went to emulsion PCR and Ion Torrent/Proton sequencing (250 or 500 flows). '
          'Marchand Table S4 summarises this as “no end repair” with direct ligation to 5′-OH/>P fragments; the SI shows that no repair precedes the first ligation, but T4 PNK and FastAP treatments occur between and after the ligations. '
          'Source inconsistency within the Birkedal SI: the P1 IVT adaptor is printed as 5′-GGCCACUACGCCUCCGCUUUCCUCUCUAUGGGCAGUCGG<b>T</b>GAU>P (p. S4, adaptor construction) and as …GGCAGUCGG<b>U</b>GAU>P (p. S4, library construction); as an in vitro transcript of the printed T7 template the RNA residue is U, so the T is read as a typographical DNA letter. The RT primer bases are not printed. '
          'GEO series <a href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE54655">GSE54655</a> (PMID 25417815) gives the same outline.')
FINAL_LIBRARY = Construct((role('P5'), role('hydrolysis-end cDNA'), role('sample index', tag='cbc', molecular_feature=feature('sample_index', 'sample_index', 'unknown')), role('P7')), name='RiboMethSeq library')
FINAL_CAPTION = 'Standard small-RNA library retaining the exact alkaline-cleavage endpoints.'
SEQ_PRIMERS = ()
READ_LENGTHS = None
SEQUENCING_INTRO = ''
SEQUENCING_ENDING = 'The Illumina implementation used 2×75 paired-end MiSeq or 50-cycle single-read HiSeq 1000 runs (Marchand Supplementary Table S2); standard kit primers enter from the ligated adapters.'
SEQUENCING_UNAVAILABLE = 'The method paper names the commercial NEBNext Small RNA kit but does not reproduce a distinct sequencing-primer oligo.'

def oligos():
    return ('NEBNext 3′ small-RNA adapter — 5′-/rApp/AGATCGGAAGAGCACACGTCT/3AmMO/-3′',)

def workflow() -> Workflow:
    return chain(molecule("total or enriched RNA", (role("RNA containing 2-prime-O-methylated riboses"),)), (
        ("Partially hydrolyse RNA under mild alkaline conditions",
         molecule("alkaline RNA fragments", (role("random fragments with gaps beside protected Nm sites"),)),
         "A 2′-O-methyl ribose protects its 3′-adjacent phosphodiester bond from cleavage. Marchand used 50 mM bicarbonate pH 9.2 at 95 °C and no gel size selection (30–200-nt fragments)."),
        ("Convert cyclic-phosphate/5′-OH ends to ligatable termini",
         molecule("end-repaired RNA fragments", (role("5-prime phosphate"), role("RNA fragment"), role("3-prime hydroxyl"))),
         "Antarctic phosphatase followed by T4 PNK creates the ends required by the standard small-RNA kit; the repaired RNA is purified on an RNeasy MinElute column."),
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

def _column_cleanup(graph, current, reviewed):
    """Marchand2016 Methods, End repair: RNeasy MinElute purification after T4 PNK."""
    column = graph.input("minelute_condition",
                         "RNeasy MinElute column with 675 µl 96% ethanol binding; 10 µl water elution",
                         MaterialKind.CONDITION)
    kept = graph.product("minelute_product", "column-purified end-repaired RNA fragments",
                         MaterialKind.MOLECULAR, state=current.state)
    waste = graph.product("minelute_waste", "column flow-through and washes", MaterialKind.WASTE)
    graph.react("minelute", "Purify end-repaired RNA on an RNeasy MinElute column",
                StepKind.CLEANUP, ((current, InputRole.TEMPLATE), (column, InputRole.CONDITION)),
                (kept, waste))
    graph.discard(waste, "physically removed in the column flow-through and washes")
    reviewed.extend(("minelute_condition", "minelute", "minelute_product", "minelute_waste"))
    return kept


def protocol_graph() -> ProtocolGraph:
    graph = ProtocolGraph(TITLE, authority=SOURCE)
    start = graph.input("input", TITLE + " source RNA", MaterialKind.BIOLOGICAL,
                        state=molecule("source RNA", (role("input RNA"),)), primary=True)
    reviewed = ["input"]
    current = start
    for identifier, label, kind in RECIPE:
        current = _apply_linear_step(graph, current, identifier, label, kind, reviewed)
        if identifier == 'phosphorylate':
            current = _column_cleanup(graph, current, reviewed)
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
