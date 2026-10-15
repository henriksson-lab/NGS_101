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


TITLE = 'GLORI — absolute m6A quantification by A-to-I conversion'
NOTES = '01_glori.html'
SOURCE = 'Liu et al., Nature Biotechnology (2023), <a href="https://doi.org/10.1038/s41587-022-01487-9">doi:10.1038/s41587-022-01487-9</a>; detailed protocol: Shen et al., Nature Protocols (2024), <a href="https://doi.org/10.1038/s41596-023-00937-1">doi:10.1038/s41596-023-00937-1</a>.'
SUMMARY = 'Glyoxal-assisted nitrite treatment converts unmethylated adenosine to inosine while m6A remains adenosine; A retention therefore reports absolute methylation stoichiometry.'
CAVEAT = "This page follows the modified-eCLIP library branch. Exact ligation adapters and RT primer come from the authors' 2024 detailed protocol; equivalence to every 2023 experimental lot is not established. The authors also validated NEBNext Small RNA and SMARTer stranded kit alternatives."
FINAL_LIBRARY = Construct((role('P5'), role('10-nt UMI', tag='umi', length_bp=10, molecular_feature=feature('umi', 'umi', 'random')), role('A-to-I-converted cDNA'), role('sample index', tag='cbc', molecular_feature=feature('sample_index', 'sample_index', 'unknown')), role('P7')), name='GLORI modified-eCLIP library')
FINAL_CAPTION = 'Modified-eCLIP branch. The 10-nt UMI is inline; converted A is read as G while m6A is retained as A.'
SEQ_PRIMERS = ()
READ_LENGTHS = None
SEQUENCING_INTRO = ''
SEQUENCING_ENDING = 'Defining modified-eCLIP libraries were sequenced paired-end on an Illumina HiSeq X Ten; reads begin with the 10-nt UMI before converted insert sequence.'
SEQUENCING_UNAVAILABLE = "The 2024 protocol's Supplementary Table 1 gives both ligation adapters and the RT primer. Complete library amplification and sequencing-primer oligos remain unreported for this branch; its final primer landing sites are not reconstructed."

def oligos():
    return ('2024 detailed protocol, 3′ RNA adaptor — 5′-rAppAGATCGGAAGAGCGTCGTG-3′', '2024 detailed protocol, RT primer — 5′-ACACGACGCTCTTCCGATCT-3′', '2024 detailed protocol, 3′ DNA adaptor — 5′-Phos-NNNNNNNNNNAGATCGGAAGAGCACACGTCTG-3′')

def workflow() -> Workflow:
    readout = BaseConversionReadout("glyoxal plus nitrite", "adenosine", "inosine", "G", "m6A")
    return chain(molecule("fragmented poly(A) RNA", (role("RNA containing A and m6A"),)), (
        ("Protect guanosines with glyoxal",
         molecule("glyoxal-protected RNA", (role("RNA"), role("protected G"), role("A and m6A exposed"))),
         "G protection suppresses unwanted guanosine deamination during the nitrite reaction."),
        ("Deaminate unmethylated A with nitrite",
         molecule("GLORI-converted RNA", (role(f"{readout.converted_base} converted to {readout.product}"), role(f"{readout.protected_base} retained as A"))),
         "The defining reaction uses 750 mM nitrite at 16 °C for 8 h followed by 4 °C overnight."),
        ("Heat-deprotect guanosines and repair RNA ends",
         molecule("deprotected converted RNA", (role("A-to-I converted RNA with intact m6A"), role("ligatable termini"))),
         "Deprotection at 95 °C restores G; T4 PNK prepares fragmented RNA for adapter ligation."),
        ("Ligate the modified-eCLIP 3-prime RNA adapter",
         molecule("GLORI adapter-ligated RNA", (role("converted RNA"),
                  Segment("published 3-prime RNA adapter", "AGATCGGAAGAGCGTCGTG", "r2"))),
         "The 2024 detailed protocol prints the pre-adenylated RNA adapter sequence in Supplementary Table 1."),
        ("Reverse-transcribe with the published RT primer",
         molecule("GLORI cDNA", (role("converted cDNA"),
                  Segment("published RT primer", "ACACGACGCTCTTCCGATCT", "r2"))),
         "The 2024 detailed protocol prints this RT primer in Supplementary Table 1."),
        ("Ligate the ten-base-UMI DNA adapter to cDNA",
         molecule("UMI-tagged GLORI cDNA", (role("converted cDNA"), role("10-nt UMI"),
                  Segment("published 3-prime cDNA adapter arm", "AGATCGGAAGAGCACACGTCTG", "r2"))),
         "The 2024 detailed protocol prints a 5-prime-phosphorylated DNA adapter with ten degenerate UMI bases."),
        ("Index-amplify and paired-end sequence",
         molecule("GLORI library", (role("P5"), role("10-nt UMI"), role("converted cDNA"), role("sample index"), role("P7"))),
         f"Unmodified A is read as {readout.read_as}; m6A remains A, enabling site-level stoichiometry."),
    ))


PCR_GEOMETRY = None
ACQUIRED_ROLES = ()
READ_USES = None

RECIPE = (('protect_g', 'Protect guanosines with glyoxal', StepKind.CHEMICAL_CONVERSION), ('deaminate', 'Deaminate unmethylated adenosine with nitrite', StepKind.CHEMICAL_CONVERSION), ('deprotect', 'Heat-deprotect guanosines', StepKind.THERMAL_TREATMENT), ('repair', 'Repair converted RNA ends', StepKind.PHOSPHORYLATION), ('ligate', 'Ligate the modified-eCLIP 3-prime RNA adapter', StepKind.LIGATION), ('rt', 'Reverse-transcribe converted RNA', StepKind.REVERSE_TRANSCRIPTION), ('ligate_cdna', 'Ligate the UMI-bearing DNA adapter to cDNA', StepKind.LIGATION))

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
