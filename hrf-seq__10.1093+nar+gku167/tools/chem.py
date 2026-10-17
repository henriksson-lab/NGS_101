from dataclasses import replace
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

HRF_RT_RANDOM = "AGACGTGTGCTCTTCCGATCT" + "N" * 8 + "S"

HRF_RT_CASSETTE = "AGACGTGTGCTCTTCCGATCTGAACCGGACCGAAGCCCG"

HRF_LIGATION_ADAPTER = "N" * 7 + "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"

HRF_PCR_FORWARD = "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCT"

HRF_INDEXES = ("GGAACT", "GGACGG", "CGTACG", "GCTACC")

def _hrf_reverse(index: str) -> str:
    return il.P7 + index + il.TRUSEQ_READ2


HRF_PCR_TEMPLATE = Construct((
    Segment("ligation-adapter forward site", HRF_PCR_FORWARD[-19:], "r1"),
    Segment("accessibility-defined cDNA", "X" * 42, placeholder=True),
    Segment("random-primer coordinate", "N" * 8 + "S", placeholder=True),
    Segment("RT-handle reverse site", revcomp(_hrf_reverse(HRF_INDEXES[0])[-21:]), "r2"),
), name="HRF-Seq PCR-ready template")

HRF_PCR_GEOMETRY = PCRGeometry.exact(
    HRF_PCR_TEMPLATE,
    PCRPrimer("PCR_forward", HRF_PCR_FORWARD, anneal_nt=19),
    PCRPrimer("PCR_REVERSE_INDEX.14_AGTTCC", _hrf_reverse(HRF_INDEXES[0]),
              anneal_nt=21),
    statement=("Supplementary Table 1 prints the random RT primer, ligation adapter, "
               "forward primer and all four indexed reverse primers."))

HRF_FINAL_LIBRARY = Construct((
    Segment("P5 / partial Read 1 PCR primer", HRF_PCR_FORWARD, "p5"),
    Segment("remaining Read 1 site", il.TRUSEQ_READ1[-10:], "r1"),
    Segment("7-nt molecular barcode", "N" * 7, "umi", placeholder=True,
            feature=feature("hrf_molecular_barcode", "umi", "random")),
    Segment("accessibility-defined cDNA", "X" * 42, placeholder=True),
    Segment("random-primer coordinate", "N" * 8 + "S", placeholder=True),
    Segment("Index 1 primer site", revcomp(il.TRUSEQ_READ2), "r2"),
    Segment("sample index", revcomp(HRF_INDEXES[0]), "cbc",
            feature=feature("hrf_sample_index", "sample_index", "whitelist",
                            whitelist="Kielpinski and Vinther Supplementary Table 1")),
    Segment("P7 reverse complement", il.P7_RC, "p7"),
), name="HRF-Seq representative indexed library")

HRF_SEQ_PRIMERS = (sp.TRUSEQ["R1"], sp.TRUSEQ["I1"], sp.TRUSEQ["R2"])


TITLE = 'HRF-Seq — hydroxyl-radical RNA accessibility'
NOTES = '01_hrf_seq.html'
SOURCE = 'Kielpinski and Vinther, Nucleic Acids Research (2014), <a href="https://doi.org/10.1093/nar/gku167">doi:10.1093/nar/gku167</a>.'
SUMMARY = 'Hydroxyl radicals cleave solvent-exposed ribose backbones; random-primed RT, a seven-base molecular barcode and paired-end sequencing recover both break and priming coordinates.'
CAVEAT = 'The defining demonstration probes purified RNAs. Hydroxyl-radical accessibility reflects backbone solvent exposure and tertiary packing, not only base pairing.'
FINAL_LIBRARY = HRF_FINAL_LIBRARY
PCR_GEOMETRY = HRF_PCR_GEOMETRY
FINAL_CAPTION = 'Representative indexed member of the published four-primer family. Read 1 begins with the seven-base molecular barcode; Read 2 reports the random-priming end.'
SEQ_PRIMERS = HRF_SEQ_PRIMERS
ACQUIRED_ROLES = ('Read 1', 'Index 1 (i7)', 'Read 2')
READ_LENGTHS = {'Read 1': 100, 'Index 1 (i7)': 6, 'Read 2': 100}
READ_USES = {'Read 1': 'extract_then_include', 'Index 1 (i7)': 'demultiplex', 'Read 2': 'include'}
SEQUENCING_INTRO = 'The source uses standard paired-end Illumina primers. Their sites, the first-read molecular barcode, sample index and second-read random-primer coordinate are derived from the published oligos.'
SEQUENCING_ENDING = ''
SEQUENCING_UNAVAILABLE = ''

def oligos():
    return ('RT_random_primer — 5′-AGACGTGTGCTCTTCCGATCTNNNNNNNNS-3′', 'RT_structure_cassette — 5′-AGACGTGTGCTCTTCCGATCTGAACCGGACCGAAGCCCG-3′', 'LIGATION_ADAPTER_RB — 5′-PHO-NNNNNNNAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT-3NHC3-3′', 'PCR_forward — 5′-AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCT-3′', *tuple((f'PCR reverse {index} — 5′-{_hrf_reverse(index)}-3′' for index in HRF_INDEXES)))

def workflow() -> Workflow:
    return chain(molecule("folded RNA pool", (role("RNA tertiary structure"),)), (
        ("Expose RNA to hydroxyl radicals alongside an untreated control",
         molecule("hydroxyl-radical-cleaved RNA", (role("intact RNA regions"),
                  role("backbone break at solvent-accessible ribose"))),
         "Backbone cleavage reports solvent accessibility rather than base identity."),
        ("Random-prime reverse transcription from primers carrying a 5-prime Illumina tail",
         molecule("break-terminated cDNA", (role("Illumina-tailed random primer"),
                  role("cDNA ending at radical-induced break"))),
         "The cDNA endpoint records the probing position; random priming records the opposite end."),
        ("Ligate a seven-base molecular-barcode adapter to each cDNA 3-prime end",
         molecule("barcoded HRF cDNA", (role("Illumina-tailed cDNA"),
                  role("7-nt molecular barcode", tag="umi", length_bp=7,
                       molecular_feature=feature("umi", "umi", "random")),
                  role("ligation adapter"))),
         "Unique barcode counting corrects PCR and fragment-length sampling bias."),
        ("PCR-add indexes, size-select and paired-end sequence",
         molecule("HRF-Seq library", (role("P5"), role("random-primer coordinate"),
                  role("accessibility-defined cDNA"), role("7-nt molecular barcode", tag="umi",
                       length_bp=7, molecular_feature=feature("umi", "umi", "random")),
                  role("sample index"), role("P7"))),
         "The paired reads recover the priming position and the hydroxyl-radical termination position."),
    ))



RECIPE = (('cleave', 'Expose RNA to hydroxyl radicals', StepKind.FRAGMENTATION), ('rt', 'Reverse-transcribe radical-induced breaks with PrimeScript', StepKind.REVERSE_TRANSCRIPTION), ('ligate', 'Ligate molecular-barcode adapter to cDNA with CircLigase', StepKind.LIGATION))

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
    graph.steps["library_pcr"] = replace(graph.steps["library_pcr"], label="Amplify with Phusion (4 plus 20 cycles)")
    graph.declare_reviewed_scope(reviewed, SOURCE)
    return graph.require_publishable()
