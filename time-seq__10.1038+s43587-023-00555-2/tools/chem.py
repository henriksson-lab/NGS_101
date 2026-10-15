"""Protocol-owned source, oligos, workflow and reaction graph."""
import illumina as il
import nextera as nx
import seqprimers as sp
from adapter_ligation import AdapterLigationGeometry, reported_adapter_end
from batch_ngs import nextera_library, seg, truseq_library
from chemdraw import Construct, MolecularState, Scene, Segment, Workflow, feature, revcomp
from circular import CircularizationGeometry
from endprep import ThreePrimeOverhangAdapter, ligation_product_scene, react_dA_tail_then_ligate, repair_and_dA_tail
from protocol_model import closed_circle, molecule, role
from pcr import PCRGeometry, PCRPrimer
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind
from accuracy_protocols import _workflow

TIME_B_PREFIX = "TGGGTGGAGGGTGG"

TIME_BARCODE = "AAAGA"                 # TS301, one published representative

TIME_TN5_REVERSE = "/5Phos/CTGTCTCTTATACACATCT"

TIME_METHYLATED_A = "T/iMe-dC/GT/iMe-dC/GG/iMe-dC/AG/iMe-dC/GT/iMe-dC/AGATGTGTATAAGAGA/iMe-dC/AG"

TIME_B_TS301 = TIME_B_PREFIX + TIME_BARCODE + nx.ME

TIME_BLOCK_B = "CTGTCTCTTATACACATCTHHHHHCCACCCTCCACCCA"

TIME_BLOCK_A = "CTGTCTCTTATACACATCTCCGAGCCCACGAGAC"

TIME_N501 = "AATGATACGGCGACCACCGAGATCTACACTAGATCGCTCGTCGGCAGCGTC"

TIME_N701 = "CAAGCAGAAGACGGCATACGAGATTCGCCTTAAAGCAGTGGTATCAACGCAGATCTGGGTGGAGGGTGG"

TIME_I5 = "TAGATCGC"

TIME_I7 = "TCGCCTTA"

TIME_CUSTOM_R2 = "AAGCAGTGGTATCAACGCAGATC" + TIME_B_PREFIX

TIME_CUSTOM_I1 = revcomp(TIME_CUSTOM_R2)


def _time_library():
    lib = Construct([
        seg("P5", il.P5, "p5"),
        seg("i5 (N501)", TIME_I5, "cbc",
            feature=feature("pool_index_i5", "sample_index", "whitelist",
                            whitelist="TIME-seq Supplementary Table 5 indexed PCR primers")),
        seg("S5", nx.S5, "s5"), seg("left mosaic end", nx.ME, "me"),
        role("bisulfite-converted captured insert"),
        seg("right mosaic end reverse complement", nx.ME_RC, "me"),
        seg("Tn5 sample barcode reverse complement", revcomp(TIME_BARCODE), "cbc",
            feature=feature("tn5_sample_barcode", "inline_barcode", "whitelist",
                            whitelist="TIME-seq Supplementary Table 5 TS301–TS364")),
        seg("B-adaptor prefix reverse complement", revcomp(TIME_B_PREFIX)),
        seg("custom sequencing handle reverse complement",
            revcomp("AAGCAGTGGTATCAACGCAGATC"), "r2"),
        seg("i7 reverse complement (N701)", revcomp(TIME_I7), "cbc",
            feature=feature("pool_index_i7", "sample_index", "whitelist",
                            whitelist="TIME-seq Supplementary Table 5 indexed PCR primers")),
        seg("P7 reverse complement", il.P7_RC, "p7"),
    ], name="TIME-seq final library")
    primers = (sp.NEXTERA["R1"],
               sp.custom("Index 1 (i7)", "TIME-seq custom Index 1", TIME_CUSTOM_I1,
                         "Griffin et al. Supplementary Table 5"),
               sp.NEXTERA["I2"],
               sp.custom("Read 2", "TIME-seq custom Read 2", TIME_CUSTOM_R2,
                         "Griffin et al. Supplementary Table 5"))
    errors = sp.verify(lib, primers, required_roles=tuple(p.role for p in primers))
    if errors:
        raise ValueError("TIME-seq final library: " + "; ".join(errors))
    return lib, primers


def workflow() -> Workflow:
    dna = molecule("sample genomic DNA", (role("genomic DNA"),))
    tagged = molecule("sample-indexed tagmentation product", (
        role("methylated A adapter"), role("genomic fragment"),
        role("5-nt Tn5 sample barcode", tag="cbc",
             molecular_feature=feature("tn5_sample_barcode", "inline_barcode", "whitelist",
                                       whitelist="TIME-seq Supplementary Table 5 TS301–TS364")),
        role("B adapter")),
        junctions=(("methylated A adapter", "genomic fragment", "Tn5 strand transfer"),
                   ("genomic fragment", "5-nt Tn5 sample barcode", "Tn5 strand transfer")))
    pooled = molecule("pooled indexed DNA", (role("many Tn5-barcoded samples"),))
    captured = molecule("hybrid-captured pool", (role("biotin-RNA bait"), role("target loci")))
    converted = molecule("bisulfite-converted captured DNA", (role("converted target loci"),))
    return _workflow(dna, (
        ("Tagment each sample with methylated A and indexed B transposomes", tagged,
         "The 5-nt B-adapter barcode identifies the sample before pooling."),
        ("Pool samples and perform methylated end repair", pooled,
         "5-methyl-dCTP replaces dCTP so adapter cytosines survive bisulfite treatment."),
        ("Hybridize biotinylated RNA baits and enrich target loci", captured,
         "Capture occurs before bisulfite conversion."),
        ("Bisulfite-convert the captured DNA", converted,
         "Unmethylated cytosines become uracils while protected methylcytosines remain."),
        ("PCR with indexed P5 and P7 primers", molecule("final TIME-seq library", (
            role("P5/i5 end"), role("converted captured insert"), role("Tn5 sample barcode"),
            role("P7/i7 end"))),
         "Pool indices and the internal Tn5 barcode support two-level demultiplexing."),
    ))


FINAL_LIBRARY, SEQ_PRIMERS = _time_library()
TITLE = 'TIME-seq — tagmentation-based indexing for methylation sequencing'
NOTES = '01_time_seq.html'
SOURCE = '<a href="https://doi.org/10.1038/s43587-023-00555-2">Griffin et al., <i>Nature Aging</i> (2024)</a>, Supplementary Information and Supplementary Table 5.'
SUMMARY = 'Bisulfite-resistant Tn5 adapters install a 5-nt sample barcode before pooling; methylated end repair, hybrid capture, bisulfite conversion and indexed PCR produce the targeted methylation library.'
CAVEAT = 'The displayed B-adapter uses published representative TS301; the paper publishes 64 alternative 5-nt sample indexes.'
READ_LENGTHS = {'Read 1': 145, 'Index 1 (i7)': 8, 'Index 2 (i5)': 8, 'Read 2': 5}
SHOW_READ_CONTENT = True
FINAL_CAPTION = 'The internal five-base Tn5 barcode identifies the pre-pooling sample; i5 and i7 identify the PCR pool. TS301 is drawn as one published representative.'
SEQUENCING_INTRO = 'Read 2 and Index 1 use the two published custom sequencing primers; Read 1 and Index 2 use the matching Nextera sites. The economical run described in the Supplementary Information uses 145–153 Read 1 cycles, 8 i7 cycles, optional 8 i5 cycles and 5 Read 2 cycles; this view uses the 145-cycle Read 1 endpoint.'
def oligos():
    return (f'<p><b>Tn5 reverse</b>: 5′-{TIME_TN5_REVERSE}-3′</p>', f'<p><b>Methylated adaptor A</b>: 5′-{TIME_METHYLATED_A}-3′</p>', f'<p><b>Indexed adaptor B, TS301</b>: 5′-{TIME_B_TS301}-3′</p>', f'<p><b>Adaptor B blocking oligo</b>: 5′-{TIME_BLOCK_B}-3′</p>', f'<p><b>Adaptor A blocking oligo</b>: 5′-{TIME_BLOCK_A}-3′</p>', f'<p><b>Indexed PCR primer N501</b>: 5′-{TIME_N501}-3′</p>', f'<p><b>Indexed PCR primer N701</b>: 5′-{TIME_N701}-3′</p>', f'<p><b>Custom Read 2 sequencing primer</b>: 5′-{TIME_CUSTOM_R2}-3′</p>', f'<p><b>Custom Index 1 sequencing primer</b>: 5′-{TIME_CUSTOM_I1}-3′</p>')

def protocol_graph():
    library = FINAL_LIBRARY
    primers = SEQ_PRIMERS
    source = ("Griffin et al., Nature Aging (2024), Supplementary Information and "
              "Supplementary Table 5; doi:10.1038/s43587-023-00555-2")
    graph = ProtocolGraph("TIME-seq", authority=source)
    batch = graph.input(
        "sample_batch", "batch of genomic-DNA samples",
        MaterialKind.BIOLOGICAL,
        state=molecule("sample genomic DNA batch", (role("independent genomic-DNA samples"),)),
        primary=True)
    sample_a = graph.product("sample_a", "representative genomic-DNA sample A",
                             MaterialKind.MOLECULAR, state=batch.state)
    sample_b = graph.product("sample_b", "representative genomic-DNA sample B",
                             MaterialKind.MOLECULAR, state=batch.state)
    graph.react("separate_samples", "Keep samples separate for Tn5 indexing",
                StepKind.SPLIT, ((batch, InputRole.TEMPLATE),), (sample_a, sample_b),
                changes_molecule=False)
    reviewed = ["sample_batch", "separate_samples"]

    def tagment(prefix: str, sample, barcode_label: str):
        adapters = graph.input(
            prefix + "_adapters",
            "methylated A and indexed B transposome adapters (" + barcode_label + ")",
            MaterialKind.ADAPTER)
        tn5 = graph.input(prefix + "_tn5", "Tn5 transposase for " + barcode_label,
                          MaterialKind.CATALYST)
        tagged = graph.product(
            prefix + "_tagged", barcode_label + " tagmented genomic DNA",
            MaterialKind.MOLECULAR,
            state=molecule("sample-indexed tagmentation product", (
                role("methylated A adapter"), role("genomic fragment"),
                role("5-nt Tn5 sample barcode", tag="cbc",
                     molecular_feature=feature(
                         "tn5_sample_barcode", "inline_barcode", "whitelist",
                         whitelist="TIME-seq Supplementary Table 5 TS301–TS364")),
                role("B adapter"))))
        losses = graph.product(prefix + "_losses", barcode_label + " non-library Tn5 products",
                               MaterialKind.WASTE)
        graph.react(
            prefix + "_tagment", "Tagment " + barcode_label + " with its indexed B transposome",
            StepKind.TAGMENTATION,
            ((sample, InputRole.TEMPLATE), (adapters, InputRole.STRUCTURAL),
             (tn5, InputRole.CATALYST)), (tagged, losses))
        graph.discard(losses, "does not retain an amplifiable A/B tagmentation configuration")
        reviewed.extend((prefix + "_adapters", prefix + "_tn5", prefix + "_tagment"))
        return tagged

    tagged_a = tagment("a", sample_a, "sample A")
    tagged_b = tagment("b", sample_b, "sample B")
    pooled = graph.product(
        "pooled_samples", "pooled pre-indexed TIME-seq DNA",
        MaterialKind.MOLECULAR,
        state=molecule("pooled indexed DNA", (role("many Tn5-barcoded samples"),)))
    graph.react("pool_samples", "Pool samples after Tn5 barcode installation",
                StepKind.POOL,
                ((tagged_a, InputRole.TEMPLATE), (tagged_b, InputRole.TEMPLATE)),
                (pooled,), changes_molecule=False)
    reviewed.append("pool_samples")

    repair = graph.input("methylated_repair_enzymes", "methylated end-repair enzyme mix",
                         MaterialKind.CATALYST)
    methylated_dctp = graph.input("methylated_dctp", "5-methyl-dCTP repair nucleotide",
                                  MaterialKind.MOLECULAR)
    repaired = graph.product("repaired_pool", "methylated-end-repaired indexed DNA pool",
                             MaterialKind.MOLECULAR, state=pooled.state)
    graph.react("methylated_end_repair", "Perform methylated end repair",
                StepKind.END_REPAIR,
                ((pooled, InputRole.TEMPLATE), (methylated_dctp, InputRole.STRUCTURAL),
                 (repair, InputRole.CATALYST)), (repaired,))
    baits = graph.input("biotin_rna_baits", "biotinylated target-locus RNA baits",
                        MaterialKind.OLIGO)
    hybrids = graph.product(
        "bait_hybrids", "DNA hybridized to biotinylated target baits",
        MaterialKind.MOLECULAR,
        state=molecule("target-bait hybrids", (role("biotin-RNA bait"), role("target loci"))))
    graph.react("bait_hybridization", "Hybridize biotinylated RNA baits",
                StepKind.HYBRIDIZATION,
                ((repaired, InputRole.TEMPLATE), (baits, InputRole.STRUCTURAL)),
                (hybrids,))
    capture = graph.input("capture_beads", "streptavidin target-capture beads",
                          MaterialKind.SELECTION)
    captured = graph.product(
        "captured_pool", "hybrid-captured target DNA", MaterialKind.MOLECULAR,
        state=molecule("hybrid-captured pool", (role("biotin-RNA bait"), role("target loci"))))
    uncaptured = graph.product("uncaptured_pool", "uncaptured off-target DNA",
                               MaterialKind.WASTE)
    graph.react("target_capture", "Enrich bait-bound target loci", StepKind.AFFINITY_CAPTURE,
                ((hybrids, InputRole.TEMPLATE), (capture, InputRole.SELECTION)),
                (captured, uncaptured))
    graph.discard(uncaptured, "does not hybridize to the target-locus bait panel")
    bisulfite = graph.input("bisulfite", "bisulfite conversion condition",
                            MaterialKind.CONDITION)
    converted = graph.product(
        "converted_pool", "bisulfite-converted captured DNA", MaterialKind.MOLECULAR,
        state=molecule("bisulfite-converted captured DNA", (role("converted target loci"),)))
    graph.react("bisulfite_conversion", "Bisulfite-convert captured DNA",
                StepKind.CHEMICAL_CONVERSION,
                ((captured, InputRole.TEMPLATE), (bisulfite, InputRole.CONDITION)),
                (converted,))
    n501 = graph.input("n501_primer", "published indexed PCR primer N501",
                       MaterialKind.OLIGO)
    n701 = graph.input("n701_primer", "published indexed PCR primer N701",
                       MaterialKind.OLIGO)
    polymerase = graph.input("pcr_polymerase", "bisulfite-library PCR polymerase",
                             MaterialKind.CATALYST)
    final = graph.product("sequencing_library", "final TIME-seq library",
                          MaterialKind.LIBRARY, state=library)
    geometry = PCRGeometry.exact(
        library, PCRPrimer("N501", TIME_N501), PCRPrimer("N701", TIME_N701),
        statement="Griffin et al. Supplementary Table 5 publishes both indexed PCR primers.")
    graph.react("indexed_pcr", "PCR with indexed P5 and P7 primers", StepKind.PCR,
                ((converted, InputRole.TEMPLATE), (n501, InputRole.STRUCTURAL),
                 (n701, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)),
                (final,), geometry=geometry)
    graph.sequence(final, SequencingPlan(
        library, tuple(primers), tuple(p.role for p in primers),
        {"Read 1": 145, "Index 1 (i7)": 8, "Index 2 (i5)": 8, "Read 2": 5},
        {"Read 1": "include", "Index 1 (i7)": "demultiplex",
         "Index 2 (i5)": "demultiplex", "Read 2": "extract"}))
    reviewed.extend((
        "methylated_repair_enzymes", "methylated_dctp", "methylated_end_repair",
        "biotin_rna_baits", "bait_hybridization", "capture_beads", "target_capture",
        "bisulfite", "bisulfite_conversion", "n501_primer", "n701_primer",
        "pcr_polymerase", "indexed_pcr", "sequencing_library",
    ))
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()
