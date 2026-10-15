from types import SimpleNamespace
from multimodal_protocols import ident, SAMPLE
from protocol_model import chain, role_library
from multimodal_protocols import Definition, InputRole, MaterialKind, ProtocolGraph, StepKind, _finish_pcr, molecule, role
DEFINITION = Definition("scONE-seq — one-tube single-cell DNA and RNA sequencing", "10.1126/sciadv.abp8901",
        "Yu et al., Science Advances (2023)", "Custom Tn5 adapters label genomic DNA while internally primed RT labels total RNA with a distinct branch barcode and UMI; both branches are amplified together.",
        "A six-base branch barcode separates DNA from RNA and an eight-base UMI identifies molecules within each branch.", "01_scone_seq.html")


def _scone_graph(d, final):
    source = f"{d.citation}; doi:{d.doi}"
    graph = ProtocolGraph(d.title, authority=source)
    cell = graph.input("cell", "lysed single cell or nucleus containing genomic DNA and total RNA",
                       MaterialKind.BIOLOGICAL, state=_scone().initial, primary=True)
    dna = graph.product("dna_branch", "genomic-DNA branch", MaterialKind.MOLECULAR,
                        state=cell.state)
    rna = graph.product("rna_branch", "total-RNA branch", MaterialKind.MOLECULAR,
                        state=cell.state)
    graph.react("separate_modalities", "Expose genomic DNA and total RNA as parallel substrates",
                StepKind.SPLIT, ((cell, InputRole.TEMPLATE),), (dna, rna),
                changes_molecule=False)
    transposome = graph.input("dna_transposome", "DNA-branch-barcode/UMI-loaded transposome",
                              MaterialKind.ADAPTER)
    tn5 = graph.input("tn5", "Tn5 transposase", MaterialKind.CATALYST)
    tagged_dna = graph.product("tagged_dna", "branch-barcoded molecularly tagged genomic DNA",
                               MaterialKind.MOLECULAR,
                               state=molecule("tagged DNA", (DNA_TYPE, DNA_UMI,
                                   role("tagmented genomic fragment"))))
    dna_waste = graph.product("dna_waste", "non-library tagmentation products",
                              MaterialKind.WASTE)
    graph.react("tagment_dna", "Tagment genomic DNA with the DNA-branch adapter",
                StepKind.TAGMENTATION,
                ((dna, InputRole.TEMPLATE), (transposome, InputRole.STRUCTURAL),
                 (tn5, InputRole.CATALYST)), (tagged_dna, dna_waste))
    graph.discard(dna_waste, "does not retain a molecularly tagged genomic insert")
    rt_primer = graph.input("rna_rt_primer", "RNA-branch-barcode/UMI internal RT primer",
                            MaterialKind.OLIGO)
    rt = graph.input("reverse_transcriptase", "reverse transcriptase", MaterialKind.CATALYST)
    cdna = graph.product("cdna", "internally primed branch-barcoded cDNA",
                         MaterialKind.MOLECULAR,
                         state=molecule("first-strand cDNA", (RNA_TYPE, RNA_UMI,
                             role("full-transcript cDNA"))))
    graph.react("reverse_transcribe", "Internally prime total RNA",
                StepKind.REVERSE_TRANSCRIPTION,
                ((rna, InputRole.TEMPLATE), (rt_primer, InputRole.STRUCTURAL),
                 (rt, InputRole.CATALYST)), (cdna,))
    tailing_enzyme = graph.input("tailing_enzyme", "poly(C)-tailing enzyme",
                                 MaterialKind.CATALYST)
    tailed = graph.product("tailed_cdna", "poly(C)-tailed cDNA with the opposite amplification end",
                           MaterialKind.MOLECULAR, state=cdna.state)
    graph.react("tail_cdna", "Poly(C)-tail the cDNA", StepKind.TAILING,
                ((cdna, InputRole.TEMPLATE), (tailing_enzyme, InputRole.CATALYST)), (tailed,))
    pooled = graph.product("pooled_templates", "co-amplifiable tagged DNA and cDNA",
                           MaterialKind.MOLECULAR,
                           state=molecule("co-amplifiable DNA and cDNA", (
                               role("common amplification end"), role("DNA or RNA insert"),
                               role("branch barcode and UMI"))))
    graph.react("pool_branches", "Pool DNA and RNA products for one-tube amplification",
                StepKind.POOL,
                ((tagged_dna, InputRole.TEMPLATE), (tailed, InputRole.TEMPLATE)), (pooled,))
    reviewed = ["cell", "separate_modalities", "dna_transposome", "tn5",
                "tagment_dna", "rna_rt_primer", "reverse_transcriptase",
                "reverse_transcribe", "tailing_enzyme", "tail_cdna", "pool_branches"]
    _finish_pcr(graph, "combined", pooled, final, d.endpoint, reviewed,
                label="scONE-seq combined DNA/RNA")
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()


def protocol_graph():
    d = DEFINITION
    final = FINAL_LIBRARY
    primers = SEQ_PRIMERS
    key = 'scone-seq'
    return _scone_graph(d, final)

SCONE_ONE_TN5 = "GTCTCGTGGGCTCGGTCATGNNNNNNNNAGATGTGTATAAGAGACAG"
SCONE_P5_TN5 = "TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG"
SCONE_MOSAIC = "CTGTCTCTTATACACATCT"
SCONE_RT_POLYT = "GTCTCGTGGGCTCGGAACTTNNNNNNNNTTTTTTTTTTTTTTTTTTTTVN"
SCONE_RT_GGG = "GTCTCGTGGGCTCGGAACTTNNNNNNNNGGGHN"
SCONE_RT_TTT = "GTCTCGTGGGCTCGGAACTTNNNNNNNNTTTVN"
SCONE_3P_ADAPTER = SCONE_RT_GGG
SCONE_AMP_PRIMER = "GATGTGTGGAGGTCTCGTGGGCTCGG"
SCONE_ILLU_R1 = "TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG"
SCONE_ILLU_R2 = "CCGAGCCCACGAGACCTGTCTCTTATACACATCT"
SCONE_ILLU_I7 = "AGATGTGTATAAGAGACAGGTCTCGTGGGCTCGG"

def _scone_oligos():
    """Exact oligos printed in scONE-seq Supplementary Table S2.

    Index windows remain symbolic because the table defines primer families, while the
    three custom Illumina sequencing primers are invariant concrete sequences.
    """
    return (
        "<b>One-Tn5 adaptor</b> — 5′-[biotin]" + SCONE_ONE_TN5 + "-3′",
        "<b>P5-Tn5 adaptor</b> — 5′-" + SCONE_P5_TN5 + "-3′",
        "<b>Mosaic</b> — 5′-[phos]" + SCONE_MOSAIC + "-3′",
        "<b>RT-polyT</b> — 5′-[biotin]" + SCONE_RT_POLYT + "-3′",
        "<b>RT-GGG</b> — 5′-[biotin]" + SCONE_RT_GGG + "-3′",
        "<b>RT-TTT</b> — 5′-[biotin]" + SCONE_RT_TTT + "-3′",
        "<b>3′ adaptor</b> — 5′-" + SCONE_3P_ADAPTER + "-3′",
        "<b>Amp-primer</b> — 5′-" + SCONE_AMP_PRIMER + "-3′",
        "<b>P7 index-primer family</b> — 5′-CAAGCAGAAGACGGCATACGAGAT[i7]AGATGTGTATAAGAGACAGGTCTCGTGGGCTCGG-3′",
        "<b>P5 index-primer family</b> — 5′-AATGATACGGCGACCACCGAGATCTACAC[i5]TCGTCGGCAGCGTC-3′",
        "<b>ONE-Illu Read 1</b> — 5′-" + SCONE_ILLU_R1 + "-3′",
        "<b>ONE-Illu Read 2</b> — 5′-" + SCONE_ILLU_R2 + "-3′",
        "<b>ONE-Illu I7 Index Primer</b> — 5′-" + SCONE_ILLU_I7 + "-3′",
    )


DNA_TYPE = ident("DNA-branch barcode", "dna_branch", "inline_barcode", length=6)
RNA_TYPE = ident("RNA-branch barcode", "rna_branch", "inline_barcode", length=6)
DNA_UMI = ident("DNA UMI", "dna_umi", "umi", length=8, tag="umi")
RNA_UMI = ident("RNA UMI", "rna_umi", "umi", length=8, tag="umi")


def _scone():
    return chain(molecule("lysed single cell or nucleus", (role("genomic DNA"),
                                                            role("total RNA"))), (
        ("Tagment genomic DNA with a custom molecularly barcoded adapter",
         molecule("tagged DNA", (DNA_TYPE, DNA_UMI, role("tagmented genomic fragment"))),
         "The six-base branch barcode distinguishes DNA; the eight-base UMI identifies an original tagged molecule."),
        ("Internally prime total RNA with a separately barcoded RT primer",
         molecule("first-strand cDNA", (RNA_TYPE, RNA_UMI, role("full-transcript cDNA"))),
         "Internal priming captures polyadenylated and non-polyadenylated RNA."),
        ("Poly(C)-tail cDNA and add the opposite amplification adapter",
         molecule("co-amplifiable DNA and cDNA", (role("common amplification end"),
                                                   role("DNA or RNA insert"),
                                                   role("branch barcode and UMI"))),
         "DNA-specific and RNA-specific identifiers permit both products to share one-tube amplification."),
        ("Amplify together and complete the indexed sequencing library",
         molecule("scONE-seq library", (role("P5 end"), SAMPLE,
                                         role("branch barcode"), role("UMI"),
                                         role("DNA or cDNA insert"), role("P7 end"))),
         "Reads are separated by the molecular branch barcode after sequencing."),
    ))

d = DEFINITION
FINAL_LIBRARY = role_library(d.title + " final library roles", (role("DNA/RNA branch barcode"), role("branch UMI"), role("insert")), index=SAMPLE)
SPEC = SimpleNamespace(
    TITLE=d.title,
    NOTES=d.notes,
    SOURCE=f'Defining source: <a href="https://doi.org/{d.doi}">{d.citation}</a>.',
    SUMMARY=d.summary,
    CAVEAT="Supplementary Table S2 publishes the custom molecular adapters, "
                       "amplification oligos and custom Illumina Read 1, Read 2 and i7 "
                       "sequencing primers. Its ONE-Illu i5 row instead directs users to "
                       "the sequencing-equipment instructions; that one primer is not "
                       "reconstructed here.",
    FINAL_LIBRARY=FINAL_LIBRARY,
    SEQ_PRIMERS=(),
    FINAL_CAPTION=d.endpoint,
    SEQUENCING_INTRO=d.endpoint,
    READ_LENGTHS=None,
    SEQUENCING_ENDING=d.endpoint,
    SEQUENCING_UNAVAILABLE="Custom Read 1, Read 2 and i7 primer sequences are "
                                       "published and listed above. The source leaves the "
                                       "ONE-Illu i5 primer to the instrument manufacturer, "
                                       "so a complete four-primer landing model remains mixed.",
    workflow=_scone,
    oligos=_scone_oligos,
)
SPEC.protocol_graph = protocol_graph
globals().update(vars(SPEC))
