"""scMT-seq joint Smart-seq2 transcriptome and scRRBS methylome workflow."""

from chemdraw import Construct, Row
from batch_ngs import nextera_library, seg
from adapter_ligation import AdapterLigationGeometry
from multimodal_spatial import finish_pcr_branch
from pcr import PCRGeometry
from protocol_extensions import analyte_split
from protocol_graph import Disposition, Terminal, InputRole, MaterialKind, ProtocolGraph, StepKind
from dataclasses import replace
import rrbs
import rt
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "smart-seq2__10.1038+nmeth.2639" / "tools"))
import smartseq2 as s2


TITLE = "scMT-seq — transcriptome and DNA methylome"
NOTES = "01_scmt-seq.html"
SOURCE = ('Defining source: <a href="https://doi.org/10.1186/s13059-016-0950-z">'
          'Hu et al., <i>Genome Biology</i> (2016)</a>.')
SUMMARY = ("Cytoplasmic RNA and nuclear genomic DNA from one cell are physically "
           "separated for Smart-seq2 transcriptomics and scRRBS methylome sequencing.")
CAVEAT = ("Hu et al., Methods: Superscript III replaces Superscript II for RNA RT; methylome enrichment uses 20 then 12 PCR cycles, followed by 180–500 bp PAGE selection. The current graph summarizes both PCR rounds and does not separately draw this recovery. "
          "The defining paper explicitly inherits Smart-seq2/Nextera and scRRBS. "
          "Their source-described chemistry is reused here, so "
          "the earlier generic random-priming methylome representation is removed. "
          "INFERRED — the inherited scRRBS paper names standard premethylated indexed "
          "Illumina adapters but does not print their exact oligos or index member; "
          "dotted canonical TruSeq arms illustrate a compatible design. "
          "One inherited sequence is not source-verified, though: Hu et al. name the "
          "Smart-seq2 protocol but print no template-switch oligo, and Smart-seq2 "
          "has two published TSO linkers &mdash; <code>AC</code> for the chosen "
          "<code>rGrG+G</code> oligo in Picelli et al. 2013 Supplementary Table 2, "
          "<code>ACAT</code> in the <i>Nature Protocols</i> form. The longer form is "
          "modelled, so the whole TSO linker is drawn dotted as an inherited "
          "inference.")


RNA, RP = nextera_library([seg("Smart-seq2 cDNA insert", "X" * 36,
                               placeholder=True)], "scMT-seq RNA library")
DNA, DP = rrbs.final_library(), tuple(rrbs.SEQ_PRIMERS)
RNA_PREAMP = s2.preamp_template()
FINAL_LIBRARIES = (
    ("Transcriptome library", RNA, RP,
     "Smart-seq2 cDNA after Nextera tagmentation and indexing.",
     "Nextera sequencing primers on the RNA branch."),
    ("DNA methylome library", DNA, DP,
     "MspI-reduced, bisulfite-converted scRRBS library.",
     "TruSeq sequencing primers on the scRRBS branch."),
)


def oligos():
    return [
        f"Smart-seq2 oligo-dT30VN: 5′-{s2.OLIGO_DT}-3′",
        f"Smart-seq2 TSO (linker inherited, not printed by this paper): "
        f"5′-{rt.SMART_HANDLE}<inf>{s2.TSO_LINKER}</inf>{rt.TSO_G_TAIL_LNA}-3′",
        f"Smart-seq2 ISPCR primer: 5′-{s2.ISPCR}-3′",
        "scRRBS premethylated indexed TruSeq adapter: see reaction geometry",
    ]


def sections():
    return [
        ("Lyse one cell and retain its nucleus",
         [Row(chunks=[("cell → cytoplasmic RNA + nuclear genomic DNA", "r1", False)])],
         "The fractions are separated before either amplification."),
        ("Separate the matched analytes",
         analyte_split("cytoplasmic RNA → Smart-seq2", "nuclear DNA → scRRBS"),
         "Cell identity is maintained by plate position."),
    ]


def protocol_graph():
    g = ProtocolGraph(TITLE, authority=SOURCE)
    cell = g.input("cell", "one intact cell", MaterialKind.BIOLOGICAL,
                   state=sections()[0][1], primary=True)
    cytoplasm = g.product("cytoplasm", "cytoplasmic RNA fraction",
                          MaterialKind.MOLECULAR, state=sections()[1][1][:1])
    nucleus = g.product("nucleus", "nuclear genomic-DNA fraction",
                        MaterialKind.MOLECULAR, state=sections()[1][1][1:])
    g.react("separate", "Lyse the cell while retaining and separating its nucleus",
            StepKind.SPLIT, ((cell, InputRole.TEMPLATE),),
            (cytoplasm, nucleus), changes_molecule=False)

    rt_primer = g.input(
        "rt_primer", "Smart-seq2 anchored oligo(dT) primer", MaterialKind.OLIGO,
        state=(seg("SMART handle", rt.SMART_HANDLE, "tso"),
               seg("oligo-dT linker", s2.DT_LINKER),
               seg("dT30VN", "T" * s2.DT_LEN + s2.DT_ANCHOR,
                   placeholder=True)))
    rt_enzyme = g.input("rt", "Superscript III reverse transcriptase (study substitution)",
                        MaterialKind.CATALYST)
    copied = g.product("copied", "first-strand Smart-seq2 cDNA",
                       MaterialKind.MOLECULAR, state=RNA_PREAMP)
    g.react("reverse_transcribe", "Reverse-transcribe cytoplasmic poly(A) RNA",
            StepKind.REVERSE_TRANSCRIPTION,
            ((cytoplasm, InputRole.TEMPLATE), (rt_primer, InputRole.STRUCTURAL),
             (rt_enzyme, InputRole.CATALYST)), (copied,))
    tso = g.input(
        "tso", "Smart-seq2 template-switch oligo", MaterialKind.OLIGO,
        state=(seg("SMART handle", rt.SMART_HANDLE, "tso"),
               seg("TSO linker", s2.TSO_LINKER, inferred=True),
               seg("rGrG+G", "GGG", "tso")))
    switched = g.product("switched", "template-switched full-length cDNA",
                         MaterialKind.MOLECULAR, state=RNA_PREAMP)
    g.react("template_switch",
            "INFERRED — Template-switch the nascent cDNA (the TSO linker is "
            "inherited from Smart-seq2; this paper prints no TSO)",
            StepKind.TEMPLATE_SWITCHING,
            ((copied, InputRole.TEMPLATE), (tso, InputRole.STRUCTURAL),
             (rt_enzyme, InputRole.CATALYST)), (switched,))
    ispcr = g.input("ispcr", "Smart-seq2 ISPCR primer", MaterialKind.OLIGO,
                    state=(seg("ISPCR", s2.ISPCR, "tso"),))
    prepol = g.input("preamp_polymerase", "Smart-seq2 cDNA PCR polymerase",
                     MaterialKind.CATALYST)
    preamp = g.product("preamp", "amplified full-length Smart-seq2 cDNA",
                       MaterialKind.MOLECULAR, state=RNA_PREAMP)
    g.react("preamp_rna", "Amplify full-length cDNA with ISPCR", StepKind.PCR,
            ((switched, InputRole.TEMPLATE), (ispcr, InputRole.STRUCTURAL),
             (prepol, InputRole.CATALYST)), (preamp,),
            geometry=s2.ispcr_geometry(RNA_PREAMP))
    transposome = g.input("rna_transposome", "Nextera transposome",
                          MaterialKind.ADAPTER)
    tn5 = g.input("rna_tn5", "Tn5 transposase", MaterialKind.CATALYST)
    rna_tagged = g.product("rna_tagged", "Nextera-tagged Smart-seq2 cDNA",
                           MaterialKind.MOLECULAR, state=RNA)
    tn5_waste = g.product("rna_tn5_waste",
                          "candidate non-library cDNA tagmentation products", MaterialKind.MOLECULAR)
    g.react("tagment_rna", "Tagment amplified cDNA with Nextera",
            StepKind.TAGMENTATION,
            ((preamp, InputRole.TEMPLATE), (transposome, InputRole.STRUCTURAL),
             (tn5, InputRole.CATALYST)), (rna_tagged, tn5_waste))
    g.mark_terminal(Terminal(tn5_waste.id, Disposition.UNRESOLVED, "Candidate tagmentation products are summarized; molecular abundance and physical recovery were not established in this pass."))

    mspi = g.input("mspi", "MspI restriction endonuclease", MaterialKind.CATALYST)
    fragments = g.product("mspi_fragments", "MspI-digested nuclear DNA",
                          MaterialKind.MOLECULAR, state=rrbs.digested_fragment_scene())
    g.react("digest", "Digest nuclear DNA at C^CGG", StepKind.DIGESTION,
            ((nucleus, InputRole.TEMPLATE), (mspi, InputRole.CATALYST)), (fragments,))
    klenow = g.input("klenow", "Klenow exo- end-repair/dA-tailing enzyme",
                     MaterialKind.CATALYST)
    prepared = g.product("a_tailed_fragments", "filled and dA-tailed MspI fragments",
                         MaterialKind.MOLECULAR, state=rrbs.a_tailed_fragment_scene())
    g.react("end_prepare", "Fill recessed ends and add one 3′-A",
            StepKind.END_REPAIR,
            ((fragments, InputRole.TEMPLATE), (klenow, InputRole.CATALYST)),
            (prepared,))
    adapters = g.input("methylated_adapters",
                       "canonical premethylated indexed TruSeq illustration", MaterialKind.ADAPTER,
                       state=rrbs.adapter_scene(), inferred=True)
    ligase = g.input("ligase", "DNA ligase", MaterialKind.CATALYST)
    ligated = g.product("ligated_fragments", "adapter-ligated genomic fragments",
                        MaterialKind.MOLECULAR)
    geometry = AdapterLigationGeometry.unavailable(
        "The inherited scRRBS source does not print the ligation-facing adapter oligos or exact end geometry.")
    g.react("adapter_ligation", "INFERRED — ligate a compatible premethylated indexed Y-adapter",
            StepKind.LIGATION,
            ((prepared, InputRole.TEMPLATE), (adapters, InputRole.STRUCTURAL),
             (ligase, InputRole.CATALYST)), (ligated,), geometry=geometry)
    bisulfite = g.input("bisulfite", "bisulfite-conversion condition",
                        MaterialKind.CONDITION)
    converted = g.product("converted", "bisulfite-converted scRRBS templates",
                          MaterialKind.MOLECULAR, state=DNA)
    g.react("bisulfite_convert", "Bisulfite-convert adapter-ligated DNA",
            StepKind.CHEMICAL_CONVERSION,
            ((ligated, InputRole.TEMPLATE), (bisulfite, InputRole.CONDITION)),
            (converted,))

    reviewed = [
        "cell", "separate", "rt_primer", "rt", "reverse_transcribe", "tso",
        "template_switch", "ispcr", "preamp_polymerase", "preamp_rna",
        "rna_transposome", "rna_tn5", "tagment_rna", "mspi", "digest",
        "klenow", "end_prepare", "methylated_adapters", "ligase",
        "adapter_ligation", "bisulfite", "bisulfite_convert",
    ]
    finish_pcr_branch(g, "rna", rna_tagged, RNA, RP, reviewed,
                      label="scMT-seq transcriptome",
                      geometry=s2.nextera_index_geometry(RNA))
    finish_pcr_branch(
        g, "dna", converted, DNA, DP, reviewed, label="scMT-seq methylome",
        geometry=PCRGeometry.unavailable(
            DNA, reason=("The inherited scRRBS paper does not sequence the library-PCR "
                         "primers. Its canonical final Illumina arms are vendor-derived, "
                         "and the bisulfite-converted pre-PCR intermediate is not "
                         "represented base-for-base.")))
    g.steps["dna_pcr"] = replace(g.steps["dna_pcr"], label="Methylome enrichment: two PCR rounds, 20 then 12 cycles (summarized)")
    g.declare_reviewed_scope(reviewed, SOURCE)
    return g.require_publishable()
