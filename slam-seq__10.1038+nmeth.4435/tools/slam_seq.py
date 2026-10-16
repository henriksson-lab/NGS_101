"""Herzog et al. SLAM-seq chemistry and QuantSeq endpoint."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
from base_conversion import BasePath
from batch_ngs import seg, truseq_library
from chemdraw import Row
from pcr import PCRGeometry
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind)



def slam_path(*, labelled: bool) -> BasePath:
    """SLAM-seq base path through IAA derivatization and RT readout."""
    return BasePath(("s4U", "IAA-s4U", "C")) if labelled else BasePath(("U", "T"))

TITLE = "SLAM-seq — thiol-linked metabolic RNA sequencing"
NOTES = "01_slam-seq.html"
SOURCE = ('Defining source: <a href="https://doi.org/10.1038/nmeth.4435">Herzog et al., <i>Nature Methods</i> (2017)</a>.')
SUMMARY = ("Metabolically label new RNA with 4-thiouridine, alkylate its thiol with iodoacetamide, and prepare a 3′ mRNA library in which reverse transcription records labeled uridines as T-to-C substitutions.")
CAVEAT = ("The paper specifies the commercial Lexogen QuantSeq FWD kit but does not print its adapter oligos. The dotted adapter skeleton is therefore an inferred single-index Illumina-compatible endpoint; the s4U chemistry and readout are published.")
ALKYLATION = {"IAA_mM": 10, "DMSO_percent": 50, "phosphate_mM": 50, "pH": 8, "temperature_C": 50, "minutes": 15}

FINAL_LIBRARY, SEQ_PRIMERS = truseq_library(
    [seg("INFERRED — QuantSeq 3′-end cDNA", "X" * 42, placeholder=True, inferred=True)],
    "SLAM-seq QuantSeq library", dual_index=False, inferred_adapters=True)
FINAL_CAPTION = ("INFERRED — a single-index Illumina-compatible QuantSeq endpoint. Within the cDNA insert, IAA-alkylated s4U is reported as C where the reference RNA has U.")
SEQUENCING_INTRO = ("The inferred QuantSeq FWD adapter skeleton places standard TruSeq Read 1, Index 1 and Read 2 sites. The authors' detailed SLAMseq protocol recommends single-read 100 (SR100) sequencing on a HiSeq 2500; Read 1 is the informative 3′-end read.")

def sections():
    return [
        ("Label newly synthesized RNA", [Row(chunks=[("nascent RNA incorporates s4U during the pulse; pre-existing RNA retains U", "w1", False)])], "The nucleotide analogue marks synthesis time without biochemical enrichment."),
        ("Alkylate the thiol", [Row(chunks=[(slam_path(labelled=True).text(), "w1", False)]), Row(chunks=[("10 mM IAA · 50% DMSO · 50 mM phosphate pH 8 · 50 °C · 15 min → DTT quench", None, False)])], "Iodoacetamide covalently carboxyamidomethylates s4U; ordinary U is unchanged."),
        ("Reverse-transcribe and build the 3′-end library", [Row(chunks=[("poly(A) RNA → oligo(dT)-anchored QuantSeq RT → second strand → indexed PCR", None, True)]), Row(chunks=[("labeled position: RNA s4U → cDNA G → sequenced C     ordinary RNA U → cDNA A → sequenced T", "w1", False)])], "INFERRED — the kit-specific adapter bases are unavailable. The defining paper establishes the 3′-end workflow and T-to-C readout."),
    ]


def protocol_graph():
    graph=ProtocolGraph(TITLE,authority=SOURCE)
    cells=graph.input("cells","cells during the metabolic pulse",MaterialKind.BIOLOGICAL,
                      primary=True)
    s4u=graph.input("s4u","4-thiouridine metabolic label",MaterialKind.MOLECULAR)
    rnap=graph.input("rna_polymerase","cellular RNA polymerases",MaterialKind.CATALYST)
    labelled=graph.product("labelled_rna","newly synthesized s4U-labelled RNA",
                           MaterialKind.MOLECULAR)
    graph.react("pulse_label","Transcribe RNA during the s4U pulse",StepKind.TRANSCRIPTION,
                ((cells,InputRole.TEMPLATE),(s4u,InputRole.STRUCTURAL),
                 (rnap,InputRole.CATALYST)),(labelled,))
    iaa=graph.input("iaa","iodoacetamide alkylation condition",MaterialKind.CONDITION)
    alkylated=graph.product("alkylated","IAA-alkylated s4U RNA",MaterialKind.MOLECULAR)
    graph.react("alkylate","Alkylate s4U thiols with iodoacetamide",
                StepKind.CHEMICAL_CONVERSION,
                ((labelled,InputRole.TEMPLATE),(iaa,InputRole.CONDITION)),(alkylated,))
    rt_primer=graph.input("oligo_dt","QuantSeq oligo-dT RT primer",MaterialKind.OLIGO,
                          inferred=True,source=SOURCE+" (kit oligo sequence undisclosed)")
    rt=graph.input("reverse_transcriptase","QuantSeq reverse transcriptase",
                   MaterialKind.CATALYST)
    cdna=graph.product("cdna","3-prime-end first-strand cDNA",MaterialKind.MOLECULAR,
                       inferred=True,source=SOURCE+" (commercial kit construct)")
    graph.react("reverse_transcribe","Reverse-transcribe from the poly(A) boundary",
                StepKind.REVERSE_TRANSCRIPTION,
                ((alkylated,InputRole.TEMPLATE),(rt_primer,InputRole.STRUCTURAL),
                 (rt,InputRole.CATALYST)),(cdna,))
    second_primer=graph.input("second_primer","QuantSeq second-strand primer/adapter",
                              MaterialKind.OLIGO,inferred=True,
                              source=SOURCE+" (kit oligo sequence undisclosed)")
    second_polymerase=graph.input("second_polymerase","second-strand polymerase",
                                  MaterialKind.CATALYST)
    duplex=graph.product("duplex","QuantSeq double-stranded library precursor",
                         MaterialKind.MOLECULAR,inferred=True,
                         source=SOURCE+" (commercial kit construct)")
    graph.react("second_strand","Synthesize the second strand and opposite adapter end",
                StepKind.EXTENSION,
                ((cdna,InputRole.TEMPLATE),(second_primer,InputRole.STRUCTURAL),
                 (second_polymerase,InputRole.CATALYST)),(duplex,))
    fwd=graph.input("pcr_forward","QuantSeq PCR forward primer",MaterialKind.OLIGO,
                    inferred=True,source=SOURCE+" (kit oligo sequence undisclosed)")
    rev=graph.input("pcr_reverse","QuantSeq indexed PCR reverse primer",MaterialKind.OLIGO,
                    inferred=True,source=SOURCE+" (kit oligo sequence undisclosed)")
    polymerase=graph.input("pcr_polymerase","library PCR polymerase",MaterialKind.CATALYST)
    library=graph.product("library",FINAL_LIBRARY.name,MaterialKind.LIBRARY,
                          state=FINAL_LIBRARY,inferred=True,
                          source=SOURCE+" (inferred Illumina-compatible kit endpoint)")
    graph.react("library_pcr","Amplify and index the 3-prime-end library",StepKind.PCR,
                ((duplex,InputRole.TEMPLATE),(fwd,InputRole.STRUCTURAL),
                 (rev,InputRole.STRUCTURAL),(polymerase,InputRole.CATALYST)),
                (library,),geometry=PCRGeometry.unavailable(
                    FINAL_LIBRARY,reason="The QuantSeq kit primers are not printed in the defining paper."))
    graph.sequence(library,SequencingPlan(library=FINAL_LIBRARY,
                                          primers=tuple(SEQ_PRIMERS)))
    graph.declare_reviewed_scope((
        "cells","s4u","rna_polymerase","pulse_label","iaa","alkylate",
        "oligo_dt","reverse_transcriptase","reverse_transcribe","second_primer",
        "second_polymerase","second_strand","pcr_forward","pcr_reverse",
        "pcr_polymerase","library_pcr","library"),SOURCE)
    return graph.require_publishable()
