"""Herzog et al. SLAM-seq chemistry and QuantSeq endpoint."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
from base_conversion import slam_path
from batch_ngs import seg, truseq_library
from chemdraw import Row

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
SEQUENCING_INTRO = ("The inferred QuantSeq FWD adapter skeleton places standard TruSeq Read 1, Index 1 and Read 2 sites. The defining experiment used single-read sequencing; Read 1 is the informative 3′-end read.")

def sections():
    return [
        ("Label newly synthesized RNA", [Row(chunks=[("nascent RNA incorporates s4U during the pulse; pre-existing RNA retains U", "w1", False)])], "The nucleotide analogue marks synthesis time without biochemical enrichment."),
        ("Alkylate the thiol", [Row(chunks=[(slam_path(labelled=True).text(), "w1", False)]), Row(chunks=[("10 mM IAA · 50% DMSO · 50 mM phosphate pH 8 · 50 °C · 15 min → DTT quench", None, False)])], "Iodoacetamide covalently carboxyamidomethylates s4U; ordinary U is unchanged."),
        ("Reverse-transcribe and build the 3′-end library", [Row(chunks=[("poly(A) RNA → oligo(dT)-anchored QuantSeq RT → second strand → indexed PCR", None, True)]), Row(chunks=[("labeled position: RNA s4U → cDNA G → sequenced C     ordinary RNA U → cDNA A → sequenced T", "w1", False)])], "INFERRED — the kit-specific adapter bases are unavailable. The defining paper establishes the 3′-end workflow and T-to-C readout."),
    ]
