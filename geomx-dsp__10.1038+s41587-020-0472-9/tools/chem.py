"""GeoMx DSP NGS readout page model."""
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))

import digital_spatial_tags as dsp
from pcr import PCRGeometry
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            SourceClaim, SourceDiscrepancy, StepKind)


TITLE = "GeoMx DSP — photocleaved spatial tags with NGS readout"
NOTES = "01_geomx_dsp.html"
SOURCE = ('Defining method: <a href="https://doi.org/10.1038/s41587-020-0472-9">'
          'Merritt et al., <i>Nature Biotechnology</i> (2020)</a>. Library preparation: '
          '<a href="https://www.brukerspatialbiology.com/wp-content/uploads/2022/06/'
          'MAN-10153-01-GeoMx-DSP-NGS-Readout-User-Manual.pdf">GeoMx DSP NGS Readout '
          'User Manual MAN-10153-01</a>.')
SUMMARY = ("Patterned UV releases target-identifying synthetic tags from a selected "
           "tissue area. A well-specific unique-dual-index PCR preserves that AOI, while "
           "a 12-nt RTS ID identifies the RNA or protein target and a 14-nt UMI counts molecules.")
CAVEAT = ("The vendor publishes the construct order and exact lengths but not the tag or "
          "Seq Code bases. Canonical TruSeq bases are shown with inferred styling solely "
          "to explain standard-primer binding; UMI, RTS-ID and UDI values remain placeholders.")
FINAL_LIBRARY, SEQ_PRIMERS = dsp.library()
READ_LENGTHS = {"Read 1": 27, "Index 1 (i7)": 8, "Index 2 (i5)": 8, "Read 2": 27}
SHOW_READ_CONTENT = True
FINAL_CAPTION = ("Vendor-published final order and lengths. Inferred canonical adapter "
                 "bases are visibly marked; the central 14-nt UMI and 12-nt target ID, "
                 "and the two AOI indexes, are structured features with undisclosed values.")
SEQUENCING_INTRO = ("MAN-10153-01 specifies paired-end 2×27 plus two 8-cycle index "
                    "reads. Read 1 sees UMI then RTS ID; Read 2 sees RTS ID then UMI. "
                    "Index 2 orientation depends on the Illumina instrument workflow. "
                    "The defining article's own runs used longer reads (MiSeq 2×75 or "
                    "NextSeq 2×37); see the source note above.")


def oligos():
    return ()


def final_rows():
    return dsp.final_rows(FINAL_LIBRARY)


def protocol_graph():
    graph=ProtocolGraph(TITLE,authority=SOURCE)
    tissue=graph.input("tissue","fixed tissue with target molecules",
                       MaterialKind.BIOLOGICAL,primary=True)
    probe=graph.input("probe","affinity probe carrying a photocleavable DSP tag",
                      MaterialKind.MOLECULAR)
    bound=graph.product("bound","probe-bound fixed tissue",MaterialKind.MOLECULAR,
                        state=tuple(dsp.probe_rows()))
    graph.react("bind_probe","Bind target-specific tagged probes in tissue",
                StepKind.GENERIC,
                ((tissue,InputRole.TEMPLATE),(probe,InputRole.STRUCTURAL)),
                (bound,),changes_molecule=False)
    uv=graph.input("patterned_uv","patterned UV cleavage condition",MaterialKind.CONDITION)
    tags=graph.product("tags","photocleaved tags collected from one AOI",
                       MaterialKind.MOLECULAR,state=dsp.photocleaved_tag())
    residual=graph.product("residual_tissue","residual tissue-bound probe material",
                           MaterialKind.WASTE)
    graph.react("photocleave","Release synthetic tags from the selected AOI",
                StepKind.CHEMICAL_CONVERSION,
                ((bound,InputRole.TEMPLATE),(uv,InputRole.CONDITION)),(tags,residual))
    graph.discard(residual,"tissue-bound probe remnants are not part of the collected library")
    # MAN-10153-01 Library Preparation steps 3-7 (Table 5): the collection plate is
    # dried (bench overnight or 65 C, 1 h, lid open) and rehydrated in nuclease-free
    # water (10 uL RNA, 80 uL protein) before 4 uL enters each PCR well.
    drydown=graph.input("drydown_rehydration",
                        "dry-down (65 °C 1 h or overnight) and rehydration in water "
                        "(10 µL RNA / 80 µL protein, 10 min RT)",MaterialKind.CONDITION)
    rehydrated=graph.product("rehydrated_tags","rehydrated AOI aspirate in one plate well",
                             MaterialKind.MOLECULAR,state=dsp.photocleaved_tag())
    graph.react("rehydrate","Dry down the collection plate and rehydrate each aspirate",
                StepKind.GENERIC,
                ((tags,InputRole.TEMPLATE),(drydown,InputRole.CONDITION)),
                (rehydrated,),changes_molecule=False)
    i5=graph.input("i5_primer","well-specific GeoMx i5 Seq Code primer",
                   MaterialKind.OLIGO,inferred=True,
                   source=SOURCE+" (Seq Code bases proprietary)")
    i7=graph.input("i7_primer","well-specific GeoMx i7 Seq Code primer",
                   MaterialKind.OLIGO,inferred=True,
                   source=SOURCE+" (Seq Code bases proprietary)")
    polymerase=graph.input("pcr_polymerase",
                           "GeoMx 5X PCR Master Mix (composition not disclosed; Table 6 "
                           "program begins with 37 °C 30 min UDG incubation, 50 °C 10 min "
                           "UDG deactivation)",MaterialKind.CATALYST)
    amplified=graph.product("amplified",FINAL_LIBRARY.name,MaterialKind.LIBRARY,
                            state=FINAL_LIBRARY,inferred=True,
                            source=SOURCE+" (canonical TruSeq arms inferred)")
    graph.react("seqcode_pcr","UDG pre-incubation, then 18-cycle PCR with the AOI-specific unique-dual-index pair",
                StepKind.PCR,
                ((rehydrated,InputRole.TEMPLATE),(i5,InputRole.STRUCTURAL),
                 (i7,InputRole.STRUCTURAL),(polymerase,InputRole.CATALYST)),
                (amplified,),geometry=PCRGeometry.unavailable(
                    FINAL_LIBRARY,reason="The vendor publishes region order and lengths but not proprietary Seq Code or SPR bases."))
    beads=graph.input("ampure","AMPure XP cleanup beads",MaterialKind.SELECTION)
    library=graph.product("library","purified pooled GeoMx DSP library",
                          MaterialKind.LIBRARY,state=FINAL_LIBRARY,inferred=True,
                          source=SOURCE+" (canonical TruSeq arms inferred)")
    cleanup_waste=graph.product("cleanup_waste","discarded cleanup fraction",
                                MaterialKind.WASTE)
    graph.react("cleanup","Pool 4 µL per well and perform two successive 1.2× AMPure XP cleanups",StepKind.CLEANUP,
                ((amplified,InputRole.TEMPLATE),(beads,InputRole.SELECTION)),
                (library,cleanup_waste),changes_molecule=False)
    graph.discard(cleanup_waste,"the rejected bead-cleanup fraction is discarded")
    # Merritt 2020 published Online Methods, "Next-generation sequencing readout for
    # photocleaved oligonucleotide counting" (publisher PDF, read 2026-10-17): same PCR
    # (4 uL aspirate, 1 uM i5/i7, 1x NSTG master mix, 37 C 30 min, 50 C 10 min, 18
    # cycles, two 1.2x AMPure XP) but MiSeq 2x75 or NextSeq 2x37 runs.
    graph.annotate(SourceDiscrepancy(
        "library","Run length: defining article versus commercial manual",
        (SourceClaim("Merritt et al. 2020, Online Methods, 'Next-generation sequencing "
                     "readout for photocleaved oligonucleotide counting'",
                     "Libraries were paired-end sequenced on a MiSeq (2 × 75) or "
                     "NextSeq (2 × 37) instrument."),
         SourceClaim("MAN-10153-01 (April 2022), sequencing run settings",
                     "Paired-end Read 1 = 27, Index 1 = 8, Index 2 = 8, Read 2 = 27.")),
        "The page draws the commercial 2×27 run with two 8-cycle index reads, because "
        "the page describes the commercial Seq Code library prepared per MAN-10153-01.",
        "The article and manual agree on the PCR program and cleanups; only the run "
        "length differs. The 26 informative tag bases fit within all three read lengths.",
        "The article does not state its index-read lengths or whether its prototype tag "
        "arms equal the commercial SPR1/SPR2 lengths."))
    graph.sequence(library,SequencingPlan(
        library=FINAL_LIBRARY,primers=tuple(SEQ_PRIMERS),read_lengths=READ_LENGTHS))
    graph.declare_reviewed_scope((
        "tissue","probe","bind_probe","patterned_uv","photocleave",
        "drydown_rehydration","rehydrate","i5_primer",
        "i7_primer","pcr_polymerase","seqcode_pcr","ampure","cleanup","library"),SOURCE)
    return graph.require_publishable()
