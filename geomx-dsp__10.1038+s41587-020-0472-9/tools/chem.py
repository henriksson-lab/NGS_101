"""GeoMx DSP NGS readout page model."""
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))

import digital_spatial_tags as dsp
from pcr import PCRGeometry
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind)


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
                    "Index 2 orientation depends on the Illumina instrument workflow.")


def oligos():
    return ()


def workflow():
    return dsp.workflow()


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
    i5=graph.input("i5_primer","well-specific GeoMx i5 Seq Code primer",
                   MaterialKind.OLIGO,inferred=True,
                   source=SOURCE+" (Seq Code bases proprietary)")
    i7=graph.input("i7_primer","well-specific GeoMx i7 Seq Code primer",
                   MaterialKind.OLIGO,inferred=True,
                   source=SOURCE+" (Seq Code bases proprietary)")
    polymerase=graph.input("pcr_polymerase","GeoMx library PCR polymerase",
                           MaterialKind.CATALYST)
    amplified=graph.product("amplified",FINAL_LIBRARY.name,MaterialKind.LIBRARY,
                            state=FINAL_LIBRARY,inferred=True,
                            source=SOURCE+" (canonical TruSeq arms inferred)")
    graph.react("seqcode_pcr","Amplify with the AOI-specific unique-dual-index pair",
                StepKind.PCR,
                ((tags,InputRole.TEMPLATE),(i5,InputRole.STRUCTURAL),
                 (i7,InputRole.STRUCTURAL),(polymerase,InputRole.CATALYST)),
                (amplified,),geometry=PCRGeometry.unavailable(
                    FINAL_LIBRARY,reason="The vendor publishes region order and lengths but not proprietary Seq Code or SPR bases."))
    beads=graph.input("ampure","AMPure XP cleanup beads",MaterialKind.SELECTION)
    library=graph.product("library","purified pooled GeoMx DSP library",
                          MaterialKind.LIBRARY,state=FINAL_LIBRARY,inferred=True,
                          source=SOURCE+" (canonical TruSeq arms inferred)")
    cleanup_waste=graph.product("cleanup_waste","discarded cleanup fraction",
                                MaterialKind.WASTE)
    graph.react("cleanup","Perform two successive AMPure XP cleanups",StepKind.CLEANUP,
                ((amplified,InputRole.TEMPLATE),(beads,InputRole.SELECTION)),
                (library,cleanup_waste),changes_molecule=False)
    graph.discard(cleanup_waste,"the rejected bead-cleanup fraction is discarded")
    graph.sequence(library,SequencingPlan(
        library=FINAL_LIBRARY,primers=tuple(SEQ_PRIMERS),read_lengths=READ_LENGTHS))
    graph.declare_reviewed_scope((
        "tissue","probe","bind_probe","patterned_uv","photocleave","i5_primer",
        "i7_primer","pcr_polymerase","seqcode_pcr","ampure","cleanup","library"),SOURCE)
    return graph.require_publishable()
