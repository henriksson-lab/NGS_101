"""Construct model for the individual dscATAC-seq and dsciATAC-seq pages."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"lib"))
import illumina as il
import nextera as nx
from chemdraw import Construct, Scene, Segment, feature
from pcr import PCRGeometry
from protocol_graph import (Disposition, InputRole, MaterialKind, ProtocolGraph,
                            SequencingPlan, StepKind, Terminal)

C1="TATGCATGAC"; C2="AGTCACTGAG"; BC_NT=7; TN5_BC_NT=6
REP_TN5_BC="AAAGAA"

def _s(n,t,g=None,**kw): return Segment(n,t,g,**kw)

def bead_published_segments():
    return [_s("BC1","N"*BC_NT,"cbc",placeholder=True,feature=feature("cell_bc1","cell_barcode","whitelist",whitelist="published bead barcode set",group="cell_id",part="BC1")),
            _s("phase block","NN",placeholder=True),_s("constant 1",C1),
            _s("BC2","N"*BC_NT,"cbc",placeholder=True,feature=feature("cell_bc2","cell_barcode","whitelist",whitelist="published bead barcode set",group="cell_id",part="BC2")),_s("constant 2",C2),
            _s("BC3","N"*BC_NT,"cbc",placeholder=True,feature=feature("cell_bc3","cell_barcode","whitelist",whitelist="published bead barcode set",group="cell_id",part="BC3")),_s("s5 primer",nx.S5,"s5")]

def tn5_r1(protocol):
    if protocol=="dscATAC-seq": return [_s("s5",nx.S5,"s5"),_s("ME",nx.ME,"me")]
    if protocol=="dsciATAC-seq": return [_s("s5",nx.S5,"s5"),
        _s("Tn5 barcode",REP_TN5_BC,"cbc",feature=feature("cell_tn5","cell_barcode","whitelist",whitelist="published Tn5 barcode set",group="cell_id",part="Tn5")),_s("ME",nx.ME,"me")]
    raise ValueError(protocol)

def bead_priming_scene():
    template=[_s("s5'",nx.S5_RC,"s5"),_s("ME'",nx.ME_RC,"me"),
              _s("insert","XXXXXXXX...",placeholder=True)]
    bead=[_s("unpublished bead 5' arm","XXXXXXXXXXXX",inferred=True,placeholder=True),
          *bead_published_segments()]
    sc=Scene(); sc.strand("tagmented strand",template,label="fragment")
    sc.anneal("released bead oligo",bead,to="tagmented strand",
              pair=("s5 primer","s5'"),label="bead primer",
              unpaired=("unpublished bead 5' arm","BC1","phase block","constant 1",
                        "BC2","constant 2","BC3"))
    sc.arrow("released bead oligo","polymerase")
    return sc

def filled_fragment(protocol):
    return Construct([*tn5_r1(protocol),_s("genomic insert","XXXXXXXX...XXXXXXXX",placeholder=True),
                      _s("ME'",nx.ME_RC,"me"),_s("s7'",nx.S7_RC,"s7")],
                     name=f"{protocol} gap-filled fragment")

def final_library(protocol):
    inner=[*bead_published_segments()]
    if protocol=="dsciATAC-seq": inner.append(_s("Tn5 barcode", "N"*TN5_BC_NT,"cbc",placeholder=True,feature=feature("cell_tn5","cell_barcode","whitelist",whitelist="published Tn5 barcode set",group="cell_id",part="Tn5")))
    elif protocol!="dscATAC-seq": raise ValueError(protocol)
    inner += [_s("ME",nx.ME,"me"),_s("genomic insert","XXXXXXXX...XXXXXXXX",placeholder=True),
              _s("ME'",nx.ME_RC,"me"),_s("s7'",nx.S7_RC,"s7"),
              _s("i7'","N"*8,"cbc",placeholder=True,feature=feature("sample_i7","sample_index","unknown")),_s("P7'",il.P7_RC,"p7")]
    return Construct([_s("unpublished bead/P5 arm","[BEAD RELEASE / P5 / R1 ARM]",inferred=True,
                         placeholder=True),*inner],name=f"{protocol} final library")

def _validate():
    if "".join(s.top for s in tn5_r1("dsciATAC-seq")) != nx.S5+REP_TN5_BC+nx.ME:
        raise ValueError("barcoded Tn5 adaptor no longer matches Supplementary Table 6")
    bead_priming_scene().rows()
_validate()

SOURCE = "Defining source: https://doi.org/10.1038/s41587-019-0147-6"
METHODS = SOURCE + " (author manuscript PMC10299900, Methods 'Droplet library preparation and sequencing')"
SUPP_NOTE = SOURCE + " (Supplementary Note, MOESM1: bead oligos are cleaved from the bead and prime by PCR)"
TABLE6 = SOURCE + " (Supplementary Table 6, MOESM8)"

# Graph targets shared by the dscATAC and dsciATAC assemblies (from partition on).
DROPLET_SCOPE = (
    "beads", "droplet_loading", "partition", "bead_release_condition", "bead_release",
    "pcr_polymerase", "gap_fill", "bead_oligo", "index_primer", "droplet_pcr",
    "break_emulsion", "ampure_1", "cleanup_beads", "second_pcr_forward",
    "second_pcr_reverse", "second_pcr", "ampure_2", "library")


def protocol_graph(protocol):
    """The defining dscATAC-seq graph; dsciATAC owns its assembly locally."""
    if protocol != "dscATAC-seq":
        raise ValueError("dsciATAC-seq assembles its graph in its own directory")
    graph = ProtocolGraph(protocol, authority=SOURCE)
    nuclei = graph.input("nuclei", "intact cells or nuclei", MaterialKind.BIOLOGICAL,
                         primary=True)
    tagged = tagment_cells(graph, nuclei, protocol=protocol)
    droplet_library(graph, tagged, protocol=protocol)
    graph.declare_reviewed_scope(("nuclei", "transposome", "tn5", "tagment",
                                  *DROPLET_SCOPE), SOURCE)
    return graph.require_publishable()


def tagment_cells(graph: ProtocolGraph, nuclei, *, protocol: str):
    """Bulk tagmentation with the kit Tn5 (37 C 30 min, ThermoMixer); cells stay intact."""
    transposome = graph.input("transposome", "Tn5 transposome (kit ATAC Tagmentation Enzyme)",
                              MaterialKind.ADAPTER, state=tuple(tn5_r1(protocol)), inferred=True,
                              source=METHODS + "; kit enzyme adaptor bases are not printed, the Nextera s5-ME "
                                     "arm is inferred from the Table 6 bead 3' s5 primer")
    tn5 = graph.input("tn5", "Tn5 transposase", MaterialKind.CATALYST)
    tagged = graph.product("tagmented_cells", "intact cells containing tagmented chromatin",
                           MaterialKind.MOLECULAR)
    graph.react("tagment", "Tagment accessible chromatin in bulk (37 C 30 min) while retaining intact cells",
                StepKind.TAGMENTATION,
                ((nuclei, InputRole.TEMPLATE), (transposome, InputRole.STRUCTURAL),
                 (tn5, InputRole.CATALYST)), (tagged,), source=METHODS)
    return tagged


def droplet_library(graph: ProtocolGraph, tagged, *, protocol: str):
    """Partition, in-droplet programme, emulsion break, cleanups and second PCR.

    Order follows the Methods PCR programme: 37 C 30 min, 85 C 10 min, 72 C 5 min,
    98 C 30 s, 8 cycles, 72 C 5 min; then emulsion break, AMPure XP, 7-9 cycle second
    PCR, AMPure XP.  Bead-release chemistry, the 85 C purpose and the 72 C purpose are
    not stated in the article; they are kept as a warned boundary and an inferred step.
    """
    sci = protocol == "dsciATAC-seq"
    beads = graph.input("beads", "super-loaded barcode beads (published 3' bead-oligo portion; 5' arm unpublished)",
                        MaterialKind.OLIGO, state=tuple(bead_published_segments()), inferred=True,
                        source=TABLE6 + "; bead linker cell is empty")
    loading = graph.input("droplet_loading",
                          "ddSEQ Single-Cell Isolator co-encapsulation with PCR mix; beads super-loaded",
                          MaterialKind.CONDITION, source=METHODS)
    cell_label = ("droplets with several Tn5-barcoded cells and one or more beads" if sci else
                  "droplets with a tagmented cell and one or more beads")
    droplets = graph.product("cell_droplets", cell_label, MaterialKind.MOLECULAR)
    cell_free = graph.product("cell_free_droplets", "droplets without a cell (beads and PCR mix only)",
                              MaterialKind.MOLECULAR)
    graph.react("partition",
                ("Co-encapsulate pooled Tn5-barcoded cells above Poisson loading with super-loaded beads"
                 if sci else "Co-encapsulate tagmented cells with super-loaded barcode beads"),
                StepKind.PARTITION,
                ((tagged, InputRole.TEMPLATE), (beads, InputRole.STRUCTURAL),
                 (loading, InputRole.CONDITION)), (droplets, cell_free), source=METHODS)
    release_condition = graph.input("bead_release_condition", "37 C 30 min then 85 C 10 min in droplets",
                                    MaterialKind.CONDITION, source=METHODS)
    released = graph.product("released_droplets", "droplets with bead oligos released from the beads",
                             MaterialKind.MOLECULAR, source=SUPP_NOTE)
    graph.react("bead_release", "Release bead oligos into droplets", StepKind.MAGIC,
                ((droplets, InputRole.TEMPLATE), (release_condition, InputRole.CONDITION)),
                (released,), source=METHODS,
                reason="The Supplementary Note says bead oligos are cleaved from the physical bead; "
                       "the article gives 37 C 30 min and 85 C 10 min first but names neither the "
                       "cleavage chemistry nor the purpose of 85 C.")
    polymerase = graph.input("pcr_polymerase", "kit droplet PCR polymerase (not named in the article)",
                             MaterialKind.CATALYST)
    filled = graph.product("gap_filled_droplets", "droplets with gap-filled tagmented fragments",
                           MaterialKind.MOLECULAR, state=filled_fragment(protocol), inferred=True,
                           source=METHODS + "; 72 C 5 min before denaturation modelled as Tn5 gap fill (purpose not stated)")
    graph.react("gap_fill", "INFERRED: fill Tn5 9-nt gaps at 72 C 5 min", StepKind.EXTENSION,
                ((released, InputRole.TEMPLATE), (polymerase, InputRole.CATALYST)), (filled,),
                source=METHODS)
    bead = graph.input("bead_oligo", "released bead oligo (published 3' portion; 5' arm unpublished)",
                       MaterialKind.OLIGO, state=tuple(bead_published_segments()), inferred=True,
                       source=SUPP_NOTE + "; 5-prime bead/P5 arm is unpublished")
    index = graph.input("index_primer",
                        "N7xx sample-index primer (P7-i7-s7; Table 6 lists N701-N729, kit 12009360 identity not stated)",
                        MaterialKind.OLIGO,
                        state=(_s("P7",il.P7,"p7"),
                               _s("i7","N"*8,"cbc",placeholder=True,feature=feature("sample_i7","sample_index","unknown")),
                               _s("s7",nx.S7,"s7")),
                        inferred=True,
                        source=TABLE6 + " N70X sheet; the Methods name only SureCell ATAC-Seq Index Kit 12009360")
    final = final_library(protocol)
    amplicons = graph.product("droplet_amplicons", "bead-barcoded s5/s7 amplicons in droplets",
                              MaterialKind.MOLECULAR, inferred=True,
                              source=SUPP_NOTE + "; bead/P5 arm unresolved")
    same_end = graph.product("same_end_amplicons",
                             "s5/s5 (bead primer only) and s7/s7 (index primer only) products",
                             MaterialKind.MOLECULAR, inferred=True,
                             source="Derived from Tn5 end combinations; not reported in the article")
    graph.react("droplet_pcr", "Bead barcoding and sample indexing PCR in droplets (8 cycles)",
                StepKind.PCR,
                ((filled, InputRole.TEMPLATE), (bead, InputRole.STRUCTURAL),
                 (index, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)),
                (amplicons, same_end), source=METHODS, geometry=PCRGeometry.unavailable(
                    final, reason="The custom Read-1 and library-primer bases are unresolved from the checked article, Supplementary Table 6 and Bio-Rad guide."))
    graph.mark_terminal(Terminal("same_end_amplicons", Disposition.UNRESOLVED,
        "Candidate products lacking one of P5/P7; no represented separation removes them, "
        "and they are not expected to form clusters."))
    broken = graph.product("pooled_products", "pooled products after breaking the emulsion",
                           MaterialKind.MOLECULAR)
    graph.react("break_emulsion", "Break the emulsion, pooling all droplets", StepKind.POOL,
                ((amplicons, InputRole.TEMPLATE), (cell_free, InputRole.TEMPLATE)), (broken,),
                source=METHODS)
    ampure = graph.input("cleanup_beads", "AMPure XP beads", MaterialKind.SELECTION, source=METHODS)
    clean1 = graph.product("cleaned_amplicons", "AMPure-purified barcoded amplicons", MaterialKind.MOLECULAR)
    waste1 = graph.product("ampure_1_waste", "material not retained by AMPure XP", MaterialKind.WASTE)
    graph.react("ampure_1", "AMPure XP cleanup", StepKind.CLEANUP,
                ((broken, InputRole.TEMPLATE), (ampure, InputRole.SELECTION)), (clean1, waste1),
                source=METHODS)
    graph.discard(waste1, "AMPure XP supernatant is removed")
    fwd = graph.input("second_pcr_forward", "second-PCR forward primer (kit; not named in the article)",
                      MaterialKind.OLIGO, inferred=True, source=METHODS + "; primers not named")
    rev = graph.input("second_pcr_reverse", "second-PCR reverse primer (kit; not named in the article)",
                      MaterialKind.OLIGO, inferred=True, source=METHODS + "; primers not named")
    amplified = graph.product("amplified_library", "re-amplified library", MaterialKind.MOLECULAR,
                              inferred=True, source=METHODS + "; bead/P5 arm unresolved")
    graph.react("second_pcr", "Re-amplify barcoded amplicons (7-9 cycles by cell input)", StepKind.PCR,
                ((clean1, InputRole.TEMPLATE), (fwd, InputRole.STRUCTURAL),
                 (rev, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)),
                (amplified,), source=METHODS, geometry=PCRGeometry.unavailable(
                    final, reason="The article does not name the second-PCR primers; it refers to the Bio-Rad User Guide."))
    library = graph.product("library", final.name, MaterialKind.LIBRARY, inferred=True,
                            source=SOURCE + " (bead/P5 arm unresolved from checked sources)")
    waste2 = graph.product("ampure_2_waste", "material not retained by AMPure XP", MaterialKind.WASTE)
    graph.react("ampure_2", "AMPure XP cleanup; Bioanalyzer HS quantification", StepKind.CLEANUP,
                ((amplified, InputRole.TEMPLATE), (ampure, InputRole.SELECTION)), (library, waste2),
                source=METHODS)
    graph.discard(waste2, "AMPure XP supernatant is removed")
    graph.sequence(library, SequencingPlan(
        limitation="The article reports read roles and cycles; the custom Read-1 and complete bead/P5 primer bases remain unresolved from the checked article, Supplementary Table 6 and Bio-Rad guide.",
        acquired_roles=("Read 1", "Index 1 (i7)", "Read 2"),
        read_lengths={"Read 1": 118, "Index 1 (i7)": 8, "Read 2": 40},
        read_uses={"Read 1": "extract_then_include", "Index 1 (i7)": "demultiplex",
                   "Read 2": "include"}))
    return library
