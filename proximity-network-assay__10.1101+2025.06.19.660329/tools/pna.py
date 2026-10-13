"""Pixelgen Proximity Network Assay (PNA) / Proxiome kit v2 chemistry."""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))

from chemdraw import Construct, Segment, feature
from circular import CircularizationGeometry
from pcr import PCRGeometry
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind
from protocol_model import molecule, role
import seqprimers as sp
from proximity_topology import ProximityNetworkTopology


TITLE = "PNA — Pixelgen Proximity Network Assay"
NOTES = "01_pna.html"
PAPER = "https://doi.org/10.1101/2025.06.19.660329"
MANUAL = "https://www.pixelgen.com/user-manual-v2/"
SOURCE = (
    'Defining source: <a href="https://doi.org/10.1101/2025.06.19.660329">'
    'Karlsson et al., <i>bioRxiv</i> (2025; updated 2026)</a>; '
    '<a href="https://www.pixelgen.com/user-manual-v2/">Pixelgen Proxiome Kit v2 manual</a>; '
    '<a href="https://github.com/PixelgenTechnologies/pixelator">Pixelator assay definitions</a>; '
    '<a href="https://patents.google.com/patent/WO2025158200A1/en">related Pixelgen patent WO2025158200A1</a>.'
)
SUMMARY = (
    "Two oligo-conjugate types are placed on each antibody target. Padlock gap-fill and "
    "rolling-circle amplification turn every bound antibody into a repeated, addressable "
    "node; a linker spanning neighboring type-1 and type-2 RCA products is then filled and "
    "ligated into a paired UMI–protein-ID sequencing molecule."
)
CAVEAT = (
    "The paper and manual omit complete oligo molecules. Pixelator nevertheless publishes "
    "the Proxiome v2 internal amplicon schema, both fixed linker-binding sequences and PID "
    "panels, so those regions are drawn base by base. The padlocks, conjugate oligos and "
    "complete PCR/flow-cell ends remain unresolved; patent examples are not assumed to be "
    "identical to the current kit."
)

UMI1 = feature("pna.umi1", "umi", "random", group="pna-link", part="type-1")
UMI2 = feature("pna.umi2", "umi", "random", group="pna-link", part="type-2")
PID1 = feature("pna.pid1", "feature_barcode", "whitelist", group="pna-link",
               part="type-1", whitelist="Pixelgen Proxiome v2 panel")
PID2 = feature("pna.pid2", "feature_barcode", "whitelist", group="pna-link",
               part="type-2", whitelist="Pixelgen Proxiome v2 panel")
I7 = feature("pna.i7", "sample_index", "whitelist", group="pna-udi", part="i7",
             whitelist="Pixelgen UM002 Appendix 3")
I5 = feature("pna.i5", "sample_index", "whitelist", group="pna-udi", part="i5",
             whitelist="Pixelgen UM002 Appendix 3")
UEI = feature("pna.uei", "umi", "random", group="pna-link", part="event")


def seg(name, token, tag=None, *, molecular_feature=None):
    return Segment(name, token, tag, placeholder=True, feature=molecular_feature)


AMPLICON_SEGMENTS = (
    Segment("type-1 UMI", "N" * 28, "umi", feature=UMI1),
    Segment("type-1 protein ID", "N" * 10, "cbc", feature=PID1),
    Segment("type-1 linker-binding site", "CAAGTGACGCTGGTCAAACACTCATGTC", "ad"),
    Segment("unique event identifier", "N" * 15, "umi", feature=UEI),
    Segment("type-2 linker-binding site", "CGCTTAGATGTCGG", "ad"),
    Segment("type-2 protein ID", "N" * 10, "cbc", feature=PID2),
    Segment("type-2 UMI", "N" * 28, "umi", feature=UMI2),
)
INTERNAL_AMPLICON = Construct(AMPLICON_SEGMENTS, name="Proxiome v2 internal amplicon")
TYPE1_LINKER_SITE = next(s.top for s in AMPLICON_SEGMENTS
                         if s.name == "type-1 linker-binding site")
TYPE2_LINKER_SITE = next(s.top for s in AMPLICON_SEGMENTS
                         if s.name == "type-2 linker-binding site")
TOPOLOGY_PANELS = ProximityNetworkTopology(
    TYPE1_LINKER_SITE, TYPE2_LINKER_SITE, len(AMPLICON_SEGMENTS[3])).panels()

# The paper datasets were made with the earlier UM001 / Proxiome v1 chemistry. Pixelator's
# v1 definition and the released FASTQs agree on the longer fixed linker-binding regions.
PAPER_V1_AMPLICON = Construct((
    Segment("type-1 UMI", "N" * 28, "umi", feature=UMI1),
    Segment("type-1 protein ID", "N" * 10, "cbc", feature=PID1),
    Segment("type-1 linker-binding site", "CAAGTGACGCTGGGCATTCAAACACTCATGTC", "ad"),
    Segment("unique event identifier", "N" * 15, "umi", feature=UEI),
    Segment("type-2 linker-binding site", "GCTTCGCTTAGATGTCGG", "ad"),
    Segment("type-2 protein ID", "N" * 10, "cbc", feature=PID2),
    Segment("type-2 UMI", "N" * 28, "umi", feature=UMI2),
), name="paper-era Proxiome v1 internal amplicon")
PREAMP_LIBRARY = Construct([
    seg("left pre-amplification handle", "[left pre-amp handle]", "r1"),
    *AMPLICON_SEGMENTS,
    seg("right pre-amplification handle", "[right pre-amp handle]", "r2"),
], name="pre-amplified proximity link")

FINAL_LIBRARY = Construct([
    seg("P5-side flow-cell end", "[P5-side end]", "p5"),
    seg("8-nt i5 sample index", "N" * 8, "cbc", molecular_feature=I5),
    seg("TruSeq Read 1 handle", "[TruSeq R1 handle]", "r1"),
    *AMPLICON_SEGMENTS,
    seg("Nextera Read 2 handle", "[Nextera R2 handle]", "r2"),
    seg("8-nt i7 sample index", "N" * 8, "cbc", molecular_feature=I7),
    seg("P7-side flow-cell end", "[P7-side end]", "p7"),
], name="PNA proximity-link library")

SEQ_PRIMERS = ()
FINAL_CAPTION = (
    "One 133-nt internal amplicon records a link between a type-1 and a type-2 antibody "
    "node. Region order, lengths and both fixed linker-binding sequences come directly "
    "from Pixelator's Proxiome v2 assay definition. PCR/flow-cell ends remain role-labelled "
    "because their complete primer molecules are not disclosed there."
)
EXTRA_MOLECULAR_STATES = ((
    "Paper dataset: Proxiome v1 internal amplicon",
    PAPER_V1_AMPLICON,
    "The manuscript's released FASTQs are from UM001 / Proxiome v1. Pixelator defines "
    "this 141-nt internal amplicon, and the FASTQs confirm its orientation: R1 (44 nt) "
    "contains UMI-1(28) + PID-1(10) + 6 nt of LBS-1; R2 (78 nt) contains, in reverse-"
    "complement view, UMI-2(28) + PID-2(10) + LBS-2(18) + UEI(15) + 7 nt of LBS-1. "
    "The current v2 construct above shortens LBS-1 by 4 nt and LBS-2 by 4 nt.",
),)
SEQUENCING_ENDING = (
    "The v2 manual specifies paired-end Illumina-compatible sequencing: Read 1 44 cycles, "
    "Read 2 78 cycles, i7 8 cycles and i5 8 cycles. Both indexes are acquired, but the "
    "recommended default demultiplexes on i7 and masks i5."
)
SEQUENCING_UNAVAILABLE = (
    "Pixelator exposes the entire 133-nt internal amplicon and therefore its read content, "
    "but not the complete PCR primers or physical flow-cell constructs. Primer landing "
    "positions at the outer handles therefore remain role-level rather than base-level."
)


def _read_window(cycles: int, *, from_right: bool = False) -> str:
    """Describe a read directly from the machine-readable internal amplicon."""
    segments = tuple(reversed(AMPLICON_SEGMENTS)) if from_right else AMPLICON_SEGMENTS
    remaining = cycles
    fields = []
    for segment in segments:
        if remaining <= 0:
            break
        n = min(remaining, len(segment))
        if n == len(segment):
            fields.append(f"{n}-nt {segment.name}")
        else:
            edge = "last" if from_right else "first"
            fields.append(f"{edge} {n} nt of {segment.name}")
        remaining -= n
    if remaining:
        raise ValueError(f"{cycles}-cycle read exceeds the internal PNA amplicon")
    prefix = "reverse-complement view: " if from_right else ""
    return prefix + " + ".join(fields)


READ_SPECS = (
    sp.read_spec("Read 1", 44, _read_window(44), "include"),
    sp.read_spec("Index 1 (i7)", 8, "i7 sample index", "demultiplex"),
    sp.read_spec("Index 2 (i5)", 8, "i5 sample index; acquired but masked by default", "technical"),
    sp.read_spec("Read 2", 78, _read_window(78, from_right=True), "include"),
)


def _state(label, *parts):
    return molecule(label, tuple(role(x) for x in parts))


def protocol_graph():
    g = ProtocolGraph(TITLE, authority=SOURCE)
    cells = g.input("cells", "intact cells carrying surface-protein targets",
                    MaterialKind.BIOLOGICAL,
                    state=_state("cell surface", "surface-protein targets"), primary=True)
    pfa = g.input("pfa", "paraformaldehyde fixation condition", MaterialKind.CONDITION)
    fixed = g.product("fixed", "PFA-fixed cells", MaterialKind.MOLECULAR, state=cells.state)
    g.react("fix", "Fix cells with PFA", StepKind.CROSS_LINKING,
            ((cells, InputRole.TEMPLATE), (pfa, InputRole.CONDITION)), (fixed,),
            changes_molecule=False)

    hashtags = g.input("hashtags", "DNA-barcoded sample-hashtag antibodies", MaterialKind.SELECTION)
    hashed = g.product("hashed", "hashtag-labelled pooled cells", MaterialKind.MOLECULAR,
                       state=_state("sample-tagged cells", "sample hashtag", "surface-protein targets"))
    free_hashtags = g.product("free_hashtags", "unbound hashtag antibodies", MaterialKind.WASTE)
    g.react("hashtag", "Bind sample hashtags and pool samples", StepKind.BINDING,
            ((fixed, InputRole.TEMPLATE), (hashtags, InputRole.SELECTION)),
            (hashed, free_hashtags))
    g.discard(free_hashtags, "removed during the hashtag-binding washes")

    antibodies = g.input("antibodies", "type-1 and type-2 DNA-barcoded antibody panel",
                         MaterialKind.SELECTION)
    stained = g.product("stained", "cells bearing type-1 and type-2 antibody oligos",
                        MaterialKind.MOLECULAR,
                        state=_state("antibody-labelled cell", "type-1 PID + UMI",
                                     "surface protein", "type-2 PID + UMI"))
    free_antibodies = g.product("free_antibodies", "unbound barcoded antibodies",
                                MaterialKind.WASTE)
    g.react("stain", "Bind the barcoded antibody panel", StepKind.BINDING,
            ((hashed, InputRole.TEMPLATE), (antibodies, InputRole.SELECTION)),
            (stained, free_antibodies))
    g.discard(free_antibodies, "removed during antibody-panel washes")

    padlocks = g.input("padlocks", "type-1 and type-2 padlock probes", MaterialKind.OLIGO)
    fill_pol = g.input("fill_polymerase", "gap-fill polymerase", MaterialKind.CATALYST)
    fill_ligase = g.input("fill_ligase", "gap-fill ligase", MaterialKind.CATALYST)
    circles = g.product("circles", "cell-bound padlock circles carrying copied PID and UMI",
                        MaterialKind.MOLECULAR,
                        state=_state("two padlock classes", "type-1 PID + UMI circle",
                                     "type-2 PID + UMI circle"))
    failed_circles = g.product("failed_circles", "unclosed padlock products", MaterialKind.WASTE)
    g.react("close_padlocks", "Hybridize, gap-fill and close both padlock classes",
            StepKind.CIRCULARIZATION,
            ((stained, InputRole.TEMPLATE), (padlocks, InputRole.STRUCTURAL),
             (fill_pol, InputRole.CATALYST), (fill_ligase, InputRole.CATALYST)),
            (circles, failed_circles), geometry=CircularizationGeometry.known(
                "gap-sealing probe closure",
                statement="The antibody oligo templates PID and UMI copying across the padlock gap before ligase seals the circle."))
    g.discard(failed_circles, "does not form a rolling-circle template")

    rca_pol = g.input("rca_polymerase", "strand-displacing RCA polymerase", MaterialKind.CATALYST)
    rcps = g.product("rcps", "localized type-1 and type-2 rolling-circle products",
                    MaterialKind.MOLECULAR,
                    state=_state("neighboring spatial nodes", "type-1 PID + UMI repeats",
                                 "type-2 PID + UMI repeats"))
    g.react("rca", "Amplify each padlock circle in situ", StepKind.EXTENSION,
            ((circles, InputRole.TEMPLATE), (rca_pol, InputRole.CATALYST)), (rcps,))

    linkers = g.input("linkers", "proximity linker oligos", MaterialKind.OLIGO)
    link_pol = g.input("link_polymerase", "proximity-link gap-fill polymerase", MaterialKind.CATALYST)
    link_ligase = g.input("link_ligase", "proximity-link ligase", MaterialKind.CATALYST)
    links = g.product("links", "type-1/linker/type-2 proximity products",
                      MaterialKind.MOLECULAR, state=INTERNAL_AMPLICON)
    unlinked = g.product("unlinked", "unlinked or same-type RCA products", MaterialKind.WASTE)
    g.react("link", "Bridge neighboring type-1 and type-2 RCPs; gap-fill and ligate",
            StepKind.PROXIMITY_LIGATION,
            ((rcps, InputRole.TEMPLATE), (linkers, InputRole.STRUCTURAL),
             (link_pol, InputRole.CATALYST), (link_ligase, InputRole.CATALYST)),
            (links, unlinked))
    g.discard(unlinked, "does not encode a retained type-1/type-2 proximity link")

    digest = g.input("digest_enzyme", "cell/RCP digestion enzymes", MaterialKind.CATALYST)
    released = g.product("released", "released proximity-link templates",
                         MaterialKind.MOLECULAR, state=INTERNAL_AMPLICON)
    debris = g.product("debris", "digested cell and RCP material", MaterialKind.WASTE)
    g.react("digest", "Digest cells and RCP scaffolds to release link templates",
            StepKind.DIGESTION,
            ((links, InputRole.TEMPLATE), (digest, InputRole.CATALYST)), (released, debris))
    g.discard(debris, "does not carry an intact proximity-link amplicon")

    pre_f = g.input("pre_f", "pre-amplification forward primer", MaterialKind.OLIGO)
    pre_r = g.input("pre_r", "pre-amplification reverse primer", MaterialKind.OLIGO)
    pre_pol = g.input("pre_pol", "pre-amplification PCR polymerase", MaterialKind.CATALYST)
    preamp = g.product("preamp", "pre-amplified proximity links", MaterialKind.MOLECULAR,
                       state=PREAMP_LIBRARY)
    g.react("preamp_pcr", "Pre-amplify released proximity links", StepKind.PCR,
            ((released, InputRole.TEMPLATE), (pre_f, InputRole.STRUCTURAL),
             (pre_r, InputRole.STRUCTURAL), (pre_pol, InputRole.CATALYST)), (preamp,),
            geometry=PCRGeometry.unavailable(PREAMP_LIBRARY,
                reason="Pixelgen names the pre-amplification primer mix but does not publish its bases."))

    idx_f = g.input("index_f", "i5 indexing PCR primer", MaterialKind.OLIGO)
    idx_r = g.input("index_r", "i7 indexing PCR primer", MaterialKind.OLIGO)
    idx_pol = g.input("index_pol", "indexing PCR polymerase", MaterialKind.CATALYST)
    library = g.product("library", "dual-indexed PNA sequencing library",
                        MaterialKind.LIBRARY, state=FINAL_LIBRARY)
    g.react("index_pcr", "Add unique dual indexes by PCR", StepKind.PCR,
            ((preamp, InputRole.TEMPLATE), (idx_f, InputRole.STRUCTURAL),
             (idx_r, InputRole.STRUCTURAL), (idx_pol, InputRole.CATALYST)), (library,),
            geometry=PCRGeometry.unavailable(FINAL_LIBRARY,
                reason="The v2 manual prints index values but not the complete primer molecules or binding sites."))
    g.sequence(library, SequencingPlan(
        acquired_roles=("Read 1", "Index 1 (i7)", "Index 2 (i5)", "Read 2"),
        read_lengths={"Read 1": 44, "Index 1 (i7)": 8, "Index 2 (i5)": 8, "Read 2": 78},
        read_uses={"Index 1 (i7)": "demultiplex", "Index 2 (i5)": "technical"},
        limitation=SEQUENCING_UNAVAILABLE))

    reviewed = [
        "cells", "pfa", "fix", "hashtags", "hashtag", "antibodies", "stain",
        "padlocks", "fill_polymerase", "fill_ligase", "close_padlocks",
        "rca_polymerase", "rca", "linkers", "link_polymerase", "link_ligase", "link",
        "digest_enzyme", "digest", "pre_f", "pre_r", "pre_pol", "preamp_pcr",
        "index_f", "index_r", "index_pol", "index_pcr", "library",
    ]
    g.declare_reviewed_scope(reviewed, SOURCE)
    return g.require_publishable()
