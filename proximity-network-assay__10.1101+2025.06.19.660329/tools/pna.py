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


TITLE = "PNA — Pixelgen Proximity Network Assay"
NOTES = "01_pna.html"
PAPER = "https://doi.org/10.1101/2025.06.19.660329"
MANUAL = "https://www.pixelgen.com/user-manual-v2/"
SOURCE = (
    'Defining source: <a href="https://doi.org/10.1101/2025.06.19.660329">'
    'Karlsson et al., <i>bioRxiv</i> (2025; updated 2026)</a>; '
    '<a href="https://www.pixelgen.com/user-manual-v2/">Pixelgen Proxiome Kit v2 manual</a>.'
)
SUMMARY = (
    "Two oligo-conjugate types are placed on each antibody target. Padlock gap-fill and "
    "rolling-circle amplification turn every bound antibody into a repeated, addressable "
    "node; a linker spanning neighboring type-1 and type-2 RCA products is then filled and "
    "ligated into a paired UMI–protein-ID sequencing molecule."
)
CAVEAT = (
    "Pixelgen publishes the molecular roles and the complete reaction order, but not the "
    "kit oligo bases or exact primer-binding boundaries. The diagrams therefore preserve "
    "the published topology and read lengths while showing proprietary regions as named "
    "placeholders."
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


def seg(name, token, tag=None, *, molecular_feature=None):
    return Segment(name, token, tag, placeholder=True, feature=molecular_feature)


FINAL_LIBRARY = Construct([
    seg("P5-side flow-cell end", "[P5-side end]", "p5"),
    seg("8-nt i5 sample index", "N" * 8, "cbc", molecular_feature=I5),
    seg("TruSeq Read 1 handle", "[TruSeq R1 handle]", "r1"),
    seg("type-1 UMI", "[UMI-1]", "umi", molecular_feature=UMI1),
    seg("type-1 protein ID", "[PID-1]", "cbc", molecular_feature=PID1),
    seg("type-1 linker-binding region", "[LBS-1]"),
    seg("proximity linker", "[linker]", "ad"),
    seg("type-2 linker-binding region", "[LBS-2]"),
    seg("type-2 protein ID", "[PID-2]", "cbc", molecular_feature=PID2),
    seg("type-2 UMI", "[UMI-2]", "umi", molecular_feature=UMI2),
    seg("Nextera Read 2 handle", "[Nextera R2 handle]", "r2"),
    seg("8-nt i7 sample index", "N" * 8, "cbc", molecular_feature=I7),
    seg("P7-side flow-cell end", "[P7-side end]", "p7"),
], name="PNA proximity-link library")

SEQ_PRIMERS = ()
FINAL_CAPTION = (
    "One amplicon records a link between a type-1 and a type-2 antibody node. The order "
    "of UMI-1, PID-1, linker, PID-2 and UMI-2 is published; kit bases and exact boundaries "
    "are proprietary."
)
SEQUENCING_ENDING = (
    "The v2 manual specifies paired-end Illumina-compatible sequencing: Read 1 44 cycles, "
    "Read 2 78 cycles, i7 8 cycles and i5 8 cycles. Both indexes are acquired, but the "
    "recommended default demultiplexes on i7 and masks i5."
)
SEQUENCING_UNAVAILABLE = (
    "The exact Proxiome v2 sequencing-primer sites and proprietary linker-library bases "
    "are not published, so primer landing positions cannot be computed base by base."
)
READ_SPECS = (
    sp.read_spec("Read 1", 44, "type-1 UMI and protein ID, followed by linker-side sequence", "include"),
    sp.read_spec("Index 1 (i7)", 8, "i7 sample index", "demultiplex"),
    sp.read_spec("Index 2 (i5)", 8, "i5 sample index; acquired but masked by default", "technical"),
    sp.read_spec("Read 2", 78, "type-2 UMI and protein ID, followed by linker-side sequence", "include"),
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
                      MaterialKind.MOLECULAR, state=FINAL_LIBRARY)
    unlinked = g.product("unlinked", "unlinked or same-type RCA products", MaterialKind.WASTE)
    g.react("link", "Bridge neighboring type-1 and type-2 RCPs; gap-fill and ligate",
            StepKind.PROXIMITY_LIGATION,
            ((rcps, InputRole.TEMPLATE), (linkers, InputRole.STRUCTURAL),
             (link_pol, InputRole.CATALYST), (link_ligase, InputRole.CATALYST)),
            (links, unlinked))
    g.discard(unlinked, "does not encode a retained type-1/type-2 proximity link")

    digest = g.input("digest_enzyme", "cell/RCP digestion enzymes", MaterialKind.CATALYST)
    released = g.product("released", "released proximity-link templates",
                         MaterialKind.MOLECULAR, state=FINAL_LIBRARY)
    debris = g.product("debris", "digested cell and RCP material", MaterialKind.WASTE)
    g.react("digest", "Digest cells and RCP scaffolds to release link templates",
            StepKind.DIGESTION,
            ((links, InputRole.TEMPLATE), (digest, InputRole.CATALYST)), (released, debris))
    g.discard(debris, "does not carry an intact proximity-link amplicon")

    pre_f = g.input("pre_f", "pre-amplification forward primer", MaterialKind.OLIGO)
    pre_r = g.input("pre_r", "pre-amplification reverse primer", MaterialKind.OLIGO)
    pre_pol = g.input("pre_pol", "pre-amplification PCR polymerase", MaterialKind.CATALYST)
    preamp = g.product("preamp", "pre-amplified proximity links", MaterialKind.MOLECULAR,
                       state=FINAL_LIBRARY)
    g.react("preamp_pcr", "Pre-amplify released proximity links", StepKind.PCR,
            ((released, InputRole.TEMPLATE), (pre_f, InputRole.STRUCTURAL),
             (pre_r, InputRole.STRUCTURAL), (pre_pol, InputRole.CATALYST)), (preamp,),
            geometry=PCRGeometry.unavailable(FINAL_LIBRARY,
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
