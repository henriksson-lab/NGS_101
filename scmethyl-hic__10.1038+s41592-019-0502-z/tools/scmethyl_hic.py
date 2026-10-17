"""scMethyl-HiC single-nucleus branch of Li et al., Nature Methods 2019."""
from __future__ import annotations

import html
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "lib"))
from storage import data_root

import seqprimers as sp
import illumina as il
from adapter_ligation import AdapterLigationGeometry
from chemdraw import Construct, Scene, Segment, duplex_rows, feature, junction_row, panel, protocol_graph_panels, revcomp
from checks import Check, Source, have
from page import head, info
from pcr import PCRGeometry
from protocol_graph import Disposition, InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind, Terminal
from protocol_model import molecule, role
from restriction import DPNII
from single_cell_hic import ContactWorkflow

PAGE_FILE = "scmethyl-hic.html"
SOURCE = "https://doi.org/10.1038/s41592-019-0502-z"
SUPPLEMENT = "Li et al. 2019 Supplementary Protocol, MOESM1, pp. 5 and 15–18"
TITLE = "scMethyl-HiC — single-nucleus contacts and DNA methylation"

# Verbatim oligos from the original Supplementary Protocol, extracted lines 213–258.
# IDT's (N1:25252525) followed by eight (N1) positions is displayed as N9 below;
# the unusual synthesis notation and /5SpC3/ remain visible in the page and note.
P5_RANDOM_PREFIX = "TTCCCTACACGACGCTCTTCCGATCT"
INLINE_4PLEX = {"AD002": "CGATGT", "AD006": "GCCAAT",
                "AD008": "ACTTGA", "AD010": "TAGCTT"}
P5_PCR_TAIL = "ACACTCTTTCCCTACACGACGCTCT"
P7_PCR_TAIL = "GTGACTGGAGTTCAGACGTGTGCTCTT"
P5_INDEXES = {"D501": "TATAGCCT", "D502": "ATAGAGGC", "D503": "CCTATCCT",
              "D504": "GGCTCTGA", "D505": "AGGCGAAG", "D506": "TAATCTTA",
              "D507": "CAGGACGT", "D508": "GTACTGAC"}
P7_INDEXES = {"D701": "CGAGTAAT", "D702": "TCTCCGGA", "D703": "AATGAGCG",
              "D704": "GGAATCTC", "D705": "TTCTGAAT", "D706": "ACGAATTC",
              "D707": "AGCTTCAG", "D708": "GCGCATTA", "D709": "CATAGCCG",
              "D710": "TTCGCGGA", "D711": "GCGCGAGA", "D712": "CTATCGCT"}
P5_PCR = {name: il.P5 + index + P5_PCR_TAIL for name, index in P5_INDEXES.items()}
P7_PCR = {name: il.P7 + index + P7_PCR_TAIL for name, index in P7_INDEXES.items()}
P5_D501 = P5_PCR["D501"]
P7_D701 = P7_PCR["D701"]
WORKFLOW = ContactWorkflow("scMethyl-HiC", DPNII, "single nucleus after proximity ligation",
                           "biotin-fill", "bisulfite-PCR", capture=False)

POOL_FEATURE = feature("four_plex_inline_index", "sample_index", "whitelist",
                       whitelist="Li 2019 Supplementary Protocol AD002/006/008/010",
                       group="single_nucleus_pool")
I5_FEATURE = feature("p5_pcr_index", "sample_index", "unknown")
I7_FEATURE = feature("p7_pcr_index", "sample_index", "unknown")


def contact_product() -> Construct:
    return WORKFLOW.junction()


def random_primer(member: str = "AD002") -> Construct:
    if member not in INLINE_4PLEX:
        raise ValueError("random-primer member is outside the published four-plex set")
    return Construct((
        Segment("P5L 5′ handle", P5_RANDOM_PREFIX, tag="r1"),
        Segment("inline pool index", INLINE_4PLEX[member], tag="cbc", feature=POOL_FEATURE),
        Segment("random nonamer", "N" * 9, placeholder=True),
    ), name=f"/5SpC3/ P5L_{member} random primer")


def final_library() -> Construct:
    """Representative published D501/D701/AD002 combination, not a claimed well assignment."""
    if not P5_D501.startswith(il.P5) or not P7_D701.startswith(il.P7):
        raise ValueError("published P5/P7 PCR primers disagree with canonical flow-cell arms")
    p5_index_start = len(il.P5)
    p5_index_end = p5_index_start + 8
    p7_index_start = len(il.P7)
    p7_index_end = p7_index_start + 8
    overlap = next(n for n in range(len(P5_RANDOM_PREFIX), 0, -1)
                   if P5_D501.endswith(P5_RANDOM_PREFIX[:n]))
    return Construct((
        Segment("P5 flow-cell arm", il.P5, tag="p5"),
        Segment("P5 PCR index, D501 example", P5_D501[p5_index_start:p5_index_end],
                tag="cbc", feature=I5_FEATURE),
        Segment("P5 PCR annealing arm", P5_D501[p5_index_end:], tag="r1"),
        Segment("copied random-primer handle", P5_RANDOM_PREFIX[overlap:], tag="r1"),
        Segment("inline pool index, AD002 example", INLINE_4PLEX["AD002"], tag="cbc",
                feature=POOL_FEATURE),
        Segment("random-primed tract", "N" * 9, placeholder=True),
        Segment("bisulfite-converted contact insert", "X" * 48, placeholder=True),
        Segment("Adaptase-added P7-side junction (length and bases unresolved)", "...",
                placeholder=True, inferred=True),
        Segment("P7 PCR annealing arm reverse complement", revcomp(P7_D701[p7_index_end:]),
                tag="r2"),
        Segment("P7 PCR index, D701 example, reverse complement",
                revcomp(P7_D701[p7_index_start:p7_index_end]), tag="cbc", feature=I7_FEATURE),
        Segment("P7 flow-cell arm reverse complement", revcomp(il.P7), tag="p7"),
    ), name="scMethyl-HiC dual-index library, representative published primers")


def protocol_graph() -> ProtocolGraph:
    g = ProtocolGraph(TITLE, authority=SOURCE)
    sample = g.input("sample", "crosslinked nuclei", MaterialKind.BIOLOGICAL,
                     state=molecule("nuclear chromatin", (role("chromosomal DNA"),)), primary=True)
    enzyme = g.input("dpnii", "DpnII", MaterialKind.CATALYST)
    cut = g.product("cut", "DpnII-digested chromatin", MaterialKind.MOLECULAR,
                    state=sample.state)
    g.react("digest", "Digest chromatin with DpnII", StepKind.DIGESTION,
            ((sample, InputRole.TEMPLATE), (enzyme, InputRole.CATALYST)), (cut,))
    biotin = g.input("biotin_dntp", "biotin-14-dATP fill-in nucleotide", MaterialKind.MOLECULAR)
    fill = g.input("fill_polymerase", "Klenow DNA polymerase I large fragment", MaterialKind.CATALYST)
    marked = g.product("marked", "biotin-filled restriction ends", MaterialKind.MOLECULAR,
                       state=cut.state)
    g.react("fill", "Fill restriction ends with biotin", StepKind.EXTENSION,
            ((cut, InputRole.TEMPLATE), (biotin, InputRole.STRUCTURAL),
             (fill, InputRole.CATALYST)), (marked,))
    ligase = g.input("contact_ligase", "T4 DNA ligase", MaterialKind.CATALYST)
    contact = g.product("contact", "proximity-ligated contacts", MaterialKind.MOLECULAR,
                        state=contact_product())
    g.react("proximity_ligate", "Join spatially adjacent restriction fragments",
            StepKind.PROXIMITY_LIGATION,
            ((marked, InputRole.TEMPLATE), (ligase, InputRole.CATALYST)), (contact,))

    sorter = g.input("sorter", "DRAQ7-stained single-nucleus FACS sorting", MaterialKind.SELECTION)
    nucleus = g.product("nucleus", "one sorted contact-bearing nucleus", MaterialKind.MOLECULAR,
                        state=contact.state)
    excluded = g.product("other_nuclei", "other nuclei outside the traced well", MaterialKind.MOLECULAR)
    g.react("sort", "Sort one nucleus into each bisulfite well", StepKind.SAMPLING,
            ((contact, InputRole.TEMPLATE), (sorter, InputRole.SELECTION)),
            (nucleus, excluded), changes_molecule=False)
    g.mark_terminal(Terminal(excluded.id, Disposition.UNRESOLVED, "Other nuclei are outside this traced well; their later disposition is not represented."))
    bisulfite = g.input("bisulfite", "EZ-96 DNA Methylation-Direct conversion",
                        MaterialKind.CONDITION,
                        source="Li 2019 article Methods (D5020); Supplementary Protocol reagent list (D5023)")
    lambda_control = g.input("lambda_control", "0.5% sonicated unmethylated lambda DNA control", MaterialKind.MOLECULAR)
    lysis = g.input("lysis", "Zymo lysis buffer containing proteinase K (supplementary route)", MaterialKind.CONDITION)
    converted = g.product("converted", "bisulfite-converted contact DNA", MaterialKind.MOLECULAR,
                          state=molecule("converted single-nucleus contact DNA",
                                         (role("C-to-T-converted contact insert"),)))
    g.react("bisulfite_convert", "Convert cytosines in the sorted nucleus",
            StepKind.CHEMICAL_CONVERSION,
            ((nucleus, InputRole.TEMPLATE), (lambda_control, InputRole.SPIKE_IN),
             (lysis, InputRole.CONDITION), (bisulfite, InputRole.CONDITION)), (converted,))

    primer_state = Scene()
    primer_state.strand("P5L_AD002", list(random_primer()), mod5="/5SpC3/")
    primers = g.input("random_primers", "published /5SpC3/ P5L indexed random-primer set",
                      MaterialKind.OLIGO, state=primer_state)
    klenow = g.input("klenow", "Klenow exo-minus", MaterialKind.CATALYST)
    primed = g.product("primed", "P5-tailed random-primed contact molecules",
                       MaterialKind.MOLECULAR, state=converted.state)
    g.react("random_prime", "Copy bisulfite DNA from P5L random primers",
            StepKind.EXTENSION,
            ((converted, InputRole.TEMPLATE), (primers, InputRole.STRUCTURAL),
             (klenow, InputRole.CATALYST)), (primed,))
    exo = g.input("exonuclease", "Exonuclease I", MaterialKind.CATALYST)
    sap = g.input("phosphatase", "shrimp alkaline phosphatase", MaterialKind.CATALYST)
    treated = g.product("enzyme_treated", "products after Exonuclease I/SAP treatment", MaterialKind.MOLECULAR)
    g.react("enzyme_cleanup", "Digest residual primers with Exonuclease I and inactivate dNTPs with SAP",
            StepKind.CHEMICAL_CONVERSION,
            ((primed, InputRole.TEMPLATE), (exo, InputRole.CATALYST), (sap, InputRole.CATALYST)), (treated,))
    spri = g.input("spri", "0.8x SPRI beads before four-plex pooling", MaterialKind.SELECTION)
    clean = g.product("clean", "AD002 reaction mixed with SPRI beads before washing", MaterialKind.MOLECULAR)
    g.react("spri_binding", "Bind treated reaction to 0.8x SPRI beads before pooling",
            StepKind.BINDING, ((treated, InputRole.TEMPLATE), (spri, InputRole.SELECTION)), (clean,))
    # One AD002 well is traced above. The other three are parallel, separately
    # processed wells; make their material inputs explicit before pooling.
    other_plates = tuple(g.input(
        "clean_" + member.lower(), f"separately primed {member} reaction mixed with 0.8x SPRI beads",
        MaterialKind.MOLECULAR, state=clean.state,
        source=SUPPLEMENT + " (parallel indexed well population)")
        for member in ("AD006", "AD008", "AD010"))
    pooled = g.product("pooled", "four indexed reaction/bead mixtures pooled before washing",
                       MaterialKind.MOLECULAR, state=clean.state)
    g.react("pool_four", "Combine AD002, AD006, AD008 and AD010 plates for four-plex",
            StepKind.POOL, ((clean, InputRole.TEMPLATE),
                            *((plate, InputRole.TEMPLATE) for plate in other_plates)), (pooled,),
            changes_molecule=False)

    washes = g.input("pool_washes", "three 80% ethanol washes and EB elution for four-plex route", MaterialKind.SELECTION)
    eluted = g.product("eluted_pool", "washed and eluted four-plex DNA", MaterialKind.MOLECULAR)
    wash_waste = g.product("pool_wash_waste", "material removed by pooled bead washes", MaterialKind.WASTE)
    g.react("wash_elute_pool", "Wash pooled beads three times and elute before Adaptase", StepKind.CLEANUP,
            ((pooled, InputRole.TEMPLATE), (washes, InputRole.SELECTION)), (eluted, wash_waste))
    g.discard(wash_waste, "removed in pooled bead washes; no exact composition assigned")
    adaptase = g.input("adaptase", "Swift Accel-NGS Adaptase Module 33096",
                       MaterialKind.CATALYST)
    p7_stub = g.input("adaptase_p7_stub", "2019 Adaptase-added P7-side sequence",
                      MaterialKind.ADAPTER, inferred=True,
                      source=SUPPLEMENT + " (exact bases unresolved in checked sources)")
    adapted = g.product("adapted", "Adaptase-treated P7-side products",
                        MaterialKind.MOLECULAR, inferred=True,
                        source=SUPPLEMENT + " (added bases and attachment chemistry unresolved)")
    g.react("adaptase_attach", "Install the P7-side Adaptase junction",
            StepKind.ADAPTER_ATTACHMENT,
            ((eluted, InputRole.TEMPLATE), (p7_stub, InputRole.STRUCTURAL),
             (adaptase, InputRole.CATALYST)), (adapted,))
    p5 = g.input("p5_pcr", "published P5L_D501–D508 PCR primer family",
                 MaterialKind.OLIGO, state=Construct((Segment("P5L_D501", P5_D501),)))
    p7 = g.input("p7_pcr", "published P7L_D701–D712 PCR primer family",
                 MaterialKind.OLIGO, state=Construct((Segment("P7L_D701", P7_D701),)))
    polymerase = g.input("pcr_polymerase", "KAPA HiFi HotStart PCR mix",
                         MaterialKind.CATALYST)
    library = g.product("amplified", final_library().name, MaterialKind.MOLECULAR,
                        state=final_library(), inferred=True,
                        source=SUPPLEMENT + " (D501/D701/AD002 are published examples, not an assigned well)")
    g.react("index_pcr", "Amplify for 17 cycles with selected P5L and P7L indexing primers",
            StepKind.PCR,
            ((adapted, InputRole.TEMPLATE), (p5, InputRole.STRUCTURAL),
             (p7, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)),
            (library,), geometry=PCRGeometry.unavailable(
                final_library(), reason="The checked 2019 sources do not resolve the Adaptase-added junction; exact P7 primer landing cannot be located."))
    final_spri = g.input("final_spri", "two rounds of 0.8x SPRI purification", MaterialKind.SELECTION)
    selected = g.product("library", final_library().name, MaterialKind.LIBRARY, state=final_library(), inferred=True)
    final_waste = g.product("final_waste", "unrecovered final SPRI material", MaterialKind.WASTE)
    g.react("final_cleanup", "Purify indexed libraries twice with 0.8x SPRI beads", StepKind.CLEANUP,
            ((library, InputRole.TEMPLATE), (final_spri, InputRole.SELECTION)), (selected, final_waste))
    g.discard(final_waste, "not recovered in final library eluate")
    g.sequence(selected, SequencingPlan(
        acquired_roles=("Read 1", "Index 1 (i7)", "Index 2 (i5)", "Read 2"),
        limitation="The 2019 Adaptase junction and selected PCR index members are unresolved; only Read 1 has a source-derived exact site in this representative shell."))
    g.declare_reviewed_scope(tuple(g.steps) + tuple(g.materials), SOURCE)
    return g.require_publishable()


def render_page() -> str:
    g = protocol_graph()
    final = final_library()
    contact = contact_product()
    rows = duplex_rows(contact, label="one contact")
    left_label = len("one contact") + 1 + len("5'- ")
    contact_rows = [*rows[:2], junction_row(contact, "locus A", "contact junction",
                    "ligation", prefix_width=left_label), *rows[2:]]
    caveat = ("INFERRED — D501, D701 and AD002 are exact published example oligos, "
              "but their assignment to a specific well and the 2019 Adaptase-added "
              "P7 junction bases are unresolved from the reviewed Li article and supplement. "
              "The historical vendor reagent composition remains under source review.")
    oligos = [("P5L_AD002 random primer", "/5SpC3/" + P5_RANDOM_PREFIX +
               INLINE_4PLEX["AD002"] + "(N1:25252525)(N1)×8"),
              ("P5L_D501 indexing PCR", P5_D501),
              ("P7L_D701 indexing PCR", P7_D701)]
    oligo_html = "".join(f'<li><strong>{html.escape(name)}</strong> <code>{html.escape(seq)}</code></li>'
                         for name, seq in oligos)
    r1 = sp.section(final, (sp.TRUSEQ["R1"],), heading="Computed Read 1 site",
                    intro="The P5L primer plus random-primer handle form this exact entry site in the representative molecule.",
                    required_roles=("Read 1",))
    return "\n".join((
        head(TITLE), '<div class="wrap" data-protocol-complete="true">',
        f'<h1>{html.escape(TITLE)}</h1>',
        '<p>Published enzyme-label difference: the main Methods and supplementary reagent list specify DpnII; the shared Biotin Fill-in step1 says Inactivate MboI. The model retains DpnII from the main Methods; the intended heat-inactivation instruction remains unresolved.</p>',
        '<p class="research-notes"><a href="01_scmethyl-hic.html">Research notes</a></p>',
        info(f'Defining source: <a href="{SOURCE}">Li et al., <i>Nature Methods</i> (2019)</a>; '
             'original Supplementary Protocol, single-cell branch and oligo table.'),
        info("DpnII proximity ligation is followed by single-nucleus sorting and bisulfite conversion. "
             "A P5-tailed random primer carries one of four six-base pool indexes; Adaptase and "
             "dual-index PCR finish the sequencing library."),
        info(caveat),
        info("SOURCE DISCREPANCY — The article's single-cell Methods (PMC6765429 "
             "extracted line 85) names Zymo D5020; the Supplementary Protocol reagent "
             "list (MOESM1 lines 202–204) names D5023, deep-well format. The graph "
             "uses their shared EZ-96 DNA Methylation-Direct kit family; the catalogue "
             "number used for each preparation is unresolved."),
        info("The article single-cell Methods describe 9ul PBS per sorting well; the detailed protocol specifies 10ul Zymo lysis buffer with proteinase K and fragmented lambda DNA. This graph follows the detailed single-cell well preparation, without asserting that these descriptions are equivalent. Enzyme digestion precedes physical SPRI recovery; four-plex reactions are pooled with beads before washing/elution. The source specifies 17 index-PCR cycles and two final 0.8x SPRI purifications."),
        '<h2>Reaction workflow</h2>', protocol_graph_panels(g, cls="long"),
        '<h2>Published oligos</h2>', '<ul>' + oligo_html + '</ul>',
        '<p>The source prints four P5L random-primer members, eight P5L PCR members and '
        'twelve P7L PCR members. The displayed three are representative published members, '
        'not a recorded assignment. The first random N uses the source’s synthesis notation.</p>',
        '<h2>Proximity-ligation product</h2>',
        panel(contact_rows, cls="small", caption="DpnII filled contact junction; ** marks the ligation."),
        '<h2>Final sequencing molecule</h2>',
        panel(duplex_rows(final, label="library"), cls="small",
              caption="INFERRED — P7-side Adaptase junction remains unresolved from the "
                      "reviewed 2019 sources; D501, D701 and AD002 are published illustrative members."),
        r1,
        '<p>Index 1, Index 2 and Read 2 are acquired, but an exact placement diagram for '
        'all four reads depends on the Adaptase junction unresolved from the checked 2019 '
        'article and supplement, and on the selected index members.</p>',
        '<p>Source boundary: two independent readers checked the article, every listed supplement, '
        'the 2021 Swift-to-IDT product mapping and the 2022 IDT Adaptase successor guide. '
        'The version-matched 2019 Swift manual remains unavailable; the successor guide describes '
        'a related stub and synthetic tail but does not establish the exact 2019 reagent sequence. '
        'The <a href="01_scmethyl-hic.html">source checklist</a> records each link and its status.</p>',
        '</div>'))


def selftest() -> None:
    c = Check()
    c.section("scMethyl-HiC source and branch behavior")
    c("four published pool indexes", tuple(INLINE_4PLEX),
      ("AD002", "AD006", "AD008", "AD010"))
    data = data_root()
    original = Source(data / "sources" / HERE.parent.name /
        "s41592-019-0502-z_41592_2019_502_MOESM1_ESM.pdf.txt",
        "Li et al. 2019 Springer Supplementary Protocol, MOESM1; tools/get_sources.py --doi 10.1038/s41592-019-0502-z")
    if have(c, original, label="source transcription of single-cell primer families"):
        oligos = re.sub(r"\s+", "", original.path.read_text(encoding="utf-8"))
        c("all P5 and P7 PCR primer members transcribed",
          (all(seq in oligos for seq in P5_PCR.values()),
           all(seq in oligos for seq in P7_PCR.values())), (True, True))
        c("all four modified random primers present in source",
          all("/5SpC3/" + P5_RANDOM_PREFIX + value in oligos
              for value in INLINE_4PLEX.values()), True)
    c("D501 and D701 index slices are derived from P5/P7 arms",
      (next(s.top for s in final_library() if s.name == "P5 PCR index, D501 example"),
       revcomp(next(s.top for s in final_library()
                    if s.name == "P7 PCR index, D701 example, reverse complement"))),
      (P5_D501[len(il.P5):len(il.P5)+8], P7_D701[len(il.P7):len(il.P7)+8]))
    c("Read 1 lands on representative library", sp.locate(final_library(), sp.TRUSEQ["R1"]) is not None, True)
    c.report()
