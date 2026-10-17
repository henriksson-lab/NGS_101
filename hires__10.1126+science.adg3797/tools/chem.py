"""HiRES joint Hi-C and RNA library, as described by the author protocol."""
from multibranch import illumina_branch
from pcr import PCRGeometry
from chemdraw import Construct, Segment
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind

TITLE = "HiRES — single-cell Hi-C and RNA-seq"
NOTES = "01_hires.html"
SOURCE = ('Defining <a href="https://doi.org/10.1126/science.adg3797">Science paper</a>; '
          '<a href="https://doi.org/10.21769/BioProtoc.4886">author protocol</a>; '
          'authors\' <a href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE223917">GEO GSE223917</a> sample records.')
SUMMARY = ("The same fixed nucleus is reverse-transcribed, MboI-digested and proximity-ligated, then sorted into its own well. "
           "A joint amplification feeds tagmentation, an RNA-P7/i5 enrichment and i7 indexing PCR.")
CAVEAT = ("The defining Science article and supplement were not available; this page follows the authors' Bio-protocol (doi:10.21769/BioProtoc.4886). "
          "The authors' GEO GSE223917 sample records repeat the tagmentation, two-stage index PCR, RNA-P7 sequence and paired 150-base NovaSeq 6000 reads, "
          "and state that RNA reads are recognized by the GAT5-RT arm at the 5′ end of Read 2 and deduplicated by UMI. "
          "The author protocol prints four HiRES primers and specifies per-well i5/i7 indexing. "
          "The TruePrep manual prints a standard finished-library shell, but its loaded adapter oligos, "
          "the HiRES RNA-P7 junction and selected index-member bases remain unresolved. "
          "The two views below are read classes from one library, not separate index-ligation rounds.")

CONTACT, CP = illumina_branch("HiRES contact reads", "DNA–DNA contact insert",
                               inferred_adapters=True)
RNA, RP = illumina_branch("HiRES transcript reads", "RNA-derived cDNA insert",
                           inferred_adapters=True)
FINAL_LIBRARIES = (
    ("Contact-read view", CONTACT, CP,
     "INFERRED — contact reads in the joint indexed library; TruePrep adapter bases are unresolved.",
     "Paired reads carrying Hi-C contacts are recovered from the indexed library."),
    ("Transcript-read view", RNA, RP,
     "INFERRED — transcript reads in the same library; TruePrep adapter bases are unresolved.",
     "Transcript reads from the same indexed wells are classified during analysis."),
)
READ_LENGTHS = {name: {"Read 1": 150, "Read 2": 150}
                for name, *_ in FINAL_LIBRARIES}

GAT5_RT = "GTAGGTGTGAGTGATGGTTGAGGTAGTATTGCGCAATGNNNNNNNNTTTTTTTTTTTTTTTVN"
GAT5_7N = "GTAGGTGTGAGTGATGGTTGAGGTAGTNNNNNNN"
GAT5 = "GTAGGTGTGAGTGATGGTTGAGGTAGT"
RNA_P7 = "GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAGGGTTGAGGTAGTATTGCGCAATG"

def oligos():
    return [f"<b>GAT5-RT</b> — 5′-{GAT5_RT}-3′",
            f"<b>GAT5-7N</b> — 5′-{GAT5_7N}-3′",
            f"<b>GAT5</b> — 5′-{GAT5}-3′",
            f"<b>RNA-P7</b> — 5′-{RNA_P7}-3′"]

def _oligo(name, seq):
    """Printed primer as fixed DNA runs and degenerate placeholder runs (N, V)."""
    import re
    parts = [Segment(f"{name} {'degenerate' if run[0] in 'NV' else 'fixed'} {i + 1}", run,
                     placeholder=run[0] in "NV")
             for i, run in enumerate(re.findall(r"[ACGT]+|[NV]+", seq))]
    return Construct(parts, name=name)

def protocol_graph():
    source = "Chen et al. 2023: https://doi.org/10.21769/BioProtoc.4886"
    g = ProtocolGraph(TITLE, authority=source)
    cells = g.input("cells", "embryo single-cell suspension (~100,000 in 200 µl PBS) or cortex nuclei (~250,000 in 500 µl storage buffer)",
                    MaterialKind.BIOLOGICAL, primary=True)
    formaldehyde = g.input("formaldehyde", "1.75% formaldehyde 10 min RT, BSA quench; Igepal lysis; 0.5% SDS 65 °C 10 min, Triton quench",
                           MaterialKind.CONDITION)
    nuclei = g.product("nuclei", "fixed, permeabilized nuclei containing chromatin and RNA", MaterialKind.MOLECULAR)
    g.react("fix_permeabilize", "Fix, lyse and SDS-permeabilize nuclei (procedure C)", StepKind.CROSS_LINKING,
            ((cells, InputRole.TEMPLATE), (formaldehyde, InputRole.CONDITION)), (nuclei,))
    rt_primer = g.input("gat5_rt", "published GAT5-RT anchored primer", MaterialKind.OLIGO,
                    state=_oligo("GAT5-RT", GAT5_RT))
    rt = g.input("rt", "SuperScript IV reverse transcriptase (10 U/µl, 55 °C 15 min)", MaterialKind.CATALYST)
    copied = g.product("copied", "nuclei containing in-situ cDNA and chromatin", MaterialKind.MOLECULAR)
    g.react("reverse_transcribe", "Reverse-transcribe RNA in situ", StepKind.REVERSE_TRANSCRIPTION,
            ((nuclei, InputRole.TEMPLATE), (rt_primer, InputRole.STRUCTURAL),
             (rt, InputRole.CATALYST)), (copied,))
    mboi = g.input("mboi", "MboI restriction nuclease (37 °C overnight; 65 °C 20 min inactivation)", MaterialKind.CATALYST)
    digested = g.product("digested", "MboI-cut chromatin with retained cDNA", MaterialKind.MOLECULAR)
    g.react("digest", "Digest chromatin with MboI", StepKind.DIGESTION,
            ((copied, InputRole.TEMPLATE), (mboi, InputRole.CATALYST)), (digested,))
    ligase = g.input("ligase", "T4 DNA ligase (Invitrogen, 16 °C 4 h, diluted to ~1 ml)", MaterialKind.CATALYST)
    joined = g.product("joined", "nuclei with proximity-ligated DNA and cDNA", MaterialKind.MOLECULAR)
    g.react("proximity_ligate", "Ligate nearby genomic ends", StepKind.PROXIMITY_LIGATION,
            ((digested, InputRole.TEMPLATE), (ligase, InputRole.CATALYST)), (joined,))
    sorter = g.input("facs", "FACSAria single-nucleus sort into dry 96-well plates (FSC/SSC gating)",
                     MaterialKind.CONDITION)
    sorted_nuclei = g.product("sorted", "one ligated nucleus per well", MaterialKind.MOLECULAR)
    unsorted = g.product("unsorted", "debris, clumps and unsorted nuclei", MaterialKind.WASTE)
    g.react("sort", "Sort single nuclei into wells (procedure F)", StepKind.PARTITION,
            ((joined, InputRole.TEMPLATE), (sorter, InputRole.CONDITION)), (sorted_nuclei, unsorted),
            changes_molecule=False)
    g.discard(unsorted, "excluded by FSC/SSC gating or left in the suspension")
    protease = g.input("qiagen_protease", "Qiagen protease in SC lysis buffer with GAT5 (50 °C 3 h, 70 °C 1 h)",
                       MaterialKind.CATALYST)
    lysed = g.product("lysed", "deproteinized single-nucleus contact DNA and cDNA", MaterialKind.MOLECULAR)
    g.react("lyse", "Lyse each nucleus with protease (procedure G1-3)", StepKind.DEPROTEINIZATION,
            ((sorted_nuclei, InputRole.TEMPLATE), (protease, InputRole.CATALYST)), (lysed,))
    gat5 = g.input("gat5", "published GAT5 primer", MaterialKind.OLIGO,
                    state=_oligo("GAT5", GAT5))
    gat5_7n = g.input("gat5_7n", "published GAT5-7N primer", MaterialKind.OLIGO,
                    state=_oligo("GAT5-7N", GAT5_7N))
    polymerase = g.input("polymerase", "Deep Vent (exo-) polymerase, ThermoPol buffer", MaterialKind.CATALYST)
    amplified = g.product("amplified", "joint single-nucleus amplicons", MaterialKind.MOLECULAR)
    g.react("preamp", "Amplify contact DNA and transcript cDNA in one well (9 stepped random-priming cycles, then 14 cycles)", StepKind.PCR,
            ((lysed, InputRole.TEMPLATE), (gat5, InputRole.STRUCTURAL),
             (gat5_7n, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)),
            (amplified,), geometry=PCRGeometry.unavailable(
                CONTACT, reason="The mixed preamplification product has no single primer-landing construct."))
    transposome = g.input("transposome", "TruePrep V2 TTE Mix V50 in TTBL", MaterialKind.ADAPTER,
                          inferred=True, source=source + " (adapter oligos unavailable)")
    tn5 = g.input("tn5", "TruePrep transposase (55 °C 10 min; 0.2% SDS stop)", MaterialKind.CATALYST)
    tagged = g.product("tagged", "tagmented joint amplicons", MaterialKind.MOLECULAR,
                       inferred=True, source=source + " (adapter oligos unavailable)")
    tagment_waste = g.product("tagment_waste", "unusable tagmentation fragments", MaterialKind.WASTE)
    g.react("tagment", "Tagment 2 µl diluted preamplification product", StepKind.TAGMENTATION,
            ((amplified, InputRole.TEMPLATE), (transposome, InputRole.STRUCTURAL),
             (tn5, InputRole.CATALYST)), (tagged, tagment_waste))
    g.discard(tagment_waste, "does not retain library-eligible adapter ends")
    rna_p7 = g.input("rna_p7", "published RNA-P7 enrichment primer", MaterialKind.OLIGO,
                    state=_oligo("RNA-P7", RNA_P7))
    i5 = g.input("i5", "well-specific i5 UD index primer", MaterialKind.OLIGO,
                 inferred=True, source=source + " (index member bases unavailable)")
    pcr_pol = g.input("pcr_polymerase", "KAPA HiFi polymerase, GC buffer", MaterialKind.CATALYST)
    enriched = g.product("enriched", "i5-extended products with RNA-P7-enriched cDNA reads",
                         MaterialKind.MOLECULAR, inferred=True,
                         source=source + " (TruePrep adapter bases unavailable)")
    g.react("enrich_pcr", "72 °C 3 min, then 5 cycles with RNA-P7 and i5 only (slight cDNA enrichment)", StepKind.PCR,
            ((tagged, InputRole.TEMPLATE), (rna_p7, InputRole.STRUCTURAL),
             (i5, InputRole.STRUCTURAL), (pcr_pol, InputRole.CATALYST)), (enriched,),
            geometry=PCRGeometry.unavailable(
                CONTACT, reason="The source supplies RNA-P7 and the i5 role, but the TruePrep landing geometry is unavailable."))
    i7 = g.input("i7", "well-specific i7 UD index primer", MaterialKind.OLIGO,
                 inferred=True, source=source + " (index member bases unavailable)")
    indexed = g.product("indexed", "per-well dual-indexed HiRES library", MaterialKind.LIBRARY,
                        inferred=True, source=source + " (TruePrep adapter bases unavailable)")
    g.react("index_pcr", "Add i7 to the same well (i5, RNA-P7 and KAPA HiFi carried over) and run 5 more cycles", StepKind.PCR,
            ((enriched, InputRole.TEMPLATE), (i5, InputRole.STRUCTURAL), (i7, InputRole.STRUCTURAL),
             (pcr_pol, InputRole.CATALYST)), (indexed,),
            geometry=PCRGeometry.unavailable(
                CONTACT, reason="The i7 UD adapter arm and TruePrep landing geometry are not reconstructed."))
    ampure = g.input("ampure", "0.6x then 0.15x AMPure XP double-sided selection", MaterialKind.SELECTION)
    joint = g.product("joint_library", "pooled, size-selected joint HiRES library", MaterialKind.LIBRARY,
                      inferred=True, source=source + " (TruePrep adapter bases unavailable)")
    size_waste = g.product("size_waste", "fragments outside the bead size window and leftover primers", MaterialKind.WASTE)
    g.react("pool_select", "Pool wells and size-select with AMPure XP", StepKind.CLEANUP,
            ((indexed, InputRole.TEMPLATE), (ampure, InputRole.SELECTION)), (joint, size_waste),
            changes_molecule=False)
    g.discard(size_waste, "removed by the double-sided AMPure selection")
    contact = g.product("contact_view", CONTACT.name, MaterialKind.LIBRARY,
                        state=CONTACT, inferred=True, source=source)
    transcript = g.product("transcript_view", RNA.name, MaterialKind.LIBRARY,
                           state=RNA, inferred=True, source=source)
    g.react("read_classes", "Classify contact and transcript reads from the joint library",
            StepKind.SPLIT, ((joint, InputRole.TEMPLATE),), (contact, transcript),
            changes_molecule=False)
    for product, library, declared in ((contact, CONTACT, CP), (transcript, RNA, RP)):
        g.sequence(product, SequencingPlan(
            library=library, primers=declared,
            acquired_roles=tuple(p.role for p in declared),
            read_uses={p.role: ("demultiplex" if "Index" in p.role else "include")
                       for p in declared}))
    g.declare_reviewed_scope((
        "cells", "formaldehyde", "fix_permeabilize", "nuclei",
        "gat5_rt", "rt", "reverse_transcribe", "mboi", "digest",
        "ligase", "proximity_ligate", "facs", "sort", "qiagen_protease", "lyse",
        "gat5", "gat5_7n", "polymerase", "preamp",
        "transposome", "tn5", "tagment", "rna_p7", "i5", "pcr_polymerase", "enrich_pcr",
        "i7", "index_pcr", "ampure", "pool_select",
        "read_classes", "contact_view", "transcript_view"), source)
    return g.require_publishable()
