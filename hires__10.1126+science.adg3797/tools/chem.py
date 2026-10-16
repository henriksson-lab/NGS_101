"""HiRES joint Hi-C and RNA library, as described by the author protocol."""
from multibranch import illumina_branch
from pcr import PCRGeometry
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind

TITLE = "HiRES — single-cell Hi-C and RNA-seq"
NOTES = "01_hires.html"
SOURCE = ('Defining <a href="https://doi.org/10.1126/science.adg3797">Science paper</a>; '
          '<a href="https://doi.org/10.21769/BioProtoc.4886">author protocol</a>.')
SUMMARY = ("The same fixed nucleus is reverse-transcribed, MboI-digested and proximity-ligated. "
           "A joint amplification then feeds tagmentation and per-well indexed PCR.")
CAVEAT = ("The author protocol prints four HiRES primers and specifies per-well i5/i7 indexing. "
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

def protocol_graph():
    source = "Chen et al. 2023: https://doi.org/10.21769/BioProtoc.4886"
    g = ProtocolGraph(TITLE, authority=source)
    nuclei = g.input("nuclei", "fixed nuclei containing chromatin and RNA",
                     MaterialKind.BIOLOGICAL, primary=True)
    rt_primer = g.input("gat5_rt", "published GAT5-RT anchored primer", MaterialKind.OLIGO)
    rt = g.input("rt", "SuperScript IV reverse transcriptase", MaterialKind.CATALYST)
    copied = g.product("copied", "nuclei containing in-situ cDNA and chromatin", MaterialKind.MOLECULAR)
    g.react("reverse_transcribe", "Reverse-transcribe RNA in situ", StepKind.REVERSE_TRANSCRIPTION,
            ((nuclei, InputRole.TEMPLATE), (rt_primer, InputRole.STRUCTURAL),
             (rt, InputRole.CATALYST)), (copied,))
    mboi = g.input("mboi", "MboI restriction nuclease", MaterialKind.CATALYST)
    digested = g.product("digested", "MboI-cut chromatin with retained cDNA", MaterialKind.MOLECULAR)
    g.react("digest", "Digest chromatin with MboI", StepKind.DIGESTION,
            ((copied, InputRole.TEMPLATE), (mboi, InputRole.CATALYST)), (digested,))
    ligase = g.input("ligase", "T4 DNA ligase", MaterialKind.CATALYST)
    joined = g.product("joined", "nuclei with proximity-ligated DNA and cDNA", MaterialKind.MOLECULAR)
    g.react("proximity_ligate", "Ligate nearby genomic ends", StepKind.PROXIMITY_LIGATION,
            ((digested, InputRole.TEMPLATE), (ligase, InputRole.CATALYST)), (joined,))
    gat5 = g.input("gat5", "published GAT5 primer", MaterialKind.OLIGO)
    gat5_7n = g.input("gat5_7n", "published GAT5-7N primer", MaterialKind.OLIGO)
    polymerase = g.input("polymerase", "Deep Vent exo-minus polymerase", MaterialKind.CATALYST)
    amplified = g.product("amplified", "joint single-nucleus amplicons", MaterialKind.MOLECULAR)
    g.react("preamp", "Amplify contact DNA and transcript cDNA in one well", StepKind.PCR,
            ((joined, InputRole.TEMPLATE), (gat5, InputRole.STRUCTURAL),
             (gat5_7n, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)),
            (amplified,), geometry=PCRGeometry.unavailable(
                CONTACT, reason="The mixed preamplification product has no single primer-landing construct."))
    transposome = g.input("transposome", "TruePrep V2 tagmentation reagent", MaterialKind.ADAPTER,
                          inferred=True, source=source + " (adapter oligos unavailable)")
    tn5 = g.input("tn5", "TruePrep transposase", MaterialKind.CATALYST)
    tagged = g.product("tagged", "tagmented joint amplicons", MaterialKind.MOLECULAR,
                       inferred=True, source=source + " (adapter oligos unavailable)")
    tagment_waste = g.product("tagment_waste", "unusable tagmentation fragments", MaterialKind.WASTE)
    g.react("tagment", "Tagment joint amplicons", StepKind.TAGMENTATION,
            ((amplified, InputRole.TEMPLATE), (transposome, InputRole.STRUCTURAL),
             (tn5, InputRole.CATALYST)), (tagged, tagment_waste))
    g.discard(tagment_waste, "does not retain library-eligible adapter ends")
    rna_p7 = g.input("rna_p7", "published RNA-P7 enrichment primer", MaterialKind.OLIGO)
    indexes = g.input("indexes", "well-specific i5/i7 dual-index primers", MaterialKind.OLIGO,
                      inferred=True, source=source + " (index member bases unavailable)")
    pcr_pol = g.input("pcr_polymerase", "KAPA HiFi polymerase", MaterialKind.CATALYST)
    joint = g.product("joint_library", "indexed joint HiRES library", MaterialKind.LIBRARY,
                      inferred=True, source=source + " (TruePrep adapter bases unavailable)")
    g.react("index_pcr", "Enrich transcripts and add per-well i5/i7 indexes", StepKind.PCR,
            ((tagged, InputRole.TEMPLATE), (rna_p7, InputRole.STRUCTURAL),
             (indexes, InputRole.STRUCTURAL), (pcr_pol, InputRole.CATALYST)), (joint,),
            geometry=PCRGeometry.unavailable(
                CONTACT, reason="The source supplies RNA-P7 and kit i5/i7 roles, but the TruePrep landing geometry is unavailable."))
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
        "nuclei", "gat5_rt", "rt", "reverse_transcribe", "mboi", "digest",
        "ligase", "proximity_ligate", "gat5", "gat5_7n", "polymerase", "preamp",
        "transposome", "tn5", "tagment", "rna_p7", "indexes", "pcr_polymerase",
        "index_pcr", "read_classes", "contact_view", "transcript_view"), source)
    return g.require_publishable()
