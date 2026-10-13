"""HiRES joint chromatin-contact and RNA profiling."""
from chemdraw import Row
from adapter_ligation import AdapterLigationGeometry
from multibranch import illumina_branch
from multimodal_spatial import finish_pcr_branch
from protocol_extensions import analyte_split, proximity_product
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, StepKind
TITLE="HiRES — single-cell Hi-C and RNA-seq"
NOTES="01_hires.html"
SOURCE='Defining source: <a href="https://doi.org/10.1126/science.adg3797">Liu et al., <i>Science</i> (2023)</a>.'
SUMMARY="Fixed nuclei undergo in situ proximity ligation and RNA capture with coordinated combinatorial indexes, producing separately amplified contact and transcript libraries from the same cell."
CAVEAT="The defining paper establishes three coordinated cell-index rounds. This page has not yet transcribed their strand-resolved adapter overhangs or the exact branch-finishing PCR-primer pairs from the supplement; those geometries remain explicitly unavailable."
CELL_PARTS=(("round 1 index",8),("round 2 index",8),("round 3 index",8))
CONTACT,CP=illumina_branch("HiRES contacts","DNA–DNA ligation junction",cell_parts=CELL_PARTS)
RNA,RP=illumina_branch("HiRES RNA","RNA-derived cDNA",cell_parts=CELL_PARTS,umi=8)
FINAL_LIBRARIES=(("Chromatin-contact library",CONTACT,CP,"Paired genomic inserts flank an in situ ligation junction.","Illumina primers on the contact branch."),("RNA library",RNA,RP,"Transcript cDNA carries the matching combinatorial cell identity.","Illumina primers on the RNA branch."))
def sections():
    return [("Digest fixed chromatin and ligate nearby DNA ends",proximity_product(left="contact locus A",right="contact locus B",bridge="Hi-C junction").rows(),"In situ ligation records spatially proximal genomic loci."),("Capture and reverse-transcribe RNA in the same nucleus",[Row(chunks=[("cellular RNA → indexed cDNA", "r1", False)])],"The RNA route remains distinct from the DNA–DNA contact junction."),("Apply coordinated combinatorial indexes",[Row(chunks=[("Hi-C product ** cell-index rounds", "cbc", False)]),Row(chunks=[("RNA cDNA ** matching cell-index rounds", "cbc", False)])],"The index combination links the two modalities to one cell."),("Amplify branch-specific libraries",analyte_split("chromatin-contact library","RNA library"),"Separate primer selections retain the shared cellular provenance."),]


def protocol_graph():
    g = ProtocolGraph(TITLE, authority=SOURCE)
    nuclei = g.input(
        "nuclei", "fixed nuclei containing chromatin and RNA",
        MaterialKind.BIOLOGICAL, state=sections()[0][1], primary=True)
    nuclease = g.input("nuclease", "chromatin restriction nuclease",
                       MaterialKind.CATALYST)
    digested = g.product(
        "digested", "restriction-digested chromatin with retained RNA",
        MaterialKind.MOLECULAR, state=nuclei.state)
    g.react("digest", "Digest fixed chromatin in situ", StepKind.DIGESTION,
            ((nuclei, InputRole.TEMPLATE), (nuclease, InputRole.CATALYST)),
            (digested,))
    fill_polymerase = g.input("fill_polymerase", "restriction-end fill-in polymerase",
                              MaterialKind.CATALYST)
    filled = g.product("filled", "filled chromatin ends with retained RNA",
                       MaterialKind.MOLECULAR, state=digested.state)
    g.react("fill", "Fill restriction-fragment ends", StepKind.EXTENSION,
            ((digested, InputRole.TEMPLATE),
             (fill_polymerase, InputRole.CATALYST)), (filled,))
    ligase = g.input("contact_ligase", "in-situ proximity ligase",
                     MaterialKind.CATALYST)
    contacts = g.product(
        "contacts", "nuclei retaining proximity-ligated DNA and RNA",
        MaterialKind.MOLECULAR, state=sections()[0][1])
    g.react("proximity_ligate", "Ligate spatially adjacent chromatin ends",
            StepKind.PROXIMITY_LIGATION,
            ((filled, InputRole.TEMPLATE), (ligase, InputRole.CATALYST)),
            (contacts,))

    rt_primer = g.input("rt_primer", "RNA-capture reverse-transcription primer",
                        MaterialKind.OLIGO, inferred=True,
                        source=SOURCE + " (exact sequence not yet transcribed)")
    rt = g.input("rt", "reverse transcriptase", MaterialKind.CATALYST)
    copied = g.product(
        "copied", "nuclei retaining proximity-ligated DNA and indexed cDNA",
        MaterialKind.MOLECULAR, state=sections()[1][1], inferred=True,
        source=SOURCE + " (exact RT-primer sequence not yet transcribed)")
    g.react("reverse_transcribe", "Capture and reverse-transcribe RNA in the same nucleus",
            StepKind.REVERSE_TRANSCRIPTION,
            ((contacts, InputRole.TEMPLATE), (rt_primer, InputRole.STRUCTURAL),
             (rt, InputRole.CATALYST)), (copied,))

    barcode_ligase = g.input("barcode_ligase", "combinatorial-index DNA ligase",
                             MaterialKind.CATALYST)
    current = copied
    reviewed = [
        "nuclei", "nuclease", "digest", "fill_polymerase", "fill",
        "contact_ligase", "proximity_ligate", "rt_primer", "rt",
        "reverse_transcribe", "barcode_ligase",
    ]
    for number in (1, 2, 3):
        adapter = g.input(
            f"barcode_{number}", f"coordinated round-{number} cell-index adapter",
            MaterialKind.ADAPTER, inferred=True,
            source=SOURCE + " (strand-resolved adapter geometry not yet transcribed)")
        indexed = g.product(
            f"indexed_{number}",
            f"contact DNA and cDNA after coordinated index round {number}",
            MaterialKind.MOLECULAR, state=sections()[2][1], inferred=True,
            source=SOURCE + " (cell-index identity known; junction geometry unresolved)")
        g.react(
            f"ligate_{number}", f"Ligate coordinated cell-index round {number}",
            StepKind.LIGATION,
            ((current, InputRole.TEMPLATE), (adapter, InputRole.STRUCTURAL),
             (barcode_ligase, InputRole.CATALYST)), (indexed,),
            geometry=AdapterLigationGeometry.unavailable(
                f"Round {number} is source-defined, but its strand-resolved overhang has not yet been transcribed from the supplement."))
        reviewed.extend((adapter.id, f"ligate_{number}"))
        current = indexed

    contact = g.product("contact_branch", "chromatin-contact finishing branch",
                        MaterialKind.MOLECULAR, state=CONTACT)
    rna = g.product("rna_branch", "transcript finishing branch",
                    MaterialKind.MOLECULAR, state=RNA)
    g.react("split", "Allocate indexed products to contact and RNA primer selections",
            StepKind.SPLIT, ((current, InputRole.TEMPLATE),), (contact, rna),
            changes_molecule=False)
    reviewed.append("split")
    finish_pcr_branch(
        g, "contact", contact, CONTACT, CP, reviewed,
        label="HiRES chromatin-contact", inferred_primers=True,
        limitation="The exact branch-finishing PCR-primer pair has not yet been transcribed from the supplement.")
    finish_pcr_branch(
        g, "rna", rna, RNA, RP, reviewed, label="HiRES RNA",
        inferred_primers=True,
        limitation="The exact branch-finishing PCR-primer pair has not yet been transcribed from the supplement.")
    g.declare_reviewed_scope(reviewed, SOURCE)
    return g.require_publishable()
