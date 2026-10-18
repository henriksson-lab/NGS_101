"""Digenome-seq: source-owned WGS boundary with explicit version limits."""
from chemdraw import Construct, Scene
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind
from protocol_model import molecule, role

TITLE = "Digenome-seq — digested-genome whole-genome sequencing"
NOTES = "01_digenome_seq.html"
SOURCE = ('<a href="https://doi.org/10.1038/nmeth.3284">Kim et al., Nature Methods (2015)</a>, published supplement; '
          'detailed source-limited route from the authors’ <a href="https://doi.org/10.1101/gr.199588.115">2016 multiplex Digenome-seq companion</a>, Discussion and Methods (PMC4772022).')
SUMMARY = "Cas9–sgRNA cleaves cell-free genomic DNA; unenriched whole-genome sequencing detects recurrent genomic read-start boundaries. The source-limited preparation below uses the explicitly identified2016 companion method."
CAVEAT = ("The defining2015 main Methods were not inspected. Its published supplement and the2016 companion are distinct sources; the later route does not establish the exact2015 kit or run. "
          "The2016 Discussion explicitly says no PCR amplification precedes WGS; its Methods names a TruSeq DNA library prep kit without a catalogue member. "
          "End preparation, physical adapter geometry, index architecture/member/length and acquired cycle counts remain unassigned. No generic TruSeq or dA-tailed endpoint is supplied.")
FINAL_LIBRARY = Construct([
    role("WGS adapter end; physical sequence/geometry unassigned"),
    role("unenriched genomic insert; Cas9 and/or shear endpoints"),
    role("opposite WGS adapter end; index architecture/width unassigned"),
], name="Digenome source-limited WGS endpoint")
SEQ_PRIMERS = ()
FINAL_CAPTION = "Source-limited2016 TruSeq-family WGS boundary; this does not establish the original2015 physical library ends or index choices."
READ_LENGTHS = None
SEQUENCING_INTRO = "The2016 companion reports HiSeqXTen WGS at30–40× depth. Defining2015 Supplementary Note1 mentions HiSeqX10 for knockout-clone WGS, a different sample scope. Acquired biological/index cycles and exact primer sites are unassigned."
SEQUENCING_ENDING = SEQUENCING_INTRO
SEQUENCING_UNAVAILABLE = "The cleavage signal is the genomic5′ boundary distribution; no literal cut-end tag or exact historical library/run geometry is inferred."


def oligos():
    return ()


def fragment_populations():
    scene = Scene()
    scene.strand("cas9_shear", [role("representative fragment with a Cas9-created end and a shear-created end")], label="Cas9/shear boundary class")
    scene.strand("shear_shear", [role("independent fragment with random shearing boundaries")], label="random-shear background class")
    return scene


def protocol_graph():
    source = "Kim etal.2016 companion PMC4772022 Discussion114 and Methods132–136; original2015 library details not established by this source"
    graph = ProtocolGraph(TITLE, authority=SOURCE)
    dna = graph.input("dna", "purified cell-free genomic DNA", MaterialKind.BIOLOGICAL,
                      state=molecule("genomic DNA", (role("genomic DNA"),)), primary=True)
    sgrna = graph.input("sgrna", "study-specific sgRNA; sequence not assigned to one target", MaterialKind.OLIGO)
    cas9 = graph.input("cas9", "Cas9 nuclease for in-vitro cleavage", MaterialKind.CATALYST)
    cut = graph.product("cut", "Cas9-digested whole genomic population", MaterialKind.MOLECULAR,
                        state=molecule("digested genomic DNA", (role("genomic DNA retaining cut and uncut regions"),)))
    graph.react("cleave", "Cleave purified genomic DNA with Cas9–sgRNA", StepKind.DIGESTION,
                ((dna,InputRole.TEMPLATE),(sgrna,InputRole.STRUCTURAL),(cas9,InputRole.CATALYST)), (cut,), source=source)
    cleanup = graph.input("cleanup_condition", "2016 companion:RNaseA50µg/ml, then DNeasy tissue-kit repurification", MaterialKind.CONDITION)
    clean = graph.product("clean", "repurified genomic DNA retaining nuclease boundaries", MaterialKind.MOLECULAR, state=cut.state)
    graph.react("cleanup", "Remove sgRNA and repurify digested DNA", StepKind.EXTRACTION,
                ((cut,InputRole.TEMPLATE),(cleanup,InputRole.CONDITION)), (clean,), source=source)
    shear = graph.input("shear_condition", "2016 companion:Covaris fragmentation to400–500bp", MaterialKind.CONDITION)
    fragments = graph.product("fragments", "independent Cas9/shear and shear/shear WGS fragments", MaterialKind.MOLECULAR, state=fragment_populations())
    graph.react("shear", "Shear the unenriched genomic DNA population", StepKind.FRAGMENTATION,
                ((clean,InputRole.TEMPLATE),(shear,InputRole.CONDITION)), (fragments,), source=source)
    kit = graph.input("truseq_library_kit", "2016 companion:TruSeq DNA library prep kit; exact member/end preparation unassigned", MaterialKind.CONDITION)
    adapters = graph.input("truseq_adapters", "source-named TruSeq-family adapters; physical oligos/index format unassigned", MaterialKind.ADAPTER)
    library = graph.product("library", FINAL_LIBRARY.name, MaterialKind.LIBRARY, state=FINAL_LIBRARY)
    graph.react("prepare_wgs_library", "Prepare the source-limited TruSeq-family WGS library", StepKind.MAGIC,
                ((fragments,InputRole.TEMPLATE),(adapters,InputRole.STRUCTURAL),(kit,InputRole.CONDITION)), (library,), source=source,
                reason="The2016 source states adapter ligation with a TruSeq DNA kit and no PCR before WGS, but does not assign an exact kit member. Individual repair, size-selection, dA-tail, adapter-end/index choices from the former generic recipe are unassigned; the defining2015 main Methods are not inspected.")
    graph.sequence(library, SequencingPlan(library=FINAL_LIBRARY, acquisition="unspecified", limitation=SEQUENCING_INTRO+" "+SEQUENCING_UNAVAILABLE))
    graph.declare_reviewed_scope(("dna","sgrna","cas9","cleave","cleanup_condition","cleanup","shear_condition","shear",
                                  "truseq_library_kit","truseq_adapters","prepare_wgs_library","library"), source)
    return graph.require_publishable()
