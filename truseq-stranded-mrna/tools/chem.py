"""TruSeq Stranded mRNA v00 owns its workflow and conditional CD endpoint."""
import illumina as il
import seqprimers as sp
from adapter_ligation import AdapterLigationGeometry
from chemdraw import Construct, Scene, Segment, feature
from endprep import EndPreparedDuplex, dA_tailed_scene
from pcr import PCRGeometry
from protocol_graph import (Disposition, Highlight, InputRole, MaterialKind,
                            ProtocolGraph, SequencingPlan, StepKind, Terminal)
from protocol_model import role

TITLE = "TruSeq Stranded mRNA"
NOTES = "01_truseq_stranded_mrna.html"
GUIDE_URL = "https://emea.support.illumina.com/content/dam/illumina-support/documents/documentation/chemistry_documentation/samplepreps_truseq/truseq-stranded-mrna-workflow/truseq-stranded-mrna-workflow-reference-1000000040498-00.pdf"
SOURCE = f'Illumina <a href="{GUIDE_URL}">TruSeq Stranded mRNA Reference Guide 1000000040498 v00</a> (October 2017).'
SUMMARY = ("Two rounds of poly(A) capture precede random-primed first-strand synthesis, "
           "dUTP-marked second-strand synthesis, dA-tailing, indexed-adapter ligation and selective PCR.")
CAVEAT = ("The workflow follows the standard HS handling branch, omitting optional inline DNA controls. "
          "INFERRED — the final duplex is a source-compatible RNA combinatorial-dual (CD) reference, "
          "using the later official Adapter Sequences v22 templates. It does not establish a particular 2017 adapter revision or D-index pair. "
          "The actual PPC primer strands, sequencing platform, physical index-2 priming workflow and acquired cycles remain unassigned model data.")
PCR_CYCLES = 15
SEQ_PRIMERS = ()
READ_SPECS = (
    sp.read_spec("Read 1", None, "Antisense to the original RNA; acquired length unassigned.", "include"),
    sp.read_spec("Read 2", None, "Sense to the original RNA if paired-end sequencing is selected; acquisition unassigned.", "include"),
    sp.read_spec("Index 1 (i7)", None, "CD reference i7 site is 8 nt; actual acquisition and member unassigned.", "extract"),
    sp.read_spec("Index 2 (i5)", None, "CD reference i5 site is 8 nt; actual acquisition, orientation and member unassigned.", "extract"),
)
SEQUENCING_ENDING = ("The v00 guide defines Read 1 as antisense and Read 2 as sense. It delegates cluster generation to the selected instrument guide. "
                     "No acquisition mask or universal index-2 primer is assigned here; source-compatible adapter sites do not prove an acquired run.")
SEQUENCING_UNAVAILABLE = "The model has not selected an instrument/run or actual physical sequencing-primer set. This is not a publication-absence claim."


def state(name, *parts):
    return Construct([role(p) for p in parts], name=name)


def cd_library(*, inferred=True):
    """Positive v22 CD templates; historical reagent assignment is conditional."""
    return Construct([
        Segment("P5", il.P5, "p5", inferred=inferred),
        role("i5 CD member", length_bp=8, inferred=inferred,
             molecular_feature=feature("sample_i5", "sample_index", "unknown")),
        Segment("Read 1 adapter arm", il.TRUSEQ_READ1, "r1", inferred=inferred),
        role("variable antisense-derived cDNA interval"),
        Segment("insert dA junction", "A"),
        Segment("i7 adapter arm", il.INDEX1_PRIMER, "r2", inferred=inferred),
        role("i7 CD member", length_bp=8, inferred=inferred,
             molecular_feature=feature("sample_i7", "sample_index", "unknown")),
        Segment("P7 reverse complement", il.P7_RC, "p7", inferred=inferred),
    ], name="conditional RNA CD library")


FINAL_LIBRARY = cd_library()
FINAL_CAPTION = ("INFERRED — RNA CD reference assembled from the official v22 adapter templates, with a variable cDNA interval. "
                 "The selected historical adapter revision and individual indexes remain unassigned; no fixed biological insert length is implied.")


def final_rows():
    unknown = tuple(s.name + "'" for s in FINAL_LIBRARY if s.is_role_token())
    sc = Scene.duplex(list(FINAL_LIBRARY), label=FINAL_LIBRARY.name, unpaired=unknown)
    sc.junction("top", "Read 1 adapter arm", "variable antisense-derived cDNA interval", "ligation boundary")
    sc.junction("top", "insert dA junction", "i7 adapter arm", "ligation boundary")
    return sc.rows()


def protocol_graph():
    g = ProtocolGraph(TITLE, authority=SOURCE)
    reviewed = []
    def reagent(mid, label, kind):
        m = g.input(mid, label, kind); reviewed.append(mid); return m
    def product(mid, label, *, drawing=None, kind=MaterialKind.MOLECULAR):
        return g.product(mid, label, kind, state=drawing)
    def react(mid, label, kind, inputs, outputs, **kwargs):
        reviewed.append(mid); return g.react(mid, label, kind, inputs, outputs, **kwargs)
    def cleanup(mid, template, label, *, library=False):
        beads = reagent(mid + "_beads", "AMPure XP beads / wash and elution", MaterialKind.SELECTION)
        clean = product(mid + "_clean", label, drawing=template.state,
                        kind=MaterialKind.LIBRARY if library else MaterialKind.MOLECULAR)
        if library: reviewed.append(clean.id)
        waste = product(mid + "_waste", "nonretained supernatant and washes; composition unassigned", kind=MaterialKind.WASTE)
        react(mid, "AMPure purification: " + label, StepKind.CLEANUP,
              ((template, InputRole.TEMPLATE), (beads, InputRole.SELECTION)), (clean, waste))
        g.discard(waste, "The source discards nonretained supernatants/washes and elutes retained material; exact excluded species and complete removal are unassigned.")
        return clean

    rna = g.input("rna", "total RNA", MaterialKind.BIOLOGICAL,
                  state=state("total RNA", "variable RNA pool"), primary=True)
    reviewed.append("rna")
    beads = reagent("polyt_beads", "oligo-dT magnetic beads: two capture rounds", MaterialKind.SELECTION)
    selected = product("selected", "bead-enriched poly(A) RNA", drawing=state("selected RNA", "variable poly(A)-containing RNA"))
    rejected = product("selection_waste", "discarded capture supernatants and washes", kind=MaterialKind.WASTE)
    react("select_polya", "Capture poly(A) RNA twice", StepKind.AFFINITY_CAPTURE,
          ((rna, InputRole.TEMPLATE), (beads, InputRole.SELECTION)), (selected, rejected))
    g.discard(rejected, "Two captures discard supernatants/washes. Their composition is not equated to all and only non-polyadenylated RNA.")
    condition = reagent("fragmentation", "FPF: divalent-cation fragmentation and priming, Elution 2-Frag-Prime 94 °C 8 min", MaterialKind.CONDITION)
    fragmented = product("fragmented", "fragmented selected RNA", drawing=state("fragmented RNA", "variable selected RNA fragments"))
    react("fragment", "Fragment and prime during the second elution", StepKind.FRAGMENTATION,
          ((selected, InputRole.TEMPLATE), (condition, InputRole.CONDITION)), (fragmented,))
    primer = reagent("rt_primer", "random hexamer first-strand primers; ordered pool unassigned", MaterialKind.OLIGO)
    rt = reagent("rt", "SuperScript II with FSA / Actinomycin D; selected listed option; 25 °C 10 min, 42 °C 15 min, 70 °C 15 min", MaterialKind.CATALYST)
    first = product("first", "first-strand cDNA / RNA mixture", drawing=state("first-strand mixture", "antisense first-strand cDNA with RNA template"))
    react("reverse_transcribe", "Random-prime first-strand cDNA synthesis", StepKind.REVERSE_TRANSCRIPTION,
          ((fragmented, InputRole.TEMPLATE), (primer, InputRole.STRUCTURAL), (rt, InputRole.CATALYST)), (first,))
    smm = reagent("smm", "SMM: RNase H, DNA polymerase I and dUTP marking / end repair; 16 °C 1 h", MaterialKind.CATALYST)
    marked = product("marked", "blunt-ended, dUTP-marked cDNA duplex", drawing=state("marked duplex", "first strand plus dUTP-marked second strand; positions unassigned"))
    degraded = product("removed_rna", "RNA-template degradation products; subsequent recovery unassigned")
    react("second_strand", "Remove RNA and synthesize the dUTP-marked blunt second strand", StepKind.EXTENSION,
          ((first, InputRole.TEMPLATE), (smm, InputRole.CATALYST)), (marked, degraded))
    g.mark_terminal(Terminal(degraded.id, Disposition.UNRESOLVED,
        "RNA removal occurs during SMM synthesis; downstream cleanup is shown, but the specific degradation products and their individual recovery are unassigned."))
    clean = cleanup("clean_cdna", marked, "cleaned blunt cDNA")
    atl = reagent("atl", "ATL A-Tailing Mix; 37 °C 30 min, then 70 °C 5 min", MaterialKind.CATALYST)
    ends = EndPreparedDuplex(Construct([role("variable marked cDNA interval")]),
                             five_prime_phosphate=True, three_prime_overhang="A", overhang_inferred=False)
    tailed = product("tailed", "dA-tailed cDNA", drawing=dA_tailed_scene(ends))
    react("tail", "Add one 3′ dA to blunt cDNA ends", StepKind.TAILING,
          ((clean, InputRole.TEMPLATE), (atl, InputRole.CATALYST)), (tailed,))
    adapter_scene = Scene()
    adapter_scene.strand("adapter", [role("indexed adapter duplex core"), Segment("3-prime T", "T")], label="source-known 3′-T adapter end")
    adapter_scene.overhang("adapter", "3-prime T", polarity="3-prime", owner="adapter")
    adapter = reagent("adapter", "RNA CD indexed adapters; compatible reference family, historical revision/member unassigned", MaterialKind.ADAPTER)
    ligase = reagent("ligase", "LIG Ligation Mix, 30 °C 10 min, followed by STL stop", MaterialKind.CATALYST)
    adapted = product("adapted", "indexed-adapter-ligated marked cDNA", drawing=state("ligated cDNA", "indexed adapter", "marked cDNA interval", "opposite indexed adapter"))
    react("ligate", "Ligate indexed adapters and stop ligation", StepKind.LIGATION,
          ((tailed, InputRole.TEMPLATE), (adapter, InputRole.STRUCTURAL), (ligase, InputRole.CATALYST)), (adapted,),
          geometry=AdapterLigationGeometry.known("cohesive", tailed.state, adapter_scene,
              statement="The v00 guide explicitly gives complementary single 3′-A substrate and 3′-T adapter ends. Exact historical adapter-core bases remain unassigned."))
    clean = cleanup("clean_ligation_1", adapted, "first-round purified ligated cDNA")
    clean = cleanup("clean_ligation_2", clean, "second-round purified ligated cDNA")
    ppc = reagent("ppc", "PPC forward primer; ordered strand unassigned", MaterialKind.OLIGO)
    ppc_reverse = reagent("ppc_reverse", "PPC reverse primer; ordered strand unassigned", MaterialKind.OLIGO)
    pmm = reagent("pmm", "PMM PCR Master Mix: polymerase does not incorporate past dUTP; 98 °C 30 s, cycles of 98 °C 10 s / 60 °C 30 s / 72 °C 30 s, 72 °C 5 min", MaterialKind.CATALYST)
    amplified = product("amplified", "PCR-enriched library", drawing=FINAL_LIBRARY)
    react("pcr", f"INFERRED — selectively enrich the CD reference in {PCR_CYCLES} PCR cycles", StepKind.PCR,
          ((clean, InputRole.TEMPLATE), (ppc, InputRole.STRUCTURAL), (ppc_reverse, InputRole.STRUCTURAL), (pmm, InputRole.CATALYST)), (amplified,),
          geometry=PCRGeometry.unavailable(clean.state,
              reason="The model has not assigned the ordered PPC strands or their exact sites on this ligated historical adapter. The v22 CD endpoint is a conditional source-compatible reference, not evidence that PCR installs its indexes."))
    library = cleanup("clean_pcr", amplified, "purified stranded RNA CD reference library", library=True)
    g.sequence(library, SequencingPlan(library=FINAL_LIBRARY, primers=(),
        limitation="RNA CD reference sites are conditional; no actual platform, physical primers or acquisition mask selected."))
    g.annotate(Highlight("second_strand", "Marking controls amplification, not immediate physical strand removal",
        "SMM removes the RNA template and produces blunt-ended cDNA; its dUTP-marked second strand is quenched during later PCR. "
        "The model assigns no enzyme-based destruction of the entire marked strand and no complete background elimination."))
    g.annotate(Highlight("adapter", "Index family is a separate choice from HS/LS handling",
        "The overview illustrates LS single-index and describes HS dual-index libraries; the purchasable components separately include RNA single, CD and UD families. "
        "This page selects HS handling with a compatible CD reference. It does not assign an individual D index, force UD pairing, or infer acquisition from adapter compatibility."))
    g.annotate(Highlight("pcr", "Ligation installs the indexes",
        "PPC anneals to the adapter ends and PCR enriches fragments bearing adapters at both ends. "
        "The source says zero-adapter fragments do not hybridize to the flow cell and one-adapter fragments cannot form clusters; "
        "their remaining abundance and molecular fate are unassigned. The two ligation cleanups are physical recovery steps, not proof that every off-product is gone."))
    g.declare_reviewed_scope(reviewed, SOURCE)
    return g.require_publishable()
