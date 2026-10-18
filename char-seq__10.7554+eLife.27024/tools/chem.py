"""ChAR-seq RNA–DNA bridge chemistry."""
from adapter_ligation import AdapterLigationGeometry
import illumina as il
import seqprimers as sp
from protocol_model import role
from chemdraw import Construct, Scene, Segment, feature, revcomp
from endprep import dA_tailed_scene, repair_and_dA_tail
from pcr import PCRGeometry
from protocol_graph import (Disposition, Terminal, InputRole, MaterialKind, ProtocolGraph,
                            SequencingPlan, StepKind)


TITLE = "ChAR-seq — chromatin-associated RNA sequencing"
NOTES = "01_char-seq.html"
SOURCE = ('Defining source: <a href="https://doi.org/10.7554/eLife.27024">'
          'Bell et al., <i>eLife</i> (2018)</a>, Figure 1—figure supplement 1 and caption.')
SUMMARY = ("A pre-adenylated, desthiobiotinylated DNA bridge first captures an RNA 3′ end. "
           "DpnII then opens the bridge's DNA-facing end and cuts nearby genomic DNA, "
           "creating a directional RNA–DNA chimera.")
CAVEAT = ("<b>Published bridge variants:</b> <a href=\"https://cdn.elifesciences.org/articles/27024/elife-27024-fig1-figsupp1-v2.jpg\">"
          "Bell et al. 2018 Figure 1—figure supplement 1</a> "
          "draws 5′-App-AAAAAAACCGGCGTCCAAGGATCTAG-3′-Sp3 paired with "
          "3′-TGGCCGCAGGTTCCTAGATCCTAG-5′, with internal desthiobiotin-TEG and no UMI. "
          "<a href=\"https://pmc.ncbi.nlm.nih.gov/articles/PMC7607860/\">Limouse et al. 2020, Methods in Molecular Biology §2.3</a> instead prints top "
          "/5rApp/AANNNAAACCGGCGTCCAAGGATCTTTAATTAAGTCGCAG/3SpC3/ and bottom "
          "/5Phos/GATCTGCGACTTAATTAAAGATCCTTGGACGCCGG/iBiodT/T, adding NNN and a PacI site. "
          "This page models the defining 2018 bridge; "
          "the figure does not resolve the exact TEG attachment bond. "
          "The 2018 Methods name a NEBNext USER hairpin and TruSeq-compatible indexing, "
          "but give no kit revision or selected index; the sequencing shell is inferred from "
          "the NEBNext family, not assigned to a particular E7335 or E7600 member.")

BRIDGE_TOP = "AAAAAAACCGGCGTCCAAGGATCTAG"
# Ordered 5′→3′; the original figure draws this strand 3′→5′.
BRIDGE_BOTTOM = "GATCCTAGATCCTTGGACGCCGGT"
# Published graphic, not a preprint exception. Both ordering strings are pinned by
# tools/selftest.py; the source guard cannot extract bases from this raster original.
SOURCE_OLIGO_EXCEPTIONS = {
    "bridge": "Visually verified in published eLife Figure 1—figure supplement 1, "
              "elife-27024-fig1-figsupp1.jpg inside the PMC5962340 supplementary ZIP; "
              "5′→3′ top/bottom transcriptions pinned in tools/selftest.py. "
              "This is an image-only transcription exemption, not a text-match pass."
}
TOP_TAIL = "AAAAAA"  # six unpaired bases; the seventh A enters the duplex
TOP_CORE = "ACCGGCGTCCAAGGATCTAG"
BOTTOM_DPNII_TAIL = "GATC"
BOTTOM_CORE = "CTAGATCCTTGGACGCCGGT"
RETAINED_BRIDGE = "ACCGGCGTCCAAG"  # upstream of the DpnII cut

# Conditional TruSeq-compatible reference only; the 2018 kit member is unassigned.
FINAL_LIBRARY = Construct([
    Segment("P5", il.TRUSEQ_P5_FULL[:-len(il.TRUSEQ_READ1)], "p5", inferred=True),
    Segment("Read 1 arm", il.TRUSEQ_READ1, "r1", inferred=True),
    Segment("RNA-derived cDNA", "R" * 38, placeholder=True),
    Segment("bridge RNA-ligation tail", TOP_TAIL),
    Segment("retained bridge through DpnII junction", RETAINED_BRIDGE, "me"),
    Segment("DpnII genomic junction", "GATC", "me"),
    Segment("nearby genomic DNA", "D" * 42, placeholder=True),
    Segment("source dA junction", "A"),
    Segment("conditional Index 1 / Read 2 arm", il.INDEX1_PRIMER, "r2", inferred=True),
    role("i7 index; length/member unassigned", tag="cbc", inferred=True,
         molecular_feature=feature("sample_index_i7", "sample_index", "unknown")),
    Segment("P7 reverse complement", il.P7_RC, "p7", inferred=True),
], name="Conditional ChAR-seq TruSeq-compatible reference")
# Only the source-positive biological acquisition is declared here. Omitting
# an index reference primer is not a claim that no index read was acquired.
SEQ_PRIMERS = (sp.TRUSEQ["R1"],)
_problems = sp.verify(FINAL_LIBRARY, SEQ_PRIMERS, required_roles=("Read 1",))
if _problems:
    raise ValueError("Conditional ChAR-seq shell: " + "; ".join(_problems))

FINAL_CAPTION = ("INFERRED — The source-defined bridge separates the RNA-derived cDNA "
                 "from nearby genomic DNA. The defining experiment used single-end 152-base "
                 "sequencing; the TruSeq-compatible shell depicts an unverified NEBNext "
                 "index-kit member, and the opposite site is not a reported second read.")
READ_LENGTHS = {"Read 1": 152}
SHOW_READ_CONTENT = False
READ_SPECS = (
    sp.read_spec('Read 1', 152, 'RNA-derived cDNA, bridge junction and nearby genomic DNA when reached', 'include'),
    sp.read_spec('Index 1 (i7)', 'unassigned', 'Selected index and acquisition/cycles unassigned', 'undetermined'),
)
READ_SPECS_NOTE = 'The defining experiment reports single-end152; the conditional outer shell does not establish actual index acquisition.'
SEQUENCING_INTRO = ("The source names a NEBNext USER-opened TruSeq-compatible library. "
                    "Its exact kit revision and index member are unverified; the paper reports "
                    "single-end 152-base sequencing. Actual index acquisition and its cycle count remain unassigned.")
SEQUENCING_UNAVAILABLE = ("No second biological read was collected in the defining experiment, "
                          "actual index acquisition/cycles remain unassigned, and the vendor hairpin adapter's exact bases are not asserted here.")


def oligos():
    return (
        f"<b>2018 bridge top</b> — 5′-App-{BRIDGE_TOP}-3′-Sp3",
        f"<b>2018 bridge bottom</b> — 5′-{BRIDGE_BOTTOM}-3′ (internal desthiobiotin-TEG; precise attachment bond unresolved)",
    )


def bridge_scene() -> Scene:
    top = [Segment("5-prime RNA-ligation tail", TOP_TAIL),
           Segment("paired bridge top", TOP_CORE, "me")]
    bottom = [Segment("5-prime DpnII-side tail", BOTTOM_DPNII_TAIL, "me"),
              Segment("paired bridge bottom", BOTTOM_CORE, "me")]
    sc = Scene()
    sc.strand("bridge top", top, label="5′-adenylated top", mod5="App",
              mod3="Sp3")
    sc.anneal("bridge bottom", bottom, to="bridge top",
              pair=("paired bridge bottom", "paired bridge top"),
              unpaired=("5-prime DpnII-side tail", "5-prime RNA-ligation tail"),
              label="bottom strand with internal desthiobiotin-TEG")
    sc.mark("bridge bottom", "paired bridge bottom", "internal desthiobiotin-TEG; exact bond unresolved")
    sc.ligation_end("bridge top", "5-prime RNA-ligation tail", kind="single-strand",
                    polarity="5-prime", owner="adapter",
                    text="5′-App end ligates an RNA 3′-OH")
    sc.labels("bridge top")
    return sc


def _rna_ligation_geometry() -> AdapterLigationGeometry:
    sc = bridge_scene()
    sc.strand("RNA substrate", [Segment("chromatin RNA", "R" * 28,
                                        placeholder=True)], label="fragmented chromatin RNA")
    sc.ligation_end("RNA substrate", "chromatin RNA", kind="single-strand",
                    polarity="3-prime", owner="substrate", text="RNA 3′-OH")
    return AdapterLigationGeometry.known(
        "single-strand", sc,
        statement="The printed 5′-App bridge top is selectively ligated to fragmented RNA 3′-OH ends by truncated T4 RNA ligase 2.")


def _rna_joined_scene(*, copied=False, opened=False, genome=True) -> Scene:
    """RNA ligation and Bst copying retain distinct RNA/DNA strand identities."""
    sc = Scene()
    core = RETAINED_BRIDGE if opened else TOP_CORE
    sc.strand("RNA–bridge", [Segment("chromatin RNA", "R" * 28, placeholder=True),
        Segment("bridge RNA-ligation tail", TOP_TAIL),
        Segment("bridge remainder", core, "me")],
        label="RNA covalently joined to bridge DNA", mod3="" if opened else "Sp3")
    sc.junction("RNA–bridge", "chromatin RNA", "bridge RNA-ligation tail", "T4 Rnl2(tr) ligation")
    bottom = [Segment("DpnII overhang" if opened else "original bottom 5-prime tail", "GATC", "me"),
              Segment("bridge bottom core", revcomp(core), "me")]
    if copied:
        bottom += [Segment("copied RNA-ligation tail", revcomp(TOP_TAIL)),
                   Segment("RNA-derived cDNA", "r" * 28, placeholder=True)]
    sc.anneal("bridge copy", bottom, to="RNA–bridge",
        pair=("bridge bottom core", "bridge remainder"),
        unpaired=(bottom[0].name,), label="Bst-extended DNA" if copied else "original DNA bridge primer")
    sc.mark("bridge copy", "bridge bottom core", "internal desthiobiotin-TEG; attachment bond unassigned")
    if opened:
        sc.overhang("bridge copy", "DpnII overhang", polarity="5-prime", owner="adapter", text="DpnII-compatible bridge end")
    if genome:
        sc.strand("genomic DNA", [role("nearby chromatin DNA; separate from RNA–bridge")], label="genomic DNA population")
    return sc


def _genomic_end(sc=None):
    sc = Scene() if sc is None else sc
    sc.strand("genomic top", [Segment("genomic DpnII overhang", "GATC", "me"),
        Segment("genomic DNA body", "D" * 30, placeholder=True)], label="DpnII-cut genomic DNA")
    sc.anneal("genomic bottom", [Segment("genomic DNA complement", "d" * 30, placeholder=True)],
        to="genomic top", pair=("genomic DNA complement", "genomic DNA body"),
        label="genomic complementary strand")
    sc.overhang("genomic top", "genomic DpnII overhang", polarity="5-prime", owner="substrate", text="genomic 5-prime GATC end")
    return sc


def _dpnii_geometry() -> AdapterLigationGeometry:
    sc = _genomic_end(_rna_joined_scene(copied=True, opened=True, genome=False))
    return AdapterLigationGeometry.known(
        "cohesive", sc,
        statement="DpnII removes the bridge's 3′ C3-blocked distal piece and produces a GATC-compatible bridge end while cutting nearby genomic DNA to the same cohesive-end class.")


def _contact_scene(*, rna=False) -> Scene:
    sc = Scene.duplex([
        Segment("retained RNA" if rna else "RNA-derived cDNA", "R" * 28, placeholder=True),
        Segment("bridge RNA-ligation tail", TOP_TAIL),
        Segment("retained desthiobiotin bridge", RETAINED_BRIDGE, "me"),
        Segment("DpnII genomic junction", "GATC", "me"),
        Segment("nearby genomic DNA", "D" * 30, placeholder=True),
    ], label="RNA/DNA strand paired to DNA; RNA persists" if rna else "cDNA–bridge–genomic DNA chimera")
    if rna:
        sc.note("top", "Mixed-polymer strand: RNA persists on the RNA side; bridge and genomic portions are DNA.")
        sc.note("bottom", "DNA strand: Bst-generated cDNA across the RNA side, with bridge/genomic DNA at the other end.")
    sc.junction("top", "retained RNA" if rna else "RNA-derived cDNA", "bridge RNA-ligation tail", "RNA-side bridge")
    sc.junction("top", "retained desthiobiotin bridge", "DpnII genomic junction",
                "DpnII cohesive ligation")
    sc.mark("top", "retained desthiobiotin bridge", "streptavidin capture")
    sc.labels("top")
    return sc


def _opened_library_state():
    sc = Scene()
    sc.strand("USER-opened library population", [role("actual USER-opened adapter/contact precursor; strand ends unassigned")],
              label="topology/adapter ends unassigned; not the completed PCR duplex")
    return sc


def _unknown_end_duplex(label):
    """An amplified duplex without substituting conditional reference arms."""
    scene = Scene()
    scene.strand("top", [role("indexed contact amplicon strand; terminal sequences unassigned")], label=label + " strand 1")
    scene.strand("bottom", [role("complementary amplicon strand; terminal sequences unassigned")], label=label + " strand 2", rev=True)
    scene.note("top", "Source-described PCR duplex; exact primer members, ends and base-level pairing geometry remain unassigned.")
    return scene


def _distal_bridge_fragment():
    """The DpnII-detached source bridge piece retains the original C3 block."""
    distal_top = TOP_CORE[len(RETAINED_BRIDGE):]
    retained_bottom = "GATC" + revcomp(RETAINED_BRIDGE)
    if not BRIDGE_BOTTOM.endswith(retained_bottom):
        raise ValueError("source bridge cannot be partitioned at the declared DpnII cut")
    distal_bottom = BRIDGE_BOTTOM[:-len(retained_bottom)]
    scene = Scene()
    scene.strand("distal top", [Segment("new DpnII 5-prime tail", distal_top[:4]),
                                Segment("distal top paired core", distal_top[4:])],
                 label="detached C3-blocked bridge piece", mod3="Sp3")
    scene.anneal("distal bottom", [Segment("original bottom 5-prime tail", distal_bottom[:4]),
                                   Segment("distal bottom paired core", distal_bottom[4:])],
                 to="distal top", pair=("distal bottom paired core", "distal top paired core"),
                 unpaired=("new DpnII 5-prime tail", "original bottom 5-prime tail"))
    scene.labels("distal top")
    return scene


def protocol_graph() -> ProtocolGraph:
    source = "Defining source: https://doi.org/10.7554/eLife.27024 Figure 1—figure supplement 1 (published image); extended protocol Steps 1–27"
    hairpin_source = source + "; NEBNext kit family named in 2018 Methods; E7335/E7600 manuals are candidate family references, not version or index assignments"
    graph = ProtocolGraph("ChAR-seq", authority=source)
    contact = graph.input("contact", "fixed chromatin-associated RNA–DNA contact",
                          MaterialKind.BIOLOGICAL, primary=True)
    fragmentation = graph.input("rna_fragmentation", "RNA-fragmentation condition",
                                MaterialKind.CONDITION)
    fragmented = graph.product("fragmented_contact", "fragmented RNA beside chromatin DNA",
                               MaterialKind.MOLECULAR)
    graph.react("fragment_rna", "Fragment chromatin-associated RNA", StepKind.FRAGMENTATION,
                ((contact, InputRole.TEMPLATE), (fragmentation, InputRole.CONDITION)),
                (fragmented,))
    bridge = graph.input("bridge", "2018 pre-annealed desthiobiotinylated ChAR-seq bridge",
                         MaterialKind.ADAPTER, state=bridge_scene())
    rna_ligase = graph.input("rna_ligase", "truncated T4 RNA ligase 2 (R55K K227Q, trKQ) in 25% PEG",
                             MaterialKind.CATALYST)
    rna_joined = graph.product("rna_joined", "RNA-ligated bridge",
                               MaterialKind.MOLECULAR, state=_rna_joined_scene())
    graph.react("rna_ligation", "Ligate RNA 3′ ends to the bridge", StepKind.LIGATION,
                ((fragmented, InputRole.TEMPLATE), (bridge, InputRole.STRUCTURAL),
                 (rna_ligase, InputRole.CATALYST)), (rna_joined,),
                geometry=_rna_ligation_geometry())
    bst = graph.input("bst", "Bst 3.0 reverse transcriptase", MaterialKind.CATALYST)
    stabilized = graph.product("stabilized_hybrid", "RNA/cDNA–bridge hybrid",
                               MaterialKind.MOLECULAR, state=_rna_joined_scene(copied=True))
    graph.react("reverse_transcription", "Extend the bridge across RNA",
                StepKind.REVERSE_TRANSCRIPTION,
                ((rna_joined, InputRole.TEMPLATE), (rna_joined, InputRole.SELF_PRIMER), (bst, InputRole.CATALYST)), (stabilized,))
    dpnii = graph.input("dpnii", "DpnII", MaterialKind.CATALYST)
    opened_bridge = graph.product("opened_bridge", "DpnII-opened RNA/cDNA bridge",
                                  MaterialKind.MOLECULAR, state=_rna_joined_scene(copied=True, opened=True, genome=False))
    genomic_end = graph.product("genomic_end", "nearby DpnII-cut genomic DNA",
                                MaterialKind.MOLECULAR, state=_genomic_end())
    detached = graph.product("detached_bridge_piece", "DpnII-detached distal bridge fragment retaining the C3 block",
                             MaterialKind.MOLECULAR, state=_distal_bridge_fragment())
    graph.react("dpnii_digest", "Open the bridge and cut nearby genomic DNA with DpnII",
                StepKind.DIGESTION,
                ((stabilized, InputRole.TEMPLATE), (dpnii, InputRole.CATALYST)),
                (opened_bridge, genomic_end, detached))
    graph.mark_terminal(Terminal(detached.id, Disposition.UNRESOLVED,
        "DpnII cleavage releases this distal bridge piece; no separate physical removal is assigned at this step."))
    dna_ligase = graph.input("dna_ligase", "T4 DNA ligase", MaterialKind.CATALYST)
    chimera = graph.product("contact_chimera", "RNA–bridge–DNA contact chimera",
                            MaterialKind.MOLECULAR, state=_contact_scene(rna=True))
    graph.react("contact_ligation", "Ligate the opened bridge to nearby genomic DNA",
                StepKind.LIGATION,
                ((opened_bridge, InputRole.TEMPLATE), (genomic_end, InputRole.STRUCTURAL),
                 (dna_ligase, InputRole.CATALYST)), (chimera,), geometry=_dpnii_geometry())
    dna_pol = graph.input("dna_polymerase", "E. coli DNA polymerase I with limiting RNase H",
                          MaterialKind.CATALYST)
    duplex = graph.product("contact_duplex", "double-stranded contact chimera",
                           MaterialKind.MOLECULAR, state=_contact_scene())
    graph.react("second_strand", "Make the second DNA strand", StepKind.EXTENSION,
                ((chimera, InputRole.TEMPLATE), (dna_pol, InputRole.CATALYST)), (duplex,))
    protk = graph.input("proteinase_k", "SDS/NaCl with proteinase K, 70 C overnight",
                        MaterialKind.CATALYST)
    reversed_dna = graph.product("reversed_dna", "crosslink-reversed, precipitated contact DNA",
                                 MaterialKind.MOLECULAR, state=_contact_scene())
    graph.react("reverse_crosslinks", "Reverse crosslinks and precipitate the DNA",
                StepKind.DEPROTEINIZATION,
                ((duplex, InputRole.TEMPLATE), (protk, InputRole.CATALYST)), (reversed_dna,))
    shear = graph.input("shearing_condition", "Covaris shearing to ~200 bp", MaterialKind.CONDITION)
    sheared = graph.product("sheared_contacts", "sheared bridge-bearing contacts",
                            MaterialKind.MOLECULAR, state=_contact_scene())
    graph.react("shear", "Shear the contact duplexes", StepKind.FRAGMENTATION,
                ((reversed_dna, InputRole.TEMPLATE), (shear, InputRole.CONDITION)), (sheared,))
    streptavidin = graph.input("streptavidin", "Dynabeads MyOne Streptavidin T1",
                               MaterialKind.SELECTION)
    captured = graph.product("captured_contacts", "purified biotin-bearing contacts",
                             MaterialKind.MOLECULAR, state=_contact_scene())
    uncaptured = graph.product("uncaptured_contacts", "non-bridge-bearing fragments",
                               MaterialKind.WASTE)
    graph.react("capture", "Enrich biotin-bearing contacts", StepKind.AFFINITY_CAPTURE,
                ((sheared, InputRole.TEMPLATE), (streptavidin, InputRole.SELECTION)),
                (captured, uncaptured))
    graph.discard(uncaptured, "fragments lacking the internal biotin bridge are washed away")
    end_repair = graph.input("end_repair_enzymes", "NEBNext End Prep enzyme mix (on bead)",
                             MaterialKind.CATALYST)
    repaired = graph.product("repaired_contacts", "end-repaired contact duplexes",
                             MaterialKind.MOLECULAR)
    graph.react("end_repair", "End-repair the captured contacts", StepKind.END_REPAIR,
                ((captured, InputRole.TEMPLATE), (end_repair, InputRole.CATALYST)),
                (repaired,))
    tailing = graph.input("tailing_polymerase", "NEBNext End Prep enzyme mix (dA tailing)",
                          MaterialKind.CATALYST)
    insert = Construct([
        Segment("RNA-derived cDNA", "R" * 28, placeholder=True),
        Segment("bridge RNA-ligation tail", TOP_TAIL),
        Segment("retained bridge", RETAINED_BRIDGE, "me"),
        Segment("DpnII genomic junction", "GATC", "me"),
        Segment("nearby genomic DNA", "D" * 30, placeholder=True),
    ], name="sheared bridge-bearing contact")
    tailed_state = repair_and_dA_tail(insert)
    tailed = graph.product("tailed_contacts", "dA-tailed contact duplexes",
                           MaterialKind.MOLECULAR,
                           state=dA_tailed_scene(tailed_state, label="dA-tailed contact duplex"))
    graph.react("da_tail", "Add one 3′ dA to each end", StepKind.TAILING,
                ((repaired, InputRole.TEMPLATE), (tailing, InputRole.CATALYST)), (tailed,))
    hairpin = graph.input("hairpin_adapter", "NEBNext USER hairpin adapter", MaterialKind.ADAPTER,
                          source=hairpin_source, inferred=True)
    hairpin_ligated = graph.product("hairpin_ligated", "hairpin-adapted contact library",
                                    MaterialKind.MOLECULAR, source=hairpin_source)
    graph.react("hairpin_ligation", "Ligate NEBNext hairpin adapters", StepKind.LIGATION,
                ((tailed, InputRole.TEMPLATE), (hairpin, InputRole.STRUCTURAL),
                 (dna_ligase, InputRole.CATALYST)), (hairpin_ligated,),
                geometry=AdapterLigationGeometry.unavailable(
                    "The exact two-ended hairpin product is topological and cannot be flattened into this compact single-junction Scene."),
                source=hairpin_source)
    user = graph.input("user_enzyme", "USER enzyme", MaterialKind.CATALYST)
    opened = graph.product("opened_library", "USER-opened adapter library",
                           MaterialKind.MOLECULAR, state=_opened_library_state(), source=hairpin_source)
    graph.react("user_opening", "Open the hairpins with USER", StepKind.DIGESTION,
                ((hairpin_ligated, InputRole.TEMPLATE), (user, InputRole.CATALYST)),
                (opened,), source=hairpin_source)
    post_user_wash = graph.input("post_user_wash", "Step22: twice750µl1×TWB; first wash warmed50°C2min", MaterialKind.SELECTION)
    washed_opened = graph.product("washed_opened_library", "washed on-bead USER-opened contacts", MaterialKind.MOLECULAR, state=opened.state)
    user_waste = graph.product("post_user_waste", "post-ligation/USER supernatant and TWB washes", MaterialKind.WASTE)
    graph.react("wash_user_opened", "Wash the USER-treated on-bead library", StepKind.CLEANUP,
                ((opened, InputRole.TEMPLATE), (post_user_wash, InputRole.SELECTION)), (washed_opened, user_waste))
    graph.discard(user_waste, "Step22 explicitly discards the supernatant and wash fractions")
    p5 = graph.input("p5_primer", "NEBNext Universal PCR primer", MaterialKind.OLIGO,
                     source=hairpin_source, inferred=True)
    p7 = graph.input("p7_primer", "NEBNext indexed P7 PCR primer", MaterialKind.OLIGO,
                     source=hairpin_source, inferred=True)
    pcr = graph.input("pcr_polymerase", "NEBNext High-Fidelity PCR master mix",
                      MaterialKind.CATALYST)
    pcr1 = graph.product("pcr1", "on-bead eight-cycle indexed amplification product",
                         MaterialKind.MOLECULAR, state=_unknown_end_duplex("eight-cycle PCR1 amplicon"),
                         source=hairpin_source)
    graph.react("index_pcr", "Amplify on bead with the Universal and Indexing primers (8 cycles)",
                StepKind.PCR,
                ((washed_opened, InputRole.TEMPLATE), (p5, InputRole.STRUCTURAL),
                 (p7, InputRole.STRUCTURAL), (pcr, InputRole.CATALYST)), (pcr1,),
                geometry=PCRGeometry.unavailable(
                    washed_opened.state,
                    reason="The actual USER-opened precursor has unassigned adapter ends. The 2018 source names NEBNext Indexing Primers for Illumina but not the exact kit revision, primer strands, or selected index; E7335/E7600 manuals establish candidate family sequences only."),
                source=hairpin_source)
    bead_separation = graph.input("pcr1_bead_separation", "Step23: magnetically pellet original streptavidin beads after PCR1", MaterialKind.SELECTION)
    pcr1_solution = graph.product("pcr1_supernatant", "soluble eight-cycle PCR1 amplicons", MaterialKind.MOLECULAR, state=pcr1.state)
    original_beads = graph.product("pcr1_beads", "original streptavidin beads and remaining bound material", MaterialKind.MOLECULAR)
    graph.react("recover_pcr1_supernatant", "Collect PCR1 supernatant after magnetic separation", StepKind.CLEANUP,
                ((pcr1, InputRole.TEMPLATE), (bead_separation, InputRole.SELECTION)), (pcr1_solution, original_beads), changes_molecule=False)
    graph.mark_terminal(Terminal(original_beads.id, Disposition.UNRESOLVED, "Step24 processes PCR supernatant; no later disposition of the original magnetic beads is assigned."))
    first_ampure = graph.input("first_ampure", "Step24: equal-volume1:1 AMPureXP;70%EtOHwash;30µl10mMTrispH8eluate", MaterialKind.SELECTION)
    cleaned1 = graph.product("pcr1_cleaned", "Step24 eluate of eight-cycle amplicons", MaterialKind.MOLECULAR, state=pcr1.state)
    saved_first_supernatant = graph.product("pcr1_saved_supernatant", "Step24 AMPure supernatant retained until completed-library QC", MaterialKind.MOLECULAR)
    first_wash = graph.product("pcr1_ethanol_wash", "Step24 70%ethanol wash supernatant", MaterialKind.WASTE)
    graph.react("first_ampure_cleanup", "Clean PCR1 supernatant with1:1 AMPure and elute30µl", StepKind.CLEANUP,
                ((pcr1_solution, InputRole.TEMPLATE), (first_ampure, InputRole.SELECTION)),
                (cleaned1, saved_first_supernatant, first_wash))
    graph.mark_terminal(Terminal(saved_first_supernatant.id, Disposition.STORED, "Step24 explicitly saves this supernatant until library validation by qPCR/Bioanalyzer, for troubleshooting."))
    graph.discard(first_wash, "Step24 removes the ethanol wash supernatant before elution")
    qpcr_aliquot = graph.product("qpcr_aliquot", "5µl Step24 eluate for cycle-determination qPCR", MaterialKind.MOLECULAR, state=cleaned1.state)
    production_aliquot = graph.product("production_aliquot", "remaining25µl Step24 eluate for production PCR2", MaterialKind.MOLECULAR, state=cleaned1.state)
    graph.react("split_pcr1_eluate", "Allocate5µl to qPCR and25µl to production amplification", StepKind.SPLIT,
                ((cleaned1, InputRole.TEMPLATE),), (qpcr_aliquot, production_aliquot), changes_molecule=False)
    qpcr_mix = graph.input("cycle_qpcr_mix", "Step25: reserved10µl Step23PCRmix plus0.5µl33×SYBRGreen", MaterialKind.CONDITION)
    qpcr_product = graph.product("qpcr_product", "qPCR-amplified analytical5µl aliquot", MaterialKind.QC, state=_unknown_end_duplex("cycle-determination qPCR product"))
    graph.react("cycle_qpcr", "Determine additional cycles to reach about one-third saturation", StepKind.PCR,
                ((qpcr_aliquot, InputRole.TEMPLATE), (p5, InputRole.STRUCTURAL), (p7, InputRole.STRUCTURAL),
                 (pcr, InputRole.CATALYST), (qpcr_mix, InputRole.CONDITION)), (qpcr_product,),
                geometry=PCRGeometry.unavailable(qpcr_aliquot.state, reason="Cycle-determination qPCR uses a separate5µl aliquot of the actual cleaned PCR1 product; its physical indexed ends remain unassigned."))
    graph.mark_terminal(Terminal(qpcr_product.id, Disposition.QC, "Analytical qPCR aliquot sets the additional-cycle count; it is not returned to production PCR2."))
    pcr2 = graph.product("pcr2_product", "off-bead production-PCR2 amplicons", MaterialKind.MOLECULAR,
                         state=_unknown_end_duplex("off-bead PCR2 product"))
    graph.react("index_pcr2", "Re-amplify the25µl production aliquot off bead for qPCR-determined cycles", StepKind.PCR,
                ((production_aliquot, InputRole.TEMPLATE), (p5, InputRole.STRUCTURAL),
                 (p7, InputRole.STRUCTURAL), (pcr, InputRole.CATALYST)), (pcr2,),
                geometry=PCRGeometry.unavailable(production_aliquot.state,
                    reason="Actual cleaned PCR1 production aliquot is the template; additional cycles are sample-specific and physical indexed primer members/end geometry remain unassigned."), source=hairpin_source)
    high_beads = graph.input("high_size_ampure", "Step26:0.5×AMPureXP (30µl beads per60µl PCRmix), RT5min", MaterialKind.SELECTION)
    upper_supernatant = graph.product("upper_size_supernatant", "Step26 first-bind supernatant, source-described0–500bp fraction", MaterialKind.MOLECULAR, state=pcr2.state)
    stored_high = graph.product("stored_high_size_beads", "Step26 first-bind high-size beads with70%EtOH added", MaterialKind.MOLECULAR)
    graph.react("high_size_selection", "Bind the high-size fraction at0.5× and retain supernatant", StepKind.SELECTION,
                ((pcr2, InputRole.TEMPLATE), (high_beads, InputRole.SELECTION)), (upper_supernatant, stored_high))
    graph.mark_terminal(Terminal(stored_high.id, Disposition.STORED, "Step26 adds70%EtOH to the larger-fragment beads and sets them aside."))
    low_beads = graph.input("low_size_ampure", "Step26:add another30µl beads (0.5×original60µl; cumulative1:1)", MaterialKind.SELECTION)
    target_bound = graph.product("polishing_bound", "Step26 cumulative1:1 bound target fraction", MaterialKind.MOLECULAR, state=pcr2.state)
    saved_low = graph.product("saved_low_size_supernatant", "Step26 second-bind low-size/adaptor-dimer supernatant", MaterialKind.MOLECULAR)
    graph.react("low_size_selection", "Bind the target fraction at cumulative1:1; save low-size supernatant", StepKind.SELECTION,
                ((upper_supernatant, InputRole.TEMPLATE), (low_beads, InputRole.SELECTION)), (target_bound, saved_low))
    graph.mark_terminal(Terminal(saved_low.id, Disposition.STORED, "Step26 explicitly saves this low-size supernatant just in case; it is not recorded as discarded."))
    ethanol = graph.input("polishing_wash", "Step26:two700µl70%EtOHwashes;dry~5min;elute50µl20mMTrispH8", MaterialKind.SELECTION)
    library = graph.product("library", "final source-described~100–500bp eluate; conditional shell is a reference", MaterialKind.LIBRARY,
                            state=FINAL_LIBRARY, source=hairpin_source, inferred=True)
    ethanol_waste = graph.product("polishing_ethanol_waste", "Step26 discarded ethanol wash supernatants", MaterialKind.WASTE)
    graph.react("polishing_elution", "Wash and elute the final~100–500bp production library", StepKind.CLEANUP,
                ((target_bound, InputRole.TEMPLATE), (ethanol, InputRole.SELECTION)), (library, ethanol_waste))
    graph.discard(ethanol_waste, "Step26 explicitly discards both ethanol wash supernatants")
    graph.sequence(library, SequencingPlan(
        library=FINAL_LIBRARY, primers=SEQ_PRIMERS, acquired_roles=("Read 1",),
        read_lengths=READ_LENGTHS, read_uses={"Read 1": "include"},
        limitation=SEQUENCING_UNAVAILABLE))
    graph.declare_reviewed_scope((
        "contact", "rna_fragmentation", "fragment_rna", "bridge", "rna_ligase",
        "rna_ligation", "bst", "reverse_transcription", "dpnii", "dpnii_digest", "detached_bridge_piece",
        "dna_ligase", "contact_ligation", "dna_polymerase", "second_strand",
        "proteinase_k", "reverse_crosslinks",
        "shearing_condition", "shear", "streptavidin", "capture",
        "end_repair_enzymes", "end_repair", "tailing_polymerase", "da_tail",
        "hairpin_adapter", "hairpin_ligation", "user_enzyme", "user_opening",
        "post_user_wash", "wash_user_opened", "washed_opened_library",
        "p5_primer", "p7_primer", "pcr_polymerase", "index_pcr", "pcr1",
        "pcr1_bead_separation", "recover_pcr1_supernatant", "first_ampure", "first_ampure_cleanup",
        "pcr1_cleaned", "pcr1_saved_supernatant", "split_pcr1_eluate", "qpcr_aliquot", "production_aliquot",
        "cycle_qpcr_mix", "cycle_qpcr", "qpcr_product", "index_pcr2", "pcr2_product",
        "high_size_ampure", "high_size_selection", "stored_high_size_beads", "low_size_ampure",
        "low_size_selection", "saved_low_size_supernatant", "polishing_wash", "polishing_elution", "library"), source)
    return graph.require_publishable()
