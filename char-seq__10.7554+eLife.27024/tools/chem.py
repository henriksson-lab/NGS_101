"""ChAR-seq RNA–DNA bridge chemistry."""
import illumina as il
from adapter_ligation import AdapterLigationGeometry
from batch_ngs import seg, truseq_library
from chemdraw import Construct, MolecularState, Scene, Segment, Workflow, feature
from endprep import dA_tailed_scene, repair_and_dA_tail
from pcr import PCRGeometry, PCRPrimer
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph,
                            SequencingPlan, StepKind)


TITLE = "ChAR-seq — chromatin-associated RNA sequencing"
NOTES = "01_char-seq.html"
SOURCE = ('Defining source: <a href="https://doi.org/10.7554/eLife.27024">'
          'Bell et al., <i>eLife</i> (2018)</a>. Exact bridge oligos are printed '
          'in the authors\' detailed protocol.')
SUMMARY = ("A pre-adenylated, biotinylated DNA bridge first captures an RNA 3′ end. "
           "DpnII then opens the bridge's DNA-facing end and cuts nearby genomic DNA, "
           "creating a directional RNA–DNA chimera with an inline three-base UMI.")
CAVEAT = ("The bridge sequence, modifications and DpnII-compatible ligation are source-exact. "
          "The paper names the NEBNext USER hairpin workflow; its exact reusable adapter and "
          "index-primer chemistry comes from the authoritative NEB oligo manuals.")

BRIDGE_TOP = "AANNNAAACCGGCGTCCAAGGATCTTTAATTAAGTCGCAG"
BRIDGE_BOTTOM = "GATCTGCGACTTAATTAAAGATCCTTGGACGCCGGTT"
TOP_TAIL = "AA"
UMI_BASES = "NNN"
TOP_SPACER = "AAA"
TOP_CORE = "CCGGCGTCCAAGGATCTTTAATTAAGTCGCA"
TOP_BLOCKED_BASE = "G"
BOTTOM_DPNII_TAIL = "GATC"
BOTTOM_CORE = "TGCGACTTAATTAAAGATCCTTGGACGCCGG"
BOTTOM_BIOTIN_END = "TT"

BRIDGE_UMI = feature("char_bridge_umi", "umi", "random",
                     note="three random bases in the bridge top strand")

FINAL_LIBRARY, SEQ_PRIMERS = truseq_library([
    seg("RNA-derived cDNA", "R" * 38, placeholder=True),
    seg("bridge RNA-ligation tail", TOP_TAIL),
    seg("bridge UMI", UMI_BASES, "umi", placeholder=True, feature=BRIDGE_UMI),
    seg("bridge spacer", TOP_SPACER),
    seg("retained bridge through DpnII junction", "CCGGCGTCCAAG", "me"),
    seg("nearby genomic DNA", "D" * 42, placeholder=True),
], "ChAR-seq RNA–DNA contact library", dual_index=False)
FINAL_CAPTION = ("A source-defined three-base UMI and bridge separate the RNA-derived cDNA "
                 "from nearby genomic DNA. The defining experiment used single-end 152-base "
                 "sequencing; the opposite TruSeq-compatible site is a library end, not a "
                 "reported second read.")
READ_LENGTHS = {"Read 1": 152}
SHOW_READ_CONTENT = True
SEQUENCING_INTRO = ("The NEBNext USER-opened hairpin library is TruSeq-compatible. The paper "
                    "reports single-end 152-base sequencing, so only Read 1 is included in the run.")
SEQUENCING_UNAVAILABLE = ("No second biological read was collected in the defining experiment, "
                          "and the vendor hairpin adapter's exact bases are not asserted here.")


def oligos():
    return (
        f"<b>Bridge top</b> — 5′-/5rApp/{BRIDGE_TOP}/3SpC3/-3′",
        f"<b>Bridge bottom</b> — 5′-/5Phos/{BRIDGE_BOTTOM[:-2]}/iBiodT/T-3′",
        f"<b>NEBNext USER hairpin adaptor</b> — 5′-/5Phos/{il.NEBNEXT_HAIRPIN[:-1]}*{il.NEBNEXT_HAIRPIN[-1]}-3′",
    )


def bridge_scene() -> Scene:
    top = [Segment("5-prime RNA-ligation tail", TOP_TAIL),
           Segment("three-base UMI", UMI_BASES, "umi", placeholder=True,
                   feature=BRIDGE_UMI),
           Segment("top spacer", TOP_SPACER),
           Segment("paired bridge top", TOP_CORE, "me"),
           Segment("3-prime blocked base", TOP_BLOCKED_BASE)]
    bottom = [Segment("5-prime DpnII-side tail", BOTTOM_DPNII_TAIL, "me"),
              Segment("paired bridge bottom", BOTTOM_CORE, "me"),
              Segment("internal-biotin terminal bases", BOTTOM_BIOTIN_END)]
    sc = Scene()
    sc.strand("bridge top", top, label="5′-adenylated top", mod5="/5rApp/",
              mod3="/3SpC3/")
    sc.anneal("bridge bottom", bottom, to="bridge top",
              pair=("paired bridge bottom", "paired bridge top"),
              unpaired=("5-prime DpnII-side tail", "internal-biotin terminal bases",
                        "5-prime RNA-ligation tail", "three-base UMI", "top spacer",
                        "3-prime blocked base"),
              label="5′-phosphorylated, internally biotinylated bottom", mod5="p")
    sc.mark("bridge top", "three-base UMI", "molecule identity")
    sc.mark("bridge bottom", "internal-biotin terminal bases", "internal /iBiodT/")
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
        statement="The printed /5rApp/ bridge top is selectively ligated to fragmented RNA 3′-OH ends by truncated T4 RNA ligase 2.")


def _rna_joined_scene() -> Scene:
    sc = Scene()
    sc.strand("RNA–bridge", [Segment("chromatin RNA", "R" * 28, placeholder=True),
                              Segment("bridge RNA-ligation tail", TOP_TAIL),
                              Segment("bridge UMI", UMI_BASES, "umi", placeholder=True,
                                      feature=BRIDGE_UMI),
                              Segment("bridge remainder", BRIDGE_TOP[8:], "me")],
              label="RNA ligated to bridge top", mod3="/3SpC3/")
    sc.junction("RNA–bridge", "chromatin RNA", "bridge RNA-ligation tail",
                "T4 Rnl2(tr) ligation")
    sc.labels("RNA–bridge")
    return sc


def _dpnii_geometry() -> AdapterLigationGeometry:
    sc = Scene()
    sc.strand("opened bridge", [Segment("bridge DpnII overhang", "GATC", "me"),
                                 Segment("RNA/cDNA and retained bridge", "X" * 30,
                                         placeholder=True)],
              label="DpnII-opened bridge")
    sc.overhang("opened bridge", "bridge DpnII overhang", polarity="5-prime",
                owner="adapter", text="DpnII-compatible 5′ overhang")
    sc.strand("genomic substrate", [Segment("genomic DpnII overhang", "GATC", "me"),
                                     Segment("nearby genomic DNA", "D" * 28,
                                             placeholder=True)],
              label="nearby DpnII-cut genomic DNA")
    sc.overhang("genomic substrate", "genomic DpnII overhang", polarity="5-prime",
                owner="substrate", text="genomic DpnII 5′ overhang")
    return AdapterLigationGeometry.known(
        "cohesive", sc,
        statement="DpnII removes the bridge's 3′ C3-blocked distal piece and produces a GATC-compatible bridge end while cutting nearby genomic DNA to the same cohesive-end class.")


def _contact_scene() -> Scene:
    sc = Scene.duplex([
        Segment("RNA-derived cDNA", "R" * 28, placeholder=True),
        Segment("three-base UMI", UMI_BASES, "umi", placeholder=True,
                feature=BRIDGE_UMI),
        Segment("retained biotin bridge", "CCGGCGTCCAAG", "me"),
        Segment("nearby genomic DNA", "D" * 30, placeholder=True),
    ], label="cDNA–bridge–genomic DNA chimera")
    sc.junction("top", "RNA-derived cDNA", "three-base UMI", "RNA-side bridge")
    sc.junction("top", "retained biotin bridge", "nearby genomic DNA",
                "DpnII cohesive ligation")
    sc.mark("top", "retained biotin bridge", "streptavidin capture")
    sc.labels("top")
    return sc


def workflow() -> Workflow:
    initial = Scene()
    initial.strand("RNA", [Segment("chromatin RNA", "R" * 28, placeholder=True)],
                   label="crosslinked fragmented RNA")
    initial.strand("DNA", [Segment("nearby genomic DNA", "D" * 34, placeholder=True)],
                   label="spatially adjacent chromatin DNA")
    w = Workflow(MolecularState("Fixed RNA–chromatin contact", tuple(initial.rows())))
    w.react("Anneal the exact two-strand bridge", bridge_scene().rows(),
            name="Pre-annealed ChAR-seq bridge",
            note="The top is /5rApp/ and /3SpC3/ blocked; the bottom is 5′ phosphorylated and internally biotinylated. The top carries the three-base UMI.")
    w.ligate_adapter("Ligate RNA 3′ ends to the bridge", _rna_ligation_geometry(),
                      _rna_joined_scene().rows(), name="RNA-ligated bridge",
                      note="The pre-adenylated bridge makes this ligation directional and ATP-independent.")
    w.react("Extend the bridge across RNA by reverse transcription",
            _rna_joined_scene().rows(), name="RNA/cDNA–bridge hybrid",
            note="Bst 3.0 extends from the bridge bottom across the DNA–RNA junction and stabilizes captured RNA.")
    w.ligate_adapter("Digest with DpnII and ligate nearby genomic DNA",
                      _dpnii_geometry(), _contact_scene().rows(),
                      name="RNA–bridge–DNA contact chimera",
                      note="DpnII both opens the blocked DNA-facing bridge end and supplies compatible nearby genomic ends.")
    w.react("Make the second strand, shear, and enrich biotin-bearing contacts",
            _contact_scene().rows(), name="Purified contact duplex",
            note="RNase H/DNA polymerase I complete cDNA; streptavidin retains the internal bridge after sonication.")
    purified_insert = Construct([
        Segment("RNA-derived cDNA", "R" * 28, placeholder=True),
        Segment("bridge UMI", UMI_BASES, "umi", placeholder=True,
                feature=BRIDGE_UMI),
        Segment("retained bridge", "CCGGCGTCCAAG", "me"),
        Segment("nearby genomic DNA", "D" * 30, placeholder=True),
    ], name="sheared bridge-bearing contact")
    end_prepared = repair_and_dA_tail(purified_insert)
    w.react("End-repair and add one 3′ dA to each end",
            dA_tailed_scene(end_prepared, label="dA-tailed contact duplex").rows(),
            name="dA-tailed contact duplex",
            note="The NEBNext end-prep module creates the ligation-ready dA-tailed state explicitly.")
    final = Scene.duplex(list(FINAL_LIBRARY), label="single-end ChAR-seq library",
                         unpaired=tuple(s.name + "'" for s in FINAL_LIBRARY if s.inferred))
    hairpin_geometry = AdapterLigationGeometry.unavailable(
        "The exact NEBNext hairpin is represented by its shared dumbbell/topology model; this compact linear Scene does not flatten both ligated hairpins.")
    w.ligate_adapter("Ligate NEBNext hairpin adapters", hairpin_geometry,
                      final.rows(), name="NEBNext hairpin-adapted contact",
                      note="The paper names the kit; the exact hairpin sequence is supplied by the authoritative NEB oligo manual and centralized in illumina.py.")
    w.react("Open the hairpins with USER and index-amplify",
            final.rows(), name="TruSeq-compatible ChAR-seq library",
            note="USER opening and NEBNext indexing create the reported TruSeq-compatible single-end library.")
    return w


def protocol_graph() -> ProtocolGraph:
    source = "Defining source: https://doi.org/10.7554/eLife.27024; exact bridge oligos from the authors' detailed protocol"
    hairpin_source = source + "; exact hairpin and index-primer sequences from NEB E7600/E7335 oligo manuals"
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
    bridge = graph.input("bridge", "exact pre-annealed biotinylated ChAR-seq bridge",
                         MaterialKind.ADAPTER, state=bridge_scene())
    rna_ligase = graph.input("rna_ligase", "truncated T4 RNA ligase 2",
                             MaterialKind.CATALYST)
    rna_joined = graph.product("rna_joined", "RNA-ligated bridge",
                               MaterialKind.MOLECULAR, state=_rna_joined_scene())
    graph.react("rna_ligation", "Ligate RNA 3′ ends to the bridge", StepKind.LIGATION,
                ((fragmented, InputRole.TEMPLATE), (bridge, InputRole.STRUCTURAL),
                 (rna_ligase, InputRole.CATALYST)), (rna_joined,),
                geometry=_rna_ligation_geometry())
    bst = graph.input("bst", "Bst 3.0 reverse transcriptase", MaterialKind.CATALYST)
    stabilized = graph.product("stabilized_hybrid", "RNA/cDNA–bridge hybrid",
                               MaterialKind.MOLECULAR, state=_rna_joined_scene())
    graph.react("reverse_transcription", "Extend the bridge across RNA",
                StepKind.REVERSE_TRANSCRIPTION,
                ((rna_joined, InputRole.TEMPLATE), (bridge, InputRole.STRUCTURAL),
                 (bst, InputRole.CATALYST)), (stabilized,))
    dpnii = graph.input("dpnii", "DpnII", MaterialKind.CATALYST)
    opened_bridge = graph.product("opened_bridge", "DpnII-opened RNA/cDNA bridge",
                                  MaterialKind.MOLECULAR)
    genomic_end = graph.product("genomic_end", "nearby DpnII-cut genomic DNA",
                                MaterialKind.MOLECULAR)
    graph.react("dpnii_digest", "Open the bridge and cut nearby genomic DNA with DpnII",
                StepKind.DIGESTION,
                ((stabilized, InputRole.TEMPLATE), (dpnii, InputRole.CATALYST)),
                (opened_bridge, genomic_end))
    dna_ligase = graph.input("dna_ligase", "DNA ligase", MaterialKind.CATALYST)
    chimera = graph.product("contact_chimera", "RNA–bridge–DNA contact chimera",
                            MaterialKind.MOLECULAR, state=_contact_scene())
    graph.react("contact_ligation", "Ligate the opened bridge to nearby genomic DNA",
                StepKind.LIGATION,
                ((opened_bridge, InputRole.TEMPLATE), (genomic_end, InputRole.STRUCTURAL),
                 (dna_ligase, InputRole.CATALYST)), (chimera,), geometry=_dpnii_geometry())
    dna_pol = graph.input("dna_polymerase", "DNA polymerase I", MaterialKind.CATALYST)
    duplex = graph.product("contact_duplex", "double-stranded contact chimera",
                           MaterialKind.MOLECULAR, state=_contact_scene())
    graph.react("second_strand", "Make the second DNA strand", StepKind.EXTENSION,
                ((chimera, InputRole.TEMPLATE), (dna_pol, InputRole.CATALYST)), (duplex,))
    shear = graph.input("shearing_condition", "sonication condition", MaterialKind.CONDITION)
    sheared = graph.product("sheared_contacts", "sheared bridge-bearing contacts",
                            MaterialKind.MOLECULAR, state=_contact_scene())
    graph.react("shear", "Shear the contact duplexes", StepKind.FRAGMENTATION,
                ((duplex, InputRole.TEMPLATE), (shear, InputRole.CONDITION)), (sheared,))
    streptavidin = graph.input("streptavidin", "streptavidin beads", MaterialKind.SELECTION)
    captured = graph.product("captured_contacts", "purified biotin-bearing contacts",
                             MaterialKind.MOLECULAR, state=_contact_scene())
    uncaptured = graph.product("uncaptured_contacts", "non-bridge-bearing fragments",
                               MaterialKind.WASTE)
    graph.react("capture", "Enrich biotin-bearing contacts", StepKind.AFFINITY_CAPTURE,
                ((sheared, InputRole.TEMPLATE), (streptavidin, InputRole.SELECTION)),
                (captured, uncaptured))
    graph.discard(uncaptured, "fragments lacking the internal biotin bridge are washed away")
    end_repair = graph.input("end_repair_enzymes", "NEBNext end-repair enzymes",
                             MaterialKind.CATALYST)
    repaired = graph.product("repaired_contacts", "end-repaired contact duplexes",
                             MaterialKind.MOLECULAR)
    graph.react("end_repair", "End-repair the captured contacts", StepKind.END_REPAIR,
                ((captured, InputRole.TEMPLATE), (end_repair, InputRole.CATALYST)),
                (repaired,))
    tailing = graph.input("tailing_polymerase", "dA-tailing polymerase", MaterialKind.CATALYST)
    insert = Construct([
        Segment("RNA-derived cDNA", "R" * 28, placeholder=True),
        Segment("bridge UMI", UMI_BASES, "umi", placeholder=True, feature=BRIDGE_UMI),
        Segment("retained bridge", "CCGGCGTCCAAG", "me"),
        Segment("nearby genomic DNA", "D" * 30, placeholder=True),
    ], name="sheared bridge-bearing contact")
    tailed_state = repair_and_dA_tail(insert)
    tailed = graph.product("tailed_contacts", "dA-tailed contact duplexes",
                           MaterialKind.MOLECULAR,
                           state=dA_tailed_scene(tailed_state, label="dA-tailed contact duplex"))
    graph.react("da_tail", "Add one 3′ dA to each end", StepKind.TAILING,
                ((repaired, InputRole.TEMPLATE), (tailing, InputRole.CATALYST)), (tailed,))
    hairpin = graph.input("hairpin_adapter", "NEBNext USER hairpin adapter", MaterialKind.ADAPTER,
                          state=il.NEBNEXT_HAIRPIN_MODEL, source=hairpin_source)
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
                           MaterialKind.MOLECULAR, source=hairpin_source)
    graph.react("user_opening", "Open the hairpins with USER", StepKind.DIGESTION,
                ((hairpin_ligated, InputRole.TEMPLATE), (user, InputRole.CATALYST)),
                (opened,), source=hairpin_source)
    p5 = graph.input("p5_primer", "NEBNext Universal PCR primer", MaterialKind.OLIGO,
                     source=hairpin_source)
    p7 = graph.input("p7_primer", "NEBNext indexed P7 PCR primer", MaterialKind.OLIGO,
                     source=hairpin_source)
    pcr = graph.input("pcr_polymerase", "index-PCR polymerase", MaterialKind.CATALYST)
    library = graph.product("library", FINAL_LIBRARY.name, MaterialKind.LIBRARY,
                            state=FINAL_LIBRARY, source=hairpin_source)
    graph.react("index_pcr", "Index-amplify the TruSeq-compatible library", StepKind.PCR,
                ((opened, InputRole.TEMPLATE), (p5, InputRole.STRUCTURAL),
                 (p7, InputRole.STRUCTURAL), (pcr, InputRole.CATALYST)), (library,),
                geometry=PCRGeometry.exact(
                    FINAL_LIBRARY,
                    PCRPrimer("NEBNext Universal", il.NEBNEXT_UNIVERSAL_PRIMER),
                    PCRPrimer("NEBNext indexed P7", il.nebnext_i7_primer("N"*8)),
                    statement="The named NEBNext E7335 primer set is sequence-defined in the authoritative vendor oligo manual."),
                source=hairpin_source)
    graph.sequence(library, SequencingPlan(
        library=FINAL_LIBRARY, primers=SEQ_PRIMERS, acquired_roles=("Read 1",),
        read_lengths=READ_LENGTHS, read_uses={"Read 1": "include"},
        limitation=SEQUENCING_UNAVAILABLE))
    graph.declare_reviewed_scope((
        "contact", "rna_fragmentation", "fragment_rna", "bridge", "rna_ligase",
        "rna_ligation", "bst", "reverse_transcription", "dpnii", "dpnii_digest",
        "dna_ligase", "contact_ligation", "dna_polymerase", "second_strand",
        "shearing_condition", "shear", "streptavidin", "capture",
        "end_repair_enzymes", "end_repair", "tailing_polymerase", "da_tail",
        "hairpin_adapter", "hairpin_ligation", "user_enzyme", "user_opening",
        "p5_primer", "p7_primer", "pcr_polymerase", "index_pcr", "library"), source)
    return graph.require_publishable()
