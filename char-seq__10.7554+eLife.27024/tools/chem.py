"""ChAR-seq RNA–DNA bridge chemistry."""
from adapter_ligation import AdapterLigationGeometry
from batch_ngs import seg, truseq_library
from chemdraw import Construct, MolecularState, Scene, Segment, Workflow
from endprep import dA_tailed_scene, repair_and_dA_tail
from pcr import PCRGeometry
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph,
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

FINAL_LIBRARY, SEQ_PRIMERS = truseq_library([
    seg("RNA-derived cDNA", "R" * 38, placeholder=True),
    seg("bridge RNA-ligation tail", TOP_TAIL),
    seg("retained bridge through DpnII junction", RETAINED_BRIDGE, "me"),
    seg("nearby genomic DNA", "D" * 42, placeholder=True),
], "ChAR-seq RNA–DNA contact library", dual_index=False, inferred_adapters=True)
FINAL_CAPTION = ("INFERRED — The source-defined bridge separates the RNA-derived cDNA "
                 "from nearby genomic DNA. The defining experiment used single-end 152-base "
                 "sequencing; the TruSeq-compatible shell depicts an unverified NEBNext "
                 "index-kit member, and the opposite site is not a reported second read.")
READ_LENGTHS = {"Read 1": 152}
SHOW_READ_CONTENT = True
SEQUENCING_INTRO = ("The source names a NEBNext USER-opened TruSeq-compatible library. "
                    "Its exact kit revision and index member are unverified; the paper reports "
                    "single-end 152-base sequencing, so only Read 1 is included in the run.")
SEQUENCING_UNAVAILABLE = ("No second biological read was collected in the defining experiment, "
                          "and the vendor hairpin adapter's exact bases are not asserted here.")


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


def _rna_joined_scene() -> Scene:
    sc = Scene()
    sc.strand("RNA–bridge", [Segment("chromatin RNA", "R" * 28, placeholder=True),
                              Segment("bridge RNA-ligation tail", TOP_TAIL),
                              Segment("bridge remainder", TOP_CORE, "me")],
              label="RNA ligated to bridge top", mod3="Sp3")
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
        Segment("bridge RNA-ligation tail", TOP_TAIL),
        Segment("retained desthiobiotin bridge", RETAINED_BRIDGE, "me"),
        Segment("nearby genomic DNA", "D" * 30, placeholder=True),
    ], label="cDNA–bridge–genomic DNA chimera")
    sc.junction("top", "RNA-derived cDNA", "bridge RNA-ligation tail", "RNA-side bridge")
    sc.junction("top", "retained desthiobiotin bridge", "nearby genomic DNA",
                "DpnII cohesive ligation")
    sc.mark("top", "retained desthiobiotin bridge", "streptavidin capture")
    sc.labels("top")
    return sc


def workflow() -> Workflow:
    initial = Scene()
    initial.strand("RNA", [Segment("chromatin RNA", "R" * 28, placeholder=True)],
                   label="crosslinked fragmented RNA")
    initial.strand("DNA", [Segment("nearby genomic DNA", "D" * 34, placeholder=True)],
                   label="spatially adjacent chromatin DNA")
    w = Workflow(MolecularState("Fixed RNA–chromatin contact", tuple(initial.rows())))
    w.react("Anneal the 2018 two-strand bridge", bridge_scene().rows(),
            name="Pre-annealed ChAR-seq bridge",
            note="The 2018 top has 5′ App and 3′ Sp3; the bottom carries an internal desthiobiotin-TEG at a bond not resolved by the figure. The 2018 bridge has no UMI.")
    w.ligate_adapter("Ligate RNA 3′ ends to the bridge", _rna_ligation_geometry(),
                      _rna_joined_scene().rows(), name="RNA-ligated bridge",
                      note="The pre-adenylated bridge makes this ligation directional and ATP-independent; T4 trKQ RNA ligase acts overnight at 23 °C in 25% PEG inside permeabilized nuclei (extended protocol step 5).")
    w.react("Extend the bridge across RNA by reverse transcription",
            _rna_joined_scene().rows(), name="RNA/cDNA–bridge hybrid",
            note="Bst 3.0 extends from the bridge bottom across the DNA–RNA junction and stabilizes captured RNA.")
    w.ligate_adapter("Digest with DpnII and ligate nearby genomic DNA",
                      _dpnii_geometry(), _contact_scene().rows(),
                      name="RNA–bridge–DNA contact chimera",
                      note="DpnII both opens the blocked DNA-facing bridge end and supplies compatible nearby genomic ends.")
    w.react("Make the second strand, shear, and enrich biotin-bearing contacts",
            _contact_scene().rows(), name="Purified contact duplex",
            note="Limiting RNase H and E. coli DNA polymerase I complete the cDNA (16 °C, step 15); crosslinks are then reversed with SDS/NaCl/proteinase K at 70 °C overnight and the DNA is ethanol precipitated before Covaris shearing to ~200 bp and capture on Dynabeads MyOne Streptavidin T1 (steps 16–19).")
    purified_insert = Construct([
        Segment("RNA-derived cDNA", "R" * 28, placeholder=True),
        Segment("bridge RNA-ligation tail", TOP_TAIL),
        Segment("retained bridge", RETAINED_BRIDGE, "me"),
        Segment("nearby genomic DNA", "D" * 30, placeholder=True),
    ], name="sheared bridge-bearing contact")
    end_prepared = repair_and_dA_tail(purified_insert)
    w.react("End-repair and add one 3′ dA to each end",
            dA_tailed_scene(end_prepared, label="dA-tailed contact duplex").rows(),
            name="dA-tailed contact duplex",
            note="The NEBNext End Prep enzyme mix acts on the bead-bound fragments (20 °C 20 min, 65 °C 30 min, step 20) and creates the ligation-ready dA-tailed state explicitly.")
    final = Scene.duplex(list(FINAL_LIBRARY), label="single-end ChAR-seq library",
                         unpaired=tuple(s.name + "'" for s in FINAL_LIBRARY if s.inferred))
    hairpin_geometry = AdapterLigationGeometry.unavailable(
        "The exact NEBNext hairpin is represented by its shared dumbbell/topology model; this compact linear Scene does not flatten both ligated hairpins.")
    w.ligate_adapter("Ligate NEBNext hairpin adapters", hairpin_geometry,
                      final.rows(), name="NEBNext hairpin-adapted contact",
                      note="Ligation is on bead with NEBNext Ultra II Ligation Master Mix and ligation enhancer (step 21); the protocol describes the adaptor as a 65-nt hairpin, but neither its exact kit revision nor the selected indexed primer is established for this 2018 experiment.")
    w.react("Open the hairpins with USER and index-amplify",
            final.rows(), name="TruSeq-compatible ChAR-seq library",
            note="USER opens the hairpin on bead (37 °C, 15 min); eight on-bead cycles with the NEBNext Universal and Indexing primers are followed by AMPure size selection, a qPCR cycle estimate and a second off-bead PCR of N further cycles (steps 21-25).")
    return w


def protocol_graph() -> ProtocolGraph:
    source = "Defining source: https://doi.org/10.7554/eLife.27024 Figure 1—figure supplement 1; bridge strings corroborated by bioRxiv 118786 supplement 1"
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
    dna_ligase = graph.input("dna_ligase", "T4 DNA ligase", MaterialKind.CATALYST)
    chimera = graph.product("contact_chimera", "RNA–bridge–DNA contact chimera",
                            MaterialKind.MOLECULAR, state=_contact_scene())
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
                           MaterialKind.MOLECULAR, source=hairpin_source)
    graph.react("user_opening", "Open the hairpins with USER", StepKind.DIGESTION,
                ((hairpin_ligated, InputRole.TEMPLATE), (user, InputRole.CATALYST)),
                (opened,), source=hairpin_source)
    p5 = graph.input("p5_primer", "NEBNext Universal PCR primer", MaterialKind.OLIGO,
                     source=hairpin_source, inferred=True)
    p7 = graph.input("p7_primer", "NEBNext indexed P7 PCR primer", MaterialKind.OLIGO,
                     source=hairpin_source, inferred=True)
    pcr = graph.input("pcr_polymerase", "NEBNext High-Fidelity PCR master mix",
                      MaterialKind.CATALYST)
    pcr1 = graph.product("pcr1", "on-bead eight-cycle indexed amplification product",
                         MaterialKind.MOLECULAR, state=FINAL_LIBRARY,
                         source=hairpin_source, inferred=True)
    graph.react("index_pcr", "Amplify on bead with the Universal and Indexing primers (8 cycles)",
                StepKind.PCR,
                ((opened, InputRole.TEMPLATE), (p5, InputRole.STRUCTURAL),
                 (p7, InputRole.STRUCTURAL), (pcr, InputRole.CATALYST)), (pcr1,),
                geometry=PCRGeometry.unavailable(
                    FINAL_LIBRARY,
                    reason="The 2018 source names NEBNext Indexing Primers for Illumina but not the exact kit revision, primer strands, or selected index; E7335/E7600 manuals establish candidate family sequences only."),
                source=hairpin_source)
    library = graph.product("library", FINAL_LIBRARY.name, MaterialKind.LIBRARY,
                            state=FINAL_LIBRARY, source=hairpin_source, inferred=True)
    graph.react("index_pcr2",
                "Size-select and re-amplify off bead for a qPCR-determined cycle number",
                StepKind.PCR,
                ((pcr1, InputRole.TEMPLATE), (p5, InputRole.STRUCTURAL),
                 (p7, InputRole.STRUCTURAL), (pcr, InputRole.CATALYST)), (library,),
                geometry=PCRGeometry.unavailable(
                    FINAL_LIBRARY,
                    reason="The second-round cycle number is set per sample from the step-25 qPCR, and the indexed primer member is still unassigned."),
                source=hairpin_source)
    graph.sequence(library, SequencingPlan(
        library=FINAL_LIBRARY, primers=SEQ_PRIMERS, acquired_roles=("Read 1",),
        read_lengths=READ_LENGTHS, read_uses={"Read 1": "include"},
        limitation=SEQUENCING_UNAVAILABLE))
    graph.declare_reviewed_scope((
        "contact", "rna_fragmentation", "fragment_rna", "bridge", "rna_ligase",
        "rna_ligation", "bst", "reverse_transcription", "dpnii", "dpnii_digest",
        "dna_ligase", "contact_ligation", "dna_polymerase", "second_strand",
        "proteinase_k", "reverse_crosslinks",
        "shearing_condition", "shear", "streptavidin", "capture",
        "end_repair_enzymes", "end_repair", "tailing_polymerase", "da_tail",
        "hairpin_adapter", "hairpin_ligation", "user_enzyme", "user_opening",
        "p5_primer", "p7_primer", "pcr_polymerase", "index_pcr", "pcr1",
        "index_pcr2", "library"), source)
    return graph.require_publishable()
