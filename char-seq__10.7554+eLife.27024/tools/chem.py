"""ChAR-seq RNA–DNA bridge chemistry."""
from adapter_ligation import AdapterLigationGeometry
from batch_ngs import seg, truseq_library
from chemdraw import Construct, MolecularState, Scene, Segment, Workflow, feature
from endprep import dA_tailed_scene, repair_and_dA_tail


TITLE = "ChAR-seq — chromatin-associated RNA sequencing"
NOTES = "01_char-seq.html"
SOURCE = ('Defining source: <a href="https://doi.org/10.7554/eLife.27024">'
          'Bell et al., <i>eLife</i> (2018)</a>. Exact bridge oligos are printed '
          'in the authors\' detailed protocol.')
SUMMARY = ("A pre-adenylated, biotinylated DNA bridge first captures an RNA 3′ end. "
           "DpnII then opens the bridge's DNA-facing end and cuts nearby genomic DNA, "
           "creating a directional RNA–DNA chimera with an inline three-base UMI.")
CAVEAT = ("The bridge sequence, modifications and DpnII-compatible ligation are source-exact. "
          "The terminal NEBNext hairpin oligo bases were not disclosed in the sources read; "
          "the final TruSeq-compatible sequencing shell is consequently marked inferred.")

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
], "ChAR-seq RNA–DNA contact library", dual_index=False, inferred_adapters=True)
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
        "The paper names the NEBNext hairpin adapter but does not print the exact adapter lot's cut-facing strand or terminal base sequence.")
    w.ligate_adapter("Ligate NEBNext hairpin adapters", hairpin_geometry,
                      final.rows(), name="NEBNext hairpin-adapted contact",
                      note="The commercial adapter identity is known, but this source does not establish its exact oligo sequence.")
    w.react("Open the hairpins with USER and index-amplify",
            final.rows(), name="TruSeq-compatible ChAR-seq library",
            note="USER opening and NEBNext indexing create the reported TruSeq-compatible single-end library.")
    return w
