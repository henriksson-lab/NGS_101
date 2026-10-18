"""Source-exact scDam&T-seq DamID and CEL-seq2 molecular states."""
from dataclasses import replace
from adapter_ligation import AdapterLigationGeometry
from chemdraw import (Construct, MolecularState, Row, Scene, Segment, Workflow,
                      feature, revcomp)
import illumina as il
import seqprimers as sp
from pcr import PCRGeometry, PCRPrimer
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph,
                            SequencingPlan, StepKind, Terminal, Disposition)

TITLE = "scDam&T-seq — 2020 detailed-protocol example"
NOTES = "01_scdam-t-seq.html"
SOURCE = ('Defining source: <a href="https://doi.org/10.1038/s41587-019-0150-y">'
          'Rooijers et al., <i>Nature Biotechnology</i> (2019)</a>; detailed protocol '
          '<a href="https://doi.org/10.1038/s41596-020-0314-8">'
          '<i>Nature Protocols</i> (2020)</a>.')
SUMMARY = ("One well-specific CEL-seq2 primer labels transcript cDNA while a different "
           "well-specific forked adapter labels DpnI-cut Dam-methylated DNA. Both carry "
           "T7 promoters, so the two branches are pooled and linearly amplified together.")
CAVEAT = ("<b>Scope.</b> This page models the explicitly printed 2020 split-barcode example, not an established 2019 experiment-to-oligo allocation. Exact first-cycle PCR geometry, T7 initiation register and run-specific sequencing-primer setup remain unassigned. Candidate other RT products and RNA/DNA carryover are inferred; no observed abundance or complete physical removal is claimed. ""<b>Published concentration inconsistency.</b> The 2020 randomhexRT mix table (PMC7779467 JATS878–884) lists a 20 µM primer stock and 13.33 mM final concentration. Dilution cannot raise that stock concentration; the intended unit remains unverified. This model assigns no corrected reaction concentration. "
          "<b>Published barcode-example discrepancy.</b> The 2020 Reagents, DamID adapter example prints TGCATATG in its top/bottom oligos and names that barcode in prose (PMC7779467 JATS284–288), but the following split-layout example repeats NNN-GATG-NNN-TCAT from the RNA primer. The model follows the printed DamID oligos, NNN-TGCA-NNN-TATG; the repeated prose example appears copied from the RNA branch and is not used as a second DamID variant. "
          "<b>Published read-label discrepancy.</b> Markodimitraki et al. 2020, Introduction/Overview (PMC7779467 JATS61), assigns transcript sequence to R1 and transcript identifiers to R2. Its data-processing instructions (JATS1102, CEL-Seq2 read processing) instead use R2 for transcript alignment, agreeing with Rooijers et al. 2019, Methods, Demultiplexing (PMC6609448 JATS113). This model follows the latter assignment: RNA R1 carries identifiers and RNA R2 supplies the transcript. The overview labels appear reversed; the precise custom run configuration remains to be checked. "
          "<b>Adapter branch shown:</b> the six-base fork and split 8-nt barcode below "
          "are the example in Markodimitraki et al. 2020, Supplementary Table 2. "
          "Rooijers et al. 2019, Methods ‘DamID adapters’ and Supplementary Table 2, "
          "also describe a four-base-fork adapter for other cells; the Dam-RING1B "
          "mESC branch uses a six-base fork. These are distinct published variants.")

T7 = "GCCGGTAATACGACTCACTATAGGG"
RA5 = il.TRUSEQ_SMALL_RNA_RA5[6:]  # exact CEL primer prefix in 2020 Reagents
UMI1, UMI2 = "NNN", "NNN"
RNA_BC1, RNA_BC2 = "GATG", "TCAT"       # protocol's printed example
DNA_BC1, DNA_BC2 = "TGCA", "TATG"       # protocol's printed example
POLYT = "T" * 25
CEL_PRIMER = T7 + RA5 + UMI1 + RNA_BC1 + UMI2 + RNA_BC2 + POLYT + "V"

DAM_FORK_TOP = "GGTGAT"
DAM_FORK_BOTTOM = "GAGCTT"
DAM_CORE_BODY = "CCGGTAATACGACTCACTATAGGGGTTCAGAGTTCTACAGTCCGACGATC"
DAM_TOP = DAM_FORK_TOP + DAM_CORE_BODY + UMI1 + DNA_BC1 + UMI2 + DNA_BC2 + "GA"
DAM_BOTTOM = ("TC" + revcomp(DNA_BC2) + UMI2 + revcomp(DNA_BC1) + UMI1
              + revcomp(DAM_CORE_BODY) + DAM_FORK_BOTTOM)

RANDOMHEX_TAIL = "GCCTTGGCACCCGAGAATTCCA"
RANDOMHEX_RT = RANDOMHEX_TAIL + "N" * 6
RNA_PCR1 = il.P5 + il.TRUSEQ_SMALL_RNA_RA5[:-5]  # source-printed RP1, 2020 Reagents
SAMPLE_INDEX = "CGTGAT"  # 2020 Reagents printed RPi example
RNA_PCR_INDEX = (il.P7 + SAMPLE_INDEX
                 + "GTGACTGGAGTTCCTTGGCACCCGAGAATTCCA")
RNA_PCR_INDEX_ANNEAL = "GTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"


def _umi(name):
    branch = "rna" if name.startswith(("CEL", "RNA")) else "damid"
    part = "part1" if "1" in name else "part2"
    return Segment(name, "NNN", "umi", feature=feature(branch + "_umi_" + part,
        "umi", "random", group=branch + "_umi", part=part))


def _barcode(name, bases, identity):
    identity = identity.lower().replace("cel_cell", "rna_cell")
    part = "part1" if "1" in name else "part2"
    return Segment(name, bases, "cbc", feature=feature(identity + "_" + part,
        "cell_barcode", "fixed", group=identity, part=part))


def cel_primer_segments():
    return [Segment("CEL T7 promoter", T7), Segment("CEL RA5", RA5, "r1"),
            _umi("CEL UMI part 1"), _barcode("CEL barcode part 1", RNA_BC1, "cel_cell"),
            _umi("CEL UMI part 2"), _barcode("CEL barcode part 2", RNA_BC2, "cel_cell"),
            Segment("CEL poly(dT)25", POLYT),
            Segment("CEL anchoring V", "V", placeholder=True)]


def cel_rt_scene():
    sc = Scene()
    sc.strand("mRNA", [Segment("mRNA body", "X" * 40, placeholder=True),
                        Segment("poly(A)-proximal B", "B", placeholder=True),
                        Segment("poly(A) tail", "A" * 25)], label="polyadenylated mRNA")
    sc.anneal("CEL-seq2 primer", cel_primer_segments(), to="mRNA",
              pair=("CEL anchoring V", "poly(A)-proximal B"),
              label="well-specific CEL-seq2 primer")
    sc.arrow("CEL-seq2 primer", "reverse transcriptase copies the transcript")
    sc.labels("CEL-seq2 primer")
    return sc


def damid_adapter_scene():
    top = [Segment("top six-base fork", DAM_FORK_TOP),
           Segment("DamID constant core", DAM_CORE_BODY), _umi("DamID UMI part 1"),
           _barcode("DamID barcode part 1", DNA_BC1, "damid_cell"),
           _umi("DamID UMI part 2"),
           _barcode("DamID barcode part 2", DNA_BC2, "damid_cell"),
           Segment("terminal GA", "GA")]
    bottom = [Segment("terminal TC", "TC"),
              _barcode("DamID barcode part 2 bottom", revcomp(DNA_BC2), "damid_cell"),
              _umi("DamID UMI part 2 bottom"),
              _barcode("DamID barcode part 1 bottom", revcomp(DNA_BC1), "damid_cell"),
              _umi("DamID UMI part 1 bottom"),
              Segment("DamID constant-core bottom", revcomp(DAM_CORE_BODY)),
              Segment("bottom six-base fork", DAM_FORK_BOTTOM)]
    sc = Scene(); sc.strand("DamID top", top, label="DamID top oligo")
    sc.anneal("DamID bottom", bottom, to="DamID top",
              pair=("terminal TC", "terminal GA"), label="DamID bottom oligo",
              mod5="/5Phos/", unpaired=("bottom six-base fork",))
    sc.mark("DamID top", "top six-base fork", "non-complementary six-base fork")
    sc.mark("DamID bottom", "bottom six-base fork", "prevents adapter concatemers")
    sc.ligation_end("DamID top", "terminal GA", kind="blunt", owner="adapter",
                    text="blunt terminal GA/TC half-GATC end")
    sc.labels("DamID top")
    return sc


def dpni_fragment_scene():
    sc = Scene.duplex([Segment("left DpnI half-site", "TC"),
                       Segment("methylated-GATC genomic fragment", "X" * 36,
                               placeholder=True),
                       Segment("right DpnI half-site", "GA")],
                      label="DpnI fragment with blunt half-GATC ends")
    sc.ligation_end("top", "left DpnI half-site", kind="blunt", owner="substrate",
                    text="blunt DpnI substrate end")
    sc.mark("top", "left DpnI half-site", "DpnI cut G^mATC leaves TC on this end")
    sc.labels("top")
    return sc


def damid_geometry():
    return AdapterLigationGeometry.known(
        "blunt", damid_adapter_scene(), dpni_fragment_scene(),
        statement="The detailed protocol prints both adapter strands and states that DpnI leaves blunt ends; terminal adapter GA plus genomic TC restores GATC.")


def ligated_damid_scene():
    top = [Segment("top six-base fork", DAM_FORK_TOP),
           Segment("DamID constant core", DAM_CORE_BODY), _umi("DamID UMI part 1"),
           _barcode("DamID barcode part 1", DNA_BC1, "damid_cell"),
           _umi("DamID UMI part 2"),
           _barcode("DamID barcode part 2", DNA_BC2, "damid_cell"),
           Segment("adapter half-site GA", "GA"), Segment("genomic half-site TC", "TC"),
           Segment("DpnI-selected genomic DNA", "X" * 40, placeholder=True)]
    bottom = [Segment("genomic bottom", "X" * 40, placeholder=True),
              Segment("genomic half-site bottom", "GA"),
              Segment("adapter half-site bottom", "TC"),
              _barcode("DamID barcode part 2 bottom", revcomp(DNA_BC2), "damid_cell"),
              _umi("DamID UMI part 2 bottom"),
              _barcode("DamID barcode part 1 bottom", revcomp(DNA_BC1), "damid_cell"),
              _umi("DamID UMI part 1 bottom"),
              Segment("DamID constant-core bottom", revcomp(DAM_CORE_BODY)),
              Segment("bottom six-base fork", DAM_FORK_BOTTOM)]
    sc = Scene(); sc.strand("top", top, label="one adapter-ligated DpnI end")
    sc.anneal("bottom", bottom, to="top",
              pair=("genomic bottom", "DpnI-selected genomic DNA"),
              label="one adapter-ligated DpnI end",
              unpaired=("bottom six-base fork",))
    sc.junction("top", "adapter half-site GA", "genomic half-site TC",
                "blunt ligation restores GATC")
    sc.mark("top", "DamID barcode part 1", "DamID-specific well barcode")
    sc.labels("top")
    return sc


# Native declarations below select the explicitly printed 2020 example, not an
# assumed allocation to a 2019 experiment. Both populations remain in one tube.
METHODS = "Markodimitraki et al. 2020 PMC7779467, Procedure RT/second strand/DpnI/ligation/IVT/library prep; printed 384-well example"
ARTICLE = "Rooijers et al. 2019 PMC6609448, Methods Robotic preparation and Raw data preprocessing"
REAGENTS = "Markodimitraki et al. 2020 PMC7779467, Reagents CEL primer/DamID adapter/Library primers, lines 274–302"


def _oligo_scene(name, sequence, mod5="", mod3=""):
    sc = Scene(); sc.strand(name, cel_primer_segments() if name == "CEL-seq2 primer" else [Segment(name, sequence, placeholder="X" in sequence)], mod5=mod5, mod3=mod3)
    return sc


class MolecularMixture(Scene):
    """Independent source Scene components, retaining semantic strands/features.

    Rendering delegates to those same components. The combined strand inventory
    namespaces names and partner links, without copying a molecular declaration.
    """
    def __init__(self, *components, unplaced_segments=()):
        super().__init__()
        self.components = tuple(part for value in components
            for part in (value.components if isinstance(value, MolecularMixture) else (value,)))
        self.unplaced_segments = tuple(unplaced_segments) + tuple(seg for value in components
            for seg in getattr(value, "unplaced_segments", ()))
        for index, scene in enumerate(self.components):
            if not isinstance(scene, Scene):
                raise TypeError("molecular mixtures retain Scene components, not rendered rows")
            for name, strand in scene.strands.items():
                key = f"population{index}:{name}"
                self.strands[key] = replace(strand, name=key,
                    partner=f"population{index}:{strand.partner}" if strand.partner else None)

    def rows(self, omit=()):
        rows = []
        for scene in self.components:
            if rows:
                rows.append(Row())
            rows.extend(scene.rows())
        return rows


def _mixed(*components):
    return MolecularMixture(*components)


def _possible_carry(name):
    sc = Scene(); sc.strand(name, [Segment(name, "[" + name + "]", placeholder=True, inferred=True)],
        label="INFERRED — possible co-carried " + name)
    sc.footer("Candidate physical carryover; no yield or complete removal is asserted.", name)
    return sc


def _transcript_cdna(duplex=False):
    parts = [*cel_primer_segments(), Segment("copied transcript", "X" * 40, placeholder=True)]
    if duplex:
        return Scene.duplex(parts, label="transcript-derived ds cDNA")
    sc = Scene(); sc.strand("transcript cDNA", parts, label="first-strand transcript cDNA")
    sc.anneal("mRNA", [Segment("mRNA body", "X" * 40, placeholder=True),
        Segment("poly(A)-proximal B", "B", placeholder=True), Segment("poly(A) tail", "A" * 25)],
        to="transcript cDNA", pair=("poly(A) tail", "CEL poly(dT)25"), label="template mRNA")
    return sc


def _mrna():
    sc = Scene(); sc.strand("mRNA", [Segment("mRNA body", "X" * 40, placeholder=True),
        Segment("poly(A)-proximal B", "B", placeholder=True), Segment("poly(A) tail", "A" * 25)], label="native polyadenylated RNA")
    return sc


def _genomic():
    sc = Scene.duplex([Segment("genomic DNA", "X" * 40, placeholder=True)], label="genomic DNA")
    sc.footer("Dam deposits N6-methyladenine; methylation is not copied to cDNA.", "top")
    return sc


def _arna(branch):
    # Source-qualified representative T7 transcripts. Payload length is schematic;
    # the exact initiation register and opposite-end adapter readthrough are not
    # claimed. The common RA5/identifier portion is retained for native RT/PCR.
    sc = Scene()
    ids = [_umi(branch + " UMI part 1"),
           _barcode(branch + " barcode part 1", RNA_BC1 if branch == "RNA" else DNA_BC1, branch + "_cell"),
           _umi(branch + " UMI part 2"),
           _barcode(branch + " barcode part 2", RNA_BC2 if branch == "RNA" else DNA_BC2, branch + "_cell")]
    parts = [Segment(branch + " RA5", RA5.replace("T", "U"), placeholder=True), *ids,
             *([Segment("copied poly(U)25", "U" * 25, placeholder=True), Segment("copied anchoring V", "V", placeholder=True)] if branch == "RNA" else [Segment("restored GAUC", "GAUC", placeholder=True)]),
             Segment(branch + " payload", "[transcript-derived aRNA]" if branch == "RNA" else "[GATC-derived aRNA]", placeholder=True)]
    sc.strand(branch + " aRNA", parts, label=branch + " aRNA (representative identifier-bearing region)")
    sc.footer("RNA copy of DNA sequence only; the ligated adapter phosphate does not become an RNA sequence modification.", branch + " aRNA")
    return sc



def _fragmented(*source_rna):
    sc = Scene()
    for branch in ("RNA", "DamID"):
        sc.strand(branch + " identifier-bearing fragment", [Segment(branch + " fragment",
            "[RA5/identifiers and variable downstream RNA fragment]", placeholder=True)], label=branch + " fragment class")
        sc.strand(branch + " other fragments", [Segment(branch + " other RNA fragments",
            "[partial identifiers or downstream-only RNA fragments]", placeholder=True)], label=branch + " other fragment classes")
    inventory = tuple(seg for component in source_rna
        for strand in component.strands.values() for seg in strand.segs if seg.feature)
    return MolecularMixture(sc, unplaced_segments=inventory)


def _library_scene(library, branch):
    sc = Scene.duplex(list(library), label=branch + " library molecules", unpaired=(branch + " insert'",))
    sc.mark("top", "RP1 sequence", "original RP1 phosphorothioate at its penultimate/final base bond; copied complement has no inherited phosphorothioate")
    sc.mark("bottom", "RPi annealing arm reverse complement'", "original RPi phosphorothioate at its penultimate/final base bond; copied complement has no inherited phosphorothioate")
    return sc

def _candidate_cdna():
    sc = Scene(); sc.strand("other cDNAs", [Segment("candidate other cDNA", "X" * 36, placeholder=True, inferred=True)], label="INFERRED — candidate other RT products")
    return sc


def _native_cdna():
    """Actual single-stranded randomhexRT products from two RNA populations."""
    sc = Scene()
    for branch, bc1, bc2 in (("RNA", RNA_BC1, RNA_BC2), ("DamID", DNA_BC1, DNA_BC2)):
        # Reverse-transcribed identifier orientation is reverse-complemented.
        parts = [Segment(branch + " randomhexRT tail", RANDOMHEX_TAIL),
                 Segment(branch + " random primer N6", "NNNNNN", placeholder=True),
                 Segment(branch + " copied payload", "[copied aRNA fragment]", placeholder=True),
                 *([Segment("copied anchoring B", "B", placeholder=True), Segment("copied poly(A)25", "A" * 25)] if branch == "RNA" else [Segment("copied GATC", "GATC")]),
                 _barcode(branch + " copied barcode part 2", revcomp(bc2), branch + "_cell"),
                 _umi(branch + " copied UMI part 2"),
                 _barcode(branch + " copied barcode part 1", revcomp(bc1), branch + "_cell"),
                 _umi(branch + " copied UMI part 1"),
                 Segment(branch + " copied RA5", revcomp(RA5))]
        sc.strand(branch + " ss cDNA", parts, label=branch + " identifier-bearing randomhexRT cDNA")
    sc.footer("Two molecular populations in the same reaction. The first PCR cycle extends RP1 on these ss templates; a duplex precursor is not supplied.", "RNA ss cDNA")
    return sc


def _endpoint(branch, bc1, bc2):
    # RP1 anneals within RA5, leaving its CGATC suffix before identifiers. RPi
    # initially anneals through randomhexRT's 21-base tail; its other bases are
    # added by PCR, not pre-existing on the ss cDNA input.
    forward_overlap = next(n for n in range(min(len(RA5), len(RNA_PCR1)), 0, -1)
                           if RA5.startswith(RNA_PCR1[-n:]))
    if forward_overlap < 15:
        raise ValueError("RP1 lacks its source-defined RA5 overlap")
    return Construct([Segment("RP1 sequence", RNA_PCR1, "p5"),
        Segment("retained RA5 suffix", RA5[forward_overlap:]),
        _umi(branch + " UMI part 1"), _barcode(branch + " barcode part 1", bc1, branch + "_cell"),
        _umi(branch + " UMI part 2"), _barcode(branch + " barcode part 2", bc2, branch + "_cell"),
        *([Segment("copied poly(dT)25", POLYT), Segment("copied anchoring V", "V", placeholder=True)] if branch == "RNA" else [Segment("restored GATC", "GATC")]),
        Segment(branch + " insert", "[transcript-derived insert]" if branch == "RNA" else "[DpnI genomic insert]", placeholder=True),
        Segment("random hexamer copy", "NNNNNN", placeholder=True),
        Segment("RPi annealing arm reverse complement", revcomp(RNA_PCR_INDEX_ANNEAL), "r2"),
        Segment("sample index", revcomp(SAMPLE_INDEX), "cbc", feature=feature("library_index", "sample_index", "fixed")),
        Segment("P7 reverse complement", revcomp(il.P7), "p7")], name=branch + " molecules in mixed indexed library")


RNA = _endpoint("RNA", RNA_BC1, RNA_BC2)
DNA = _endpoint("DamID", DNA_BC1, DNA_BC2)
# Exact sequencing-primer identity/custom run setup is not established in this
# capped pass. Source read roles stay explicit; no standard TruSeq primers added.
RP = DP = ()
FINAL_LIBRARIES = (
    ("Transcriptome library", RNA, RP, "RNA-derived molecules in the same physical library.", ""),
    ("DamID library", DNA, DP, "DamID-derived molecules in the same physical library.", ""),
)
SEQUENCING_UNAVAILABLE = {name: "Small RNA library oligos are source-known; exact sequencing-primer/run setup is unassigned in this pass. No generic TruSeq sequencing shell is substituted." for name, *_ in FINAL_LIBRARIES}
READ_SPECS = {
    "Transcriptome library": [sp.read_spec("Read 1", 75, "split UMI/cell barcode; identifiers extracted computationally", "extract"), sp.read_spec("Read 2", 75, "transcript sequence", "include")],
    "DamID library": [sp.read_spec("Read 1", 75, "DamID UMI/cell barcode then genomic sequence; identifiers extracted before mapping", "extract_then_include"), sp.read_spec("Read 2", 75, "paired genomic fragment", "undetermined")],
}


def oligos():
    return [f"<b>CEL-seq2 example primer</b> — 5′-{CEL_PRIMER}-3′",
            f"<b>DamID example top</b> — 5′-{DAM_TOP}-3′",
            f"<b>DamID example bottom</b> — 5′-/5Phos/{DAM_BOTTOM}-3′",
            f"<b>randomhexRT</b> — 5′-{RANDOMHEX_RT}-3′",
            f"<b>RP1</b> — 5′-{RNA_PCR1[:-1]}*{RNA_PCR1[-1]}-3′",
            f"<b>RPi printed CGTGAT example</b> — 5′-{RNA_PCR_INDEX[:-1]}*{RNA_PCR_INDEX[-1]}-3′"]


def _cleanup(g, current, ident, label, reagent):
    beads = g.input(ident + "_beads", reagent, MaterialKind.SELECTION, source=METHODS)
    kept = g.product(ident + "_kept", label, MaterialKind.MOLECULAR, state=current.state, source=METHODS)
    waste = g.product(ident + "_waste", "bead-cleanup aqueous phase and ethanol washes", MaterialKind.WASTE, source=METHODS)
    g.react(ident, "Bead cleanup of both molecular populations together", StepKind.CLEANUP,
            ((current, InputRole.TEMPLATE), (beads, InputRole.SELECTION)), (kept, waste), changes_molecule=False, source=METHODS)
    g.discard(waste, "physically removed with aqueous phase/ethanol; no assertion that all off-products are eliminated")
    return kept


def protocol_graph():
    g = ProtocolGraph(TITLE, authority=SOURCE)
    genomic = _genomic()
    cell = g.input("cell", "one Dam-marked cell: separate poly(A) RNA and genomic DNA populations", MaterialKind.BIOLOGICAL,
        state=_mixed(_mrna(), genomic), primary=True, source=ARTICLE)
    cel = g.input("cel_primer", "CEL-seq2 split-barcode example", MaterialKind.OLIGO,
        state=_oligo_scene("CEL-seq2 primer", CEL_PRIMER), source=REAGENTS)
    rt = g.input("cel_rt", "Superscript II", MaterialKind.CATALYST, source=METHODS)
    first = g.product("first_strand", "transcript cDNA/RNA hybrids plus independent gDNA", MaterialKind.MOLECULAR,
        state=_mixed(_transcript_cdna(), genomic), source=METHODS)
    g.react("cel_rt_reaction", "Reverse-transcribe poly(A) RNA; genomic DNA remains in well", StepKind.REVERSE_TRANSCRIPTION,
        ((cell, InputRole.TEMPLATE), (cel, InputRole.STRUCTURAL), (rt, InputRole.CATALYST)), (first,), source=METHODS)
    rna_hydrolysates = _possible_carry("original RNA hydrolysates and residual RNA")
    ss = g.input("second_strand_mix", "DNA polymerase I, E. coli ligase, RNase H", MaterialKind.CATALYST, source=METHODS)
    ds = g.product("ds_templates", "ds transcript cDNA plus independent gDNA", MaterialKind.MOLECULAR,
        state=_mixed(_transcript_cdna(True), genomic, rna_hydrolysates), source=METHODS)
    g.react("second_strand", "RNase H-assisted second-strand synthesis (16 °C 2 h)", StepKind.EXTENSION,
        ((first, InputRole.TEMPLATE), (ss, InputRole.CATALYST)), (ds,), source=METHODS)
    protease = g.input("proteinase_k", "Proteinase K", MaterialKind.CATALYST, source=METHODS)
    deprot = g.product("deproteinated", "deproteinated mixed DNA templates", MaterialKind.MOLECULAR, state=ds.state, source=METHODS)
    g.react("deproteinize", "Proteinase K (50 °C overnight; 80 °C inactivation)", StepKind.DEPROTEINIZATION,
        ((ds, InputRole.TEMPLATE), (protease, InputRole.CATALYST)), (deprot,), source=METHODS)
    dpni = g.input("dpni", "DpnI", MaterialKind.CATALYST, source=METHODS)
    digested = g.product("digested", "ds transcript cDNA plus DpnI-cut genomic fragments", MaterialKind.MOLECULAR,
        state=_mixed(_transcript_cdna(True), dpni_fragment_scene(), rna_hydrolysates), source=METHODS)
    g.react("digest", "DpnI cuts methylated genomic DNA; transcript cDNA retained", StepKind.DIGESTION,
        ((deprot, InputRole.TEMPLATE), (dpni, InputRole.CATALYST)), (digested,), source=METHODS)
    adapter = g.input("dam_adapter", "annealed forked DamID adapter, bottom /5Phos/", MaterialKind.ADAPTER, state=damid_adapter_scene(), source=REAGENTS)
    ligase = g.input("dam_ligase", "T4 DNA ligase", MaterialKind.CATALYST, source=METHODS)
    tagged = g.product("tagged", "ds transcript cDNA and barcoded DamID fragments (same well)", MaterialKind.MOLECULAR,
        state=_mixed(_transcript_cdna(True), ligated_damid_scene(), rna_hydrolysates), source=METHODS)
    g.react("ligate_damid", "Ligate forked adapter to DpnI ends (16 °C overnight)", StepKind.LIGATION,
        ((digested, InputRole.TEMPLATE), (adapter, InputRole.STRUCTURAL), (ligase, InputRole.CATALYST)), (tagged,), geometry=damid_geometry(), source=METHODS)
    other_wells = g.input("other_wells", "other uniquely indexed wells; exact plate/member allocation unassigned",
        MaterialKind.MOLECULAR, state=_mixed(_oligo_scene("other-well transcript cDNA", "X" * 36),
        _oligo_scene("other-well DamID DNA", "X" * 36)), source=METHODS)
    pooled = g.product("pooled_wells", "pooled wells; independent transcript and DamID template populations",
        MaterialKind.MOLECULAR, state=_mixed(tagged.state, other_wells.state), source=METHODS + "; representative member sequences shown for the pooled families")
    g.react("pool_wells", "Pool distinct indexed wells before IVT (no covalent linkage)", StepKind.POOL,
        ((tagged, InputRole.TEMPLATE), (other_wells, InputRole.TEMPLATE)), (pooled,), changes_molecule=False, source=METHODS)
    cleaned = _cleanup(g, pooled, "pre_ivt_cleanup", "pooled RNA/DamID templates",
        "bead purification; exact reagent/dilution unassigned in this capped 2020 reading")
    t7 = g.input("t7", "MEGAscript T7 enzyme mix", MaterialKind.CATALYST, source=METHODS)
    rna_arna, damid_arna = _arna("RNA"), _arna("DamID")
    arna = g.product("arna", "RNA-derived and DamID-derived aRNA populations", MaterialKind.MOLECULAR,
        state=_mixed(rna_arna, damid_arna, cleaned.state), source=METHODS)
    g.react("ivt", "Joint T7 IVT (37 °C overnight); DNA templates remain in mixture", StepKind.TRANSCRIPTION,
        ((cleaned, InputRole.TEMPLATE), (t7, InputRole.CATALYST)), (arna,), source=METHODS)
    cleaned_arna = _cleanup(g, arna, "arna_cleanup", "cleaned amplified RNA", "bead purification; exact reagent choice unassigned in this capped 2020 reading")
    condition = g.input("fragment_buffer", "Tris-acetate/KOAc/MgOAc fragmentation buffer; EDTA stop", MaterialKind.CONDITION, source=METHODS)
    fragments = g.product("fragments", "fragmented aRNA including identifier-bearing and other fragments", MaterialKind.MOLECULAR, state=_mixed(_fragmented(rna_arna, damid_arna), cleaned.state), source=METHODS)
    g.react("fragment", "Fragment aRNA (94 °C 2 min), EDTA stop on ice", StepKind.FRAGMENTATION,
        ((cleaned_arna, InputRole.TEMPLATE), (condition, InputRole.CONDITION)), (fragments,), source=METHODS)
    postfrag = _cleanup(g, fragments, "fragment_cleanup", "12 µL fragmented aRNA eluate", "0.8 volume undiluted AMPure XP; two 80% ethanol washes")
    aliquot = g.product("rt_aliquot", "5 µL fragmented aRNA allocated to RT", MaterialKind.MOLECULAR, state=postfrag.state, source=METHODS)
    rest = g.product("remaining_arna", "remaining fragmented aRNA eluate (QC/remainder allocation unresolved)", MaterialKind.MOLECULAR, state=postfrag.state, source=METHODS)
    g.react("allocate_rt", "Allocate 5 µL from 12 µL fragmented aRNA eluate", StepKind.SPLIT,
        ((postfrag, InputRole.TEMPLATE),), (aliquot, rest), changes_molecule=False, source=METHODS)
    g.mark_terminal(Terminal(rest.id, Disposition.UNRESOLVED, "not consumed as if the complete eluate entered RT; detailed QC/remainder allocation deferred"))
    hexprimer = g.input("randomhex", "randomhexRT", MaterialKind.OLIGO, state=_oligo_scene("randomhexRT", RANDOMHEX_RT), source=REAGENTS)
    rt2 = g.input("rt2", "Superscript II; no correction assigned to published concentration inconsistency", MaterialKind.CATALYST, source=METHODS)
    cdna_state = _mixed(_native_cdna(), aliquot.state, _candidate_cdna())
    cdna = g.product("cdna", "two single-stranded randomhexRT cDNA populations, same tube", MaterialKind.MOLECULAR, state=cdna_state, source=METHODS)
    g.react("randomhex_rt", "Reverse-transcribe fragmented aRNA; retain original RNA and candidate other cDNAs", StepKind.REVERSE_TRANSCRIPTION,
        ((aliquot, InputRole.TEMPLATE), (hexprimer, InputRole.STRUCTURAL), (rt2, InputRole.CATALYST)), (cdna,), source=METHODS)
    p1 = g.input("rp1", "RP1 (terminal phosphorothioate)", MaterialKind.OLIGO, state=_oligo_scene("RP1", RNA_PCR1, mod3="terminal phosphorothioate"), source=REAGENTS)
    pi = g.input("rpi", "RPi CGTGAT example (terminal phosphorothioate)", MaterialKind.OLIGO, state=_oligo_scene("RPi", RNA_PCR_INDEX, mod3="terminal phosphorothioate"), source=REAGENTS)
    enzyme = g.input("pcr_mix", "NEBNext High-Fidelity PCR Master Mix", MaterialKind.CATALYST, source=METHODS)
    product = g.product("indexed", "indexed mixed RNA/DamID molecular populations", MaterialKind.MOLECULAR,
        state=_mixed(_library_scene(RNA, "RNA"), _library_scene(DNA, "DamID"), cdna_state), source=METHODS)
    g.react("pcr", "Terminal indexing PCR on actual ss cDNA; first cycle creates opposite strand", StepKind.PCR,
        ((cdna, InputRole.TEMPLATE), (p1, InputRole.STRUCTURAL), (pi, InputRole.STRUCTURAL), (enzyme, InputRole.CATALYST)), (product,),
        geometry=PCRGeometry.unavailable(cdna_state, reason="Actual incoming state is two ss cDNA populations. RP1 primes copied RA5; RPi initially matches only the 21-nt randomhexRT tail. Exact multi-population ss first-cycle geometry is not represented by this model; no pre-existing full RPi arm or duplex is asserted."), source=METHODS)
    final = _cleanup(g, product, "final_cleanup_1", "indexed library after first purification", "0.8 volume undiluted AMPure XP; two 80% ethanol washes")
    final = _cleanup(g, final, "final_cleanup_2", "indexed library after repeat purification", "repeat 0.8 volume undiluted AMPure XP purification; two 80% ethanol washes")
    # Material classes in one physical library, not a tube split.
    rna = g.product("rna_library", "RNA-derived molecules in mixed physical library", MaterialKind.LIBRARY, state=final.state, source=METHODS)
    dna = g.product("damid_library", "DamID-derived molecules in mixed physical library", MaterialKind.LIBRARY, state=final.state, source=METHODS)
    g.react("readout_classes", "Track co-resident molecular readout classes (no physical separation)", StepKind.GENERIC,
        ((final, InputRole.TEMPLATE),), (rna, dna), changes_molecule=False, source=METHODS)
    for material, library, heading in ((rna, RNA, "Transcriptome library"), (dna, DNA, "DamID library")):
        g.sequence(material, SequencingPlan(library=library, acquired_roles=("Read 1", "Read 2"),
            limitation=SEQUENCING_UNAVAILABLE[heading] + " NextSeq 500 paired-end 75 nt, 2020 Procedure. Barcode/UMI extraction and mapping filters are computational; no physical molecule removal is asserted."))
    # Curator-written inspected scope. New graph nodes do not auto-enrol as reviewed.
    g.declare_reviewed_scope(("cel_primer", "dam_adapter", "randomhex", "rp1", "rpi"), REAGENTS)
    g.declare_reviewed_scope(("cell",), ARTICLE)
    g.declare_reviewed_scope((
        "cel_rt", "second_strand_mix", "proteinase_k", "dpni", "dam_ligase", "other_wells",
        "pre_ivt_cleanup_beads", "t7", "arna_cleanup_beads", "fragment_buffer",
        "fragment_cleanup_beads", "rt2", "pcr_mix", "final_cleanup_1_beads", "final_cleanup_2_beads",
        "cel_rt_reaction", "second_strand", "deproteinize", "digest", "ligate_damid", "pool_wells",
        "pre_ivt_cleanup", "ivt", "arna_cleanup", "fragment", "fragment_cleanup", "allocate_rt",
        "randomhex_rt", "pcr", "final_cleanup_1", "final_cleanup_2", "readout_classes",
        "rna_library", "damid_library"), METHODS)
    return g.require_publishable()
