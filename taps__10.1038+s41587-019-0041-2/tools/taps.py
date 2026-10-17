"""Liu et al. 2019 TAPS chemistry."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))

import illumina as il
import seqprimers as sp
from base_conversion import BasePath
from adapter_ligation import AdapterLigationGeometry
from chemdraw import Construct, Row, Scene, Segment, complement_segments, feature, revcomp
from endprep import dA_tailed_scene, repair_and_dA_tail
from pcr import PCRGeometry
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            SourceClaim, SourceDiscrepancy, StepKind)



def taps_path(base: str) -> BasePath:
    """TAPS chemistry and PCR readout for one cytosine state."""
    if base in {"5mC", "5hmC"}:
        return BasePath((base, "5caC", "DHU", "T"))
    if base == "C":
        return BasePath(("C", "C"))
    raise ValueError("TAPS input must be C, 5mC or 5hmC")

TITLE = "TAPS — bisulfite-free methylation sequencing"
NOTES = "01_taps.html"
SOURCE = ('Defining source: <a href="https://doi.org/10.1038/s41587-019-0041-2">Liu et al., <i>Nature Biotechnology</i> (2019)</a>; '
          'library procedure from published Supplementary Note 7 and the authors\' accepted-manuscript Online Methods.')
SUMMARY = "TET oxidation and pyridine–borane reduction make 5mC and 5hmC read as T, while unmodified C remains C."
CAVEAT = ("INFERRED regions: Supplementary Note 7 prints only the two ligation-adapter oligos. The P5 end, the "
          "i7-side extension and P7 come from NEBNext Multiplex Oligos for Illumina index primers (E7335S named "
          "as an example) whose sequences and chosen index are not printed in the inspected sources; those bases "
          "are drawn from the standard TruSeq layout and marked inferred. The six-base Index 6 sequence printed in "
          "Supplementary Note 2 belongs to a synthetic validation molecule, not the genomic library.")
INDEX6 = "GCCAAT"  # Synthetic validation molecule only, not the genomic library index.

# Supplementary Note 7, "Preparation of sequencing library with KAPA Hyper kit", step 1:
# Illumina multiplexing adapters without phosphorothioate bonds, IDT HPLC-purified.
ADAPTER_TOP = "ACACTCTTTCCCTACACGACGCTCTTCCGATCT"
ADAPTER_BOTTOM = "GATCGGAAGAGCACACGTCT"   # 5'-/5Phos/
SUPP_NOTE_7 = "Liu et al. 2019 Supplementary Note 7 (published MOESM1)"
METHODS = "Liu et al. 2019 accepted manuscript (Oxford ORA), Methods 'Whole-genome TAPS'"
NEB_LIMIT = ("Supplementary Note 7 names NEBNext Multiplex Oligos for Illumina (e.g. Index Set 1, E7335S) "
             "but prints neither primer sequence nor the index used; the vendor original was not read.")


def _stem_split() -> tuple[str, str]:
    stem = revcomp(ADAPTER_BOTTOM[:12])           # top-strand bases paired by the bottom oligo
    if not ADAPTER_TOP[:-1].endswith(stem):
        raise ValueError("TAPS adapter oligos do not form the expected forked stem")
    return ADAPTER_TOP[:-1 - len(stem)], stem


R1_FORK, R1_STEM = _stem_split()
R2_STEM, R2_FORK = ADAPTER_BOTTOM[:12], ADAPTER_BOTTOM[12:]
INSERT = Segment("TAPS-treated genomic insert", "X" * 36, placeholder=True)


def adapter_scene() -> Scene:
    """Annealed Supplementary Note 7 adapter: 12-bp stem, 3'-T, unpaired forks."""
    sc = Scene()
    sc.strand("adapter top", [Segment("Read 1 fork arm", R1_FORK, "r1"),
                              Segment("adapter stem", R1_STEM, "r1"),
                              Segment("3-prime T", ADAPTER_TOP[-1], "r1")],
              label="adapter top oligo")
    sc.anneal("adapter bottom", [Segment("adapter stem complement", R2_STEM, "r2"),
                                 Segment("Index/Read 2 fork arm", R2_FORK, "r2")],
              to="adapter top", pair=("adapter stem complement", "adapter stem"),
              label="adapter bottom oligo", mod5="p", unpaired=("Index/Read 2 fork arm",))
    sc.ligation_end("adapter top", "3-prime T", kind="cohesive", polarity="3-prime",
                    owner="adapter", text="3′-T ligation overhang")
    return sc


def ligated_scene() -> Scene:
    """One insert carrying the forked adapter at both ends (before PCR adds P5/P7)."""
    top = [Segment("Read 1 fork arm", R1_FORK, "r1"), Segment("adapter stem", R1_STEM, "r1"),
           Segment("adapter 3-prime T", "T", "r1"), INSERT, Segment("insert 3-prime dA", "A"),
           Segment("adapter stem complement", R2_STEM, "r2"),
           Segment("Index/Read 2 fork arm", R2_FORK, "r2")]
    bottom = complement_segments(top)
    bottom[0] = Segment("bottom Read 1 fork arm", R1_FORK, "r1")
    bottom[-1] = Segment("bottom Index/Read 2 fork arm", R2_FORK, "r2")
    sc = Scene(); sc.strand("top", top, label="adapter-ligated molecule")
    sc.anneal("bottom", bottom, to="top", pair=("adapter stem'", "adapter stem"),
              label="adapter-ligated molecule",
              unpaired=("bottom Read 1 fork arm", "bottom Index/Read 2 fork arm"))
    sc.junction("top", "adapter 3-prime T", INSERT.name, "ligation")
    sc.junction("top", "insert 3-prime dA", "adapter stem complement", "ligation")
    return sc


def library():
    lib = Construct([
        Segment("P5 (NEBNext universal primer)", il.P5[:-4], "p5", inferred=True),
        Segment("adapter top oligo / Read 1 site", ADAPTER_TOP, "r1"),
        INSERT,
        Segment("dA junction / Read 2 site start", "A"),
        Segment("adapter bottom oligo", ADAPTER_BOTTOM, "r2"),
        Segment("index-primer extension", il.INDEX1_PRIMER[len(ADAPTER_BOTTOM):], "r2", inferred=True),
        Segment("sample i7 (chosen NEBNext index)", "[i7]", "cbc", placeholder=True, inferred=True,
                feature=feature("sample_i7", "sample_index", "unknown",
                                note="NEBNext Multiplex Oligos kit named, chosen index not reported")),
        Segment("P7 reverse complement", il.P7_RC, "p7", inferred=True),
    ], name="TAPS TruSeq-layout library")
    primers = (sp.TRUSEQ["R1"], sp.TRUSEQ["I1"], sp.TRUSEQ["R2"])
    problems = sp.verify(lib, primers, required_roles=tuple(p.role for p in primers))
    if problems: raise ValueError("invalid TAPS library: " + "; ".join(problems))
    return lib, primers


FINAL_LIBRARY, SEQ_PRIMERS = library()
FINAL_CAPTION = ("Single-index library after TAPS conversion and 5-cycle Uracil+ PCR; P5, the i7-side "
                 "extension and P7 are INFERRED from the named NEBNext primers.")
SEQUENCING_INTRO = ("Paired-end 80-bp reads on a NextSeq 500 (accepted-manuscript Methods); i7 identifies "
                    "the library, but the index-read cycle count is not stated.")
READ_LENGTHS = {"Read 1": 80, "Read 2": 80}
READ_USES = {"Read 1": "include", "Read 2": "include", "Index 1 (i7)": "demultiplex"}


def _conversion_rows(stage: int) -> tuple[Row, ...]:
    return tuple(Row(chunks=[(f"{base}: {taps_path(base).after(stage)}", None if base == "C" else "w1", False)])
                 for base in ("5mC", "5hmC", "C"))


def protocol_graph() -> ProtocolGraph:
    g = ProtocolGraph(TITLE, authority=SUPP_NOTE_7)
    raw = g.input("genomic_dna", "100 ng genomic DNA mixed with control spike-ins",
                  MaterialKind.BIOLOGICAL, state=Scene.duplex([INSERT]), primary=True)
    spikes = g.input("spike_ins", "0.5% CpG-methylated lambda DNA; Methods also add 0.025% unmodified 2-kb control and 0.25% N5mCNN/N5hmCNN oligos",
                     MaterialKind.MOLECULAR, source=METHODS)
    covaris = g.input("sonication", "Covaris M220 sonication in EB buffer to ~200 bp", MaterialKind.CONDITION)
    sheared = g.product("sheared", "sheared genomic DNA plus spike-ins", MaterialKind.MOLECULAR,
                        state=Scene.duplex([INSERT]))
    g.react("fragment", "Sonicate to ~200 bp", StepKind.FRAGMENTATION,
            ((raw, InputRole.TEMPLATE), (spikes, InputRole.SPIKE_IN), (covaris, InputRole.CONDITION)), (sheared,))
    sizer = g.input("size_beads", "0.55X–1X AMPure XP double-sided selection", MaterialKind.SELECTION)
    sized = g.product("sized", "200–400 bp fragments", MaterialKind.MOLECULAR, state=Scene.duplex([INSERT]))
    off_size = g.product("off_size", "fragments outside 200–400 bp", MaterialKind.WASTE)
    g.react("size_select", "Select 200–400 bp fragments", StepKind.SELECTION,
            ((sheared, InputRole.TEMPLATE), (sizer, InputRole.SELECTION)), (sized, off_size))
    g.discard(off_size, "removed by double-sided AMPure XP selection")

    tailed = repair_and_dA_tail(Construct([INSERT]))
    er = g.input("end_repair_enzymes", "KAPA Hyper Prep Kit (KK8504) end-repair and A-tailing module",
                 MaterialKind.CATALYST)
    tailed_m = g.product("da_tailed_dna", "end-repaired, 5′-phosphorylated, 3′-dA-tailed DNA",
                         MaterialKind.MOLECULAR, state=dA_tailed_scene(tailed))
    g.react("end_repair_and_tailing", "End-repair and add 3′ dA", StepKind.TAILING,
            ((sized, InputRole.TEMPLATE), (er, InputRole.CATALYST)), (tailed_m,))
    adapter = g.input("truseq_adapter", "15 µM annealed Illumina multiplexing adapter, no phosphorothioates (IDT, HPLC)",
                      MaterialKind.ADAPTER, state=adapter_scene())
    ligase = g.input("dna_ligase", "KAPA Hyper Prep ligation module; ligation extended to 30 min", MaterialKind.CATALYST)
    ligated = g.product("ligated_library", "adapter-ligated DNA (forked ends)", MaterialKind.MOLECULAR,
                        state=ligated_scene())
    g.react("adapter_ligation", "Ligate the forked multiplexing adapter", StepKind.LIGATION,
            ((tailed_m, InputRole.TEMPLATE), (adapter, InputRole.STRUCTURAL), (ligase, InputRole.CATALYST)),
            (ligated,), geometry=AdapterLigationGeometry.known(
                "cohesive", dA_tailed_scene(tailed), adapter_scene(),
                statement="Supplementary Note 7 prints both adapter oligos; their annealed 12-bp stem leaves a 3′-T that pairs with the insert 3′-dA."))
    lig_beads = g.input("post_ligation_cleanup", "KAPA post-ligation bead clean-up; elute in 20 µl water",
                        MaterialKind.SELECTION)
    ligated_clean = g.product("ligated_clean", "purified adapter-ligated DNA", MaterialKind.MOLECULAR,
                              state=ligated_scene())
    lig_waste = g.product("ligation_waste", "free adapter and ligation reagents", MaterialKind.WASTE)
    g.react("post_ligation_cleanup_step", "Clean up the ligation", StepKind.CLEANUP,
            ((ligated, InputRole.TEMPLATE), (lig_beads, InputRole.SELECTION)), (ligated_clean, lig_waste))
    g.discard(lig_waste, "removed by the kit's post-ligation bead clean-up")

    tet = g.input("tet", "4 µM mTet1CD, Fe(II), α-KG, ascorbate, ATP, DTT; 37 °C 80 min, then 0.8 U proteinase K 50 °C 1 h; two rounds",
                  MaterialKind.CATALYST)
    oxidized = g.product("oxidized_library", "5caC-containing adapter-ligated DNA", MaterialKind.MOLECULAR,
                         state=_conversion_rows(1))
    g.react("tet_oxidation", "Oxidize 5mC/5hmC to 5caC with mTet1CD (two rounds)", StepKind.CHEMICAL_CONVERSION,
            ((ligated_clean, InputRole.TEMPLATE), (tet, InputRole.CATALYST)), (oxidized,))
    ox_clean = g.input("oxidation_cleanup", "Bio-Spin P-30 column then 1.8X AMPure XP after each round; elute 35 µl",
                       MaterialKind.SELECTION)
    oxidized_clean = g.product("oxidized_clean", "purified oxidized DNA", MaterialKind.MOLECULAR,
                               state=_conversion_rows(1))
    ox_waste = g.product("oxidation_waste", "enzyme digest, cofactors and buffer", MaterialKind.WASTE)
    g.react("oxidation_cleanup_step", "Purify oxidized DNA", StepKind.CLEANUP,
            ((oxidized, InputRole.TEMPLATE), (ox_clean, InputRole.SELECTION)), (oxidized_clean, ox_waste))
    g.discard(ox_waste, "removed by gel column and bead purification")

    borane = g.input("pyridine_borane", "1 M pyridine borane, 600 mM sodium acetate pH 4.3; 37 °C 16 h, 850 rpm",
                     MaterialKind.CONDITION)
    reduced = g.product("converted_library", "DHU-containing adapter-ligated DNA", MaterialKind.MOLECULAR,
                        state=_conversion_rows(2))
    g.react("pyridine_borane_reduction", "Reduce 5caC to DHU with pyridine borane", StepKind.CHEMICAL_CONVERSION,
            ((oxidized_clean, InputRole.TEMPLATE), (borane, InputRole.CONDITION)), (reduced,))
    zymo = g.input("zymo_cleanup", "Zymo-Spin IC column with PB buffer", MaterialKind.SELECTION)
    reduced_clean = g.product("converted_clean", "purified TAPS-converted DNA", MaterialKind.MOLECULAR,
                              state=_conversion_rows(2))
    red_waste = g.product("reduction_waste", "pyridine borane and acetate buffer", MaterialKind.WASTE)
    g.react("reduction_cleanup", "Purify converted DNA", StepKind.CLEANUP,
            ((reduced, InputRole.TEMPLATE), (zymo, InputRole.SELECTION)), (reduced_clean, red_waste))
    g.discard(red_waste, "removed by spin-column purification")

    fwd = g.input("library_forward_primer", "NEBNext universal PCR primer (sequence not printed)",
                  MaterialKind.OLIGO, inferred=True)
    rev = g.input("library_reverse_primer", "NEBNext index primer, e.g. Index Set 1 E7335S (sequence and index not printed)",
                  MaterialKind.OLIGO, inferred=True)
    pol = g.input("pcr_polymerase", "KAPA HiFi HotStart Uracil+ ReadyMix, 5 cycles", MaterialKind.CATALYST)
    amplified = g.product("amplified", "amplified library (DHU copied as T)", MaterialKind.MOLECULAR,
                          state=FINAL_LIBRARY)
    g.react("library_pcr", "INFERRED — NEBNext P5/i7/P7 ends: amplify 5 cycles; DHU is copied as T", StepKind.PCR,
            ((reduced_clean, InputRole.TEMPLATE), (fwd, InputRole.STRUCTURAL), (rev, InputRole.STRUCTURAL),
             (pol, InputRole.CATALYST)), (amplified,),
            geometry=PCRGeometry.unavailable(FINAL_LIBRARY, reason=NEB_LIMIT))
    final_beads = g.input("final_cleanup", "1X AMPure XP", MaterialKind.SELECTION)
    final = g.product("taps_library", "amplified TAPS sequencing library", MaterialKind.LIBRARY,
                      state=FINAL_LIBRARY)
    pcr_waste = g.product("pcr_waste", "primers, nucleotides and polymerase", MaterialKind.WASTE)
    g.react("final_cleanup_step", "Purify the library", StepKind.CLEANUP,
            ((amplified, InputRole.TEMPLATE), (final_beads, InputRole.SELECTION)), (final, pcr_waste))
    g.discard(pcr_waste, "removed by 1X AMPure XP purification")
    g.sequence(final, SequencingPlan(FINAL_LIBRARY, tuple(SEQ_PRIMERS), tuple(p.role for p in SEQ_PRIMERS),
                                     read_lengths=READ_LENGTHS, read_uses=READ_USES))
    g.annotate(SourceDiscrepancy(
        "end_repair_and_tailing", "Library-prep kit name",
        (SourceClaim("Accepted manuscript Methods, 'Whole-genome TAPS'",
                     "End-repair/A-tailing and ligation of Illumina Multiplexing adapters were prepared with the KAPA HyperPlus kit."),
         SourceClaim("Published Supplementary Note 7, Reagents and 'Preparation of sequencing library with KAPA Hyper kit'",
                     "Lists the KAPA Hyper Prep Kit (KK8504) and uses it for end-repair, A-tailing and ligation after Covaris sonication.")),
        "The graph labels the module KAPA Hyper Prep KK8504 and models sonication, not enzymatic fragmentation.",
        "Supplementary Note 7 is the published step-by-step protocol, and both sources fragment by Covaris sonication.",
        "The typeset published Methods were not inspected, so whether the published article repeats 'HyperPlus' is unverified; "
        "the end-repair, A-tailing and ligation reagent formulations of the two kits were not compared."))
    g.declare_reviewed_scope([
        "genomic_dna", "spike_ins", "sonication", "fragment", "size_beads", "size_select",
        "end_repair_enzymes", "end_repair_and_tailing", "truseq_adapter", "dna_ligase",
        "adapter_ligation", "post_ligation_cleanup", "post_ligation_cleanup_step", "tet",
        "tet_oxidation", "oxidation_cleanup", "oxidation_cleanup_step", "pyridine_borane",
        "pyridine_borane_reduction", "zymo_cleanup", "reduction_cleanup",
        "library_forward_primer", "library_reverse_primer", "pcr_polymerase", "library_pcr",
        "final_cleanup", "final_cleanup_step", "taps_library"], SUPP_NOTE_7)
    return g.require_publishable()
