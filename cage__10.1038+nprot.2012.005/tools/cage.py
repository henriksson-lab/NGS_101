"""Takahashi et al. cap-trapper CAGE model."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))

import seqprimers as sp
from adapter_ligation import AdapterLigationGeometry
from captrap import CapTrappedHybrid
from chemdraw import (Construct, MolecularState, Scene, Segment, Workflow,
                      feature, oligo, revcomp)
from pcr import PCRGeometry, PCRPrimer
from protocol_graph import (DetailView, InputRole, MaterialKind, ProtocolGraph,
                            SequencingPlan, SourceClaim, SourceDiscrepancy, StepKind)

TITLE = "CAGE — cap analysis of gene expression"
NOTES = "01_cage.html"
SOURCE = ('Defining protocol: <a href="https://doi.org/10.1038/nprot.2012.005">'
          'Takahashi et al., <i>Nature Protocols</i> (2012)</a>.')
SUMMARY = ("Select cDNAs that reach a biotinylated RNA cap, ligate a barcoded 5′ linker, "
           "make the second strand, release a 27-nt 5′ tag with EcoP15I, and ligate the "
           "opposite Illumina linker.")

CAVEAT = ("The barcoded 5′ linker is drawn with its AGA representative of the 15 "
          "Table 1 barcodes. Table 1 labels two lower strands '5′SOL-lower-AGT'; the "
          "first (CTGCTGAGTCTG…) is complementary to the ACT upper linker, so its label "
          "appears to be a typographical error. Neither row is drawn here.")

RT_ECOP = "AAGGTCTATCAGCAG" + "N" * 15
LINKER_CORE = "CCACCGACAGGTTCAGAGTTCTACAG"
BARCODE = "AGA"
ECOP_SITE = "CAGCAG"
FWD = "AATGATACGGCGACCACCGACAGGTTCAGAGTTC"
REV = "CAAGCAGAAGACGGCATACGA"
SEQUENCING = "CGGCGACCACCGACAGGTTCAGAGTTCTACAG"
THREE_UPPER = "NN" + revcomp(REV)
FIVE_UPPER_N6 = LINKER_CORE + BARCODE + ECOP_SITE + "N" * 6
FIVE_UPPER_GN5 = LINKER_CORE + BARCODE + ECOP_SITE + "G" + "N" * 5
FIVE_LOWER = revcomp(LINKER_CORE + BARCODE + ECOP_SITE)
SECOND_STRAND_PRIMER = LINKER_CORE
THREE_LOWER = REV
CAPTURED = CapTrappedHybrid(cdna_reaches_cap=True, cap_biotinylated=True)


FWD_TAIL = FWD[:FWD.index(LINKER_CORE[:8])]
TAG_PAIRED_NT = 25   # 27-nt tag window minus the two 3′-linker NN positions


def _tag_parts() -> list:
    """Tag window plus 3′ linker, top strand.

    Table 1 and the stated 96-bp product fix the arithmetic: the 3′ upper linker's
    5′ NN occupies the last two positions of the 27-nt EcoP15I tag window, opposite
    the cDNA-derived two-base overhang.
    """
    return [Segment("cap-derived tag bases", "X" * TAG_PAIRED_NT, placeholder=True),
            Segment("3′-linker NN (tag positions 26–27)", "N" * 2, placeholder=True),
            Segment("P7-side linker / reverse-primer complement", revcomp(REV), "p7")]


def _barcode_segments() -> list:
    return [Segment("3-nt sample barcode", BARCODE, "cbc",
                    feature=feature("sample_inline", "sample_index", "whitelist",
                                    whitelist="published CAGE linker set")),
            Segment("EcoP15I recognition site", ECOP_SITE)]


def pre_pcr_library() -> Construct:
    return Construct([Segment("5′ linker handle", LINKER_CORE, "r1"),
                      *_barcode_segments(), *_tag_parts()],
                     name="dual-linker CAGE tag before PCR")


def final_library() -> tuple[Construct, tuple]:
    # The forward PCR primer adds the P5-side tail; the custom sequencing primer spans
    # the last six tail bases and the 26-base 5′ linker handle.
    if not FWD.endswith(LINKER_CORE[:len(FWD) - len(FWD_TAIL)]):
        raise ValueError("CAGE forward primer no longer overlaps the 5′ linker handle")
    lib = Construct([
        Segment("P5-side PCR tail", FWD_TAIL, "p5"),
        Segment("5′ linker handle (custom sequencing-primer site)", LINKER_CORE, "r1"),
        *_barcode_segments(), *_tag_parts(),
    ], name="96-bp CAGE tag library")
    if len(lib) != 96:
        raise ValueError("CAGE final library no longer matches the published 96 bp")
    primer = sp.custom("Read 1", "CAGE custom sequencing primer", SEQUENCING,
                       "Takahashi et al. 2012, Table 1")
    problems = sp.verify(lib, (primer,), required_roles=("Read 1",))
    if problems:
        raise ValueError("invalid CAGE endpoint: " + "; ".join(problems))
    return lib, (primer,)


FINAL_LIBRARY, SEQ_PRIMERS = final_library()
FINAL_CAPTION = ("The 96-bp single-read product. The custom primer reads the three-base "
                 "sample barcode and EcoP15I site before the 27-nt transcript 5′ tag.")
SEQUENCING_INTRO = ("This protocol explicitly uses its Table 1 custom CAGE sequencing "
                    "primer, not the standard Illumina Read 1 primer.")


def oligos():
    return (
        oligo("RT-N15-EcoP", (Segment("RT-N15-EcoP", RT_ECOP),)),
        oligo("5′ SOL N6, AGA representative",
              (Segment("SOL N6", FIVE_UPPER_N6),), three="/3Phos/-3′"),
        oligo("5′ SOL GN5, AGA representative",
              (Segment("SOL GN5", FIVE_UPPER_GN5),), three="/3Phos/-3′"),
        oligo("5′ SOL lower, AGA", (Segment("SOL lower", FIVE_LOWER),),
              mods="/5Phos/", three="/3AmMO/-3′"),
        oligo("biotinylated second-strand primer",
              (Segment("second-strand primer", SECOND_STRAND_PRIMER),), mods="/5Bio/"),
        oligo("3′ upper linker", (Segment("upper linker", THREE_UPPER),)),
        oligo("3′ lower linker", (Segment("lower linker", THREE_LOWER),)),
        oligo("PCR forward primer", (Segment("PCR forward", FWD),)),
        oligo("PCR reverse primer", (Segment("PCR reverse", REV),)),
        oligo("custom sequencing primer", (Segment("sequencing primer", SEQUENCING),)),
    )


def _rt_scene() -> Scene:
    sc = Scene()
    sc.strand("RNA", [Segment("5-prime cap", "m7G", placeholder=True),
                       Segment("RNA body", "R" * 38, placeholder=True)],
              label="capped RNA")
    sc.anneal("RT primer", [Segment("EcoP15I RT handle", RT_ECOP[:15], "r1"),
                             Segment("random N15", "N" * 15, placeholder=True)],
              to="RNA", pair=("random N15", "RNA body"),
              unpaired=("EcoP15I RT handle",), label="RT-N15-EcoP")
    sc.arrow("RT primer", "reverse transcription toward the cap")
    sc.labels("RT primer")
    return sc


def _captured_scene() -> Scene:
    sc = Scene()
    sc.strand("RNA", [Segment("biotinylated cap", "m7G", "w1", placeholder=True),
                       Segment("cap-protected RNA", "R" * 34, placeholder=True)],
              label="oxidized, biotinylated capped RNA")
    sc.anneal("cDNA", [Segment("cap-complete cDNA", "Y" * 34, placeholder=True),
                        Segment("RT EcoP15I handle", RT_ECOP[:15], "r1")],
              to="RNA", pair=("cap-complete cDNA", "cap-protected RNA"),
              unpaired=("RT EcoP15I handle",), label="cap-complete first strand")
    sc.mark("RNA", "biotinylated cap", "streptavidin-selected cap")
    sc.labels("cDNA")
    return sc


def _five_prime_geometry() -> AdapterLigationGeometry:
    sc = Scene()
    upper = [Segment("upper paired handle", LINKER_CORE),
             Segment("three-base sample barcode", BARCODE, "cbc",
                     feature=feature("sample_inline", "sample_index", "whitelist",
                                     whitelist="published CAGE linker set")),
             Segment("upper EcoP15I site", ECOP_SITE),
             Segment("N6 splint tail", "N" * 6, placeholder=True)]
    sc.strand("upper linker", upper, label="5′ SOL-N6-AGA", mod3="p")
    sc.anneal("lower linker", [Segment("lower ligating strand", FIVE_LOWER)],
              to="upper linker", pair=("lower ligating strand", "upper paired handle"),
              label="5′ SOL-lower-AGA", mod5="p", mod3="/3AmMO/")
    sc.strand("cDNA substrate", [Segment("captured cDNA", "C" * 28, placeholder=True)],
              label="cap-selected single-strand cDNA")
    sc.ligation_end("lower linker", "lower ligating strand", kind="splinted",
                    polarity="5-prime", owner="adapter",
                    text="5′-phosphate aligned by the N6/GN5 splint tail")
    sc.ligation_end("cDNA substrate", "captured cDNA", kind="splinted",
                    polarity="3-prime", owner="substrate", text="cDNA 3′-OH")
    sc.mark("upper linker", "N6 splint tail", "randomly anneals at the cDNA terminus")
    sc.labels("upper linker")
    return AdapterLigationGeometry.known(
        "splinted", sc,
        statement="Table 1 prints a 3′-phosphate-blocked N6/GN5 upper strand and a 5′-phosphorylated, 3′-amino-blocked lower strand; the random upper tail splints the lower 5′ phosphate beside the cDNA 3′-OH.")


def _five_prime_product() -> Scene:
    sc = Scene()
    sc.strand("ligated cDNA", [Segment("captured cDNA", "C" * 28, placeholder=True),
                                Segment("lower-linker-derived handle", FIVE_LOWER, "r1")],
              label="5′-linker-ligated cDNA")
    sc.junction("ligated cDNA", "captured cDNA", "lower-linker-derived handle",
                "splinted ssDNA ligation")
    sc.labels("ligated cDNA")
    return sc


def _second_strand_scene() -> Scene:
    sc = Scene.duplex([
        Segment("custom sequencing handle", LINKER_CORE, "r1"),
        Segment("three-base sample barcode", BARCODE, "cbc",
                feature=feature("sample_inline", "sample_index", "whitelist",
                                whitelist="published CAGE linker set")),
        Segment("EcoP15I site", ECOP_SITE),
        Segment("cap-derived cDNA", "X" * 36, placeholder=True),
        Segment("RT-primer EcoP15I site", RT_ECOP[:15], "r1"),
    ], label="second-strand CAGE duplex")
    sc.mark("top", "custom sequencing handle", "5′-biotinylated second-strand primer")
    sc.labels("top")
    return sc


def _three_prime_geometry() -> AdapterLigationGeometry:
    sc = Scene()
    sc.strand("3′ upper linker", [Segment("two-base 5-prime protrusion", "NN",
                                           placeholder=True),
                                   Segment("upper paired arm", revcomp(REV), "p7")],
              label="3′ upper linker")
    sc.anneal("3′ lower linker", [Segment("lower paired arm", THREE_LOWER, "p7")],
              to="3′ upper linker", pair=("lower paired arm", "upper paired arm"),
              unpaired=("two-base 5-prime protrusion",), label="3′ lower linker")
    sc.overhang("3′ upper linker", "two-base 5-prime protrusion", polarity="5-prime",
                owner="adapter", text="published 5′-NN protrusion")
    sc.strand("EcoP15I tag", [Segment("two-base EcoP15I end", "NN", placeholder=True),
                               Segment("paired tag bases", "X" * TAG_PAIRED_NT,
                                       placeholder=True)],
              label="EcoP15I-cut tag, overhang strand (27 nt of tag)")
    sc.overhang("EcoP15I tag", "two-base EcoP15I end", polarity="5-prime",
                owner="substrate", text="two-base EcoP15I end")
    return AdapterLigationGeometry.known(
        "cohesive", sc,
        statement="Table 1 writes NN at the 5′ end of the 3′ upper linker; that two-base protrusion is the cohesive end used at the EcoP15I-cut tag boundary.")


def _final_scene() -> Scene:
    sc = Scene.duplex(list(FINAL_LIBRARY), label="96-bp CAGE tag")
    sc.junction("top", "EcoP15I recognition site", "cap-derived tag bases",
                "5′ linker / cap-tag boundary")
    sc.junction("top", "cap-derived tag bases", "3′-linker NN (tag positions 26–27)",
                "3′ linker ligation")
    sc.labels("top")
    return sc


def workflow() -> Workflow:
    w = Workflow(MolecularState("Capped RNA plus RT-N15-EcoP", tuple(_rt_scene().rows())))
    w.react("Reverse-transcribe toward the RNA cap, oxidize and biotinylate RNA diols, then cap-trap",
            _captured_scene().rows(), name="Cap-complete RNA/cDNA hybrid",
            note="RNase I removes unprotected RNA; only cDNA that protects the biotinylated cap remains attached for streptavidin capture.")
    w.ligate_adapter("Ligate the barcoded 5′ SOL linker to captured cDNA",
                      _five_prime_geometry(), _five_prime_product().rows(),
                      name="5′-linker-ligated cDNA",
                      note="N6 and GN5 upper linkers are mixed; their random tail supplies the transient splint while only the lower strand is covalently joined.")
    w.react("Prime with biotinylated 2nd SOL and synthesize the second strand",
            _second_strand_scene().rows(), name="Double-stranded CAGE construct",
            note="The exact second-strand primer recreates the custom sequencing handle upstream of the sample barcode and EcoP15I site.")
    w.react("Dephosphorylate the 5′ linker and cut 27 nt into cDNA with EcoP15I",
            _second_strand_scene().rows(), name="27-nt EcoP15I cap tag",
            note="Antarctic phosphatase suppresses linker dimers before EcoP15I creates the two-base tag end.")
    w.ligate_adapter("Ligate the exact 3′ linker and PCR-amplify",
                      _three_prime_geometry(), _final_scene().rows(),
                      name="96-bp CAGE tag library",
                      note="The upper strand's leading NN forms the source-written cohesive end; PCR forward/reverse primers complete the flow-cell-compatible product.")
    return w


def protocol_graph():
    graph=ProtocolGraph(TITLE,authority=SOURCE)
    rna=graph.input("capped_rna","5 µg total RNA (capped and uncapped)",MaterialKind.MOLECULAR,
                    primary=True)
    rt_primer=graph.input("rt_primer","RT-N15-EcoP primer",MaterialKind.OLIGO,state=RT_ECOP)
    rt=graph.input("reverse_transcriptase","PrimeScript reverse transcriptase, sorbitol/trehalose",
                   MaterialKind.CATALYST)
    hybrid=graph.product("hybrid","first-strand RNA/cDNA hybrids",
                         MaterialKind.MOLECULAR,state=tuple(_rt_scene().rows()))
    graph.react("reverse_transcribe","Reverse-transcribe toward the RNA cap",
                StepKind.REVERSE_TRANSCRIPTION,
                ((rna,InputRole.TEMPLATE),(rt_primer,InputRole.STRUCTURAL),
                 (rt,InputRole.CATALYST)),(hybrid,))
    cap_chemistry=graph.input("cap_chemistry","NaIO4 diol oxidation, then long-arm biotin hydrazide",
                              MaterialKind.CONDITION)
    biotinylated=graph.product("biotinylated","hybrids biotinylated at cap and RNA 3′ diols",
                               MaterialKind.MOLECULAR)
    graph.react("biotinylate_cap","Oxidize RNA diols and biotinylate cap and 3′ ends",
                StepKind.CHEMICAL_CONVERSION,
                ((hybrid,InputRole.TEMPLATE),(cap_chemistry,InputRole.CONDITION)),
                (biotinylated,))
    rnase=graph.input("rnase_i","RNase ONE (RNase I)",MaterialKind.CATALYST)
    trimmed=graph.product("rnase_trimmed","hybrids with unprotected single-stranded RNA removed",
                          MaterialKind.MOLECULAR)
    ss_rna=graph.product("ss_rna_fragments","cleaved unprotected RNA, including biotinylated 3′ ends and caps of incomplete cDNAs",
                         MaterialKind.WASTE)
    graph.react("rnase_digest","Digest single-stranded RNA not protected by cDNA",
                StepKind.DIGESTION,((biotinylated,InputRole.TEMPLATE),(rnase,InputRole.CATALYST)),
                (trimmed,ss_rna))
    graph.discard(ss_rna,"removed by RNAClean XP purification before capture")
    streptavidin=graph.input("streptavidin","tRNA-coated MPG Streptavidin beads; four wash buffers",
                            MaterialKind.SELECTION)
    captured=graph.product("captured","bead-bound cap-complete RNA/cDNA hybrids",
                           MaterialKind.MOLECULAR,state=tuple(_captured_scene().rows()))
    incomplete=graph.product("incomplete","unbound non-cap-complete cDNA",MaterialKind.WASTE)
    graph.react("cap_trap","Capture cDNA that protects a biotinylated RNA cap",
                StepKind.AFFINITY_CAPTURE,
                ((trimmed,InputRole.TEMPLATE),(streptavidin,InputRole.SELECTION)),
                (captured,incomplete))
    graph.discard(incomplete,"left in the cleared solution and washes")
    naoh=graph.input("naoh_release","50 mM NaOH, then Tris neutralization; AMPure XP",
                     MaterialKind.CONDITION)
    released=graph.product("ss_cdna","released cap-selected single-strand cDNA",
                           MaterialKind.MOLECULAR)
    bead_rna=graph.product("bead_rna","beads with biotinylated RNA",MaterialKind.WASTE)
    graph.react("release_cdna","Release cap-selected cDNA from the beads with alkali",
                StepKind.STRAND_SEPARATION,
                ((captured,InputRole.TEMPLATE),(naoh,InputRole.CONDITION)),(released,bead_rna))
    graph.discard(bead_rna,"beads are separated magnetically; only the supernatant is kept")
    linker5=graph.input("five_prime_linker","barcoded 5′ SOL linker, N6:GN5 mixture",
                        MaterialKind.ADAPTER,state=tuple(_five_prime_geometry().rows()))
    graph.annotate(SourceDiscrepancy(
        "five_prime_linker","Published N6:GN5 linker ratio differs",
        (SourceClaim("Takahashi et al. 2012, Introduction, 'Providing a priming site for the second strand cDNA synthesis'",
                     "The GN5 linker is used at a ratio of 5:1 with the random N6."),
         SourceClaim("Takahashi et al. 2012, Box 1, barcoded 5′ linker preparation",
                     "N6 and GN5 linkers are mixed at a final ratio N6:GN5 = 1:4 (3 µg + 3 µg vs 12 µg + 12 µg).")),
        "Draw both linker species without choosing a ratio.",
        "The ratio changes linker composition, not the drawn junction geometry.",
        "Which ratio was used for the reported libraries is not stated."))
    ligase5=graph.input("five_prime_ligase","DNA Ligation Kit Mighty Mix (PEG), 16 h",
                        MaterialKind.CATALYST)
    linked=graph.product("linked","5′-linker-ligated cDNA",MaterialKind.MOLECULAR,
                         state=tuple(_five_prime_product().rows()))
    graph.react("ligate_five_prime","Ligate the barcoded SOL linker to captured cDNA",
                StepKind.LIGATION,
                ((released,InputRole.TEMPLATE),(linker5,InputRole.STRUCTURAL),
                 (ligase5,InputRole.CATALYST)),(linked,),geometry=_five_prime_geometry())
    ampure5=graph.input("linker_ampure","two AMPure XP purifications; barcoded samples pooled",
                        MaterialKind.SELECTION)
    pooled=graph.product("pooled","purified, pooled 5′-linker-ligated cDNA",MaterialKind.MOLECULAR,
                         state=linked.state)
    free_linker=graph.product("free_five_linker","excess 5′ linker",MaterialKind.WASTE)
    graph.react("purify_pool","Purify twice to remove excess 5′ linker and pool barcodes",
                StepKind.CLEANUP,((linked,InputRole.TEMPLATE),(ampure5,InputRole.SELECTION)),
                (pooled,free_linker))
    graph.discard(free_linker,"removed with the AMPure supernatants")
    second_primer=graph.input("second_primer","5′-biotinylated 2nd SOL primer",
                              MaterialKind.OLIGO,state=SECOND_STRAND_PRIMER)
    polymerase=graph.input("second_polymerase","TaKaRa LA Taq; 94/42/68 °C",
                           MaterialKind.CATALYST)
    duplex=graph.product("duplex","double-stranded CAGE construct",MaterialKind.MOLECULAR,
                         state=tuple(_second_strand_scene().rows()))
    graph.react("second_strand","Synthesize the second CAGE strand",StepKind.EXTENSION,
                ((pooled,InputRole.TEMPLATE),(second_primer,InputRole.STRUCTURAL),
                 (polymerase,InputRole.CATALYST)),(duplex,))
    phosphatase=graph.input("phosphatase","Antarctic phosphatase; then AMPure XP",
                            MaterialKind.CATALYST)
    dephosphorylated=graph.product("dephosphorylated","5′-linker-dephosphorylated duplex",
                                   MaterialKind.MOLECULAR,state=duplex.state)
    graph.react("dephosphorylate","Dephosphorylate the 5′ linker",
                StepKind.DEPHOSPHORYLATION,
                ((duplex,InputRole.TEMPLATE),(phosphatase,InputRole.CATALYST)),
                (dephosphorylated,))
    ecop=graph.input("ecop15i","EcoP15I with sinefungin, 37 °C 3 h; MgCl2, 65 °C inactivation",
                     MaterialKind.CATALYST)
    tag=graph.product("cap_tag","27-nt EcoP15I cap tag",MaterialKind.MOLECULAR,
                      state=tuple(_three_prime_geometry().rows()))
    long_piece=graph.product("long_piece","distal non-biotinylated EcoP15I cleavage product",
                             MaterialKind.MOLECULAR)
    graph.react("digest","Cut 27 nt into cap-derived cDNA",StepKind.DIGESTION,
                ((dephosphorylated,InputRole.TEMPLATE),(ecop,InputRole.CATALYST)),
                (tag,long_piece))
    graph.annotate(DetailView(
        "cap_tag","27-nt tag window",
        "The paper states a 27-nt tag and a 96-bp final product. With the Table 1 "
        "linkers, that length places the 3′ linker's 5′ NN at tag positions 26–27, "
        "opposite the cDNA-derived two-base overhang; this placement is derived."))
    linker3=graph.input("three_prime_linker","two-base-overhang 3′ linker",
                        MaterialKind.ADAPTER,state=tuple(_three_prime_geometry().rows()))
    ligase3=graph.input("three_prime_ligase","T4 DNA ligase, 16 °C 16 h",
                        MaterialKind.CATALYST)
    ligated=graph.product("ligated","dual-linker CAGE tag",MaterialKind.MOLECULAR,
                          state=pre_pcr_library())
    graph.react("ligate_three_prime","Ligate the exact 3′ linker",
                StepKind.LIGATION,
                ((tag,InputRole.TEMPLATE),(linker3,InputRole.STRUCTURAL),
                 (ligase3,InputRole.CATALYST)),(ligated,),geometry=_three_prime_geometry())
    mpg2=graph.input("tag_streptavidin","tRNA-coated MPG Streptavidin beads; four wash buffers and water",
                     MaterialKind.SELECTION)
    bead_tags=graph.product("bead_tags","bead-bound biotinylated CAGE tags (PCR template)",
                            MaterialKind.MOLECULAR,state=pre_pcr_library())
    unbound=graph.product("unbound_linkers","3′-linker dimers, excess 3′ linker and unbiotinylated cleavage products",
                          MaterialKind.WASTE)
    graph.react("remove_three_linker","Capture biotinylated tags; wash away 3′ linkers",
                StepKind.AFFINITY_CAPTURE,
                ((ligated,InputRole.TEMPLATE),(long_piece,InputRole.TEMPLATE),
                 (mpg2,InputRole.SELECTION)),(bead_tags,unbound))
    graph.discard(unbound,"left in solution and washes")
    forward=graph.input("pcr_forward","published CAGE forward PCR primer",
                        MaterialKind.OLIGO,state=FWD)
    reverse=graph.input("pcr_reverse","published CAGE reverse PCR primer",
                        MaterialKind.OLIGO,state=REV)
    pcr_polymerase=graph.input("pcr_polymerase","Phusion High-Fidelity DNA polymerase",
                               MaterialKind.CATALYST)
    amplified=graph.product("amplified",FINAL_LIBRARY.name,MaterialKind.MOLECULAR,
                            state=FINAL_LIBRARY)
    graph.react("library_pcr","Test cycle number, then bulk-PCR the 96-bp CAGE tags",StepKind.PCR,
                ((bead_tags,InputRole.TEMPLATE),(forward,InputRole.STRUCTURAL),
                 (reverse,InputRole.STRUCTURAL),(pcr_polymerase,InputRole.CATALYST)),
                (amplified,),geometry=PCRGeometry.exact(
                    pre_pcr_library(),PCRPrimer("PCR Forward",FWD,anneal_nt=len(FWD)-len(FWD_TAIL)),
                    PCRPrimer("PCR Reverse",REV),
                    statement="Table 1 primers: the forward primer's 3′ 21 nt match the 5′ linker "
                              "handle and its 5′ tail completes the P5 side; the reverse primer "
                              "equals the 3′ lower linker."))
    exo=graph.input("exonuclease_i","Exonuclease I, then QIAquick column",MaterialKind.CATALYST)
    library=graph.product("library",FINAL_LIBRARY.name,MaterialKind.LIBRARY,
                          state=FINAL_LIBRARY)
    primers_waste=graph.product("degraded_primers","degraded single-stranded PCR primers",
                                MaterialKind.WASTE)
    graph.react("exo_cleanup","Degrade excess PCR primers and purify",StepKind.DIGESTION,
                ((amplified,InputRole.TEMPLATE),(exo,InputRole.CATALYST)),(library,primers_waste))
    graph.discard(primers_waste,"removed by QIAquick purification")
    graph.sequence(library,SequencingPlan(library=FINAL_LIBRARY,
                                          primers=tuple(SEQ_PRIMERS)))
    graph.declare_reviewed_scope((
        "capped_rna","rt_primer","reverse_transcriptase","reverse_transcribe",
        "cap_chemistry","biotinylate_cap","rnase_i","rnase_digest","streptavidin","cap_trap",
        "naoh_release","release_cdna","five_prime_linker","five_prime_ligase","ligate_five_prime",
        "linker_ampure","purify_pool","second_primer","second_polymerase","second_strand",
        "phosphatase","dephosphorylate","ecop15i","digest","three_prime_linker",
        "three_prime_ligase","ligate_three_prime","tag_streptavidin","remove_three_linker",
        "pcr_forward","pcr_reverse","pcr_polymerase","library_pcr","exonuclease_i",
        "exo_cleanup","library"),SOURCE)
    return graph.require_publishable()
