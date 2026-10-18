"""Protocol-owned source, oligos, molecular workflow and reaction graph."""
import illumina as il
import seqprimers as sp
from adapter_ligation import AdapterLigationGeometry
from chemdraw import Construct, Scene, Segment, complement_segments, feature, revcomp
from pcr import PCRGeometry, PCRPrimer
from protocol_graph import (AnnotationSlot, Highlight, InputRole, MaterialKind,
                            ProtocolGraph, SequencingPlan, StepKind)
from protocol_model import chain, molecule, role

SITE_A1_F = "GTTGACATGCTGGATTGAGACTTCCTACACTCTTTCCCTACACGACGCTCTTCCGATCT"

SITE_A1_R = "GATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGGAAGTCTCAATCCAGCATGTCAAC"

SITE_A2_R = "ACTACTAATACGACT"

SITE_A2_F7 = "GTCGTATTAGTAGTANNNNNNNAGATCGGAAGAGCACACGTCTGAACTCC"

SITE_A2_F6 = "GTCGTATTAGTAGTANNNNNNAGATCGGAAGAGCACACGTCTGAACTCC"

SITE_A2_F5 = "GTCGTATTAGTAGTANNNNNAGATCGGAAGAGCACACGTCTGAACTCC"

SITE_RECOVERY_F = "GGAGTTCAGACGTGTGCTC"

SITE_RECOVERY_R = "GTTGACATGCTGGATTGAGACTTC"

SITE_INDEX_F = "AATGATACGGCGACCACCGAGATCTACACNNNNNNNNACACTCTTTCCCTACACGACG"

SITE_INDEX_R = "CAAGCAGAAGACGGCATACGAGATTANNNNNNNNGACTGGAGTTCAGACGTGTGCTC"


def workflow():
    return chain(molecule("purified HMW genomic DNA", (role("intact genomic DNA"),)), (
        ("Cleave naked genomic DNA with Cas9–sgRNA in vitro",
         molecule("Cas9-digested genome", (role("genomic DNA"), role("programmable Cas9 cleavage ends"))),
         "A concentration series of each sgRNP exposes biochemically accessible on- and off-target sites; this is not a measurement of cleavage inside cells."),
        ("Add one 3-prime dA to each Cas9 cleavage end",
         tailed_cut_scene(),
         "The Cas9-created ends are shown in their cohesive ligation-ready state."),
        ("Ligate biotinylated Adapter 1 to the dA-tailed cleavage ends",
         molecule("cut-end-labelled genome", (role("5-prime-biotin Adapter 1"), role("Cas9-cleavage-adjacent genomic DNA")),
                  junctions=(("5-prime-biotin Adapter 1", "Cas9-cleavage-adjacent genomic DNA", "T/A ligation"),)),
         "Adapter 1 is installed before fragmentation, so its biotin identifies the Cas9-created end."),
        ("Fragment with dsDNA Fragmentase to about 200–800 bp",
         molecule("labelled genomic fragments", (role("biotinylated Cas9 end"), role("cut-to-random genomic fragment"), role("Fragmentase-created distal end"))),
         "The new distal end is physically distinct from the already labelled cleavage end."),
        ("End-repair and ligate the pooled N5/N6/N7 Adapter 2",
         molecule("dual-adapter SITE-Seq fragment", (role("biotinylated Adapter 1"), role("cut-to-random genomic DNA"), role("5/6/7-nt diversity spacer + Adapter 2")),
                  junctions=(("biotinylated Adapter 1", "cut-to-random genomic DNA", "cleavage-end ligation"),
                             ("cut-to-random genomic DNA", "5/6/7-nt diversity spacer + Adapter 2", "distal ligation"))),
         "Pooling three adapter lengths creates base diversity; the paper does not define this tract as a molecule-counting UMI."),
        ("Capture biotinylated fragments on streptavidin beads",
         molecule("enriched SITE-Seq fragments", (role("Adapter 1-labelled Cas9 end"), role("adjacent genomic DNA"), role("Adapter 2 distal end"))),
         "Only molecules retaining the Adapter 1-labelled cleavage end survive the affinity selection."),
        ("Recovery PCR, then dual-index PCR",
         molecule("SITE-Seq sequencing library", (role("P5"), role("i5 sample index"), role("Read 1 site"), role("DNA beginning at a Cas9 cut"), role("5/6/7-nt diversity spacer"), role("Index 1 site"), role("i7 sample index"), role("P7"))),
         "The two published PCR rounds first recover captured molecules and then install both flow-cell ends and eight-base sample indexes."),
    ))


def _site_seq_library() -> tuple[Construct, tuple]:
    """Indexed product, with the N7 Adapter 2 branch shown explicitly.

    Every constant arm is sliced out of a published oligo rather than retyped: the P5
    end comes from Index Forward, the distal end from Adapter 2 N7 Forward extended by
    the Index Reverse tail (reverse-complemented).  That tail is what makes the i7 index
    sit two bases closer to the Read-2 arm than in a stock TruSeq library.
    """
    if SITE_INDEX_F[:len(il.P5)] != il.P5 or not SITE_A1_F.endswith(il.TRUSEQ_READ1):
        raise ValueError("SITE-Seq P5/Read 1 arms no longer match the published oligos")
    index_r_top = revcomp(SITE_INDEX_R)              # the strand the Index Reverse primer writes
    a2_distal = SITE_A2_F7[22:]                      # AGATCGGAAGAGC… Read 2 / Index 1 arm
    overlap = len(a2_distal) - a2_distal.index(index_r_top[:19])
    tail = index_r_top[overlap:index_r_top.index("N")]   # bases added 5′ of the i7 index
    after_index = index_r_top[index_r_top.rindex("N") + 1:-len(il.P7_RC)]
    if not tail or not after_index or not index_r_top.endswith(il.P7_RC):
        raise ValueError("Index Reverse no longer resolves into arm/index/P7 parts")
    lib = Construct([
        Segment("P5", il.P5, "p5"),
        Segment("i5 sample index", "I" * 8, "cbc", placeholder=True,
                feature=feature("sample_index_i5", "sample_index", "unknown")),
        Segment("Read 1 arm", il.TRUSEQ_READ1, "r1"),
        role("genomic DNA beginning at the in-vitro Cas9 cleavage end"),
        Segment("Adapter 2 constant arm", SITE_A2_F7[:15]),
        Segment("Adapter 2 length-diversity spacer (N7 branch)", "N" * 7, "cbc", placeholder=True),
        Segment("Index 1 + Read 2 arm", a2_distal, "r2"),
        Segment("Index Reverse tail", tail, "r2"),
        Segment("i7 sample index reverse complement", "I" * 8, "cbc", placeholder=True,
                feature=feature("sample_index_i7", "sample_index", "unknown")),
        Segment("Index Reverse spacer", after_index, "p7"),
        Segment("P7 reverse complement", il.P7_RC, "p7"),
    ], name="SITE-Seq library")
    primers = (sp.TRUSEQ["R1"],
               sp.mismatching(sp.TRUSEQ["I1"], (
                   "the published Index Reverse primer places the i7 index directly after "
                   "…GAACTCCAGTC, so the stock TruSeq Index 1 primer's 3′ AC has no site")),
               sp.TRUSEQ["I2"])
    problems = sp.verify(lib, primers, required_roles=tuple(p.role for p in primers))
    if problems:
        raise ValueError("SITE-Seq final library: " + "; ".join(problems))
    return lib, primers


TITLE = 'SITE-Seq — in-vitro Cas9 cleavage-site sequencing'
NOTES = '01_site_seq.html'
SOURCE = '<a href="https://doi.org/10.1038/nmeth.4284">Cameron et al., Nature Methods (2017)</a>; exact chemistry from the published supplementary protocol and Supplementary Table 26.'
SUMMARY = 'Cas9–sgRNA cleaves purified genomic DNA in vitro; a biotinylated adapter labels those ends before enzymatic fragmentation, allowing streptavidin enrichment and Read 1 entry directly beside each candidate cut site.'
CAVEAT = 'The i7 end is drawn exactly as the published Index Reverse primer writes it, which places the sample index two bases closer to the Read 2 arm than a stock TruSeq library; the stock Index 1 sequencing primer therefore has no exact site here.'
FINAL_CAPTION = 'N7 Adapter 2 branch; N5 and N6 differ only in diversity-spacer length. Adapter 1 defines the Cas9-cleavage end, while Adapter 2 occupies the later Fragmentase end. Every constant arm is sliced from a published oligo: the distal arms come from Adapter 2 N7 Forward and from the reverse complement of Index Reverse.'
READ_LENGTHS = {'Read 1': 150}
SEQUENCING_INTRO = 'The supplementary protocol recommends 2–3 million 150-base single-end reads per sample (section O). Both published index primers carry eight variable bases. Read 1 enters the genome from the Adapter 1 arm at the Cas9 cleavage end; the i5 index is read off the P5 arm, while the i7 index read is shown below without an exact stock-primer site for the reason given above.'
SEQUENCING_ENDING = ''
SEQUENCING_UNAVAILABLE = ''
FINAL_LIBRARY, SEQ_PRIMERS = _site_seq_library()
def oligos():
    return (f'Adapter 1 Forward — 5′-/5Biosg/{SITE_A1_F[:-1]}*{SITE_A1_F[-1]}-3′', f'Adapter 1 Reverse — 5′-{SITE_A1_R}-3′', f'Adapter 2 N7 Forward — 5′-/5Phos/{SITE_A2_F7}-3′', f'Adapter 2 N6 Forward — 5′-/5Phos/{SITE_A2_F6}-3′', f'Adapter 2 N5 Forward — 5′-/5Phos/{SITE_A2_F5}-3′', f'Adapter 2 Reverse — 5′-{SITE_A2_R[:-1]}*{SITE_A2_R[-1]}-3′', f'Recovery PCR Forward — 5′-{SITE_RECOVERY_F}-3′', f'Recovery PCR Reverse — 5′-{SITE_RECOVERY_R}-3′', f'Index Forward — 5′-{SITE_INDEX_F}-3′', f'Index Reverse — 5′-{SITE_INDEX_R}-3′')

def adapter1_scene():
    sc=Scene();sc.strand("A1F",[Segment("Adapter1 paired core",SITE_A1_F[:-1]),Segment("A1 terminalT",SITE_A1_F[-1])],mod5="biotin")
    sc.anneal("A1R",[Segment("Adapter1 reverse",SITE_A1_R)],to="A1F",pair=("Adapter1 reverse","Adapter1 paired core"))
    sc.ligation_end("A1F","A1 terminalT",kind="cohesive",polarity="3-prime",owner="adapter",text="3′T; terminal phosphorothioate")
    return sc


def adapter2_scene(forward):
    sc=Scene();sc.strand("A2F",[Segment("A2 paired core",forward[:14]),Segment("A2 unpaired tail",forward[14:])],mod5="p")
    sc.anneal("A2R",[Segment("A2 reverse paired core",SITE_A2_R[:-1]),Segment("A2 terminalT",SITE_A2_R[-1])],
        to="A2F",pair=("A2 reverse paired core","A2 paired core"),unpaired=("A2 unpaired tail",))
    sc.ligation_end("A2R","A2 terminalT",kind="cohesive",polarity="3-prime",owner="adapter",text="3′T; terminal phosphorothioate")
    return sc


def tailed_cut_scene():
    from endprep import EndPreparedDuplex,dA_tailed_scene
    return dA_tailed_scene(EndPreparedDuplex(Construct([role("Cas9-cleavage-adjacent genomic DNA")]),
        five_prime_phosphate=True,three_prime_overhang="A",overhang_inferred=False))


def recovery_template():
    # The N7 branch is representative of the physical N5/N6/N7 adapter pool.
    return Construct([Segment("Adapter1 forward",SITE_A1_F,"r1"),
        role("cut-to-random genomic DNA"),Segment("distal dA junction","A",inferred=True),
        Segment("Adapter2 N7 forward",SITE_A2_F7,"r2")],name="N7 adapter-ligated recovery-PCR template")


def protocol_graph():
    source = SOURCE
    graph = ProtocolGraph(TITLE, authority=source)
    dna = graph.input("dna", "purified high-molecular-weight genomic DNA",
                      MaterialKind.BIOLOGICAL,
                      state=molecule("purified HMW genomic DNA", (role("intact genomic DNA"),)),
                      primary=True)
    sgrna = graph.input("sgrna", "programmable sgRNA", MaterialKind.OLIGO)
    cas9 = graph.input("cas9", "Cas9 nuclease", MaterialKind.CATALYST)
    cut = graph.product("cut", "Cas9-cleaved genomic DNA", MaterialKind.MOLECULAR,
                        state=dna.state)
    graph.react("cleave", "Cleave naked genomic DNA with Cas9–sgRNA (37 °C, 4 h, RNP 0.25–1024 nM)",
                StepKind.GENOME_EDITING,
                ((dna, InputRole.TEMPLATE), (sgrna, InputRole.STRUCTURAL),
                 (cas9, InputRole.CATALYST)), (cut,))
    protease = graph.input("protease", "Proteinase K / RNase A termination mix",
                           MaterialKind.CATALYST)
    quenched = graph.product("quenched", "deproteinized Cas9-cleaved genomic DNA",
                             MaterialKind.MOLECULAR, state=cut.state)
    graph.react("terminate", "Terminate the cleavage reaction (37 °C 20 min, 55 °C 20 min)",
                StepKind.DEPROTEINIZATION,
                ((cut, InputRole.TEMPLATE), (protease, InputRole.CATALYST)), (quenched,))
    tailing = graph.input("tailing_enzyme", "Klenow Fragment (NEBNext dA-tailing module)",
                          MaterialKind.CATALYST)
    tailed = graph.product("tailed", "dA-tailed Cas9 cleavage ends",
                           MaterialKind.MOLECULAR, state=tailed_cut_scene())
    graph.react("tail", "Add one 3-prime dA at Cas9 cleavage ends (37 °C, 30 min)",
                StepKind.TAILING,
                ((quenched, InputRole.TEMPLATE), (tailing, InputRole.CATALYST)), (tailed,))
    adapter1 = graph.input("adapter1", "biotinylated 3′-T Adapter 1 duplex (Forward/Reverse annealed 95 °C → RT)",
                           MaterialKind.ADAPTER, state=adapter1_scene())
    ligase1 = graph.input("ligase1", "NEB Quick Ligase", MaterialKind.CATALYST)
    labelled = graph.product("labelled", "Adapter-1-labelled Cas9 ends",
                             MaterialKind.MOLECULAR, state=tailed.state)
    graph.react("ligate_adapter1", "Ligate biotinylated Adapter 1 at cleavage ends (20 °C 30 min, then 16 °C overnight)",
                StepKind.LIGATION,
                ((tailed, InputRole.TEMPLATE), (adapter1, InputRole.STRUCTURAL),
                 (ligase1, InputRole.CATALYST)), (labelled,),
                geometry=AdapterLigationGeometry.known(
                    "cohesive",adapter1_scene(),
                    statement="Adapter 1 labels the dA-tailed Cas9 end before fragmentation; the biotinylated Forward strand ends …CCGATC*T, so the duplex presents a 3′-T overhang (derived by pairing the two printed strands)."))
    spri2 = graph.input("spri2", "0.5× SPRIselect reagent", MaterialKind.SELECTION)
    cleaned = graph.product("cleaned", "adapter-1-ligated DNA after SPRI 2",
                            MaterialKind.MOLECULAR, state=labelled.state)
    free_adapter = graph.product("free_adapter", "residual free Adapter 1",
                                 MaterialKind.WASTE)
    graph.react("spri_2", "SPRI 2 cleanup: remove residual Adapter 1", StepKind.CLEANUP,
                ((labelled, InputRole.TEMPLATE), (spri2, InputRole.SELECTION)),
                (cleaned, free_adapter))
    graph.discard(free_adapter, "removed on the 0.5× SPRI supernatant so it cannot form adapter dimer")
    fragmentase = graph.input("fragmentase", "NEBNext dsDNA Fragmentase", MaterialKind.CATALYST)
    fragments = graph.product("fragments", "Adapter-1-labelled 200–800-bp fragments",
                              MaterialKind.MOLECULAR, state=cleaned.state)
    graph.react("fragment", "Fragment labelled genomic DNA (37 °C, 1 h, EDTA quench)",
                StepKind.FRAGMENTATION,
                ((cleaned, InputRole.TEMPLATE), (fragmentase, InputRole.CATALYST)),
                (fragments,))
    repair = graph.input("distal_repair", "NEBNext Ultra End Repair / dA-tailing enzyme mix",
                         MaterialKind.CATALYST)
    repaired = graph.product("repaired", "end-repaired distal fragments", MaterialKind.MOLECULAR,
                             state=fragments.state)
    graph.react("repair_distal", "End-repair Fragmentase-created distal ends (20 °C 30 min, 65 °C 30 min)",
                StepKind.END_REPAIR,
                ((fragments, InputRole.TEMPLATE), (repair, InputRole.CATALYST)), (repaired,))
    adapter2 = graph.input("adapter2", "pooled N5/N6/N7 Adapter 2 duplex (annealed 95 °C → RT)",
                           MaterialKind.ADAPTER,
                           state=tuple(adapter2_scene(x) for x in (SITE_A2_F7,SITE_A2_F6,SITE_A2_F5)))
    ligase2 = graph.input("ligase2", "NEBNext Ultra Blunt/TA Ligase Master Mix + ligation enhancer",
                          MaterialKind.CATALYST)
    dual = graph.product("dual", "dual-adapter SITE-Seq fragments", MaterialKind.MOLECULAR,
                         state=recovery_template())
    graph.react("ligate_adapter2", "Ligate the length-diverse distal Adapter 2",
                StepKind.LIGATION,
                ((repaired, InputRole.TEMPLATE), (adapter2, InputRole.STRUCTURAL),
                 (ligase2, InputRole.CATALYST)), (dual,),
                geometry=AdapterLigationGeometry.known(
                    "cohesive",*(adapter2_scene(x) for x in (SITE_A2_F7,SITE_A2_F6,SITE_A2_F5)),
                    statement="The 15-base Adapter 2 Reverse strand pairs with the first 14 bases of each Forward strand and leaves its terminal T unpaired, so the duplex presents a 3′-T overhang (derived by pairing the printed strands). The protocol prepares the distal ends with the NEBNext Ultra End Repair / dA-tailing module and seals with Blunt/TA Ligase, consistent with a T/A junction; the end-repair step is not separately described as dA-tailing."))
    streptavidin = graph.input("streptavidin", "Dynabeads M-280 streptavidin",
                               MaterialKind.SELECTION)
    captured = graph.product("captured", "biotin-enriched SITE-Seq fragments",
                             MaterialKind.MOLECULAR, state=dual.state)
    rejected = graph.product("rejected", "Adapter-1-negative fragments", MaterialKind.WASTE)
    graph.react("capture", "Capture Adapter-1-labelled fragments",
                StepKind.AFFINITY_CAPTURE,
                ((dual, InputRole.TEMPLATE), (streptavidin, InputRole.SELECTION)),
                (captured, rejected))
    graph.discard(rejected, "does not retain the cleavage-end Adapter 1")
    rec_f = graph.input("recovery_forward", "recovery PCR forward primer", MaterialKind.OLIGO,
                        state=SITE_RECOVERY_F)
    rec_r = graph.input("recovery_reverse", "recovery PCR reverse primer", MaterialKind.OLIGO,
                        state=SITE_RECOVERY_R)
    polymerase = graph.input("polymerase", "Q5 Hot Start High-Fidelity 2× master mix",
                             MaterialKind.CATALYST)
    recovered = graph.product("recovered", "recovery-amplified SITE-Seq fragments",
                              MaterialKind.MOLECULAR, state=recovery_template())
    graph.react("recovery_pcr", "On-bead recovery PCR (12 cycles, 61 °C anneal)", StepKind.PCR,
                ((captured, InputRole.TEMPLATE), (rec_f, InputRole.STRUCTURAL),
                 (rec_r, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)),
                (recovered,), geometry=PCRGeometry.exact(captured.state,
                    PCRPrimer("published Recovery Reverse (construct-left primer)",SITE_RECOVERY_R),
                    PCRPrimer("published Recovery Forward",SITE_RECOVERY_F),
                    statement="Published oligo names preserved: Recovery Reverse occupies the construct-left site and Recovery Forward the construct-right site on this displayed top orientation; sites are computed on the actual adapter-ligated N7 branch."))
    reviewed = ["dna", "sgrna", "cas9", "cleave", "protease", "terminate",
                "tailing_enzyme", "tail",
                "adapter1", "ligase1", "ligate_adapter1", "spri2", "spri_2",
                "cleaned", "free_adapter", "fragmentase", "fragment",
                "distal_repair", "repair_distal", "adapter2", "ligase2", "ligate_adapter2",
                "streptavidin", "capture", "recovery_forward", "recovery_reverse",
                "polymerase", "recovery_pcr"]
    index_f=graph.input("index_forward","published Index Forward",MaterialKind.OLIGO,state=SITE_INDEX_F)
    index_r=graph.input("index_reverse","published Index Reverse",MaterialKind.OLIGO,state=SITE_INDEX_R)
    index_pol=graph.input("index_polymerase","Q5 2× Master Mix",MaterialKind.CATALYST)
    library=graph.product("index_library","SITE-Seq indexed sequencing library",MaterialKind.LIBRARY,state=FINAL_LIBRARY)
    graph.react("index_pcr","Indexing PCR on diluted recovery product (12 cycles,60°C anneal)",StepKind.PCR,
        ((recovered,InputRole.TEMPLATE),(index_f,InputRole.STRUCTURAL),(index_r,InputRole.STRUCTURAL),
         (index_pol,InputRole.CATALYST)),(library,),geometry=PCRGeometry.exact(recovered.state,
            PCRPrimer("Index Forward",SITE_INDEX_F),PCRPrimer("Index Reverse",SITE_INDEX_R),
            statement="Both index-primer3′ sites are located on the recovery amplicon; the Index Reverse5′ tail creates the noncanonical i7 end displayed in the final product."))
    graph.sequence(library,SequencingPlan(library=FINAL_LIBRARY,primers=SEQ_PRIMERS,
        acquired_roles=("Read 1",),read_lengths={"Read 1":150},
        limitation="Supplementary protocol recommends150-base single-end reads. Eight variable bases in each index primer do not establish acquired index cycles or the missing exact index-read primer site."))
    reviewed.extend(("index_forward","index_reverse","index_polymerase","index_pcr","index_library"))
    graph.annotate(Highlight(
        "index_pcr", "The published Index Reverse primer shifts the i7 index",
        "Index Reverse (Supplementary Protocol oligo table; Supplementary Table 26 row IR) reads "
        "CAAGCAGAAGACGGCATACGAGAT-TA-[8-base index]-GACTGGAGTTCAGACGTGTGCTC. Its 3′ 19 bases prime on the "
        "Adapter 2 arm, so the finished top strand runs …GAACTCCAGTC-[i7 index]-TA-P7′ instead of the stock "
        "TruSeq …GAACTCCAGTCAC-[i7 index]-P7′. The i7 index therefore sits two bases closer to the Read 2 arm "
        "and the stock TruSeq Index 1 sequencing primer has no exact site on this library; the protocol does not "
        "print a custom index-read primer, so the acquired index cycles remain unverified.",
        slot=AnnotationSlot.AFTER_WORKFLOW))
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()
