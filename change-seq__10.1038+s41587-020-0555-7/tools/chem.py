"""Protocol-owned source, oligos, molecular workflow and reaction graph."""
import illumina as il
import seqprimers as sp
from adapter_ligation import reported_adapter_end
from chemdraw import Construct, Scene, Segment, feature, revcomp
from circular import CircularizationGeometry
from pcr import PCRGeometry, PCRPrimer
from protocol_graph import Disposition, Terminal, InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind
from protocol_model import chain, closed_circle, molecule, role

CHANGE_TRANSPOSON_TOP = "ACGUAGATGTGTATAAGAGACAG"

CHANGE_TRANSPOSON_BOTTOM = "CTGTCTCTTATACACATCTACGT"


def change_library() -> tuple[Construct, tuple]:
    """NEBNext dual-index product recovered outward from a Cas9-opened DNA circle."""
    # Supplementary Protocol names E7600S Set1 and its dual-eight-cycle reads.
    # Variable windows identify the selected kit members without assigning a study well.
    lib=Construct([
        Segment("P5",il.P5,"p5"),
        Segment("i5","I"*8,"cbc",placeholder=True,feature=feature("sample_index_i5","sample_index","unknown")),
        Segment("Read 1 arm",il.TRUSEQ_READ1,"r1"),
        role("genomic DNA beginning exactly at an in-vitro Cas9 cleavage end"),
        Segment("dA junction","A"),Segment("Index 1 / Read 2 arm",il.INDEX1_PRIMER,"r2"),
        Segment("i7 reverse complement","I"*8,"cbc",placeholder=True,feature=feature("sample_index_i7","sample_index","unknown")),
        Segment("P7 reverse complement",il.P7_RC,"p7")],name="CHANGE-seq E7600S library")
    primers=(sp.TRUSEQ["R1"],sp.TRUSEQ["I1"],sp.TRUSEQ["I2"],sp.TRUSEQ["R2"])
    errors=sp.verify(lib,primers,required_roles=tuple(p.role for p in primers))
    if errors: raise ValueError("CHANGE-seq library: "+"; ".join(errors))
    return lib,primers


def custom_transposon():
    sc=Scene()
    sc.strand("oCRL225",[Segment("dU-containing transposon",CHANGE_TRANSPOSON_TOP,
              note="U is deoxyuridine (/ideoxyU/), not an RNA residue")],mod5="p")
    sc.anneal("oCRL226",[Segment("opposite transposon strand",CHANGE_TRANSPOSON_BOTTOM)],
              to="oCRL225",pair=("opposite transposon strand","dU-containing transposon"),mod5="p")
    return sc


def tailed_cut_product():
    from endprep import EndPreparedDuplex,dA_tailed_scene
    return dA_tailed_scene(EndPreparedDuplex(
        Construct([role("Cas9-to-Cas9 genomic fragment")],name="Cas9-linearized circle"),
        five_prime_phosphate=True,three_prime_overhang="A",overhang_inferred=False))


def opened_adapter_template():
    return Construct([Segment("USER-opened Read1 arm",il.NEBNEXT_ARM_READ1,"r1"),
                      role("Cas9-linearized genomic insert"),Segment("dA junction","A"),
                      Segment("USER-opened truncated Read2 arm",il.NEBNEXT_ARM_READ2,"r2")],
                     name="USER-opened NEBNext-adapted Cas9 fragment")


def workflow():
    return chain(molecule("purified high-molecular-weight genomic DNA", (role("genomic DNA"),)), (
        ("Tagment with the custom dU-bearing Tn5 transposome",
         molecule("approximately 400-bp tagmented fragments", (role("palindromic dU/Tn5 end"), role("genomic DNA"), role("palindromic dU/Tn5 end"))),
         "Tn5 installs the same 23-base duplex at both fragment ends; its 19-base mosaic end performs strand transfer and its four-base palindrome is the future ligation end."),
        ("Gap-repair the staggered transposition products",
         molecule("fully duplex tagmented fragments", (role("copied custom end"), role("duplex genomic insert"), role("copied custom end"))),
         "Kapa HiFi HotStart Uracil+ Ready Mix and Taq DNA ligase (72 °C, 30 min) fill and seal the Tn5 gaps before the engineered ends are opened."),
        ("Open dU with USER and phosphorylate with T4 PNK",
         molecule("ligation-ready linear fragments", (role("four-base cohesive end"), role("duplex genomic insert"), role("four-base cohesive end"))),
         "A single USER/T4 PNK reaction (37 °C, 1 h) cleaves at the internal deoxyuridine and supplies ligatable 5-prime phosphates."),
        ("Circularize by dilute intramolecular ligation",
         closed_circle("covalently closed genomic circles", (role("genomic insert"), role("custom circle-closing junction")), closure_label="T4 DNA ligase closure"),
         "The complementary four-base ends close each fragment; T4 DNA ligase at 5 ng/ul for 16 h at 16 C favours intramolecular ligation."),
        ("Digest every molecule that still has a free end",
         closed_circle("exonuclease-purified DNA circles", (role("closed genomic circle"),), closure_label="no free DNA end"),
         "Exonuclease I, lambda exonuclease and Plasmid-Safe DNase remove residual linear or incompletely closed DNA."),
        ("Cut circles with Cas9–sgRNA",
         molecule("target-containing circles linearized at Cas9 sites", (role("Cas9 cleavage end"), role("off-target or on-target genomic DNA"), role("Cas9 cleavage end"))),
         "Only a circle containing a sequence cleaved by the programmed RNP acquires the free ends needed for library construction."),
        ("Add one 3-prime dA to each Cas9-cleaved end",
         tailed_cut_product(),
         "Proteinase K first strips the RNP; the Kapa A-tailing enzyme then adds the dA before the hairpin adapter is introduced."),
        ("Ligate the NEBNext hairpin adapter",
         molecule("hairpin-adapter-ligated cleavage products", (role("NEBNext hairpin"), role("Cas9-to-Cas9 genomic fragment"), role("NEBNext hairpin"))),
         "Kapa DNA ligase joins the NEBNext hairpin adapter to the dA-tailed ends, preserving each newly exposed Cas9 end as a sequencing boundary."),
        ("Open the hairpins with USER and PCR with dual indexes",
         molecule("dual-index CHANGE-seq library", (role("P5 + i5 + Read 1 arm"), role("genomic DNA starting at a Cas9 cut"), role("Read 2 arm + i7 + P7"))),
         "USER makes the commercial adapter PCR-addressable; twenty Kapa HiFi cycles with NEBNext i5 and i7 primers complete the Illumina library."),
    ))


def oligos():
    i5 = il.nebnext_i5_primer(il.NEBNEXT_I5_SET1["i501"])
    i7 = il.nebnext_i7_primer(il.NEBNEXT_I7_SET1["i701"])
    return (
        "oCRL225 — 5′-/5Phos/ACG/ideoxyU/AGATGTGTATAAGAGACAG-3′",
        f"oCRL226 — 5′-/5Phos/{CHANGE_TRANSPOSON_BOTTOM}-3′",
        "NEBNext Adapter for Illumina — 5′-/5Phos/GATCGGAAGAGCACACGTCTGAACTCCAGTC[dU]ACACTCTTTCCCTACACGACGCTCTTCCGATC*T-3′",
        f"NEBNext i501 primer example — 5′-{i5}-3′",
        f"NEBNext i701 primer example — 5′-{i7}-3′",
    )


TITLE = 'CHANGE-seq — circularization for high-throughput analysis of nuclease genome-wide effects'
NOTES = '01_change_seq.html'
SOURCE = '<a href="https://doi.org/10.1038/s41587-020-0555-7">Lazzarotto et al., Nature Biotechnology (2020)</a>; exact reaction order and custom transposon from the paper’s Supplementary Protocol.'
SUMMARY = 'A dU-bearing custom Tn5 transposome converts genomic DNA into exonuclease-resistant circles; Cas9 selectively reopens circles at on- and off-target sites, which become the ends of a dual-index Illumina library.'
CAVEAT = 'INFERRED — supplied Set1 primer architectures and USER-opened arms are assigned from canonical vendor components of the named kit; study-specific primer members remain unassigned. The two custom Tn5 oligos are paper-specific. The downstream NEBNext adapter and example Set 1 index primers are canonical vendor sequences for the commercial products named by the protocol.'
FINAL_CAPTION = 'One strand of the symmetric Cas9-linearized product is shown. Read 1 and Read 2 both enter genomic DNA from newly exposed nuclease-cleavage ends; the custom Tn5 circle-closing junction is internal to the selected insert rather than a sequencing adapter.'
READ_LENGTHS = {'Read 1': 150, 'Index 1 (i7)': 8, 'Index 2 (i5)': 8, 'Read 2': 150}
SEQUENCING_INTRO = 'The Supplementary Protocol sequences the pooled library on a MiSeq 600-cycle v3 kit with 150-bp paired-end reads and 8-bp dual-index reads; the article also reports NextSeq 550 runs with 150-bp paired-end reads. Standard TruSeq primers bind the NEBNext-derived final construct; the paired genomic reads begin at the two Cas9-cleaved ends and point inward.'
SEQUENCING_ENDING = ''
SEQUENCING_UNAVAILABLE = ''
FINAL_LIBRARY, SEQ_PRIMERS = change_library()

def protocol_graph():
    source = SOURCE
    graph = ProtocolGraph(TITLE, authority=source)
    dna = graph.input("dna", "purified high-molecular-weight genomic DNA",
                      MaterialKind.BIOLOGICAL,
                      state=molecule("purified HMW genomic DNA", (role("genomic DNA"),)),
                      primary=True)
    adapter = graph.input("tn5_adapter", "custom dU-bearing Tn5 transposon",
                          MaterialKind.ADAPTER,state=custom_transposon())
    tn5 = graph.input("tn5", "Tn5 transposase", MaterialKind.CATALYST)
    tagged = graph.product("tagged", "custom-end tagmented genomic fragments",
                           MaterialKind.MOLECULAR, state=molecule("custom-end tagmented fragments",
                             (role("transferred custom end"),role("genomic DNA"),role("opposite custom end"))))
    tn5_waste = graph.product("tn5_waste", "Tn5 off-products and unusable fragments",
                              MaterialKind.MOLECULAR)
    graph.react("tagment", "Tagment with the custom dU-bearing transposome",
                StepKind.TAGMENTATION,
                ((dna, InputRole.TEMPLATE), (adapter, InputRole.STRUCTURAL),
                 (tn5, InputRole.CATALYST)), (tagged, tn5_waste))
    graph.mark_terminal(Terminal(tn5_waste.id,Disposition.UNRESOLVED,
        "Unusable end classes are candidates; absence of a recoverable circle does not establish a physical separation at tagmentation."))
    repair = graph.input("repair_mix", "Kapa HiFi HotStart Uracil+ Ready Mix with Taq DNA ligase",
                         MaterialKind.CATALYST)
    repaired = graph.product("repaired", "fully duplex custom-end fragments",
                             MaterialKind.MOLECULAR, state=molecule("gap-repaired custom-end duplex",(role("copied custom end"),role("duplex genomic insert"),role("opposite copied end"))))
    graph.react("gap_repair", "Fill and seal the staggered transposition products",
                StepKind.END_REPAIR,
                ((tagged, InputRole.TEMPLATE), (repair, InputRole.CATALYST)), (repaired,))
    user = graph.input("user_enzyme", "USER enzyme", MaterialKind.CATALYST)
    opened = graph.product("opened", "USER-opened cohesive-end fragments",
                           MaterialKind.MOLECULAR, state=molecule("USER-opened genomic fragments",(role("ACGT cohesive end"),role("genomic insert"),role("opposite cohesive end"))))
    graph.react("open_du", "Open the internal dU sites", StepKind.DIGESTION,
                ((repaired, InputRole.TEMPLATE), (user, InputRole.CATALYST)), (opened,))
    pnk = graph.input("pnk", "T4 polynucleotide kinase", MaterialKind.CATALYST)
    phosphorylated = graph.product("phosphorylated", "phosphorylated cohesive-end fragments",
                                   MaterialKind.MOLECULAR, state=molecule("PNK-treated USER-opened fragments",(role("5′ phosphorylated cohesive end"),role("genomic insert"),role("opposite phosphorylated end"))))
    graph.react("phosphorylate", "Phosphorylate the opened ligation ends",
                StepKind.PHOSPHORYLATION,
                ((opened, InputRole.TEMPLATE), (pnk, InputRole.CATALYST)), (phosphorylated,))
    circle_ligase = graph.input("circle_ligase", "T4 DNA ligase (circularization at 5 ng/ul)",
                               MaterialKind.CATALYST)
    circles = graph.product("circles", "covalently closed genomic circles",
                            MaterialKind.MOLECULAR, state=closed_circle("covalently closed genomic circle",(role("genomic insert"),role("custom closure junction")),closure_label="T4 ligase closure"))
    graph.react("circularize", "Circularize by dilute intramolecular ligation",
                StepKind.CIRCULARIZATION,
                ((phosphorylated, InputRole.TEMPLATE),
                 (circle_ligase, InputRole.CATALYST)), (circles,),
                geometry=CircularizationGeometry.known(
                    "double-strand intramolecular ligation",
                    statement="CHANGE-seq dilute-ligates the USER-opened ACGT cohesive ends of tagmented genomic fragments into covalently closed circles."))
    exonuclease = graph.input("exonuclease", "Exonuclease I, lambda exonuclease and Plasmid-Safe ATP-dependent DNase",
                              MaterialKind.CATALYST)
    purified = graph.product("purified_circles", "exonuclease-purified genomic circles",
                             MaterialKind.MOLECULAR, state=circles.state)
    graph.react("remove_linear", "Digest molecules that retain free ends",
                StepKind.DIGESTION,
                ((circles, InputRole.TEMPLATE), (exonuclease, InputRole.CATALYST)),
                (purified,))
    sgrna = graph.input("sgrna", "programmed sgRNA", MaterialKind.OLIGO)
    cas9 = graph.input("cas9", "Cas9 nuclease", MaterialKind.CATALYST)
    linear = graph.product("linear", "Cas9-linearized target-containing circles",
                           MaterialKind.MOLECULAR, state=molecule("Cas9-opened genomic circle",(role("Cas9 cut end"),role("genomic DNA including custom closure junction"),role("opposite Cas9 cut end"))))
    graph.react("cas9_cut", "Reopen target-containing circles with Cas9–sgRNA",
                StepKind.GENOME_EDITING,
                ((purified, InputRole.TEMPLATE), (sgrna, InputRole.STRUCTURAL),
                 (cas9, InputRole.CATALYST)), (linear,))
    protk = graph.input("proteinase_k", "proteinase K", MaterialKind.CATALYST)
    deproteinized = graph.product("deproteinized", "RNP-free Cas9-linearized DNA",
                                  MaterialKind.MOLECULAR, state=linear.state)
    graph.react("strip_rnp", "Strip the Cas9 RNP with proteinase K", StepKind.DEPROTEINIZATION,
                ((linear, InputRole.TEMPLATE), (protk, InputRole.CATALYST)), (deproteinized,))
    tailing = graph.input("tailing_enzyme", "Kapa A-tailing enzyme", MaterialKind.CATALYST)
    tailed = graph.product("tailed", "dA-tailed Cas9-linearized products",
                           MaterialKind.MOLECULAR, state=tailed_cut_product())
    graph.react("tail", "Add one 3-prime dA to Cas9-cleaved ends", StepKind.TAILING,
                ((deproteinized, InputRole.TEMPLATE), (tailing, InputRole.CATALYST)), (tailed,))
    hairpin = graph.input("nebnext_adapter", "NEBNext hairpin adapter", MaterialKind.ADAPTER)
    ligase = graph.input("adapter_ligase", "Kapa DNA ligase (HTP PCR-free kit)", MaterialKind.CATALYST)
    ligated = graph.product("ligated", "hairpin-adapter-ligated cleavage products",
                            MaterialKind.MOLECULAR, state=molecule("hairpin-adapted cut fragment",(role("NEBNext hairpin adapter"),role("Cas9-cut genomic fragment"),role("opposite NEBNext hairpin adapter"))))
    graph.react("ligate_adapter", "Ligate the NEBNext hairpin adapter", StepKind.LIGATION,
                ((tailed, InputRole.TEMPLATE), (hairpin, InputRole.STRUCTURAL),
                 (ligase, InputRole.CATALYST)), (ligated,),
                geometry=reported_adapter_end(
                    "NEBNext hairpin adapter", kind="cohesive", polarity="3-prime",
                    terminal_sequence="T",
                    statement="The NEBNext hairpin adapter presents the T end used with dA-tailed Cas9 products."))
    opener = graph.input("hairpin_user", "USER enzyme supplied with NEBNext Dual Index Primers Set 1",
                         MaterialKind.CATALYST)
    addressable = graph.product("addressable", "USER-opened PCR-addressable molecules",
                                MaterialKind.MOLECULAR, state=opened_adapter_template())
    graph.react("open_hairpin", "Open adapter hairpins with USER", StepKind.DIGESTION,
                ((ligated, InputRole.TEMPLATE), (opener, InputRole.CATALYST)),
                (addressable,))
    reviewed = ["dna", "tn5_adapter", "tn5", "tagment", "repair_mix", "gap_repair",
                "user_enzyme", "open_du", "pnk", "phosphorylate", "circle_ligase",
                "circularize", "exonuclease", "remove_linear", "sgrna", "cas9",
                "cas9_cut", "proteinase_k", "strip_rnp", "tailing_enzyme", "tail", "nebnext_adapter",
                "adapter_ligase", "ligate_adapter", "hairpin_user", "open_hairpin"]
    forward_seq=il.nebnext_i5_primer("N"*8)
    reverse_seq=il.nebnext_i7_primer("N"*8)
    forward=graph.input("index_forward","NEBNext Set1 i5 primer architecture",MaterialKind.OLIGO,
        state=Construct([Segment("i5 primer",forward_seq)]),inferred=True,
        source="Supplementary Protocol E7600S kit assignment; canonical vendor primer architecture")
    reverse=graph.input("index_reverse","NEBNext Set1 i7 primer architecture",MaterialKind.OLIGO,
        state=Construct([Segment("i7 primer",reverse_seq)]),inferred=True,
        source="Supplementary Protocol E7600S kit assignment; canonical vendor primer architecture")
    polymerase=graph.input("index_polymerase","2X Kapa HiFi HotStart Ready Mix",MaterialKind.CATALYST)
    library=graph.product("index_library","CHANGE-seq dual-index sequencing library",MaterialKind.LIBRARY,state=FINAL_LIBRARY)
    graph.react("index_pcr","INFERRED — PCR with kit-qualified NEBNext Set1 primer architectures",StepKind.PCR,
        ((addressable,InputRole.TEMPLATE),(forward,InputRole.STRUCTURAL),
         (reverse,InputRole.STRUCTURAL),(polymerase,InputRole.CATALYST)),(library,),
        geometry=PCRGeometry.exact(addressable.state,PCRPrimer("Set1 i5",forward_seq),
            PCRPrimer("Set1 i7",reverse_seq),statement="Kit-qualified primer architecture; actual USER-opened adapter template, not the final indexed library. Exact physical index members remain unassigned."))
    graph.sequence(library,SequencingPlan(library=FINAL_LIBRARY,primers=SEQ_PRIMERS,
        acquired_roles=tuple(p.role for p in SEQ_PRIMERS),read_lengths=READ_LENGTHS))
    reviewed.extend(("index_forward","index_reverse","index_polymerase","index_pcr","index_library"))
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()
