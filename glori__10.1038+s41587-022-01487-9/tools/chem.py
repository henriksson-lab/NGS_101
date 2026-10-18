"""Protocol-owned source, oligos, molecular workflow and reaction graph."""
from adapter_ligation import AdapterLigationGeometry
from chemdraw import Construct, Scene, Segment, Workflow, feature, revcomp
import illumina as il
from pcr import PCRGeometry
from protocol_graph import DetailView, Disposition, InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind, Terminal
from protocol_model import chain, molecule, role
import seqprimers as sp
from rna_end_readout import BaseConversionReadout


TITLE = 'GLORI — absolute m6A quantification by A-to-I conversion'
NOTES = '01_glori.html'
SOURCE = 'Liu et al., Nature Biotechnology (2023), <a href="https://doi.org/10.1038/s41587-022-01487-9">doi:10.1038/s41587-022-01487-9</a>; detailed protocol: Shen et al., Nature Protocols (2024), <a href="https://doi.org/10.1038/s41596-023-00937-1">doi:10.1038/s41596-023-00937-1</a>.'
SUMMARY = 'Glyoxal-assisted nitrite treatment converts unmethylated adenosine to inosine while m6A resists conversion; RNA-sense A retention therefore reports methylation stoichiometry.'
CAVEAT = ("INFERRED — outer library arms and exact sample-index member remain unassigned. "
          "The graph follows only Liu2023's modified-eCLIP branch. Its printed linker, RT primer "
          "and UMI adapter are explicit reagents; selected P5/indexed-P7 primer bases and landing "
          "remain unresolved in this model. The 2024 detailed protocol is a linked follow-up "
          "source, not a replacement recipe. NEBNext Small RNA and SMARTer branches are not drawn.")
# Role order derived 2026-10-17 from the 2024 Supplementary Table 1 oligos: the RNA
# adaptor is the reverse complement of the TruSeq Read 1 primer 3' end (the RT primer
# is that 3' end), and the DNA adaptor arm is the reverse complement of the TruSeq
# Read 2 primer 3' end.  Hence the cDNA strand reads Read1 site - insert - N10 UMI -
# Read2 site.  P5 on the Read 1 side is the standard Illumina convention, inferred:
# the library PCR primers are not printed in the inspected sources.
FINAL_LIBRARY = Construct((role('P5', inferred=True), role('Read 1 site from RNA adaptor / RT primer', tag='r1'), role('cDNA copied from converted RNA (C opposite inosine)'), role('10-nt UMI', tag='umi', length_bp=10, molecular_feature=feature('umi', 'umi', 'random')), role('Read 2 site from 3-prime DNA adaptor', tag='r2'), role('sample index', tag='cbc', molecular_feature=feature('sample_index', 'sample_index', 'unknown')), role('P7', inferred=True)), name='GLORI modified-eCLIP library')
FINAL_CAPTION = 'INFERRED — modified-eCLIP endpoint role order; role order derived from the published adaptors. The 10-nt UMI lies between the cDNA 3-prime end and the Read 2 site; converted A is read as G in the RNA sense while m6A yields A base calls.'
SEQ_PRIMERS = ()
READ_LENGTHS = {"Read 1": 150, "Read 2": 150}
READ_SPECS = (
    sp.read_spec("Read 1", 150, "paired-end partner; not used for GLORI site quantification", "not_used"),
    sp.read_spec("Read 2", 150, "RNA-sense sequence with initial 10-base UMI; deduplicate then trim UMI", "extract_then_include"),
)
SEQUENCING_INTRO = ''
SEQUENCING_ENDING = 'Liu 2023 Methods: the eCLIP-branch libraries were sequenced on Illumina HiSeq X Ten as paired-end 2 × 150 base reads; the 2024 Nature Protocols Step 121 gives the same paired-end 150-bp setting. Only read R2 was analysed, because R2 is sense to the original RNA; the 10-bp UMI at the 5′ end of R2 is used for duplicate removal and then trimmed. This matches the model, where Read 2 is primed from the UMI adapter arm. GEO GSE210563 additionally lists NovaSeq 6000 runs, and the deposited SRA runs (PRJNA866163) are single-end 150-base spots, consistent with only R2 being deposited. The sample-index read length is not stated.'
SEQUENCING_UNAVAILABLE = ("The selected library PCR/sequence-primer members and exact final landing "
                          "remain unassigned. The known 2×150 run and RNA-sense R2/UMI handling are "
                          "declared independently of that unresolved geometry.")

def oligos():
    return ("Liu 2023 Methods, 3′ linker — 5′rAPP-AGATCGGAAGAGCGTCGTG-3SpC3 (2024 Supplementary Table 1 prints rAppAGATCGGAAGAGCGTCGTG; 2024 Reagents state the 5′ adenylation and 3′ C3 chain)",
            'Liu 2023 Methods, reverse transcription primer — ACACGACGCTCTTCCGATCT',
            "Liu 2023 Methods, 5′ adapter (UMI) — 5′Phos-NNNNNNNNNNAGATCGGAAGAGCACACGTCTG-3SpC3 (2024 '3′ DNA adaptor'; Reagents state the 5′ phosphate and 3′ C3 spacer)")

def workflow() -> Workflow:
    readout = BaseConversionReadout("glyoxal plus nitrite", "adenosine", "inosine", "G", "m6A")
    return chain(molecule("fragmented poly(A) RNA", (role("RNA containing A and m6A"),)), (
        ("Protect guanosines with glyoxal",
         molecule("glyoxal-protected RNA", (role("RNA"), role("protected G"), role("A and m6A exposed"))),
         "G protection suppresses unwanted guanosine deamination during the nitrite reaction."),
        ("Deaminate unmethylated A with nitrite",
         molecule("GLORI-converted RNA population", (role(f"{readout.converted_base} converted to {readout.product}"), role(f"{readout.protected_base} remains chemically unconverted"))),
         "Liu 2023 Methods: 750 mM NaNO2, 40 mM MES pH 6.0 and glyoxal for 8 h at 16 °C (an earlier repository note's extra 4 °C overnight incubation is not in the Methods or 2024 Step 23)."),
        ("Heat-deprotect guanosines and repair RNA ends",
         molecule("deprotected converted RNA", (role("A-to-I converted RNA with intact m6A"), role("ligatable termini"))),
         "Liu 2023 Methods: deprotection buffer (triethylammonium acetate/formamide) at 95 °C for 10 min; T4 PNK then removes 3′ phosphates before linker ligation."),
        ("Ligate the modified-eCLIP 3-prime RNA adapter",
         molecule("GLORI adapter-ligated RNA", (role("converted RNA"),
                  Segment("published 3-prime RNA adapter", "AGATCGGAAGAGCGTCGTG", "r1"))),
         "Liu 2023 Methods: 5′-adenylated, 3′-C3-blocked linker ligated by T4 RNA ligase 2 truncated KQ; excess linker removed with 5′ deadenylase and RecJf."),
        ("Reverse-transcribe with the published RT primer",
         molecule("GLORI cDNA", (Segment("published RT primer", "ACACGACGCTCTTCCGATCT", "r1"),
                  role("cDNA copying converted RNA: C opposite inosine"))),
         "Liu 2023 Methods: SuperScript III with this primer; excess primer digested by exonuclease I."),
        ("Ligate the ten-base-UMI DNA adapter to cDNA",
         molecule("UMI-tagged GLORI cDNA", (Segment("published RT primer", "ACACGACGCTCTTCCGATCT", "r1"),
                  role("cDNA copying converted RNA: C opposite inosine"), UMI,
                  Segment("published 3-prime cDNA adapter arm", "AGATCGGAAGAGCACACGTCTG", "r2"))),
         "Liu 2023 Methods: 5′-phosphorylated, 3′-C3-blocked adapter with ten degenerate UMI bases, ligated to the cDNA 3′ end on silane beads with T4 RNA ligase 1 (high concentration)."),
        ("INFERRED — amplify with P5 and indexed P7 primers (bases not printed) and sequence",
         molecule("GLORI library", tuple(FINAL_LIBRARY)),
         f"Unmodified A is read as {readout.read_as} in the RNA sense; m6A remains A, enabling site-level stoichiometry. The 2023 Methods name P5 and Illumina-indexed P7 primers without bases; 2 × 150 paired-end HiSeq X Ten, R2 (RNA-sense, UMI first) analysed."),
    ))



M23 = "Liu2023 published Methods, s41587-022-01487-9.pdf.txt lines2614–2659, 2685–2697"
LINKER = Segment("3′ RNA linker", revcomp(il.TRUSEQ_READ1[-19:]), "r1", note="5′ rApp; 3′ SpC3")
RT_PRIMER = Construct([Segment("RT 5′ extra A", "A", "r1"),
                       Segment("RT annealing region", il.TRUSEQ_READ1[-19:], "r1")], name="published RT primer")
UMI = Segment("10-nt UMI", "N"*10, "umi", placeholder=True, length_bp=10,
              feature=feature("umi", "umi", "random"))
DNA_ADAPTER_ARM = Segment("3′ DNA-adapter arm", il.TRIM_SEEN_IN_READ1[:22], "r2", note="3′ SpC3 block")


def _adapter(name, parts, *, mod5):
    sc = Scene(); sc.strand("adapter", list(parts), label=name, mod5=mod5, mod3="SpC3 (blocked)")
    return sc


def _ligation_end(name, substrate_parts, adapter_parts, *, mod5, statement):
    sc = _adapter(name, adapter_parts, mod5=mod5)
    sc.strand("acceptor", list(substrate_parts), label="3′ OH acceptor", mod3="OH")
    sc.ligation_end("adapter", adapter_parts[0].name, kind="single-strand", polarity="5-prime",
                    owner="adapter", text="5′ " + mod5 + " donor; acceptor 3′ OH")
    return AdapterLigationGeometry.known("single-strand", sc, statement=statement)


def rt_priming_scene():
    sc = Scene()
    sc.strand("RNA", [role("converted RNA fragment"), LINKER], label="linker-ligated RNA", mod3="SpC3")
    sc.anneal("primer", list(RT_PRIMER), to="RNA", pair=("RT annealing region", LINKER.name),
              label="published RT primer (extra 5′ A)")
    sc.arrow("primer", "SuperScript III copies the converted RNA")
    return sc


def _digest_excess(graph, identifier, template, enzyme_label, product_label):
    enzyme = graph.input(identifier+"_enzyme", enzyme_label, MaterialKind.CATALYST, source=M23)
    retained = graph.product(identifier+"_product", product_label, MaterialKind.MOLECULAR, state=template.state)
    hydrolysates = graph.product(identifier+"_hydrolysates", "digested excess reagent products", MaterialKind.MOLECULAR)
    graph.react(identifier, "Digest excess reagent: " + enzyme_label, StepKind.DIGESTION,
                ((template, InputRole.TEMPLATE), (enzyme, InputRole.CATALYST)), (retained, hydrolysates))
    graph.mark_terminal(Terminal(hydrolysates.id, Disposition.UNRESOLVED,
        "The source digests excess reagent before precipitation/bead purification. This model "
        "does not resolve every hydrolysis product's subsequent physical partition."))
    return retained


def _purify(graph, identifier, template, label):
    reagent = graph.input(identifier+"_condition", label, MaterialKind.SELECTION, source=M23)
    retained = graph.product(identifier+"_product", label+" retained material", MaterialKind.MOLECULAR, state=template.state)
    removed = graph.product(identifier+"_excluded", "removed purification fraction", MaterialKind.WASTE)
    graph.react(identifier, label, StepKind.CLEANUP,
                ((template, InputRole.TEMPLATE), (reagent, InputRole.SELECTION)), (retained, removed), changes_molecule=False)
    graph.discard(removed, "physically removed during the named source purification")
    return retained


def protocol_graph() -> ProtocolGraph:
    graph = ProtocolGraph(TITLE, authority=SOURCE)
    rna = graph.input("input", "oligo(dT)-purified mRNA with m6A spike-in probes", MaterialKind.BIOLOGICAL,
                      state=molecule("source RNA", (role("mRNA containing A, G and m6A"),)), primary=True)
    frag = graph.input("fragment_condition", "NEB E6150S fragment buffer, 94 °C 4 min", MaterialKind.CONDITION, source=M23)
    fragments = graph.product("fragment_product", "fragmented mRNA", MaterialKind.MOLECULAR,
                              state=molecule("mRNA fragments", (role("RNA containing A, G and m6A"),)))
    graph.react("fragment", "Fragment the mRNA", StepKind.FRAGMENTATION,
                ((rna, InputRole.TEMPLATE), (frag, InputRole.CONDITION)), (fragments,))
    glyoxal = graph.input("protect_g_condition", "200 ng RNA, 6 µL 8.8 M glyoxal, 20 µL DMSO in 40 µL; "
                         "50 °C 30 min, then 10 µL saturated H₃BO₃, 50 °C 30 min", MaterialKind.CONDITION, source=M23)
    protected = graph.product("protect_g_product", "glyoxal-protected RNA", MaterialKind.MOLECULAR,
                              state=molecule("protected RNA", (role("glyoxal-protected G"), role("A and m6A in RNA"))))
    graph.react("protect_g", "Protect guanosines with glyoxal", StepKind.CHEMICAL_CONVERSION,
                ((fragments, InputRole.TEMPLATE), (glyoxal, InputRole.CONDITION)), (protected,))
    nitrite = graph.input("deaminate_condition", "750 mM NaNO₂, 40 mM MES pH6.0, glyoxal; 8 h 16 °C; ethanol precipitation",
                          MaterialKind.CONDITION, source=M23)
    converted = graph.product("deaminate_product", "A-to-I-converted RNA population", MaterialKind.MOLECULAR,
                              state=molecule("converted RNA", (role("inosine formed from unmethylated A"),
                                            role("m6A remains chemically unconverted"), role("protected G"))))
    graph.react("deaminate", "Deaminate unmethylated adenosines to inosine", StepKind.CHEMICAL_CONVERSION,
                ((protected, InputRole.TEMPLATE), (nitrite, InputRole.CONDITION)), (converted,))
    heat = graph.input("deprotect_condition", "triethylammonium acetate pH8.6/formamide; 95 °C 10 min; "
                       "ethanol precipitation and RNA Clean & Concentrator", MaterialKind.CONDITION, source=M23)
    deprotected = graph.product("deprotect_product", "deprotected converted RNA", MaterialKind.MOLECULAR,
                                state=molecule("converted RNA after deprotection", (role("RNA with inosine, m6A and restored G"),)))
    graph.react("deprotect", "Heat-deprotect guanosines", StepKind.THERMAL_TREATMENT,
                ((converted, InputRole.TEMPLATE), (heat, InputRole.CONDITION)), (deprotected,))
    pnk = graph.input("repair_enzyme", "T4 PNK M0201S, 3′ dephosphorylation", MaterialKind.CATALYST, source=M23)
    repaired = graph.product("repair_product", "3′-OH converted RNA", MaterialKind.MOLECULAR,
                             state=molecule("RNA acceptor", (role("converted RNA fragment"), role("3′ OH"))))
    graph.react("repair", "Remove RNA 3′ phosphates", StepKind.DEPHOSPHORYLATION,
                ((deprotected, InputRole.TEMPLATE), (pnk, InputRole.CATALYST)), (repaired,))
    repaired = _purify(graph, "repair_precipitation", repaired, "Ethanol precipitate end-repaired RNA")
    linker = graph.input("ligate_adapter", "published 5′ rApp, 3′ SpC3 RNA linker", MaterialKind.ADAPTER,
                         state=_adapter("GLORI RNA linker", (LINKER,), mod5="rApp"), source=M23)
    rnl2 = graph.input("ligate_ligase", "T4 RNA ligase2 truncated KQ M0373S", MaterialKind.CATALYST, source=M23)
    ligated_parts = (role("converted RNA fragment"), LINKER)
    tagged_state = molecule("linker-ligated converted RNA", ligated_parts,
                            junctions=(("converted RNA fragment", LINKER.name, "RNA–linker ligation"),))
    tagged = graph.product("ligate_product", "linker-tagged converted RNA", MaterialKind.MOLECULAR, state=tagged_state)
    graph.react("ligate", "Ligate preadenylated linker to RNA 3′ OH", StepKind.LIGATION,
                ((repaired, InputRole.TEMPLATE), (linker, InputRole.STRUCTURAL), (rnl2, InputRole.CATALYST)), (tagged,),
                geometry=_ligation_end("RNA linker", (role("converted RNA fragment"),), (LINKER,), mod5="rApp",
                    statement="Liu2023 Methods: preadenylated 3′-C3-blocked linker, PNK-treated RNA and Rnl2tr KQ. "
                              "Random RNA bases do not change this source-declared single-strand end chemistry."))
    tagged = _digest_excess(graph, "linker_cleanup", tagged, "5′ deadenylase M0331S + RecJf M024S", "RNA after excess-linker digestion")
    tagged = _purify(graph, "linker_precipitation", tagged, "Ethanol precipitate linker-ligated RNA")
    primer = graph.input("rt_primer", "ACACGACGCTCTTCCGATCT RT primer, 1 µL 2 µM", MaterialKind.OLIGO,
                          state=molecule("published RT primer", tuple(RT_PRIMER)), source=M23)
    rt = graph.input("rt_enzyme", "SuperScript III 18080093, 20 µL FS-buffer RT reaction", MaterialKind.CATALYST, source=M23)
    cdna_parts = (*tuple(RT_PRIMER), role("cDNA copied from converted RNA: C opposite inosine"))
    cdna_state = molecule("GLORI single-strand cDNA", cdna_parts)
    cdna = graph.product("rt_product", "cDNA copying converted RNA", MaterialKind.MOLECULAR, state=cdna_state)
    graph.react("rt", "Reverse-transcribe the converted RNA with the printed primer", StepKind.REVERSE_TRANSCRIPTION,
                ((tagged, InputRole.TEMPLATE), (primer, InputRole.STRUCTURAL), (rt, InputRole.CATALYST)), (cdna,))
    cdna = _digest_excess(graph, "primer_cleanup", cdna, "exonuclease I M0293S", "cDNA after excess-primer digestion")
    cdna = _purify(graph, "cdna_purification", cdna, "MyOne Silane bead purification before DNA-adapter ligation")
    dna_adapter = graph.input("ligate_cdna_adapter", "5′-phosphorylated N10 DNA adapter, 3′ SpC3", MaterialKind.ADAPTER,
                             state=_adapter("GLORI cDNA adapter", (UMI, DNA_ADAPTER_ARM), mod5="phosphate"), source=M23)
    rnl1 = graph.input("ligate_cdna_ligase", "T4 RNA ligase1 high concentration M0437M, ligation on silane beads",
                       MaterialKind.CATALYST, source=M23)
    tagged_cdna_parts = (*cdna_parts, UMI, DNA_ADAPTER_ARM)
    tagged_cdna_state = molecule("UMI-adapter-tagged single-strand cDNA", tagged_cdna_parts,
                        junctions=((cdna_parts[-1].name, UMI.name, "cDNA–adapter ligation"),))
    adapted = graph.product("ligate_cdna_product", "UMI-tagged GLORI cDNA", MaterialKind.MOLECULAR, state=tagged_cdna_state)
    graph.react("ligate_cdna", "Ligate the UMI DNA adapter to cDNA 3′ OH", StepKind.LIGATION,
                ((cdna, InputRole.TEMPLATE), (dna_adapter, InputRole.STRUCTURAL), (rnl1, InputRole.CATALYST)), (adapted,),
                geometry=_ligation_end("cDNA adapter", cdna_parts, (UMI, DNA_ADAPTER_ARM), mod5="phosphate",
                    statement="Liu2023 Methods: 5′-phosphorylated N10 adapter with 3′ C3 block ligated to cDNA on beads with Rnl1."))
    adapted = _purify(graph, "adapted_purification", adapted, "MyOne Silane bead purification after DNA-adapter ligation")
    forward = graph.input("library_forward", "P5 primer (selected 2023 study member unassigned)", MaterialKind.OLIGO, source=M23)
    reverse = graph.input("library_reverse", "Illumina-indexed P7 primer (selected 2023 member/index length unassigned)", MaterialKind.OLIGO, source=M23)
    pol = graph.input("library_polymerase", "library PCR polymerase identity unassigned in this declaration", MaterialKind.CATALYST, source=M23)
    amplified = graph.product("library_amplified", "GLORI library PCR product", MaterialKind.MOLECULAR, state=FINAL_LIBRARY)
    graph.react("library_pcr", "Amplify with the named P5/indexed-P7 primers", StepKind.PCR,
                ((adapted, InputRole.TEMPLATE), (forward, InputRole.STRUCTURAL), (reverse, InputRole.STRUCTURAL),
                 (pol, InputRole.CATALYST)), (amplified,), geometry=PCRGeometry.unavailable(adapted.state,
                    reason="Selected library primer members and exact outer amplicon ends are unassigned; "
                           "the actual single-strand UMI-tagged cDNA is the template, not the expected final library."))
    cleaned = _purify(graph, "library_cleanup", amplified, "1.8× XP bead cleanup")
    gel = graph.input("library_gel", "8% TBE gel, recover 165–250-bp products (2023 setting)", MaterialKind.SELECTION, source=M23)
    library = graph.product("library_library", "GLORI modified-eCLIP library", MaterialKind.LIBRARY, state=FINAL_LIBRARY)
    off = graph.product("library_gel_excluded", "PCR products outside 165–250-bp band", MaterialKind.WASTE)
    graph.react("library_size_select", "Gel select the 2023 165–250-bp library", StepKind.SELECTION,
                ((cleaned, InputRole.TEMPLATE), (gel, InputRole.SELECTION)), (library, off), changes_molecule=False)
    graph.discard(off, "excluded when the 165–250-bp gel band is excised")
    graph.sequence(library, SequencingPlan(acquired_roles=("Read 1", "Read 2"), read_lengths=READ_LENGTHS,
                    read_uses={"Read 1": "not_used", "Read 2": "extract_then_include"}, limitation=SEQUENCING_ENDING))
    graph.declare_reviewed_scope(("input", "fragment_condition", "fragment", "protect_g_condition", "protect_g",
        "deaminate_condition", "deaminate", "deprotect_condition", "deprotect", "repair_enzyme", "repair",
        "repair_precipitation_condition", "repair_precipitation", "ligate_adapter", "ligate_ligase", "ligate",
        "linker_cleanup_enzyme", "linker_cleanup", "linker_precipitation_condition", "linker_precipitation", "rt_primer", "rt_enzyme", "rt",
        "primer_cleanup_enzyme", "primer_cleanup", "cdna_purification_condition", "cdna_purification", "ligate_cdna_adapter", "ligate_cdna_ligase", "ligate_cdna",
        "adapted_purification_condition", "adapted_purification", "library_forward", "library_reverse", "library_polymerase", "library_pcr",
        "library_cleanup_condition", "library_cleanup", "library_gel", "library_size_select", "library_library"), M23)
    graph.annotate(DetailView("rt", "Printed RT-primer landing",
        "The 19-base linker pairs with the RT primer's 19-base 3′ region; the primer's "
        "extra 5′ A remains unpaired. This establishes the declared local pairing, not conversion yield.",
        state=rt_priming_scene()))
    return graph.require_publishable()
