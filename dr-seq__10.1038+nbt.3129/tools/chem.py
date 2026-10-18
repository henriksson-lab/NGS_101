"""DR-seq co-amplification and source-backed branch intermediates."""
from __future__ import annotations

from adapter_ligation import AdapterLigationGeometry
from chemdraw import (Construct, MolecularState, Row, Scene, Segment, Workflow,
                      complement_segments, duplex_rows, feature)
from pcr import PCRGeometry
from protocol_graph import (Highlight, InputRole, MaterialKind, ProtocolGraph,
                            SourceClaim, SourceDiscrepancy, StepKind, SequencingPlan, Disposition, Terminal)

import seqprimers as sp
from protocol_model import molecule,role
seg=Segment

TITLE = "DR-seq — genome and transcriptome without physical separation"
NOTES = "01_dr-seq.html"
SOURCE = ('Defining source: <a href="https://doi.org/10.1038/nbt.3129">'
          'Dey et al., <i>Nature Biotechnology</i> (2015)</a>.')
SUMMARY = ("A cell-barcoded oligo(dT) first makes cDNA in the same lysate as genomic DNA. "
           "Ad-2 random primers then quasilinearly copy both analytes before the pool is "
           "split into genome and transcriptome finishing reactions.")
CAVEAT = ("The checked DR-seq paper prints P1, P2 and P3, but not the complete Ad-1x "
          "or Ad-2 oligos. Its cited CEL-Seq source has not been reviewed, so their full "
          "working sequences are inferred. RNA library finishing follows CEL-Seq and DNA "
          "library finishing names historical NEBNext Ultra; their native outer adapter/index geometry remains unassigned, with no canonical shell substituted. "
          "The article gives Ad-2 eight random 3′ nucleotides, while the supplement says "
          "its 3′ end is GGG or TTT; see the source discrepancy below.")

AD1_SPACER = "CGATTGAGGCCGG"
T7_PROMOTER = "TAATACGACTCACTATAGGG"
RA5 = "GTTCAGAGTTCTACAGTCCGACGATC"
CELL_BARCODE_NT = 8
POLYT_NT = 24
AD2_HANDLE = "GTGAGTGATGGTTGAGGTAGTGTGGAG"
AD2_RANDOM_NT = 8          # article line 77: 27-nt defined sequence + 8 random nt
AD2_TERMINAL = "NNN"       # supplement Fig. 3 legend / Note: 3' end GGG or TTT (placeholder)
AD2_VARIABLE_NT = AD2_RANDOM_NT - len(AD2_TERMINAL)
P1 = "CGATTGAGGCCGGTAATAC"
P2 = AD2_HANDLE
P3 = "GTGAGCTGGAGTTGAGGTAGTGTGGAG"


def ad1_segments():
    return [seg("Ad-1x spacer", AD1_SPACER, inferred=True),
            seg("T7 promoter", T7_PROMOTER, "t7", inferred=True),
            seg("RA5", RA5, "r1", inferred=True),
            seg("cell barcode", "B" * CELL_BARCODE_NT, "cbc", placeholder=True,
                feature=feature("cell_barcode", "cell_barcode", "fixed"), inferred=True),
            seg("oligo(dT)24", "T" * POLYT_NT, inferred=True),
            seg("V anchor", "V", placeholder=True, inferred=True)]


def ad2_segments():
    return [seg("Ad-2 common handle", AD2_HANDLE, "tso", inferred=True),
            seg("Ad-2 random N5", "N" * AD2_VARIABLE_NT, placeholder=True,
                inferred=True),
            seg("Ad-2 3' GGG or TTT end", AD2_TERMINAL, placeholder=True,
                inferred=True)]


def p3_mismatch_positions() -> list[int]:
    """1-based positions where printed P3 differs from printed P2 (Ad-2 handle)."""
    return [i + 1 for i, (a, b) in enumerate(zip(P2, P3)) if a != b]


def rt_scene() -> Scene:
    mrna = [seg("transcript", "X" * 34, placeholder=True),
            seg("poly(A)", "A" * POLYT_NT)]
    sc = Scene(); sc.strand("mRNA", mrna, label="poly(A) RNA")
    sc.anneal("Ad-1x", ad1_segments(), to="mRNA",
              pair=("oligo(dT)24", "poly(A)"), label="cell-barcoded Ad-1x",
              unpaired=("V anchor",))
    sc.arrow("Ad-1x", "reverse transcription")
    sc.mark("Ad-1x", "cell barcode", "8-nt cell identity")
    sc.labels("Ad-1x")
    return sc


def ad2_priming_scene(template_name: str) -> Scene:
    template = [seg(template_name, "X" * 40, placeholder=True),
                seg("CCC or AAA priming site", "NNN", placeholder=True),
                seg("random priming site", "N" * AD2_VARIABLE_NT, placeholder=True)]
    sc = Scene(); sc.strand("template", template, label=template_name)
    sc.anneal("Ad-2", ad2_segments(), to="template",
              pair=("Ad-2 random N5", "random priming site"), label="Ad-2 random primer")
    sc.arrow("Ad-2", "strand-displacing extension")
    sc.labels("Ad-2")
    return sc


def gdna_amplicon() -> Construct:
    left = [*ad2_segments(), seg("genomic copy", "X" * 40, placeholder=True)]
    return Construct([*left, *complement_segments(ad2_segments(), suffix=" right complement")],
                     name="Ad-2/Ad-2 genomic full amplicon")


def cdna_amplicon() -> Construct:
    # The productive transcript-derived species is Ad-2 at its random-primed end and
    # retains Ad-1x at the poly(A)-primed end.
    return Construct([*ad2_segments(), seg("cDNA", "X" * 38, placeholder=True),
                      seg("copied poly(A)", "A" * POLYT_NT),
                      seg("cell barcode complement", "b" * CELL_BARCODE_NT, "cbc",
                          placeholder=True,
                          feature=feature("cell_barcode", "cell_barcode", "fixed")),
                      *complement_segments(ad1_segments()[:3], suffix=" complement")],
                     name="Ad-2/Ad-1x transcript amplicon")


GDNA_AMPLICON = gdna_amplicon()
CDNA_AMPLICON = cdna_amplicon()


def p3_amplicon() -> Construct:
    # Methods lines 402-404 list one primer, 5'-biotin P3, so both ends of the
    # re-amplified genomic products carry P3; the internal genomic copy is symbolic.
    p3 = [seg("5'-biotin P3", P3, "tso")]
    return Construct([*p3, seg("genomic copy", "X" * 40, placeholder=True),
                      *complement_segments(p3, suffix=" right complement")],
                     name="biotin-P3/P3 genomic amplicon")


P3_AMPLICON = p3_amplicon()


def arna_product() -> Construct:
    return Construct([
        seg("T7 initial G-rich end", "GGG", "t7",inferred=True), seg("RA5", RA5, "r1",inferred=True),
        seg("cell barcode", "B" * CELL_BARCODE_NT, "cbc", placeholder=True,
            feature=feature("cell_barcode", "cell_barcode", "fixed"),inferred=True),
        seg("poly(U) and antisense transcript", "UUU" + "X" * 38, placeholder=True),
        *complement_segments(ad2_segments(), suffix=" complement"),
    ], name="cell-barcoded antisense RNA")


ARNA = arna_product()


def arna_rows():
    sc = Scene(); sc.strand("aRNA", list(ARNA), label="cell-barcoded aRNA")
    sc.labels("aRNA")
    return sc.rows()

# The exact post-sonication DNA insert and exact post-IVT RNA insert are retained, while
# outer Illumina completion is marked inferred because those oligos are not printed here.
DNA=Construct([role("NEBNext Ultra genome library; native outer adapter/index architecture unassigned",inferred=True,
                    molecular_feature=feature("genome_library_cell_index","sample_index","unknown",note="Source uses different index primers for each single cell"))],name="DR-seq genome endpoint (bounded NEBNext route)")
RNA=Construct([role("CEL-Seq outer left end; unassigned",inferred=True),
    role("Ad-1x cell barcode",length_bp=CELL_BARCODE_NT,inferred=True,molecular_feature=feature("cell_barcode","cell_barcode","fixed")),
    role("RNA-derived insert; variable extent"),role("CEL-Seq outer right end; unassigned",inferred=True)],name="DR-seq transcriptome endpoint (bounded CEL-Seq route)")
DP=RP=()
READ_SPECS={"Genome library":(sp.read_spec("Read 1","50 or100","genomic fragment; exact native primer/index geometry unassigned","include"),sp.read_spec("Read 2","50 or100","opposite genomic fragment end","include")),
 "Transcriptome library":(sp.read_spec("Read 1",100,"CEL-Seq-derived cell identity/transcript readout; exact allocation unassigned","undetermined"),sp.read_spec("Read 2",100,"CEL-Seq-derived paired readout; exact allocation unassigned","undetermined"))}
SEQUENCING_UNAVAILABLE="Defining paper specifies paired HiSeq2500 runs; versioned NEBNext/CEL-Seq adapter, index and sequencing-primer architecture remain unassigned."

FINAL_LIBRARIES = (
    ("Genome library", DNA, DP,
     "NEBNext Ultra library from sonicated, Ad-2-depleted genomic amplicons; outer shell is inferred because its lot-specific oligos are not printed.",
     "Physical NEBNext adapter/index geometry remains unassigned; no canonical TruSeq shell is substituted."),
    ("Transcriptome library", RNA, RP,
     "Directional CEL-Seq-derived library retaining the 8-nt Ad-1x cell barcode; outer small-RNA completion is shown as inferred.",
     "Native CEL-Seq primer/index geometry remains unassigned; paired100-cycle readout is source-declared."),
)


def oligos():
    return (
        f"<b>INFERRED — Ad-1x working model</b> — 5′-<inf>{AD1_SPACER}{T7_PROMOTER}{RA5}[8-nt cell barcode]{'T' * POLYT_NT}V</inf>-3′",
        f"<b>INFERRED — Ad-2 working model</b> — 5′-<inf>{AD2_HANDLE}[5 random nt]GGG or TTT</inf>-3′",
        f"<b>P1</b> — 5′-{P1}-3′",
        f"<b>P2</b> — 5′-{P2}-3′",
        f"<b>5′-biotin P3</b> — 5′-biotin-{P3}-3′",
    )


def workflow() -> Workflow:
    start = Scene.duplex([seg("genomic DNA", "X" * 42, placeholder=True)],
                         label="single-cell genomic DNA")
    wf = Workflow(MolecularState("One lysed cell: genomic DNA plus poly(A) RNA",
                                 tuple(start.rows())))
    wf.react("INFERRED — Reverse-transcribe RNA with the working Ad-1x oligo",
             rt_scene().rows(), name="Ad-1x-tagged cDNA beside genomic DNA",
             note="The 5′ tail installs the T7 promoter, RA5 and cell barcode; genomic DNA remains in the tube.")
    wf.react("INFERRED — Random-prime genomic DNA and cDNA with working Ad-2 oligo",
             [*ad2_priming_scene("genomic DNA").rows(),
              *ad2_priming_scene("Ad-1x-tagged cDNA").rows()],
             name="First Ad-2 semi-amplicons")
    wf.react("Run seven MALBAC-like quasilinear cycles",
             [*duplex_rows(GDNA_AMPLICON, label="genome-derived Ad-2/Ad-2 species"),
              *duplex_rows(CDNA_AMPLICON, label="transcript-derived Ad-2/Ad-1x species"),
              Row(chunks=[("full Ad-2/Ad-2 amplicons self-anneal through complementary ends and leave the quasilinear pool", None, False)])],
             name="Co-amplified genome and transcript products")
    wf.react("Split the common amplification pool into DNA and RNA aliquots",
             [Row(chunks=[("DNA aliquot → P2 PCR → 5′-biotin P3 PCR → sonicate → streptavidin depletion", None, False)]),
              Row(chunks=[("RNA aliquot → P1 second strand → T7 IVT", None, False)])],
             name="Modality-specific finishing reactions")
    wf.react("Select the RNA branch with P1 and amplify it by T7 transcription",
             [*duplex_rows(CDNA_AMPLICON, label="P1-completed T7 template"),
              Row(chunks=[("T7 RNA polymerase →", "t7", False)]),
              *arna_rows()],
             name="Antisense RNA copies from transcript-derived templates",
             note="Genomic Ad-2/Ad-2 products lack the Ad-1x T7 end and do not make this aRNA.")
    wf.react("Finish the DNA and RNA aliquots as separate Illumina libraries",
             [*duplex_rows(DNA, label="genome library"),
              *duplex_rows(RNA, label="transcriptome library")],
             name="DR-seq genome and transcriptome libraries",
             note="The source names NEBNext Ultra for DNA and the CEL-Seq workflow for RNA; the unprinted outer oligos are marked inferred.")
    return wf


def common_population(*,p3=False):
    scene=Scene()
    scene.strand("genome population",tuple(P3_AMPLICON if p3 else GDNA_AMPLICON),label="genome-derived representative amplicons")
    scene.strand("transcript population",tuple(CDNA_AMPLICON),label="transcript-derived material in the same tube; recovery/amplification not resolved")
    return scene


def cell_population(*,copied=False):
    sc=Scene.duplex([role("genomic DNA population")],unpaired=("genomic DNA population'",))
    sc.strand("RNA or cDNA",[*ad1_segments(),role("first-strand cDNA")] if copied else [role("poly(A) mRNA")],label="independent first-strand cDNA" if copied else "independent mRNA")
    return sc


def protocol_graph():
    """Typed provenance for DR-seq's common amplification and two finishes."""
    g = ProtocolGraph(TITLE, authority=SOURCE)
    cell = g.input(
        "cell_lysate", "one lysed cell containing genomic DNA and poly(A) RNA",
        MaterialKind.BIOLOGICAL,
        state=cell_population(),primary=True)

    ad1 = g.input("ad1", "cell-barcoded Ad-1x oligo(dT)24V primer",
                  MaterialKind.OLIGO, state=molecule("INFERRED Ad-1x",tuple(ad1_segments())))
    rt = g.input("rt", "ArrayScript reverse transcriptase (MessageAmp II), 42 °C then 80 °C 10 min",
                 MaterialKind.CATALYST)
    copied = g.product(
        "copied", "Ad-1x-tagged cDNA beside genomic DNA",
        MaterialKind.MOLECULAR,
        state=cell_population(copied=True))
    g.react(
        "reverse_transcribe", "INFERRED Ad-1x model — reverse-transcribe RNA",
        StepKind.REVERSE_TRANSCRIPTION,
        ((cell, InputRole.TEMPLATE), (ad1, InputRole.STRUCTURAL),
         (rt, InputRole.CATALYST)), (copied,))

    ad2 = g.input("ad2", "Ad-2 common-handle random-N8 primer",
                  MaterialKind.OLIGO)
    strand_polymerase = g.input(
        "strand_polymerase",
        "Bst large fragment + Pyrophage 3173 exo- (fresh mix before each cycle)",
        MaterialKind.CATALYST)
    primed = g.product(
        "ad2_primed", "Ad-2-primed genomic DNA and cDNA copies",
        MaterialKind.MOLECULAR,
        state=[*ad2_priming_scene("genomic DNA").rows(),
               *ad2_priming_scene("Ad-1x-tagged cDNA").rows()])
    g.react(
        "random_prime", "INFERRED Ad-2 model — random-prime genomic DNA and cDNA",
        StepKind.EXTENSION,
        ((copied, InputRole.TEMPLATE), (ad2, InputRole.STRUCTURAL),
         (strand_polymerase, InputRole.CATALYST)), (primed,))
    coamplified = g.product(
        "coamplified", "quasilinearly co-amplified genome and transcript products",
        MaterialKind.MOLECULAR,
        state=common_population())
    g.react(
        "quasilinear_cycles", "INFERRED representative products — seven quasilinear cycles",
        StepKind.EXTENSION,
        ((primed, InputRole.TEMPLATE), (ad2, InputRole.STRUCTURAL),
         (strand_polymerase, InputRole.CATALYST)), (coamplified,))

    dna = g.product("dna_branch", "genome finishing aliquot",
                    MaterialKind.MOLECULAR, state=coamplified.state)
    rna = g.product("rna_branch", "transcriptome finishing aliquot",
                    MaterialKind.MOLECULAR, state=coamplified.state)
    g.react("split", "Split the common amplification pool into DNA and RNA aliquots",
            StepKind.SPLIT, ((coamplified, InputRole.TEMPLATE),), (dna, rna),
            changes_molecule=False)

    p2_forward = g.input("p2_forward", "P2 genome-amplification primer",
                         MaterialKind.OLIGO, state=Construct([seg("P2", P2)], name="P2"))
    p2_polymerase = g.input("p2_polymerase", "Deep VentR (exo-) DNA polymerase",
                            MaterialKind.CATALYST)
    p2_product = g.product(
        "p2_product", "P2-amplified genomic products", MaterialKind.MOLECULAR,
        state=common_population())
    g.react(
        "p2_pcr", "Amplify the genome aliquot with P2 (21 cycles), column clean-up", StepKind.PCR,
        ((dna, InputRole.TEMPLATE), (p2_forward, InputRole.STRUCTURAL),
         (p2_polymerase, InputRole.CATALYST)), (p2_product,),
        geometry=PCRGeometry.unavailable(
            dna.state,symmetric_species=True,
            reason="The actual aliquot retains both analytes. The source gives one P2 sequence for the repeated Ad-2 handles; this representative construct does not assert a unique oriented primer pair."))

    p3 = g.input("p3_biotin", "5-prime-biotin P3 primer", MaterialKind.OLIGO,
                 state=Construct([seg("5'-biotin P3", P3)], name="5'-biotin P3"))
    p3_product = g.product(
        "p3_product", "biotin-P3/P3 genomic amplicons",
        MaterialKind.MOLECULAR, state=common_population(p3=True))
    g.react(
        "p3_pcr", "Re-amplify 80 ng with 5-prime-biotin P3 (4 cycles 46 °C, 9 cycles 59 °C)",
        StepKind.PCR,
        ((p2_product, InputRole.TEMPLATE), (p3, InputRole.STRUCTURAL),
         (p2_polymerase, InputRole.CATALYST)), (p3_product,),
        geometry=PCRGeometry.unavailable(
            p2_product.state,symmetric_species=True,
            reason="The Methods list only biotin-P3, which differs from the Ad-2 handle at positions "
                   + ", ".join(map(str, p3_mismatch_positions()))
                   + "; the priming geometry of the mismatched first cycles is not modelled."))
    shear = g.input("shearing", "sonication condition", MaterialKind.CONDITION)
    fragmented = g.product("fragmented_dna", "sonicated genome amplicons",
                           MaterialKind.MOLECULAR, state=molecule("sonicated population",(role("internal fragments and biotin-P3 end fragments; lengths/end positions unassigned"),)))
    g.react("sonicate", "Sonicate (Bioruptor, low, 15 cycles) to ~300 bp and column-purify",
            StepKind.FRAGMENTATION,
            ((p3_product, InputRole.TEMPLATE), (shear, InputRole.CONDITION)),
            (fragmented,))
    streptavidin = g.input("streptavidin", "streptavidin depletion matrix",
                           MaterialKind.SELECTION)
    depleted = g.product("depleted_dna", "Ad-2-depleted sonicated DNA",
                         MaterialKind.MOLECULAR, state=molecule("retained supernatant+wash DNA",(role("sonicated Ad-2-depleted DNA; transcript carryover unquantified"),)))
    bound = g.product("bound_ad2", "streptavidin-bound Ad-2 end fragments",
                      MaterialKind.WASTE)
    g.react(
        "deplete_ad2", "Deplete biotin-labelled Ad-2 ends with streptavidin",
        StepKind.SELECTION,
        ((fragmented, InputRole.TEMPLATE),
         (streptavidin, InputRole.SELECTION)), (depleted, bound),
        changes_molecule=False)
    g.discard(bound, "contains the biotin-labelled Ad-2 ends removed before library prep")
    neb_kit=g.input("neb_kit","NEBNext Ultra DNA Library Prep Kit with distinct index primers per cell",MaterialKind.CONDITION)
    dna_library=g.product("dna_library","indexed genome sequencing library",MaterialKind.LIBRARY,state=DNA)
    g.react("dna_library_prep","Prepare indexed DNA library using NEBNext Ultra",StepKind.MAGIC,
        ((depleted,InputRole.TEMPLATE),(neb_kit,InputRole.CONDITION)),(dna_library,),
        reason="The original names the historical kit; its lot-specific end preparation/adapter/index PCR sequence is not reconstructed, and no later kit revision or canonical duplex is substituted.")
    g.sequence(dna_library,SequencingPlan(library=DNA,acquired_roles=("Read 1","Read 2"),limitation="HiSeq2500 paired50 or100 DNA runs; no particular alternative or physical primer/index architecture is assigned."))

    p1 = g.input("p1", "P1 transcript-branch primer", MaterialKind.OLIGO,
                 state=Construct([seg("P1", P1)], name="P1"))
    second_polymerase = g.input("second_polymerase", "second-strand polymerase",
                                MaterialKind.CATALYST)
    t7_template = g.product(
        "t7_template", "P1-completed cell-barcoded T7 template",
        MaterialKind.MOLECULAR, state=common_population())
    g.react(
        "second_strand", "P1 second strand (one PCR cycle); pool non-overlapping barcodes; column-purify",
        StepKind.EXTENSION,
        ((rna, InputRole.TEMPLATE), (p1, InputRole.STRUCTURAL),
         (second_polymerase, InputRole.CATALYST)), (t7_template,))
    t7_polymerase = g.input("t7_polymerase", "T7 RNA polymerase",
                            MaterialKind.CATALYST)
    arna = g.product("arna", "cell-barcoded antisense RNA copies",
                     MaterialKind.MOLECULAR, state=molecule("single-stranded aRNA",tuple(ARNA)))
    ivt_residual=g.product("ivt_residual","DNA templates and nontranscribed carryover outside the recovered aRNA endpoint",MaterialKind.MOLECULAR)
    g.react("ivt", "T7 IVT 37 °C 13 h (MessageAmp II), aRNA column clean-up without fragmentation",
            StepKind.TRANSCRIPTION,
            ((t7_template, InputRole.TEMPLATE),
             (t7_polymerase, InputRole.CATALYST)), (arna,ivt_residual))
    g.mark_terminal(Terminal(ivt_residual.id,Disposition.UNRESOLVED,"Lack of genomic-DNA transcription does not establish physical destruction or complete column exclusion. Component recovery through source aRNA purification remains unassigned."))
    cel_adapter = g.input(
        "cel_adapter", "CEL-Seq small-RNA finishing adapter",
        MaterialKind.ADAPTER, inferred=True,
        source=SOURCE + " (workflow referenced; outer oligos not printed)")
    rna_ligated = g.product(
        "rna_ligated", "3-prime-adapter-ligated antisense RNA",
        MaterialKind.MOLECULAR, state=molecule("adapter-ligated aRNA",(role("aRNA with Ad-1x cell identity"),role("3′adapter; species/geometry unassigned"))), inferred=True,
        source=SOURCE + " (CEL-Seq finishing delegated)")
    g.react(
        "adapt_rna", "Attach the referenced CEL-Seq small-RNA finishing adapter",
        StepKind.LIGATION,
        ((arna, InputRole.TEMPLATE), (cel_adapter, InputRole.STRUCTURAL)),
        (rna_ligated,), geometry=AdapterLigationGeometry.unavailable(
            "The defining paper refers to the CEL-Seq finishing workflow without printing its outer adapter oligos."))
    rna_rt_primer = g.input(
        "rna_rt_primer", "CEL-Seq library RT primer", MaterialKind.OLIGO, inferred=True,
        source=SOURCE + " (main text line 89-90 names reverse transcription; primer not printed)")
    library_rt_enzyme=g.input("rna_library_rt_enzyme","historical CEL-Seq library reverse transcriptase; exact component unassigned",MaterialKind.CATALYST,
        source=SOURCE+" cites historical CEL-Seq finishing; no component-specific enzyme inheritance established")
    rna_adapted = g.product(
        "rna_adapted", "reverse-transcribed CEL-Seq library template",
        MaterialKind.MOLECULAR, state=molecule("single-strand CEL-Seq library cDNA",(role("reverse-transcribed adapter-ligated aRNA; native primer/end species unassigned"),)), inferred=True,
        source=SOURCE + " (CEL-Seq finishing delegated)")
    g.react(
        "rna_library_rt", "Reverse-transcribe the adapter-ligated aRNA",
        StepKind.REVERSE_TRANSCRIPTION,
        ((rna_ligated, InputRole.TEMPLATE), (rna_rt_primer, InputRole.STRUCTURAL),
         (library_rt_enzyme, InputRole.CATALYST)), (rna_adapted,))

    reviewed = [
        "cell_lysate", "ad1", "rt", "reverse_transcribe", "ad2",
        "strand_polymerase", "random_prime", "quasilinear_cycles", "split",
        "p2_forward", "p2_polymerase", "p2_pcr", "p3_biotin",
        "p3_pcr", "shearing", "sonicate", "streptavidin",
        "deplete_ad2", "neb_kit", "dna_library_prep", "dna_library", "p1", "second_polymerase",
        "second_strand", "t7_polymerase", "ivt", "cel_adapter", "adapt_rna",
        "rna_ligated", "rna_rt_primer", "rna_library_rt_enzyme", "rna_library_rt",
    ]
    rna_forward=g.input("rna_forward","CEL-Seq finishing PCR primer; physical species unassigned",MaterialKind.OLIGO)
    rna_reverse=g.input("rna_reverse","CEL-Seq opposite PCR primer; physical species unassigned",MaterialKind.OLIGO)
    rna_pol=g.input("rna_polymerase","CEL-Seq library PCR reagent; versioned recipe not reconstructed",MaterialKind.CATALYST)
    rna_library=g.product("rna_library","cell-barcoded transcriptome sequencing library",MaterialKind.LIBRARY,state=RNA)
    g.react("rna_library_pcr","PCR-complete the CEL-Seq-derived library",StepKind.PCR,
        ((rna_adapted,InputRole.TEMPLATE),(rna_forward,InputRole.STRUCTURAL),(rna_reverse,InputRole.STRUCTURAL),(rna_pol,InputRole.CATALYST)),(rna_library,),
        geometry=PCRGeometry.unavailable(rna_adapted.state,reason="Actual input is single-stranded library cDNA; the cited historical CEL-Seq primer/reaction recipe remains a review target, not a finished-duplex PCR template."))
    g.sequence(rna_library,SequencingPlan(library=RNA,acquired_roles=("Read 1","Read 2"),read_lengths={"Read 1":100,"Read 2":100},limitation="Published paired100 RNA run retained; actual CEL-Seq primer/index geometry and exact read allocations remain unassigned."))
    reviewed.extend(("rna_forward","rna_reverse","rna_polymerase","rna_library_pcr","rna_library"))
    g.annotate(SourceDiscrepancy(
        "ad2", "Ad-2 3′ end: eight random nucleotides or GGG/TTT?",
        (SourceClaim("Article main text p.2 (PMC PDF text line 77)",
                     "Ad-2 has a defined 27-nt 5′ sequence followed by 8 random nucleotides"),
         SourceClaim("Supplementary Fig. 3 legend and Supplementary Note (publisher supplement text lines 268, 1404)",
                     "the 3′ end of Ad-2 has the sequence GGG or TTT (citing Zong et al. MALBAC)")),
        "Draw the 8-nt variable region as 5 random nt plus a named terminal GGG-or-TTT placeholder triplet; total length 8 nt.",
        "Both statements are printed in the published article and its supplement; the supplement also reports GGG/TTT as the dominant observed 3-mer binding sites (~70% of reads).",
        "The complete Ad-2 primer mix is not printed in DR-seq; whether the 8 variable positions include the fixed triplet is taken from the supplement, not from an order sheet."))
    g.annotate(Highlight(
        "p3_pcr", "Biotin-P3 is not identical to P2",
        "Printed P3 differs from printed P2 (the Ad-2 handle) at positions "
        + ", ".join(map(str, p3_mismatch_positions()))
        + "; the first four cycles anneal at 46 °C. The Methods do not explain the difference, so the page keeps both printed sequences."))
    g.annotate(Highlight("reverse_transcribe","RNA-template disposition remains unassigned",
        "The state tracks the co-amplification substrates, genomic DNA and first-strand cDNA. RNA-template recovery after initial RT is not reconstructed, and no physical removal/destruction is declared."))
    g.annotate(Highlight("rna_library_rt","Separate historical library-RT component and RNA fate",
        "ArrayScript belongs to the original cellular RT. Its use is not inferred for the later CEL-Seq library RT; that catalyst is unassigned. The library RT produces ssDNA, while the RNA template's subsequent disposition remains unassigned, not a presumed physical discard."))
    g.declare_reviewed_scope(reviewed, SOURCE)
    return g.require_publishable()
