"""Protocol-owned reaction graph for cropseq-multi."""
from types import SimpleNamespace
from chemdraw import Scene, Construct, Segment, feature
import seqprimers as sp
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind
from protocol_model import chain, molecule, role
from pcr import PCRGeometry

# The defining v2 Methods direct-capture section selects GEM-X 5′ v3 and a
# dual-index 28/10/10/90 run. Exact bead barcode/UMI boundaries remain unassigned
# until that kit's molecular declaration is reconstructed; run cycles alone do
# not fix those two component lengths.
CELL = role("GEM-X 5′ v3 droplet barcode", molecular_feature=feature("cell_barcode", "cell_barcode", "unknown"))
UMI = role("GEM-X 5′ v3 UMI", tag="umi", molecular_feature=feature("umi", "umi", "unknown"))
GUIDE = role("internal guide barcode (iBAR)", length_bp=12,
             molecular_feature=feature("guide_ibar", "guide_barcode", "unknown"))
I7 = role("i7 sample index", length_bp=10,
          molecular_feature=feature("sample_i7", "sample_index", "unknown"))
I5 = role("i5 sample index", length_bp=10,
          molecular_feature=feature("sample_i5", "sample_index", "unknown"))
DIRECT_LENGTHS = {"Read 1": 28, "Index 1 (i7)": 10, "Index 2 (i5)": 10, "Read 2": 90}
DIRECT_ROLES = tuple(DIRECT_LENGTHS)

TITLE = "CROPseq-multi — multiplexed perturbation cassette"
DOI = "10.1101/2024.03.17.585235"
CITATION = "Walton et al., bioRxiv (2024)"
SUMMARY = ("A reverse-oriented tRNA–sgRNA array in the lentiviral 3′ LTR expresses multiple "
           "guides while short internal barcodes support mRNA, direct-guide and optical readouts.")
FINAL_CAPTION = ("Direct-capture branch only: the 10x 5′ sgRNA library carries the cell "
                 "barcode, UMI and sample index around the 12-base internal guide barcode "
                 "(iBAR) and its linked spacer. The in situ mRNA and CSMv2 T7-IVT branches "
                 "instead end in rolling-circle amplicons whose iBARs are read by in situ "
                 "sequencing by synthesis; cell identity there is spatial, so those "
                 "branches carry no cell barcode, UMI or sample index.")
NOTES = "01_cropseq_multi.html"
SOURCE = f'Defining source: <a href="https://doi.org/{DOI}">{CITATION}</a>.'
CAVEAT = ("INFERRED — the direct-capture endpoint records the GEM-X 5′ v3 readout roles; "
          "exact bead, RT and dial-out junction geometry is not assembled here. The two "
          "10-base sample indexes apply only to the direct-capture libraries.")
# Scoped to the direct-capture terminal (graph material "direct_library"); the two
# in situ terminals are role-level RCA products with their own SBS limitation.
FINAL_LIBRARY = Construct([
    role("P5 / kit-specific read end", inferred=True), I5, CELL, UMI,
    GUIDE, role("linked guide spacer"), I7, role("P7 / kit-specific read end", inferred=True),
], name="CROPseq-multi direct-capture sgRNA endpoint (role order inferred)")
SEQ_PRIMERS = ()
SEQUENCING_INTRO = FINAL_CAPTION
READ_LENGTHS = DIRECT_LENGTHS
READ_SPECS = (
    sp.read_spec("Read 1", 28, "GEM-X droplet barcode and UMI; exact component boundary not reconstructed", "extract"),
    sp.read_spec("Index 1 (i7)", 10, "sample index", "demultiplex"),
    sp.read_spec("Index 2 (i5)", 10, "sample index", "demultiplex"),
    sp.read_spec("Read 2", 90, "guide spacer and iBAR; scaffold-specific landing not reconstructed", "extract"),
)
SEQUENCING_ENDING = FINAL_CAPTION
SEQUENCING_UNAVAILABLE = ("Target panels, pool barcodes and kit revisions vary; the defining "
                          "source establishes the read roles, while this page does not invent "
                          "a universal primer sequence. For the published 5-prime "
                          "direct-capture demonstration the sgRNA and gene-expression "
                          "libraries were pooled 1:10 and read on a NovaSeq X with 28 "
                          "cycles for read 1, 10 cycles for each index and 90 cycles for "
                          "read 2; the plasmid- and gDNA-screen libraries are single-step "
                          "PCR amplicons read with custom primers on MiSeq 150- or "
                          "300-cycle kits instead.")
# Exact modified primers: defining v2 supplement, sheet oligos rows 29–33.
DIRECT_PRIMERS = {
    "oRW1220": "AAGCAGTGGTATCAACGCAGAGTACcaagttgtgtacggaccagc*c",
    "oRW1221": "AAGCAGTGGTATCAACGCAGAGTACcaagttgataacggactagc*c",
    "oRW1223": "GATCTACACTCTTTCCctacacgacgctcttccgat",
    "oRW1226": "GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTgccttatttgaacttgctaagcagtc",
    "oRW1227": "GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTgccttatttcaacttgctatgctgtt",
}
DIRECT_PRIMER_SOURCE = "Defining v2 supplementary workbook, sheet oligos, rows 29–33"

def population_state(label, *populations):
    """Each supplied tuple is one molecule/population, not a shared covalent chain."""
    sc=Scene()
    for index,parts in enumerate(populations):
        sc.strand(f"population {index+1}",parts,label=f"{label}: population {index+1}")
    return sc


def oligos():
    return tuple(f"<b>{name}</b> — 5′-{modified}-3′" for name, modified in DIRECT_PRIMERS.items())

def _direct_primer(name):
    modified = DIRECT_PRIMERS[name]
    return molecule(name, (Segment(name, modified.replace("*", "").upper(),
                                   note="One terminal phosphorothioate linkage (*)" if "*" in modified else ""),))

_definition = SimpleNamespace(title=TITLE, doi=DOI, citation=CITATION, endpoint=FINAL_CAPTION)

def _cropseq_multi():
    return chain(molecule("CROPseq-multi 3-prime LTR cassette", (
        role("U6 promoter, reverse orientation"), role("tRNA"), role("sgRNA 1"),
        role("tRNA"), role("sgRNA 2"), GUIDE)), (
        ("Duplicate the engineered 3′ LTR during lentiviral integration",
         molecule("integrated perturbation cassette", (role("5-prime LTR copy"),
             role("selection-gene transcript"), role("3-prime LTR copy"))),
         "Keeping the Pol III cassette on the minus strand prevents tRNA processing from cleaving lentiviral genomic RNA or the Pol II selection transcript."),
        ("Transcribe and tRNA-process the multiplexed Pol III guide array",
         population_state('two mature perturbation guides',(role("sgRNA 1 with iBAR"),),(role("sgRNA 2 with iBAR"),)),
         "RNases P and Z release separate guides; each 12-base iBAR lies only 19 bp from its linked spacer."),
        ("Capture perturbation identity by mRNA, direct guide capture or T7-IVT",
         molecule("CROPseq-multi readout", (GUIDE, role("linked guide spacer"))),
         "The same iBAR is transcribed antisense in polyadenylated mRNA and within the Pol III sgRNA; CSMv2 additionally includes a T7 promoter. Only the direct-capture branch adds a cell barcode, UMI and sample index; the in situ branches assign cells by position."),
    ))

workflow = _cropseq_multi


DIRECT_READS = ("Illumina NovaSeq X, 28 cycles read 1, 10 cycles for each index and "
                "90 cycles read 2, with the gene-expression and sgRNA libraries pooled "
                "10:1 (Methods, 'Direct-capture single-cell RNA sequencing with 10X "
                "Genomics 5' capture').")


def _direct_capture_branch(graph, template, final, reviewed, d):
    """5-prime direct-capture branch exactly as the defining Methods section describes it.

    The source follows 10x Genomics CG000735 (GEM-X Single Cell 5' v3, CRISPR Screening)
    but substitutes two custom scaffold-specific RT primers and its own sgRNA dial-out
    PCR, because the kit's CRISPR Poly-dT Primer Mix B is not complementary to both CSM
    sgRNA scaffolds.  The kit steps themselves are not disclosed here, so they stay
    warned boundaries rather than invented chemistry.
    """
    rt1 = graph.input("direct_rt_primer_1",
                      "oRW1220 sgRNA 1 5-prime capture RT primer "
                      "(SMART handle + scaffold 1, 3-prime phosphorothioate)",
                      MaterialKind.OLIGO, state=_direct_primer("oRW1220"), source=DIRECT_PRIMER_SOURCE)
    rt2 = graph.input("direct_rt_primer_2",
                      "oRW1221 sgRNA 2 5-prime capture RT primer "
                      "(SMART handle + scaffold 2, 3-prime phosphorothioate)",
                      MaterialKind.OLIGO, state=_direct_primer("oRW1221"), source=DIRECT_PRIMER_SOURCE)
    rt_mix = graph.input("direct_rt_mix",
                         "10x GEM-X 5-prime v3 RT mastermix (RT Reagent E, Poly-dT RT "
                         "Primer B, Reducing Agent B, RT Enzyme E)",
                         MaterialKind.CATALYST)
    cdna_state = molecule("CROPseq-multi direct-capture cDNA",
                          (CELL, UMI, GUIDE, role("linked guide spacer")))
    cdna = graph.product("direct_cdna",
                         "barcoded 5-prime cDNA (mRNA and sgRNA capture products)",
                         MaterialKind.MOLECULAR, state=cdna_state)
    graph.react("direct_rt",
                "Reverse-transcribe in GEMs with the two custom sgRNA capture primers "
                "added to the CG000735 RT mastermix",
                StepKind.REVERSE_TRANSCRIPTION,
                ((template, InputRole.TEMPLATE), (rt1, InputRole.STRUCTURAL),
                 (rt2, InputRole.STRUCTURAL), (rt_mix, InputRole.CATALYST)), (cdna,))
    amplified = graph.product("direct_amplified_cdna", "amplified barcoded cDNA",
                              MaterialKind.MOLECULAR, state=cdna_state)
    graph.react("direct_cdna_amplification",
                "Amplify the barcoded cDNA as in 10x CG000735",
                StepKind.MAGIC, ((cdna, InputRole.TEMPLATE),), (amplified,),
                reason=("The defining source amplifies cDNA with the commercial kit and "
                        "does not publish its primer sequences or conditions."))
    small = graph.product("direct_sgrna_fraction",
                          "small SPRIselect cDNA fraction (sgRNA capture products)",
                          MaterialKind.MOLECULAR, state=cdna_state)
    large = graph.product("direct_gex_fraction",
                          "large SPRIselect cDNA fraction (gene-expression cDNA)",
                          MaterialKind.MOLECULAR, state=cdna_state)
    graph.react("direct_spri_split",
                "Size-separate the amplified cDNA by SPRIselect into the sgRNA and "
                "gene-expression fractions",
                StepKind.SPLIT, ((amplified, InputRole.TEMPLATE),), (small, large),
                changes_molecule=False)
    gex_library = graph.product("direct_gex_library",
                                "gene-expression library (10x CG000733)",
                                MaterialKind.LIBRARY, state=cdna_state)
    graph.react("direct_gex_library_prep",
                "Build the gene-expression library from the large fraction following "
                "the standard 10x CG000733 protocol",
                StepKind.MAGIC, ((large, InputRole.TEMPLATE),), (gex_library,),
                reason=("The gene-expression library is made with the unmodified "
                        "commercial protocol, which this source does not reproduce."))
    graph.sequence(gex_library, SequencingPlan(
        read_lengths=DIRECT_LENGTHS, acquired_roles=DIRECT_ROLES,
        limitation=("Gene-expression partner library built by the unmodified 10x "
                    "CG000733 protocol; sequenced together with the sgRNA library. "
                    + DIRECT_READS)))
    forward = graph.input("direct_dialout_forward",
                          "oRW1223 universal sgRNA dial-out forward primer",
                          MaterialKind.OLIGO, state=_direct_primer("oRW1223"), source=DIRECT_PRIMER_SOURCE)
    reverse1 = graph.input("direct_dialout_reverse_1",
                           "oRW1226 sgRNA 1 dial-out reverse primer", MaterialKind.OLIGO,
                           state=_direct_primer("oRW1226"), source=DIRECT_PRIMER_SOURCE)
    reverse2 = graph.input("direct_dialout_reverse_2",
                           "oRW1227 sgRNA 2 dial-out reverse primer", MaterialKind.OLIGO,
                           state=_direct_primer("oRW1227"), source=DIRECT_PRIMER_SOURCE)
    polymerase = graph.input("direct_polymerase",
                             "NEBNext Ultra II Q5 master mix", MaterialKind.CATALYST)
    amplicon = graph.product("direct_sgrna_amplicon",
                             "sgRNA dial-out amplicons (both sgRNAs in one reaction)",
                             MaterialKind.MOLECULAR, state=cdna_state)
    graph.react("direct_dialout_pcr",
                "Dial out both sgRNA species from the small fraction in one 11-cycle PCR",
                StepKind.PCR,
                ((small, InputRole.TEMPLATE), (forward, InputRole.STRUCTURAL),
                 (reverse1, InputRole.STRUCTURAL), (reverse2, InputRole.STRUCTURAL),
                 (polymerase, InputRole.CATALYST)), (amplicon,),
                geometry=PCRGeometry.unavailable(
                    cdna_state,
                    reason=("Spacers and iBARs are library members, so the amplicon is a "
                            "family rather than one sequence; the two reverse primers "
                            "address the two sgRNA scaffolds.")))
    index_p5 = graph.input("direct_index_primer_p5",
                           "10x dual-index sample-indexing primer, P5 side (CG000735)",
                           MaterialKind.OLIGO)
    index_p7 = graph.input("direct_index_primer_p7",
                           "10x dual-index sample-indexing primer, P7 side (CG000735)",
                           MaterialKind.OLIGO)
    index_mix = graph.input("direct_index_mix", "sample-index PCR mix (CG000735; reagent details not reconstructed)",
                            MaterialKind.CATALYST)
    library = graph.product("direct_library", "CROPseq-multi sgRNA sequencing library",
                            MaterialKind.LIBRARY, state=final)
    graph.react("direct_index_pcr",
                "Index the sgRNA amplicons for sequencing as described in 10x CG000735",
                StepKind.PCR,
                ((amplicon, InputRole.TEMPLATE), (index_p5, InputRole.STRUCTURAL),
                 (index_p7, InputRole.STRUCTURAL), (index_mix, InputRole.CATALYST)),
                (library,),
                geometry=PCRGeometry.unavailable(
                    amplicon.state,
                    reason=("The sample-indexing PCR is the commercial one; the selected "
                            "index member is not fixed by the defining source.")))
    graph.sequence(library, SequencingPlan(
        read_lengths=DIRECT_LENGTHS, acquired_roles=DIRECT_ROLES,
        limitation=("direct-guide branch: " + d.endpoint + " " + DIRECT_READS)))
    reviewed.extend((
        "direct_rt_primer_1", "direct_rt_primer_2", "direct_rt_mix", "direct_rt",
        "direct_cdna_amplification", "direct_spri_split", "direct_gex_library_prep",
        "direct_gex_library", "direct_dialout_forward", "direct_dialout_reverse_1",
        "direct_dialout_reverse_2", "direct_polymerase", "direct_dialout_pcr",
        "direct_index_primer_p5", "direct_index_primer_p7", "direct_index_mix", "direct_index_pcr",
        "direct_library"))
    return library


IN_SITU_METHODS = ("Methods, 'In situ amplification - mRNA protocol', 'In situ amplification - "
                   "T7 IVT protocol' and 'Sequencing by synthesis'; oligo workbook sheet 'oligos'")
IN_SITU_READOUT = ("In situ sequencing by synthesis with MiSeq Nano v2 incorporation and cleavage "
                   "reagents, imaged on a Nikon Ti2 at 10X; serial decoding strips and re-hybridizes "
                   "SBS primers (oRW981/oRW995/oRW1113 for iBAR 1, oRW938 for iBAR 2). Cell identity "
                   "is spatial, so no cell barcode, UMI or sample index is read.")


def _in_situ_branch(graph, prefix, template, label, reviewed):
    """Optical pooled-screen readout: in situ RT, padlock gap-fill ligation, RCA and SBS.

    The defining Methods read these branches by in situ sequencing, not by PCR and
    short-read sequencing.  Primers and reagents are named as printed; padlock geometry
    and the rolling-circle product are not assembled base by base, so that span is an
    explicit unmodelled step.
    """
    src = graph.authority + " (" + IN_SITU_METHODS + ")"
    rt_primers = graph.input(prefix + "_rt_primers",
                             "biotinylated in situ RT primers: oRW893 (tRNA_G, spacer/iBAR 1) and "
                             "oRW887, oRW965 or oRW967 (tRNA_P, _Q or _A, spacer/iBAR 2), 1 µM each",
                             MaterialKind.OLIGO, source=src)
    rt = graph.input(prefix + "_rt_enzyme",
                     "RevertAid H minus reverse transcriptase with RiboLock, 16 h at 37 °C",
                     MaterialKind.CATALYST, source=src)
    cdna = graph.product(prefix + "_cdna", label + " in situ cDNA", MaterialKind.MOLECULAR,
                         state=molecule(label + " in situ cDNA",
                                        (GUIDE, role("linked guide spacer"))))
    graph.react(prefix + "_rt", "Reverse-transcribe in situ from the biotinylated tRNA-anchored primers",
                StepKind.REVERSE_TRANSCRIPTION,
                ((template, InputRole.TEMPLATE), (rt_primers, InputRole.STRUCTURAL),
                 (rt, InputRole.CATALYST)), (cdna,))
    padlocks = graph.input(prefix + "_padlocks",
                           "5′-phosphorylated padlock probes oRW963 (spacer/iBAR 1) and oRW968, "
                           "oRW964 or oRW966 (spacer/iBAR 2), 0.1 µM each",
                           MaterialKind.OLIGO, source=src)
    enzymes = graph.input(prefix + "_gapfill_rca_enzymes",
                          "RNase H, TaqIT polymerase and Ampligase (gap fill and ligation); "
                          "Phi29 DNA polymerase (RCA)",
                          MaterialKind.CATALYST, source=src)
    rolonies = graph.product(prefix + "_rca_product", label + " rolling-circle amplicons in situ",
                             MaterialKind.LIBRARY,
                             state=molecule(label + " RCA product",
                                            (GUIDE, role("linked guide spacer"),
                                             role("padlock SBS-primer site"))))
    graph.react(prefix + "_padlock_rca",
                "Streptavidin and post-fixation, padlock gap fill and ligation, then RCA (not modelled base by base)",
                StepKind.MAGIC,
                ((cdna, InputRole.TEMPLATE), (padlocks, InputRole.STRUCTURAL),
                 (enzymes, InputRole.CATALYST)), (rolonies,),
                reason=("The Methods give streptavidin binding and 3% formaldehyde/0.1% glutaraldehyde "
                        "post-fixation, gap fill and ligation at 37 °C then 45 °C, and 16 h Phi29 RCA at "
                        "30 °C. The padlock arm landing sites and the gap-filled span are not assembled "
                        "in this role-level graph."))
    graph.sequence(rolonies, SequencingPlan(limitation=label + ": " + IN_SITU_READOUT))
    reviewed.extend((rt_primers.id, rt.id, prefix + "_rt", padlocks.id, enzymes.id,
                     prefix + "_padlock_rca", rolonies.id))
    return rolonies


def _cropseq_multi_graph(d, final):
    source = f"{d.citation}; doi:{d.doi}"
    graph = ProtocolGraph(d.title, authority=source)
    cassette = graph.input("cassette", "CROPseq-multi lentiviral 3-prime-LTR cassette",
                           MaterialKind.MOLECULAR, state=_cropseq_multi().initial, primary=True)
    cells = graph.input("target_cells", "target cells", MaterialKind.BIOLOGICAL)
    integrase = graph.input("integrase", "lentiviral integration machinery",
                            MaterialKind.CATALYST)
    integrated = graph.product("integrated_cassette", "integrated duplicated-LTR perturbation cassette",
                               MaterialKind.MOLECULAR,
                               state=molecule("integrated perturbation cassette", (
                                   role("5-prime LTR copy"), role("selection-gene transcript"),
                                   role("3-prime LTR copy"))))
    graph.react("integrate", "Integrate and duplicate the engineered 3-prime LTR",
                StepKind.INTEGRATION,
                ((cassette, InputRole.TEMPLATE), (cells, InputRole.ANALYTE),
                 (integrase, InputRole.CATALYST)), (integrated,))
    rna_polymerase = graph.input("rna_polymerase", "Pol III transcription machinery",
                                 MaterialKind.CATALYST)
    guides = graph.product("guides", "tRNA-processed guides carrying linked iBARs",
                           MaterialKind.MOLECULAR,
                           state=population_state('processed perturbation guides',(role("sgRNA 1 with iBAR"),),(role("sgRNA 2 with iBAR"),)))
    graph.react("transcribe_guides", "Transcribe and tRNA-process the multiplexed guide array",
                StepKind.TRANSCRIPTION,
                ((integrated, InputRole.TEMPLATE), (rna_polymerase, InputRole.CATALYST)),
                (guides,))
    mrna = graph.product("mrna_branch", "poly(A)-mRNA perturbation-identity branch",
                         MaterialKind.MOLECULAR,
                         state=molecule("polyadenylated Pol II identity transcript", (GUIDE, role("antisense guide-cassette transcript"), role("poly(A)"))))
    direct = graph.product("direct_branch", "direct-guide-capture branch",
                           MaterialKind.MOLECULAR, state=guides.state)
    t7 = graph.product("t7_branch", "CSMv2 T7-IVT branch", MaterialKind.MOLECULAR,
                       state=molecule("CSMv2 DNA template", (role("T7 promoter"), GUIDE, role("linked spacer"))))
    graph.react("split_readouts", "Select the mRNA, direct-guide or CSMv2 T7-IVT readout route",
                StepKind.MAGIC, ((integrated, InputRole.TEMPLATE), (guides, InputRole.TEMPLATE)),
                (mrna, direct, t7),
                reason="These are alternative source-defined assay routes, not a physical split of mature "
                       "sgRNA into mRNA and DNA. Pol II identity-transcript production and preparation "
                       "of the CSMv2 DNA template are not reconstructed base by base.")
    reviewed = ["cassette", "target_cells", "integrase", "integrate", "rna_polymerase",
                "transcribe_guides", "split_readouts"]
    _in_situ_branch(graph, "mrna", mrna, "CROPseq-multi in situ mRNA", reviewed)
    _direct_capture_branch(graph, direct, final, reviewed, d)
    t7_polymerase = graph.input("t7_rna_polymerase",
                                "T7 RNA Polymerase Mix (NEB E2040L), 0.2X, 37 °C overnight",
                                MaterialKind.CATALYST,
                                source=source + " (Methods, 'In situ amplification - T7 IVT protocol')")
    t7_rna = graph.product("t7_rna", "T7-amplified CSMv2 guide RNA", MaterialKind.MOLECULAR,
                           state=molecule("T7 guide-barcode RNA", (GUIDE, role("linked guide spacer"))))
    graph.react("t7_transcription", "Transcribe the CSMv2 guide-barcode region in situ by T7 IVT, then post-fix",
                StepKind.TRANSCRIPTION,
                ((t7, InputRole.TEMPLATE), (t7_polymerase, InputRole.CATALYST)), (t7_rna,))
    reviewed.extend(("t7_rna_polymerase", "t7_transcription"))
    _in_situ_branch(graph, "t7", t7_rna, "CROPseq-multi in situ T7-IVT", reviewed)
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()

def protocol_graph():
    """Assemble the source-scoped material graph here in this protocol directory."""
    return _cropseq_multi_graph(_definition, FINAL_LIBRARY)
