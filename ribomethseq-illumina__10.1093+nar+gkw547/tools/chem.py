"""Protocol-owned source, oligos, molecular workflow and reaction graph."""
from adapter_ligation import AdapterLigationGeometry
from chemdraw import Construct, Workflow, feature
from pcr import PCRGeometry
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind
from protocol_model import chain, molecule, role



TITLE = 'RiboMethSeq (Illumina) — 2′-O-methylation by protected cleavage'
NOTES = '01_ribomethseq.html'
SOURCE = 'Marchand et al., Nucleic Acids Research (2016), <a href="https://doi.org/10.1093/nar/gkw547">doi:10.1093/nar/gkw547</a>; original Ion Torrent method (not drawn): Birkedal et al., Angewandte Chemie (2015), <a href="https://doi.org/10.1002/anie.201408362">doi:10.1002/anie.201408362</a>.'
SUMMARY = 'Partial alkaline hydrolysis avoids phosphodiester bonds immediately 3′ of 2′-O-methylated riboses; depletion of mapped fragment ends identifies and quantifies Nm sites.'
CAVEAT = ('This page documents the 2016 Marchand Illumina implementation (NEBNext Small RNA library), which is its defining source. The original 2015 Birkedal RiboMethSeq protocol is a distinct earlier Ion Torrent implementation and is not drawn; nothing below is transferred onto the drawn construct. '
          'The Birkedal Angewandte article itself has still not been obtained; its Wiley Supporting Information (anie_201408362_sm_miscellaneous_information.pdf, obtained 2026-10-17) was read. '
          'According to that SI (Materials and Methods, pp. S3–S5; Figure S4), RNA at 1 µg/µl was mixed with an equal volume of 100 mM Na₂CO₃ pH 9.0 and heated at 90 °C for 1–10 min; 20–40-nt fragments were excised from a 10% denaturing gel and eluted near neutral pH to keep their 2′,3′-cyclic phosphates. '
          'Fragments (5′-OH, 2′,3′-cyclic phosphate) were cloned with a kinase-dead Arabidopsis tRNA ligase (AtRNL S701A, “AtRNLΔPNK”): first a T4 PNK-phosphorylated 3′ RNA adaptor corresponding to the Ion Torrent A sequence with a 4-nt barcode (CR30 5′-CGCACUGAGUCGGAG or CR35 5′-CGAGCUGAGUCGGAG; 30 °C 1 h), then one T4 PNK reaction to phosphorylate the ligation product and remove cyclic phosphates of unligated fragments, gel purification, ligation of an in vitro transcribed, HDV-ribozyme-cleaved 5′ adaptor (P1 IVT, 5′-triphosphate, 43 nt, 2′,3′-cyclic phosphate end), gel purification, FastAP removal of the junction 2′-phosphates, and reverse transcription with a short Ion A primer carrying barcode #1 (CGCA); libraries went to emulsion PCR and Ion Torrent/Proton sequencing (250 or 500 flows). '
          'Marchand Table S4 summarises this as “no end repair” with direct ligation to 5′-OH/>P fragments; the SI shows that no repair precedes the first ligation, but T4 PNK and FastAP treatments occur between and after the ligations. '
          'Source inconsistency within the Birkedal SI: the P1 IVT adaptor is printed as 5′-GGCCACUACGCCUCCGCUUUCCUCUCUAUGGGCAGUCGG<b>T</b>GAU>P (p. S4, adaptor construction) and as …GGCAGUCGG<b>U</b>GAU>P (p. S4, library construction); as an in vitro transcript of the printed T7 template the RNA residue is U, so the T is read as a typographical DNA letter. The RT primer bases are not printed. '
          'GEO series <a href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE54655">GSE54655</a> (PMID 25417815) gives the same outline.')
FINAL_LIBRARY = Construct((role('P5'), role('hydrolysis-end cDNA'), role('sample index', tag='cbc', molecular_feature=feature('sample_index', 'sample_index', 'unknown')), role('P7')), name='RiboMethSeq library')
FINAL_CAPTION = 'Standard small-RNA library retaining the exact alkaline-cleavage endpoints.'
SEQ_PRIMERS = ()
READ_LENGTHS = None
SEQUENCING_INTRO = ''
SEQUENCING_ENDING = 'The Illumina implementation used 2×75 paired-end MiSeq or 50-cycle single-read HiSeq 1000 runs (Marchand Supplementary Table S2); standard kit primers enter from the ligated adapters.'
SEQUENCING_UNAVAILABLE = 'The selected E7330S adapter/primer members, kit revision and exact primer landing remain unassigned in this model. Separate published SR50 and PE2×75 run declarations are retained.'

def oligos():
    # The previous generic adapter string had no version-qualified local source.
    # Keep the named E7330S reagents unassigned until the matching vendor original is read.
    return ()

def workflow() -> Workflow:
    return chain(molecule("total or enriched RNA", (role("RNA containing 2-prime-O-methylated riboses"),)), (
        ("Partially hydrolyse RNA under mild alkaline conditions",
         molecule("alkaline RNA fragments", (role("random fragments with gaps beside protected Nm sites"),)),
         "A 2′-O-methyl ribose protects its 3′-adjacent phosphodiester bond from cleavage. Marchand used 50 mM bicarbonate pH 9.2 at 95 °C and no gel size selection (30–200-nt fragments)."),
        ("Convert cyclic-phosphate/5′-OH ends to ligatable termini",
         molecule("end-repaired RNA fragments", (role("5-prime phosphate"), role("RNA fragment"), role("3-prime hydroxyl"))),
         "Antarctic phosphatase followed by T4 PNK creates the ends required by the standard small-RNA kit; the repaired RNA is purified on an RNeasy MinElute column."),
        ("Ligate 3′ and 5′ small-RNA adapters",
         molecule("dual-adapter RNA fragments", (role("NEBNext 5-prime small-RNA adapter"), role("hydrolysis fragment"), role("NEBNext 3-prime small-RNA adapter"))),
         "Direct ligation preserves the exact hydrolysis endpoints used for methylation scoring."),
        ("Reverse-transcribe, index by PCR and sequence fragment ends",
         molecule("RiboMethSeq library", tuple(FINAL_LIBRARY)),
         "Paired-end MiSeq can recover both ends; SR50 HiSeq analysis can use the less biased 5′-end profile."),
    ))



M16 = "Marchand2016 PMC5027498 Methods RNA fragmentation/End repair/Library preparation/Sequencing, lines62–80"
RUN_SOURCE = "Marchand2016 Supplementary Table S2, File008.doc.txt Table S2 run-type key and rows1–4"


def _column_cleanup(graph, current):
    column = graph.input("minelute_condition", "RNeasy MinElute; 675 µL 96% ethanol binding; 10 µL water elution",
                         MaterialKind.SELECTION, source=M16)
    kept = graph.product("minelute_product", "column-purified 5′-P/3′-OH RNA", MaterialKind.MOLECULAR, state=current.state)
    waste = graph.product("minelute_waste", "column flow-through and washes", MaterialKind.WASTE)
    graph.react("minelute", "Purify repaired RNA on RNeasy MinElute", StepKind.CLEANUP,
                ((current, InputRole.TEMPLATE), (column, InputRole.SELECTION)), (kept, waste), changes_molecule=False)
    graph.discard(waste, "physically removed in the column flow-through and washes")
    return kept


def protocol_graph() -> ProtocolGraph:
    graph = ProtocolGraph(TITLE, authority=SOURCE)
    rna = graph.input("input", "1–250 ng total or purified RNA (source study sample)", MaterialKind.BIOLOGICAL,
                      state=molecule("source RNA", (role("RNA with 2′-O-methylated sites"),)), primary=True)
    hydrolysis = graph.input("hydrolyse_condition", "50 mM bicarbonate pH9.2, 95 °C 4–10 min (drawn alkaline branch)",
                             MaterialKind.CONDITION, source=M16)
    fragments = graph.product("hydrolyse_product", "alkaline RNA fragments, 30–200 nt distribution", MaterialKind.MOLECULAR,
                              state=molecule("hydrolysis fragments", (role("5′ OH"), role("RNA fragment"), role("phosphate/cyclic-phosphate 3′ end"))))
    graph.react("hydrolyse", "Partially hydrolyse RNA under alkaline conditions", StepKind.CHEMICAL_CONVERSION,
                ((rna, InputRole.TEMPLATE), (hydrolysis, InputRole.CONDITION)), (fragments,))
    precipitation = graph.input("hydrolysis_precip_condition", "ethanol/NaOAc pH5.2/glycoblue precipitation; 80% ethanol wash",
                                 MaterialKind.SELECTION, source=M16)
    recovered = graph.product("hydrolysis_precip_product", "recovered fragmented RNA", MaterialKind.MOLECULAR, state=fragments.state)
    supernatant = graph.product("hydrolysis_precip_waste", "discarded precipitation liquid and wash", MaterialKind.WASTE)
    graph.react("hydrolysis_precip", "Recover RNA fragments by ethanol precipitation", StepKind.CLEANUP,
                ((fragments, InputRole.TEMPLATE), (precipitation, InputRole.SELECTION)), (recovered, supernatant), changes_molecule=False)
    graph.discard(supernatant, "physically removed during precipitation/washing; no gel size fractionation is applied")
    phosphatase = graph.input("dephosphorylate_enzyme", "5 U Antarctic phosphatase; 37 °C 30 min; inactivate 70 °C 5 min",
                              MaterialKind.CATALYST, source=M16)
    dephosphorylated = graph.product("dephosphorylate_product", "3′-dephosphorylated RNA", MaterialKind.MOLECULAR,
                                     state=molecule("RNA after phosphatase", (role("5′ OH"), role("RNA fragment"), role("3′ OH"))))
    graph.react("dephosphorylate", "Remove RNA 3′ phosphates without gel purification", StepKind.DEPHOSPHORYLATION,
                ((recovered, InputRole.TEMPLATE), (phosphatase, InputRole.CATALYST)), (dephosphorylated,))
    kinase = graph.input("phosphorylate_enzyme", "T4 PNK, 1 mM ATP; 37 °C 1 h", MaterialKind.CATALYST, source=M16)
    repaired = graph.product("phosphorylate_product", "5′-phosphorylated/3′-OH RNA", MaterialKind.MOLECULAR,
                             state=molecule("ligation-ready RNA fragment", (role("5′ phosphate"), role("RNA fragment"), role("3′ OH"))))
    graph.react("phosphorylate", "Phosphorylate RNA 5′ ends", StepKind.PHOSPHORYLATION,
                ((dephosphorylated, InputRole.TEMPLATE), (kinase, InputRole.CATALYST)), (repaired,))
    purified = _column_cleanup(graph, repaired)
    # The defining original identifies E7330S and direct RNA-end ligation; kit-specific
    # adapter members, modified ends and enzymes are not supplied from a generic recipe.
    adapter3 = graph.input("ligate_3_adapter", "NEBNext E7330S 3′ RNA-adapter reagent (member/end state unassigned)",
                           MaterialKind.ADAPTER, source=M16)
    ligase3 = graph.input("ligate_3_ligase", "E7330S 3′-adapter ligation enzyme component (unassigned)", MaterialKind.CATALYST, source=M16)
    rna3 = graph.product("ligate_3_product", "3′-adapter-ligated hydrolysis RNA", MaterialKind.MOLECULAR,
                         state=molecule("3′-tagged RNA", (role("5′ phosphate"), role("hydrolysis fragment"), role("E7330S 3′ adapter")),
                               junctions=(("hydrolysis fragment", "E7330S 3′ adapter", "direct 3′ RNA-end ligation"),)))
    graph.react("ligate_3", "Ligate the kit 3′ adapter directly to RNA ends", StepKind.LIGATION,
                ((purified, InputRole.TEMPLATE), (adapter3, InputRole.STRUCTURAL), (ligase3, InputRole.CATALYST)), (rna3,),
                geometry=AdapterLigationGeometry.unavailable("Marchand names the E7330S route and direct RNA-end ligation; "
                    "this declaration has not assigned the matching vendor adapter/end chemistry or kit revision."))
    adapter5 = graph.input("ligate_5_adapter", "NEBNext E7330S 5′ RNA-adapter reagent (member/end state unassigned)",
                           MaterialKind.ADAPTER, source=M16)
    ligase5 = graph.input("ligate_5_ligase", "E7330S 5′-adapter ligation enzyme component (unassigned)", MaterialKind.CATALYST, source=M16)
    rna5 = graph.product("ligate_5_product", "dual-adapter RNA preserving hydrolysis endpoints", MaterialKind.MOLECULAR,
                         state=molecule("dual-tagged RNA", (role("E7330S 5′ adapter"), role("hydrolysis fragment"), role("E7330S 3′ adapter")),
                               junctions=(("E7330S 5′ adapter", "hydrolysis fragment", "direct 5′ RNA-end ligation"),
                                          ("hydrolysis fragment", "E7330S 3′ adapter", "direct 3′ RNA-end ligation"))))
    graph.react("ligate_5", "Ligate the kit 5′ adapter directly to the RNA 5′ end", StepKind.LIGATION,
                ((rna3, InputRole.TEMPLATE), (adapter5, InputRole.STRUCTURAL), (ligase5, InputRole.CATALYST)), (rna5,),
                geometry=AdapterLigationGeometry.unavailable("The E7330S vendor oligo/end-state original remains a "
                    "review target; no Ion Torrent/Birkedal or generic adapter is substituted."))
    rt_primer = graph.input("rt_primer", "E7330S reverse-transcription primer (member unassigned)", MaterialKind.OLIGO, source=M16)
    rt = graph.input("rt_enzyme", "E7330S reverse-transcription enzyme (unassigned)", MaterialKind.CATALYST, source=M16)
    cdna = graph.product("rt_product", "hydrolysis-fragment cDNA", MaterialKind.MOLECULAR,
                         state=molecule("single-strand cDNA", (role("kit RT-primer-derived end"), role("antisense cDNA preserving RNA endpoints"))))
    graph.react("rt", "Reverse-transcribe the dual-adapter RNA per the selected kit", StepKind.REVERSE_TRANSCRIPTION,
                ((rna5, InputRole.TEMPLATE), (rt_primer, InputRole.STRUCTURAL), (rt, InputRole.CATALYST)), (cdna,))
    forward = graph.input("library_forward", "E7330S forward library primer (unassigned)", MaterialKind.OLIGO, source=M16)
    reverse = graph.input("library_reverse", "E7330S indexing primer (member/index length unassigned)", MaterialKind.OLIGO, source=M16)
    pol = graph.input("library_polymerase", "E7330S PCR enzyme component (unassigned)", MaterialKind.CATALYST, source=M16)
    library = graph.product("library_pre_run", "RiboMethSeq Illumina library family", MaterialKind.MOLECULAR, state=FINAL_LIBRARY)
    graph.react("library_pcr", "Amplify the endpoint-preserving library using the selected kit", StepKind.PCR,
                ((cdna, InputRole.TEMPLATE), (forward, InputRole.STRUCTURAL), (reverse, InputRole.STRUCTURAL),
                 (pol, InputRole.CATALYST)), (library,), geometry=PCRGeometry.unavailable(cdna.state,
                    reason="The selected E7330S vendor primer members and exact adapter/cDNA geometry remain "
                           "unassigned. RNA end preservation is the source readout contract, not an exact PCR proof."))
    run_choice = graph.input("run_assignment", "Published run setting / sample-library assignment (Table S2)", MaterialKind.CONDITION,
                             source=RUN_SOURCE)
    pe = graph.product("library_pe75", "MiSeq PE2×75 library/run variant", MaterialKind.LIBRARY, state=FINAL_LIBRARY)
    sr = graph.product("library_sr50", "HiSeq1000 SR50 library/run variant", MaterialKind.LIBRARY, state=FINAL_LIBRARY)
    graph.react("run_variants", "Represent the study's alternative PE75 and SR50 library/run assignments", StepKind.SPLIT,
                ((library, InputRole.TEMPLATE), (run_choice, InputRole.CONDITION)), (pe, sr), changes_molecule=False,
                reason="Table S2 assigns different study libraries to different run settings; this branch declaration "
                       "does not claim that every physical sample was sequenced in both modes.")
    graph.sequence(pe, SequencingPlan(acquired_roles=("Read 1", "Read 2"), read_lengths={"Read 1":75,"Read 2":75},
                    limitation="Marchand Table S2: MiSeq PE2×75. Paired ends support RNA-end counting; "
                               "exact primer landing and index cycles are unassigned."))
    graph.sequence(sr, SequencingPlan(acquired_roles=("Read 1",), read_lengths={"Read 1":50},
                    limitation="Marchand Methods and Table S2: HiSeq1000 SR50; source 5′-end profile analysis. "
                               "Exact primer landing and index cycles are unassigned."))
    graph.declare_reviewed_scope(("input", "hydrolyse_condition", "hydrolyse", "hydrolysis_precip_condition", "hydrolysis_precip",
        "dephosphorylate_enzyme", "dephosphorylate", "phosphorylate_enzyme", "phosphorylate", "minelute_condition", "minelute",
        "ligate_3_adapter", "ligate_3_ligase", "ligate_3", "ligate_5_adapter", "ligate_5_ligase", "ligate_5", "rt_primer", "rt_enzyme", "rt",
        "library_forward", "library_reverse", "library_polymerase", "library_pcr", "run_assignment", "run_variants", "library_pe75", "library_sr50"), M16+"; "+RUN_SOURCE)
    return graph.require_publishable()
