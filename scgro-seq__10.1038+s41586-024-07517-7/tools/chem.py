"""scGRO-seq clickable nascent-RNA chemistry."""
from chemdraw import Row
from multibranch import illumina_branch
from protocol_extensions import clicked_nascent_rna
from pcr import PCRGeometry
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph,
                            SequencingPlan, SourceClaim, SourceDiscrepancy, StepKind)
TITLE="scGRO-seq — single-cell nascent RNA by click chemistry"
NOTES="01_scgro-seq.html"
SOURCE='Defining source: <a href="https://doi.org/10.1038/s41586-024-07517-7">Mahat et al., <i>Nature</i> (2024)</a>.'
SUMMARY="Nuclear run-on in permeabilized cells or nuclei incorporates chain-terminating 3′-O-propargyl nucleotides; single cells are sorted into 96 urea-filled wells and CuAAC clicks a 5′-azide barcode/UMI oligo onto the alkyne. ZnCl2 fragmentation, template-switching RT through the triazole, denaturing-PAGE size selection, six-cycle pre-amplification, Cas9 (DASH) depletion of empty-adaptor amplicons and a final indexed PCR give the library."
FINAL_LIBRARY,SEQ_PRIMERS=illumina_branch("scGRO-seq","nascent-RNA cDNA",cell=10,umi=5,
                                          inferred_adapters=True)
FINAL_CAPTION="INFERRED — the 10-nt cell barcode and 5-nt UMI represent the third 43-nt azide adaptor in Supplementary Table 7; the complete final adapter junction is not strand-resolved here."
SEQUENCING_INTRO="The paper reports a 64/43/8/8 NovaSeq SP100 run and trims Read 1 motif CCCCTGTCTCTTATACACAT and Read 2 motif AGATCGGAAGAGCGTCGTGT before extracting cell barcodes and UMIs from Read 1. The displayed primer binding remains an inferred final-library model."
CAVEAT=("Supplementary Table 7 gives three alternative 5′-azide adaptor designs with "
        "different barcode/UMI dimensions (10/6, 12/4 and 10/5). The final-library "
        "drawing uses the third design as a representative; its Illumina shell remains inferred.")
READ_LENGTHS={"Read 1":64,"Read 2":43,"Index 1 (i7)":8,"Index 2 (i5)":8}

def oligos():
    return [
        "<b>5′-AzScBc design 1, 43 nt</b> — 5′-[azideT]JJJJJJJJJJNNNNNNGATCGTCGGACTGTAGAACTCTGAAC-3′ (J = cell barcode; N = UMI)",
        "<b>5′-AzScBc design 2, 86 nt hairpin</b> — 5′-[azideT]NNNNJJJJJJJJJJJJGATCGTCGGACTGTAGAACTCTGCGGCCGTGCTCGTTTTCGAGCACGGCCGCAGAGTTCTACAGTCCGA-3′",
        "<b>5′-AzScBc design 3, 43 nt</b> — 5′-[azideT]ACAGGNNNNNJJJJJJJJJJAGATCGGAAGAGCGTCGTGTAG[SpC3]-3′",
    ]
def sections():
    return [("Run on engaged polymerases with 3-prime O-propargyl NTPs",[Row(chunks=[("nascent RNA → RNA—O-propargyl-NMP (chain terminated)", "w1", False)])],"The terminal alkyne marks the active RNA 3-prime end."),("Click the azide barcode oligo to the alkyne",clicked_nascent_rna().rows(),"Copper-catalysed azide–alkyne cycloaddition creates a triazole linkage that RT can traverse."),("Reverse-transcribe and add sequencing adapters",[Row(chunks=[("clicked RNA → barcoded cDNA → PCR/index adapters", "r1", False)])],"The clicked oligo supplies molecule and cell identity before amplification."),]

def protocol_graph():
    source="Defining source: https://doi.org/10.1038/s41586-024-07517-7"
    graph=ProtocolGraph("scGRO-seq",authority=source)
    cells=graph.input("cultured_cells","V6.5 mouse ES cells",MaterialKind.BIOLOGICAL,primary=True)
    perm_buffer=graph.input("permeabilization_buffer",
                            "douncing or permeabilization buffer (Triton X-100, or Tween-20/NP-40)",
                            MaterialKind.CONDITION)
    complexes=graph.product("engaged_complexes",
                            "permeabilized cells or isolated nuclei with engaged transcription complexes",
                            MaterialKind.BIOLOGICAL)
    graph.react("permeabilize","Isolate nuclei or permeabilize cells and store in storage buffer",
                StepKind.PERMEABILIZATION,
                ((cells,InputRole.TEMPLATE),(perm_buffer,InputRole.CONDITION)),(complexes,))
    nucleotide=graph.input("propargyl_nucleotide",
                           "50 µM each 3′-(O-propargyl)-ATP/CTP/GTP/UTP (Jena NU-945 to NU-948)",
                           MaterialKind.MOLECULAR)
    run_on_buffer=graph.input("run_on_buffer",
                              "nuclear run-on buffer (400 mM KCl, 0.05% Sarkosyl, 5 min at 37 °C)",
                              MaterialKind.CONDITION)
    polymerase=graph.input("engaged_polymerase","engaged RNA polymerase",MaterialKind.CATALYST)
    labelled=graph.product("labelled_rna","alkyne-terminated nascent RNA in nuclei",
                           MaterialKind.MOLECULAR)
    graph.react("run_on","Run on engaged polymerases with 3′-O-propargyl NTPs (chain terminating)",
                StepKind.RUN_ON,
                ((complexes,InputRole.TEMPLATE),(nucleotide,InputRole.STRUCTURAL),
                 (run_on_buffer,InputRole.CONDITION),(polymerase,InputRole.CATALYST)),(labelled,))
    sorted_cell=graph.product("sorted_cell","one run-on cell or nucleus per barcoded well, lysed in 8 M urea",
                              MaterialKind.MOLECULAR)
    unsorted=graph.product("unsorted_suspension","cells and nuclei left in the sorted suspension",
                           MaterialKind.WASTE)
    graph.react("sort_cells","FACS-sort single cells or nuclei into 96 wells containing 2.5 µl 8 M urea",
                StepKind.PARTITION,((labelled,InputRole.TEMPLATE),),(sorted_cell,unsorted))
    graph.discard(unsorted,"material not deposited into a barcoded well is not carried forward")
    click_oligo=graph.input("click_oligo","5′-AzScBc DNA, one cell barcode per well (design 3)",
                            MaterialKind.ADAPTER,
                            source=source+" (Supplementary Table 7, design 3)")
    copper=graph.input("cuaac_catalyst",
                       "CuAAC catalyst: 1 mM CuSO4, 800 mM sodium ascorbate, 5 mM BTTAA ligand",
                       MaterialKind.CATALYST)
    cuaac_condition=graph.input("cuaac_condition",
                                "15% PEG 8000, 2.66 M urea, 2 h at 50 °C",MaterialKind.CONDITION)
    clicked=graph.product("clicked_rna","barcoded clicked nascent RNA",MaterialKind.MOLECULAR,
                          state=clicked_nascent_rna(barcode_length=10,umi_length=5))
    unclicked=graph.product("unclicked_azide","unreacted 30 nM 5′-AzScBc DNA",MaterialKind.MOLECULAR)
    graph.react("click_attachment","CuAAC: click the 5′-azide barcode oligo onto the 3′ alkyne",
                StepKind.ADAPTER_ATTACHMENT,
                ((sorted_cell,InputRole.TEMPLATE),(click_oligo,InputRole.STRUCTURAL),
                 (cuaac_condition,InputRole.CONDITION),(copper,InputRole.CATALYST)),
                (clicked,unclicked))
    cleanup_reagents=graph.input("cuaac_cleanup",
                                 "5 mM EDTA quench, TRIzol (PEG 8000) and 3 kDa Amicon filtration",
                                 MaterialKind.SELECTION)
    pooled=graph.product("pooled_clicked_rna","pooled, reagent-free clicked RNA from 96 wells",
                         MaterialKind.MOLECULAR,state=clicked.state)
    carried_azide=graph.product("carried_unclicked_azide","unclicked 5′-AzScBc DNA carried into RT",
                                MaterialKind.MOLECULAR)
    spent=graph.product("spent_cuaac_reagents","removed PEG 8000, sodium ascorbate, CuSO4 and BTTAA",
                        MaterialKind.WASTE)
    graph.react("quench_and_clean","Quench with EDTA, pool 96 wells and remove the CuAAC reagents",
                StepKind.CLEANUP,
                ((clicked,InputRole.TEMPLATE),(unclicked,InputRole.TEMPLATE),
                 (cleanup_reagents,InputRole.SELECTION)),
                (pooled,carried_azide,spent),changes_molecule=False)
    graph.discard(spent,"small-molecule reagents pass the 3 kDa membrane or partition into TRIzol")
    zinc=graph.input("zinc_chloride","10 mM ZnCl2, 5 min at 65 °C",MaterialKind.CONDITION)
    fragments=graph.product("fragmented_clicked_rna","fragmented clicked RNA–DNA conjugates",
                            MaterialKind.MOLECULAR)
    graph.react("fragment_clicked_rna","Fragment the clicked RNA with ZnCl2",StepKind.FRAGMENTATION,
                ((pooled,InputRole.TEMPLATE),(zinc,InputRole.CONDITION)),(fragments,))
    rt=graph.input("reverse_transcriptase",
                   "RNase H− M-MuLV reverse transcriptase (Thermo EP0751), 16% PEG 8000, 2 h at 50 °C",
                   MaterialKind.CATALYST)
    rt_primer=graph.input("rt_primer","1 µM RT primer (omitted for the hairpin-forming 86-nt design)",
                          MaterialKind.OLIGO,inferred=True,
                          source=source+" (ordering sequence not printed)")
    tso=graph.input("template_switch_oligo","1 µM template-switching oligo",MaterialKind.OLIGO,
                    inferred=True,source=source+" (ordering sequence not printed)")
    cdna=graph.product("cdna","barcoded nascent-RNA cDNA",MaterialKind.MOLECULAR)
    empty_cdna=graph.product("empty_cdna",
                             "empty cDNA from unclicked 5′-AzScBc DNA, non-templated CCC and TSO",
                             MaterialKind.MOLECULAR)
    graph.react("reverse_transcription","Reverse-transcribe through the triazole link with template switching",
                StepKind.REVERSE_TRANSCRIPTION,
                ((fragments,InputRole.TEMPLATE),(carried_azide,InputRole.TEMPLATE),
                 (rt_primer,InputRole.STRUCTURAL),(tso,InputRole.STRUCTURAL),
                 (rt,InputRole.CATALYST)),(cdna,empty_cdna))
    gel=graph.input("denaturing_page","10% denaturing PAGE size selection",MaterialKind.SELECTION)
    selected=graph.product("selected_cdna","gel-purified barcoded cDNA",MaterialKind.MOLECULAR)
    residual_empty=graph.product("residual_empty_cdna",
                                 "empty cDNA co-migrating with the selected size range",
                                 MaterialKind.MOLECULAR)
    excised=graph.product("excised_short_species",
                          "unclicked 5′-AzScBc DNA and empty cDNA left in the gel",MaterialKind.WASTE)
    graph.react("size_select","Size-select the cDNA away from unclicked adaptor and empty cDNA",
                StepKind.SELECTION,
                ((cdna,InputRole.TEMPLATE),(empty_cdna,InputRole.TEMPLATE),
                 (gel,InputRole.SELECTION)),(selected,residual_empty,excised),
                changes_molecule=False)
    graph.discard(excised,"species outside the excised gel band are physically removed")
    q5=graph.input("q5_master_mix","NEBNext Ultra II Q5 High-Fidelity 2× master mix (NEB M0544)",
                   MaterialKind.CATALYST)
    fwd_primer=graph.input("index_pcr_forward","0.5 µM P5-side indexed PCR primer (unique dual index)",
                           MaterialKind.OLIGO,inferred=True,
                           source=source+" (ordering sequence not printed)")
    rev_primer=graph.input("index_pcr_reverse","0.5 µM P7-side indexed PCR primer (unique dual index)",
                           MaterialKind.OLIGO,inferred=True,
                           source=source+" (ordering sequence not printed)")
    pre_library=graph.product("pre_library","pre-amplified barcoded cDNA",MaterialKind.LIBRARY)
    empty_amplicon=graph.product("empty_amplicon","pre-amplified empty-adaptor amplicon (ACAGG target)",
                                 MaterialKind.LIBRARY)
    graph.react("pre_amplify","Pre-amplify to dsDNA: 6 cycles, 98/60/72 °C (Supplementary Table 8)",
                StepKind.PCR,
                ((selected,InputRole.TEMPLATE),(residual_empty,InputRole.TEMPLATE),
                 (fwd_primer,InputRole.STRUCTURAL),(rev_primer,InputRole.STRUCTURAL),
                 (q5,InputRole.CATALYST)),(pre_library,empty_amplicon),
                geometry=PCRGeometry.unavailable(
                    FINAL_LIBRARY,reason="The pre-amplified scGRO-seq amplicon and empty-adaptor junctions are not strand-resolved."))
    cas9=graph.input("cas9","6.6 µM S. pyogenes Cas9 (NEB M0386T), 1 h at 37 °C",MaterialKind.CATALYST)
    guide=graph.input("dash_guide","20 µM gRNA against the ACAGG-containing empty-adaptor amplicon",
                      MaterialKind.OLIGO,inferred=True,source=source+" (guide ordering sequence unresolved)")
    dash_library=graph.product("dash_library","barcoded cDNA surviving DASH",MaterialKind.LIBRARY)
    cut_empty=graph.product("cut_empty","Cas9-cleaved empty-adaptor amplicons",MaterialKind.WASTE)
    graph.react("dash","DASH: cut empty-adaptor amplicons with the Cas9–gRNA complex",StepKind.DIGESTION,
                ((pre_library,InputRole.TEMPLATE),(empty_amplicon,InputRole.TEMPLATE),
                 (guide,InputRole.STRUCTURAL),(cas9,InputRole.CATALYST)),(dash_library,cut_empty))
    graph.discard(cut_empty,"Cas9 cleavage prevents amplification of the empty-adaptor artifact; the cut DNA is not physically removed")
    library=graph.product("library",FINAL_LIBRARY.name,MaterialKind.LIBRARY,state=FINAL_LIBRARY)
    graph.react("library_pcr","Final two-step indexed PCR of DASHed molecules (Supplementary Table 9)",
                StepKind.PCR,
                ((dash_library,InputRole.TEMPLATE),(fwd_primer,InputRole.STRUCTURAL),
                 (rev_primer,InputRole.STRUCTURAL),(q5,InputRole.CATALYST)),(library,),
                geometry=PCRGeometry.unavailable(
                    FINAL_LIBRARY,reason="The compact branch model records the final construct but not separate PCR-primer landing objects."))
    graph.sequence(library,SequencingPlan(
        library=FINAL_LIBRARY,primers=SEQ_PRIMERS,
        acquired_roles=tuple(p.role for p in SEQ_PRIMERS),
        read_uses={p.role:("demultiplex" if "Index" in p.role else "include") for p in SEQ_PRIMERS}))
    graph.annotate(SourceDiscrepancy(
        "library_pcr","Cycle programs printed for the two PCRs",
        (SourceClaim("Methods, 'RT through the triazole link and pre-amplification'",
                     "The purified cDNA was PCR amplified for 6 cycles (Supplementary Table 8)."),
         SourceClaim("Supplementary Tables 8 and 9",
                     "Both tables print the same 6-cycle 98 °C/60 °C/72 °C program, and the printed rows carry the anneal and extension times in separate lines.")),
        "The model states 6 cycles for the pre-amplification and labels the final PCR as the two-step program of Supplementary Table 9 without asserting a different cycle count.",
        "The Methods call the final program 'two-step' while Table 9 reproduces the pre-amplification table, so the final cycle number is not independently printed.",
        "Whether the final library PCR uses more than 6 cycles is not stated."))
    graph.declare_reviewed_scope(tuple(graph.materials)+tuple(graph.steps),source)
    return graph.require_publishable()
