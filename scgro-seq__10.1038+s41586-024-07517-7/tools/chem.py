"""scGRO-seq source-owned clickedRNA, nativeSS RT and actual mixed PCR inputs."""
from chemdraw import Construct,Scene,Segment,feature,complement_segments,revcomp
from protocol_model import role
from pcr import PCRGeometry
from protocol_graph import InputRole,MaterialKind,ProtocolGraph,SequencingPlan,SourceClaim,SourceDiscrepancy,StepKind,Terminal,Disposition
import seqprimers as sp
TITLE='scGRO-seq — single-cell nascentRNA byclick chemistry'
NOTES='01_scgro-seq.html'
SOURCE='<a href="https://doi.org/10.1038/s41586-024-07517-7">Mahat etal.,Nature2024</a>; publishedMethods andSupplementaryTables7–9.'
SUMMARY='Chain-terminatednascentRNA is clicked toazideDNA in single-wells, physicallypooled/cleaned/fragmented and reverse-transcribed throughtriazole. ActualSS cDNA feeds prePCR; DASHcleavage fragments remain molecularclasses throughfinalPCR ratherthanbeingphysicallydiscarded.'
CAVEAT='SelectedpublishedTable7design3: physical5′azideT plusACAGG,N5UMI,J10cellbarcode,sourcearm andphysical3′SpC3 block. Designs1 and2 are separate alternatives (10/6 and12/4 barcode/UMI); nohairpin primer inheritance into design3. SourceRTprimer,TSO,indexprimerstrands andcompletefirst-cycle/nativeendgeometry remain unassignedwithinboundedpass. Chemicalazide/triazole/SpC3 residues belongtooriginaltemplate, notnewDNAcopies. Source2×runonbuffer concentrations are notfinal concentrations. PublishedTable9sixcycles isselectedwhileMethodstwo-step wording remainsqualified.'
CELL=feature('scgro.cell','cell_barcode','whitelist',whitelist='sourceperwell5AzScBc set3')
UMI=feature('scgro.umi','umi','random')
I5=feature('scgro.i5','sample_index','unknown',note='unique dualindex; physicalwidthunassigned, sourceacquired8cycles')
I7=feature('scgro.i7','sample_index','unknown',note='unique dualindex; physicalwidthunassigned, sourceacquired8cycles')
ARM='AGATCGGAAGAGCGTCGTGTAG'
FINAL_LIBRARY=None
SEQ_PRIMERS=()
SHOW_READ_CONTENT=False
SEQUENCING_ENDING='PublishedNovaSeqSP100:Read1=64,Read2=43,i7=8,i5=8. Reads trimR1CCCCTGTCTCTTATACACAT/R2AGATCGGAAGAGCGTCGTGT beforeR1cellbarcode/UMI extraction; exactnativeprimerboundariesunassigned.'
READ_SPECS=tuple(sp.read_spec(r,n,'cellbarcode/UMI+nascentRNA' if r=='Read 1' else 'nascentRNA' if r=='Read 2' else 'sampledualindex','extract_then_include' if r=='Read 1' else 'include' if r=='Read 2' else 'demultiplex') for r,n in [('Read 1',64),('Read 2',43),('Index 1 (i7)',8),('Index 2 (i5)',8)])
def azide_parts():return [Segment('physicalazideT','T'),Segment('DASHdesignACAGG','ACAGG'),Segment('N5UMI','N'*5,feature=UMI),Segment('J10cellbarcode','N'*10,feature=CELL),Segment('sourceDNAarm',ARM),role('physical3′SpC3; nonDNAblock')]
def copied_ids():return [role('copiedcellbarcode; nativeorientation/positionbounded',length_bp=10,molecular_feature=CELL),role('copiedUMI',length_bp=5,molecular_feature=UMI)]
def native_state(stage):
 sc=Scene()
 if stage in ('raw','runon'):
  sc.strand('RNA',[role('engagednascentRNAwithphysicalU'),role('terminal3′OpropargylNMP;chainterminated' if stage=='runon' else 'originalnascentRNA3end')],label='originalRNA')
  sc.strand('DNA',[role('originalgenomicDNA')],label='independentDNA');sc.strand('protein',[role('independentnuclearproteins')],label='independentprotein');sc.strand('otherRNA',[role('independentotherRNA')],label='otherRNA');return sc
 if stage=='unclicked':sc.strand('azide',azide_parts(),label='originaldesign3 unclickedDNA;3′SpC3 retained');return sc
 if stage in ('clicked','fragmented'):
  sc.strand('RNA-DNA',[role('originalnascentRNA;physicalU'+(';Znfragmentedextentunassigned' if stage=='fragmented' else '')),role('sourceRNA3′O–triazole–DNA5′azideTlinkage'),*azide_parts()],label='originalclickedRNA–DNA, not a copiedDNA strand')
  sc.strand('unclicked',azide_parts(),label='possibleoriginalunclickedDNA carriedthroughcleanup');input_residuals(sc);return sc
 empty=stage in ('emptySS','emptyamp','cut')
 parts=[role('nativeRTprimer/firstcopyend'),role('copiedsourceDNAarm; fullnativeprimerextentunassigned'),*copied_ids(),role('copiedACAGG/azideTbasecontext; triazolechemicalbondnotcopied'),*([] if empty else [role('antisensenascentcDNA copiedthroughtriazole')]),role('sourceTSO-derivedDNA footprint; switchjunctionextentunassigned; SpC3notcopied')]
 if stage in ('amplified','emptyamp','dashed','indexed'):parts += [role('installedi5primerlineage',molecular_feature=I5),role('installedi7primerlineage',molecular_feature=I7)]
 sc.strand('daughter',parts,label=('Cas9cleavedemptyDNAfragments; breakextentunassigned' if stage=='cut' else 'actual '+stage+' DNA'))
 sc.strand('residuals',[role('independentoriginalclickedRNA/unclickedDNA/RToligos/earliercopies orprocessingproducts; recoveryunassigned')],label='copying/nonamplificationnotphysicaldisposal')
 if stage in ('dashed','indexed'):sc.strand('cleavedclasses',[role('possibleCas9cleavedDNA fragments anduncleavedemptyamplicons; nointerveningphysicalseparation')],label='independentDASHclasses')
 input_residuals(sc)
 return sc

def input_residuals(sc):
 sc.strand('possibleoriginalgenome',[role('possible originalgenomicDNA; classspecificphysicalrecoveryunassigned')],label='independent carriedgenomicDNA')
 sc.strand('possibleotherRNA',[role('possible originalotherRNA oritsprocessingproducts; physicalrecoveryunassigned')],label='independent carriedotherRNA')
 sc.strand('possibleprotein',[role('possible originalproteins orprocessingproducts; classspecificphysicalrecoveryunassigned')],label='independent proteinclass')

def pcr_input(*states):
 sc=Scene()
 for i,state in enumerate(states):
  for key,strand in state.strands.items():sc.strand(str(i)+key,strand.segs,label=strand.label)
 return sc

def oligos():
    return ['<b>5′-AzScBc design 1, 43 nt</b> — 5′-[azideT]JJJJJJJJJJNNNNNNGATCGTCGGACTGTAGAACTCTGAAC-3′ (J = cell barcode; N = UMI)', '<b>5′-AzScBc design 2, 86 nt hairpin</b> — 5′-[azideT]NNNNJJJJJJJJJJJJGATCGTCGGACTGTAGAACTCTGCGGCCGTGCTCGTTTTCGAGCACGGCCGCAGAGTTCTACAGTCCGA-3′', '<b>5′-AzScBc design 3, 43 nt</b> — 5′-[azideT]ACAGGNNNNNJJJJJJJJJJAGATCGGAAGAGCGTCGTGTAG[SpC3]-3′']
REVIEWED=('cultured_cells', 'permeabilization_buffer', 'engaged_complexes', 'permeabilize', 'propargyl_nucleotide', 'run_on_buffer', 'engaged_polymerase', 'labelled_rna', 'run_on', 'runon_wash_condition', 'washed_runon', 'runon_wash_waste', 'runon_washes', 'sorted_cell', 'unsorted_suspension', 'sort_cells', 'click_oligo', 'cuaac_catalyst', 'cuaac_condition', 'clicked_rna', 'unclicked_azide', 'click_attachment', 'cuaac_cleanup', 'pooled_clicked_rna', 'carried_unclicked_azide', 'spent_cuaac_reagents', 'quench_and_clean', 'zinc_chloride', 'fragmented_clicked_rna', 'fragment_clicked_rna', 'reverse_transcriptase', 'rt_primer', 'template_switch_oligo', 'cdna', 'empty_cdna', 'reverse_transcription', 'denaturing_page', 'selected_cdna', 'residual_empty_cdna', 'excised_short_species', 'size_select', 'q5_master_mix', 'index_pcr_forward', 'index_pcr_reverse', 'pre_library', 'empty_amplicon', 'pre_amplify', 'cas9', 'dash_guide', 'dash_library', 'cut_empty', 'dash', 'library', 'library_pcr')
REVIEWED += ('preamp_clean_rule','preamp_cleanup')
def protocol_graph():
    source = 'Defining source: https://doi.org/10.1038/s41586-024-07517-7'
    graph = ProtocolGraph('scGRO-seq', authority=source)
    cells = graph.input('cultured_cells', 'V6.5 mouse ES cells', MaterialKind.BIOLOGICAL, primary=True, state=native_state('raw'))
    perm_buffer = graph.input('permeabilization_buffer', 'douncing or permeabilization buffer (Triton X-100, or Tween-20/NP-40)', MaterialKind.CONDITION)
    complexes = graph.product('engaged_complexes', 'permeabilized cells or isolated nuclei with engaged transcription complexes', MaterialKind.BIOLOGICAL, state=native_state('raw'))
    graph.react('permeabilize', 'Isolate nuclei or permeabilize cells and store in storage buffer', StepKind.PERMEABILIZATION, ((cells, InputRole.TEMPLATE), (perm_buffer, InputRole.CONDITION)), (complexes,))
    nucleotide = graph.input('propargyl_nucleotide', '50µM each in2×buffer;25µMfinal each 3′-(O-propargyl)-ATP/CTP/GTP/UTP (Jena NU-945 to NU-948)', MaterialKind.MOLECULAR)
    run_on_buffer = graph.input('run_on_buffer', 'nuclear run-on buffer (2×buffer400mMKCl/0.05%Sarkosyl;final200mMKCl/0.025%Sarkosyl, 5 min at 37 °C)', MaterialKind.CONDITION)
    polymerase = graph.input('engaged_polymerase', 'engaged RNA polymerase', MaterialKind.CATALYST)
    labelled = graph.product('labelled_rna', 'alkyne-terminated nascent RNA in nuclei', MaterialKind.MOLECULAR, state=native_state('runon'))
    graph.react('run_on', 'Run on engaged polymerases with 3′-O-propargyl NTPs (chain terminating)', StepKind.RUN_ON, ((complexes, InputRole.TEMPLATE), (nucleotide, InputRole.STRUCTURAL), (run_on_buffer, InputRole.CONDITION), (polymerase, InputRole.CATALYST)), (labelled,))
    wash_rule = graph.input('runon_wash_condition', '500g2min4°Cpellet;three150µlresuspensionbufferwashes;35µmmesh', MaterialKind.SELECTION)
    washed = graph.product('washed_runon', 'washedlabellednuclei/cells', MaterialKind.MOLECULAR, state=labelled.state)
    wash_waste = graph.product('runon_wash_waste', 'aspiratedrunonsupernatant/washes', MaterialKind.WASTE)
    graph.react('runon_washes', 'Pelletandwash excessunincorporatednucleotides beforebarcoding', StepKind.CLEANUP, ((labelled, InputRole.TEMPLATE), (wash_rule, InputRole.SELECTION)), (washed, wash_waste))
    graph.discard(wash_waste, 'Actualsourceaspiratedsupernatant/washes; classpuritynotassumed.')
    sorted_cell = graph.product('sorted_cell', 'one run-on cell or nucleus per barcoded well, lysed in 8 M urea', MaterialKind.MOLECULAR, state=native_state('runon'))
    unsorted = graph.product('unsorted_suspension', 'cells and nuclei left in the sorted suspension', MaterialKind.MOLECULAR)
    graph.react('sort_cells', 'FACS-sort single cells or nuclei into 96 wells containing 2.5 µl 8 M urea', StepKind.PARTITION, ((washed, InputRole.TEMPLATE),), (sorted_cell, unsorted))
    graph.mark_terminal(Terminal(unsorted.id, Disposition.UNRESOLVED, 'Unallocatedsourcecells notassumedphysicallydestroyed.'))
    click_oligo = graph.input('click_oligo', 'source5′azideT/3′SpC3 DNA, one cellbarcode perwell (design3)', MaterialKind.ADAPTER, source=source + ' (Supplementary Table 7, design 3)', state=native_state('unclicked'))
    copper = graph.input('cuaac_catalyst', 'CuAAC catalyst: 1 mM CuSO4, 800 mM sodium ascorbate, 5 mM BTTAA ligand', MaterialKind.CATALYST)
    cuaac_condition = graph.input('cuaac_condition', '15% PEG 8000, 2.66 M urea, 2 h at 50 °C', MaterialKind.CONDITION)
    clicked = graph.product('clicked_rna', 'barcoded clicked nascent RNA', MaterialKind.MOLECULAR, state=native_state('clicked'))
    unclicked = graph.product('unclicked_azide', 'unreacted 30 nM 5′-AzScBc DNA', MaterialKind.MOLECULAR, state=native_state('unclicked'))
    graph.react('click_attachment', 'CuAAC: click the 5′-azide barcode oligo onto the 3′ alkyne', StepKind.ADAPTER_ATTACHMENT, ((sorted_cell, InputRole.TEMPLATE), (click_oligo, InputRole.STRUCTURAL), (cuaac_condition, InputRole.CONDITION), (copper, InputRole.CATALYST)), (clicked, unclicked))
    cleanup_reagents = graph.input('cuaac_cleanup', '5 mM EDTA quench, TRIzol (PEG 8000) and 3 kDa Amicon filtration', MaterialKind.SELECTION)
    pooled = graph.product('pooled_clicked_rna', 'pooled, reagent-free clicked RNA from 96 wells', MaterialKind.MOLECULAR, state=native_state('clicked'))
    carried_azide = graph.product('carried_unclicked_azide', 'unclicked 5′-AzScBc DNA carried into RT', MaterialKind.MOLECULAR, state=native_state('unclicked'))
    spent = graph.product('spent_cuaac_reagents', 'removed PEG 8000, sodium ascorbate, CuSO4 and BTTAA', MaterialKind.WASTE)
    graph.react('quench_and_clean', 'Quench with EDTA, pool 96 wells and remove the CuAAC reagents', StepKind.CLEANUP, ((clicked, InputRole.TEMPLATE), (unclicked, InputRole.TEMPLATE), (cleanup_reagents, InputRole.SELECTION)), (pooled, carried_azide, spent), changes_molecule=False)
    graph.discard(spent, 'small-molecule reagents pass the 3 kDa membrane or partition into TRIzol')
    zinc = graph.input('zinc_chloride', '10 mM ZnCl2, 5 min at 65 °C', MaterialKind.CONDITION)
    fragments = graph.product('fragmented_clicked_rna', 'fragmented clicked RNA–DNA conjugates', MaterialKind.MOLECULAR, state=native_state('fragmented'))
    graph.react('fragment_clicked_rna', 'Fragment the clicked RNA with ZnCl2', StepKind.FRAGMENTATION, ((pooled, InputRole.TEMPLATE), (zinc, InputRole.CONDITION)), (fragments,))
    rt = graph.input('reverse_transcriptase', 'RNase H− M-MuLV reverse transcriptase (Thermo EP0751), 16% PEG 8000, 2 h at 50 °C', MaterialKind.CATALYST)
    rt_primer = graph.input('rt_primer', '1 µM RT primer (omitted for the hairpin-forming 86-nt design)', MaterialKind.OLIGO, inferred=True, source=source + ' (native orderingsequence unassigned in inspectedscope)')
    tso = graph.input('template_switch_oligo', '1 µM template-switching oligo', MaterialKind.OLIGO, inferred=True, source=source + ' (native orderingsequence unassigned in inspectedscope)')
    cdna = graph.product('cdna', 'barcoded nascent-RNA cDNA', MaterialKind.MOLECULAR, state=native_state('SS'))
    empty_cdna = graph.product('empty_cdna', 'empty cDNA from unclicked 5′-AzScBc DNA, non-templated CCC and TSO', MaterialKind.MOLECULAR, state=native_state('emptySS'))
    graph.react('reverse_transcription', 'Reverse-transcribe through the triazole link with template switching', StepKind.REVERSE_TRANSCRIPTION, ((fragments, InputRole.TEMPLATE), (carried_azide, InputRole.TEMPLATE), (rt_primer, InputRole.STRUCTURAL), (tso, InputRole.STRUCTURAL), (rt, InputRole.CATALYST)), (cdna, empty_cdna))
    gel = graph.input('denaturing_page', '10% denaturing PAGE size selection', MaterialKind.SELECTION)
    selected = graph.product('selected_cdna', 'gel-purified barcoded cDNA', MaterialKind.MOLECULAR, state=native_state('SS'))
    residual_empty = graph.product('residual_empty_cdna', 'empty cDNA co-migrating with the selected size range', MaterialKind.MOLECULAR, state=native_state('emptySS'))
    excised = graph.product('excised_short_species', 'unclicked 5′-AzScBc DNA and empty cDNA left in the gel', MaterialKind.WASTE)
    graph.react('size_select', 'Size-select the cDNA away from unclicked adaptor and empty cDNA', StepKind.SELECTION, ((cdna, InputRole.TEMPLATE), (empty_cdna, InputRole.TEMPLATE), (gel, InputRole.SELECTION)), (selected, residual_empty, excised), changes_molecule=False)
    graph.discard(excised, 'species outside the excised gel band are physically removed')
    q5 = graph.input('q5_master_mix', 'NEBNext Ultra II Q5 High-Fidelity 2× master mix (NEB M0544)', MaterialKind.CATALYST)
    fwd_primer = graph.input('index_pcr_forward', '0.5 µM P5-side indexed PCR primer (unique dual index)', MaterialKind.OLIGO, inferred=True, source=source + ' (native orderingsequence unassigned in inspectedscope)')
    rev_primer = graph.input('index_pcr_reverse', '0.5 µM P7-side indexed PCR primer (unique dual index)', MaterialKind.OLIGO, inferred=True, source=source + ' (native orderingsequence unassigned in inspectedscope)')
    pre_library = graph.product('pre_library', 'pre-amplified barcoded cDNA', MaterialKind.LIBRARY, state=native_state('amplified'))
    empty_amplicon = graph.product('empty_amplicon', 'pre-amplified empty-adaptor amplicon (ACAGG target)', MaterialKind.LIBRARY, state=native_state('emptyamp'))
    graph.react('pre_amplify', 'Pre-amplify to dsDNA: 6 cycles, 98/60/72 °C (Supplementary Table 8)', StepKind.PCR, ((selected, InputRole.TEMPLATE), (residual_empty, InputRole.TEMPLATE), (fwd_primer, InputRole.STRUCTURAL), (rev_primer, InputRole.STRUCTURAL), (q5, InputRole.CATALYST)), (pre_library, empty_amplicon), geometry=PCRGeometry.unavailable(pcr_input(selected.state, residual_empty.state), reason='The pre-amplified scGRO-seq amplicon and empty-adaptor junctions are not strand-resolved.'))
    clean_rule=graph.input('preamp_clean_rule','NativeprePCR cleanup referenced ascleanedPCR reaction; reagent/ratio unassigned',MaterialKind.SELECTION)
    clean_pre=graph.product('cleaned_pre_library','recoveredpreampnascentcopies class; samephysicalpoolasemptyclass',MaterialKind.LIBRARY,state=pre_library.state)
    clean_empty=graph.product('cleaned_empty_amplicon','co-recoveredemptyamplicons; samephysicalpool, notpurifiedbranch',MaterialKind.LIBRARY,state=empty_amplicon.state)
    preamp_waste=graph.product('preamp_cleanup_waste','actualnonretainedcleanupfraction; molecularpartitionunassigned',MaterialKind.WASTE)
    graph.react('preamp_cleanup','CleanmixedprePCR products beforeDASH; nativeconditionsunassigned',StepKind.CLEANUP,((pre_library,InputRole.TEMPLATE),(empty_amplicon,InputRole.TEMPLATE),(clean_rule,InputRole.SELECTION)),(clean_pre,clean_empty,preamp_waste));graph.discard(preamp_waste,'SourcecleanedPCR reaction; physicalnonretainedfraction only, no classpurityclaim.')
    cas9 = graph.input('cas9', '6.6 µM S. pyogenes Cas9 (NEB M0386T), 1 h at 37 °C', MaterialKind.CATALYST)
    guide = graph.input('dash_guide', '20 µM gRNA against the ACAGG-containing empty-adaptor amplicon', MaterialKind.OLIGO, inferred=True, source=source + ' (guide ordering sequence unresolved)')
    dash_library = graph.product('dash_library', 'barcoded cDNA surviving DASH', MaterialKind.LIBRARY, state=native_state('dashed'))
    cut_empty = graph.product('cut_empty', 'Cas9-cleaved empty-adaptor amplicons', MaterialKind.MOLECULAR, state=native_state('cut'))
    graph.react('dash', 'DASH: cut empty-adaptor amplicons with the Cas9–gRNA complex', StepKind.DIGESTION, ((clean_pre, InputRole.TEMPLATE), (clean_empty, InputRole.TEMPLATE), (guide, InputRole.STRUCTURAL), (cas9, InputRole.CATALYST)), (dash_library, cut_empty))
    library = graph.product('library', 'native scGRO-seq library', MaterialKind.LIBRARY, state=native_state('indexed'))
    graph.react('library_pcr', 'Final indexedPCR ofmixed DASHedmolecules:Table9sixcycles, sourceMethodstwo-step terminology', StepKind.PCR, ((dash_library, InputRole.TEMPLATE), (cut_empty, InputRole.TEMPLATE), (fwd_primer, InputRole.STRUCTURAL), (rev_primer, InputRole.STRUCTURAL), (q5, InputRole.CATALYST)), (library,), geometry=PCRGeometry.unavailable(pcr_input(dash_library.state, cut_empty.state), reason='Actual mixed post-DASH population; first-cycle/fullnativeprimer geometry unassigned.'))
    graph.sequence(library, SequencingPlan(read_lengths={'Read 1': 64, 'Read 2': 43, 'Index 1 (i7)': 8, 'Index 2 (i5)': 8}, acquired_roles=('Read 1', 'Read 2', 'Index 1 (i7)', 'Index 2 (i5)'), limitation=SEQUENCING_ENDING))
    graph.annotate(SourceDiscrepancy('library_pcr', 'Cycle programs printed for the two PCRs', (SourceClaim("Methods, 'RT through the triazole link and pre-amplification'", 'The purified cDNA was PCR amplified for 6 cycles (Supplementary Table 8).'), SourceClaim('Supplementary Tables 8 and 9', 'Both tables print the same 6-cycle 98 °C/60 °C/72 °C program, and the printed rows carry the anneal and extension times in separate lines.')), 'The model states 6 cycles for the pre-amplification and labels the final PCR as the two-step program of Supplementary Table 9 without asserting a different cycle count.', "The Methods call the final program 'two-step' while Table 9 reproduces the pre-amplification table, so Table9 explicitlyprints6cycles; wholeprogramdescribedastwo-step inMethods.", 'ThisgraphusesTable9sixcycles; additionaluninspectedconditionsarenotassigned.'))
    graph.declare_reviewed_scope(REVIEWED, source)
    return graph.require_publishable()
