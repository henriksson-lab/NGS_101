"""Published QuRIE-seq branches; source statements do not imply exact products."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
import illumina as il
from chemdraw import Construct, Scene, Segment, feature, protocol_graph_panels
from pcr import PCRGeometry
from protocol_graph import (ProtocolGraph, MaterialKind as K, InputRole as R, StepKind as S,
    SequencingPlan, Terminal, Disposition, SourceDiscrepancy, SourceClaim, AnnotationSlot)
from protocol_views import render_annotations
from page import head, info

PAGE_FILE='qurie.html'
SOURCE='https://doi.org/10.1016/j.crmeth.2021.100070; published STAR Methods and Supplementary Table S3'
PE2_N6_AS_WRITTEN=il.PE_READ2_ARM+'NNNNN'
PE2=il.P5+il.PE_READ2_PRIMER[1:23]
PE1_TAIL=il.TRUSEQ_READ1[3:19]
CUSTOM_READ1=il.PE_READ2_ARM[2:]
CUSTOM_INDEX=il.INDEX2_PRIMER_RC[:30]
CUSTOM_READ2=il.TRUSEQ_READ1[3:]
RNA_IVT_C=37
PROTEIN_IVT_C=50

def pe1():
    return Construct([Segment('P7',il.P7),
        Segment('six-base library index','N'*6,placeholder=True,
                feature=feature('qurie.sample','sample_index','unknown')),
        Segment('PE1 annealing arm',PE1_TAIL)],name='Source PE1 index-primer pattern')

def random_primer():
    return Construct([Segment('PE2-N6 arm',il.PE_READ2_ARM),
        Segment('five Ns printed in Table S3','N'*5,placeholder=True,inferred=True,note='Actual manufactured tract unresolved: five printed Ns versus randomhexamer/name and38bp column.')],
        name='PE2-N6 as printed; source name/length discrepancy unresolved')

def _one(g,key,label,kind,template,reagent,role,extra=()):
    reagent_kind={R.CATALYST:K.CATALYST,R.CONDITION:K.CONDITION,R.STRUCTURAL:K.OLIGO}[role]
    r=g.input(key+'_reagent',reagent,reagent_kind)
    out=g.product(key+'_product',label,K.MOLECULAR)
    g.react(key,label,kind,((template,R.TEMPLATE),(r,role),*extra),(out,))
    return out

def _cleanup(g,key,template,reagent,extra=()):
    selection=g.input(key+'_beads',reagent,K.SELECTION)
    kept=g.product(key+'_kept',key.replace('_',' ')+' retained material',K.MOLECULAR)
    waste=g.product(key+'_waste',key.replace('_',' ')+' unretained material/reagents',K.WASTE)
    g.react(key,key.replace('_',' '),S.CLEANUP,
            ((template,R.TEMPLATE),(selection,R.SELECTION),*extra),(kept,waste))
    g.discard(waste,'Removed during this source-declared purification')
    return kept

def _inventory(g,key,label,items):
    """Co-present component views of a source-declared supply, never fresh stocks."""
    stock=g.input(key,label,K.MOLECULAR)
    parts=tuple(g.product(key+'_'+suffix,part_label,kind) for suffix,part_label,kind in items)
    g.react(key+'_components','Expose inventory: '+label,S.MAGIC,((stock,R.TEMPLATE),),parts,
        changes_molecule=False,reason='Component inventory views of the same actual source supply. This is not a physical separation or extra additions. Exact formulation, occupancy and reaction allocations are unassigned unless explicitly labelled.')
    return parts

def _branch_finish(g,prefix,template,p1,p2):
    # Source explicitly consumes a qPCR aliquot to determine the final PCR cycles.
    reserve=g.product(prefix+'_pcr_input','Actual9.5uL retained template for '+prefix+' library PCR',K.MOLECULAR)
    qc=g.product(prefix+'_qpcr','Actual0.5uL qPCR calibration aliquot for '+prefix+' PCR cycle selection',K.QC)
    g.react(prefix+'_qpcr_split','Reserve aliquot for qPCR cycle determination',S.SPLIT,
            ((template,R.TEMPLATE),),(reserve,qc))
    g.mark_terminal(Terminal(qc.id,Disposition.QC,'Actual0.5uL analytic template; source separate qPCR uses10uL2xKAPA/EvaGreen/PE1-PE2. Detailed analytic assay continuation is outside this graph; the aliquot is not returned to finalPCR.'))
    enzyme,nts,buffer=_inventory(g,prefix+'_kapa','Actual branch KAPA HiFi HotStart PCR mix; final reaction12.5uL2x',(
        ('enzyme','Polymerase component of same branch KAPA supply',K.CATALYST),
        ('nt','Nucleotide component of same branch KAPA supply',K.MOLECULAR),
        ('buffer','Buffer component of same branch KAPA supply',K.CONDITION)))
    amplified=g.product(prefix+'_amplified',prefix+' PE1/PE2-amplified library material',K.MOLECULAR)
    g.react(prefix+'_pcr','Amplify '+prefix+' branch with PE1/PE2',S.PCR,
        ((reserve,R.TEMPLATE),(p1,R.STRUCTURAL),(p2,R.STRUCTURAL),(enzyme,R.CATALYST),(nts,R.STRUCTURAL),(buffer,R.CONDITION)),(amplified,),
        geometry=PCRGeometry.unavailable(Scene(),reason='Full branch templates/bead-promoter and RT-primer geometry remain unassigned in this capped reading; individual source oligos do not determine the final molecule.'))
    ratio='0.7× AMPure XP' if prefix=='rna' else '1.2× AMPure XP'
    selection=g.input(prefix+'_final_beads',ratio,K.SELECTION)
    lib=g.product(prefix+'_library',prefix+' QuRIE-seq library',K.LIBRARY)
    waste=g.product(prefix+'_final_waste','Unretained '+prefix+' PCR/cleanup material',K.WASTE)
    g.react(prefix+'_final_cleanup','Purify '+prefix+' library',S.CLEANUP,
            ((amplified,R.TEMPLATE),(selection,R.SELECTION)),(lib,waste))
    g.discard(waste,'Removed with final library bead cleanup')
    assay=g.product(prefix+'_final_qc','Analytical material for Qubit/BioAnalyzer; amounts unassigned',K.QC)
    retained_lib=g.product(prefix+'_retained_library','Library retained after analytical QC selection',K.LIBRARY)
    g.react(prefix+'_qc_split','Reserve final-library analytical material',S.SPLIT,((lib,R.TEMPLATE),),(assay,retained_lib),changes_molecule=False)
    g.mark_terminal(Terminal(assay.id,Disposition.QC,'Source final Qubit and BioAnalyzer characterization; aliquot sizes and detailed assay reactions outside this preparation model.'))
    contribution=g.product(prefix+'_pool_contribution','Selected '+prefix+' contribution for source joint sequencing pool',K.LIBRARY)
    remainder=g.product(prefix+'_library_remainder','Untransferred '+prefix+' library stock; amount unassigned',K.LIBRARY)
    g.react(prefix+'_select_pool','Select '+prefix+' library contribution',S.SPLIT,((retained_lib,R.TEMPLATE),),(contribution,remainder),changes_molecule=False)
    g.mark_terminal(Terminal(remainder.id,Disposition.UNRESOLVED,'Source joint sequencing does not imply whole-stock consumption. Contribution and residual amounts remain unassigned.'))
    return contribution

def protocol_graph():
    g=ProtocolGraph('QuRIE-seq — BJAB RNA and intracellular-protein branches',authority=SOURCE,representation='published_instance')
    cells=g.input('cells','BJAB cells after the source experimental stimulation',K.BIOLOGICAL,primary=True)
    dsp,spdp=_inventory(g,'fix_reagent','Actual source2xDSP/SPDP fixative mixture,5mM each in PBS',(
        ('dsp','DSP component of actual fixative supply',K.MOLECULAR),('spdp','SPDP component of actual fixative supply',K.MOLECULAR)))
    fixed=g.product('fix_product','Source reversibly fixed BJAB cells',K.MOLECULAR)
    fix_condition=g.input('fix_conditions','Source1hroomtemperature fixation',K.CONDITION)
    g.react('fix','DSP/SPDP reversible fixation',S.CROSS_LINKING,((cells,R.TEMPLATE),(dsp,R.STRUCTURAL),(spdp,R.STRUCTURAL),(fix_condition,R.CONDITION)),(fixed,))
    permeable=_one(g,'permeabilize','Tris quench and Triton permeabilization (combined source treatment)',S.PERMEABILIZATION,
                    fixed,'100mM Tris-HCl,150mM NaCl,0.1% TritonX100',R.CONDITION)
    conjugate=g.input('antibody_tags','DNA-conjugated antibodies prepared with DBCO-S-S-NHS and 5′ azide oligos',K.OLIGO)
    stained=g.product('stained','Cells immunostained with intracellular/surface antibody-DNA panel',K.MOLECULAR)
    block_ds,block_rnase,block_buffer=_inventory(g,'blocking_supply','Actual source prepared blocking mixture',(
        ('ds','Source0.2mg/ml dextran sulfate component',K.MOLECULAR),('rnase','Source0.5U/ml RNasin component',K.CATALYST),('buffer','Source0.5xprotein-free blocking buffer inPBS',K.CONDITION)))
    g.react('stain','Block and stain cells',S.BINDING,((permeable,R.TEMPLATE),(conjugate,R.STRUCTURAL),(block_ds,R.STRUCTURAL),(block_rnase,R.CATALYST),(block_buffer,R.CONDITION)),(stained,))
    washed=_cleanup(g,'cell_wash',stained,'Two blocking-buffer washes after staining; recover singlets')
    suspension_bsa,suspension_rnase,suspension_buffer=_inventory(g,'cell_suspension_supply','Actual post-wash cell suspension supply',(
        ('bsa','Source5mg/mlBSA component; carryover unassigned',K.MOLECULAR),('rnase','Source0.5U/mlRNasin component',K.CATALYST),('buffer','SourcePBS suspensionbuffer',K.CONDITION)))
    optiprep=g.input('optiprep','Source15percentv/vOptiprep just before encapsulation',K.MOLECULAR)
    prepared=g.product('prepared_cells','Source washed/singlet-selected cell suspension for encapsulation',K.MOLECULAR)
    g.react('prepare_cells','Resuspend and singlet-select cells before encapsulation',S.MAGIC,
        ((washed,R.TEMPLATE),(suspension_bsa,R.STRUCTURAL),(suspension_rnase,R.CATALYST),(suspension_buffer,R.CONDITION),(optiprep,R.STRUCTURAL)),(prepared,),changes_molecule=False,
        reason='Source suspension and singlet selection interface; exact retained cell numbers and carryover composition unassigned, no claim that all cells survive.')
    bead=g.input('bead_primers','Source inDrop-based primer-loaded gel-bead supply; whole primer/matrix structure unassembled',K.MOLECULAR)
    rt,dtt,rnaseout,rt_buffer=_inventory(g,'superscript','Actual sourceRTinlet mixture:1.3xRTpremix/MgCl2/DTT/RNaseOUT/SuperScriptIII/water',(
        ('enzyme','SuperScriptIII component of sourceRTinlet supply',K.CATALYST),
        ('dtt','DTT component of same sourceRTsupply',K.MOLECULAR),
        ('rnaseout','RNaseOUT component of same sourceRTsupply',K.CATALYST),
        ('buffer','RTpremix/MgCl2/buffer classes of same supply; other components unassigned',K.CONDITION)))
    oil=g.input('droplet_oil','Actual HFE7500 oil/RANsurfactant source continuous phase',K.CONDITION)
    droplets=g.product('droplets','Aggregate source cell/bead/RT droplet collection; occupancy/multiplicity unassigned',K.MOLECULAR)
    continuous=g.product('continuous_phase','Continuous emulsion phase and unencapsulated material; composition/amounts unassigned',K.MOLECULAR)
    g.react('partition','Co-encapsulate cells, loaded beads and RTsupply before reduction',S.PARTITION,
        ((prepared,R.TEMPLATE),(bead,R.STRUCTURAL),(dtt,R.STRUCTURAL),(rt,R.CATALYST),(rnaseout,R.CATALYST),(rt_buffer,R.CONDITION),(oil,R.CONDITION)),(droplets,continuous),changes_molecule=False,
        reason='Source2nl droplet collection interface. Component ports are co-present views of the same RTinlet supply, not independent additions or duplicated stocks; cell/bead occupancy and exact arrays unassigned.')
    reduced=g.product('reduce_product','Source DTT-treated material within droplets; exact crosslink products unassembled',K.MOLECULAR)
    reduction_condition=g.input('reduce_conditions','DTT already present in droplets; precise reduction kinetics unassigned',K.CONDITION)
    g.react('reduce','Reduce source disulfide crosslinks after emulsification',S.CHEMICAL_CONVERSION,((droplets,R.TEMPLATE),(dtt,R.STRUCTURAL),(reduction_condition,R.CONDITION)),(reduced,),reason='DTT is already co-present from the sourceRTsupply; not freshly added to sealed droplets. Exact reduced protein/linker products remain unassembled.')
    uv=g.input('uv','8minUV exposure at10mW/cm² after emulsification',K.CONDITION)
    primers=g.product('released_primers','Released source poly(dT)cell-barcoding primer classes; wholebody unassembled',K.OLIGO)
    photo_material=g.product('photoreleased_material','Co-present droplet cell/bead material after primer photorelease; exact matrix products unassigned',K.MOLECULAR)
    g.react('release_primers','Photorelease bead primers within source droplets',S.CHEMICAL_CONVERSION,
            ((reduced,R.TEMPLATE),(uv,R.CONDITION)),(primers,photo_material))
    copied=g.product('barcoded_mix','Barcoded RNA-derived cDNA and antibody-derived DNA products',K.MOLECULAR)
    g.react('barcode','Prime poly(A)/poly(dA) and copy in droplets:50°C2h,70°C15min',S.REVERSE_TRANSCRIPTION,
            ((photo_material,R.TEMPLATE),(primers,R.STRUCTURAL),(rt,R.CATALYST),(rnaseout,R.CATALYST),(rt_buffer,R.CONDITION)),(copied,))
    recovered=_cleanup(g,'recover',copied,'PFO emulsion break; pellet debris; separate aqueous phase and Spin-X purification',extra=((continuous,R.TEMPLATE),))
    exo=_one(g,'exo','Exonuclease I treatment:37°C30min',S.DIGESTION,recovered,'Exonuclease I in FastDigest buffer',R.CATALYST)
    fractionator=g.input('fractionator','0.6× AMPure XP size separation',K.SELECTION)
    rna=g.product('rna_fraction','Long RNA-derived material retained on beads',K.MOLECULAR)
    protein=g.product('protein_fraction','Short antibody-derived material retained in supernatant',K.MOLECULAR)
    g.react('fractionate','Retain both RNA and protein size fractions',S.SPLIT,
            ((exo,R.TEMPLATE),(fractionator,R.SELECTION)),(rna,protein))
    p1=g.input('pe1','PE1 with six-base sample/library index',K.OLIGO,state=pe1())
    p2=g.input('pe2','PE2 primer',K.OLIGO,state=Construct([Segment('PE2',PE2)]))
    sss_buffer=g.input('rna_sss_buffer','Actualsource2uL10xSSSbuffer plus17uLsample; exact buffer composition unassembled',K.CONDITION)
    ds=_one(g,'rna_sss','RNA-branch second-strand synthesis:16°C1.5h,65°C20min',S.EXTENSION,rna,'Actualsource1uLNEBNextSSSenzyme; molecular formulation unassembled',R.CATALYST,extra=((sss_buffer,R.CONDITION),))
    rna_t7,rna_nt,rna_buf=_inventory(g,'rna_ivt_reagent','Actual60uL RNAbranchHiScribeT7IVTmixture',(
        ('enzyme','Source10xdilutedT7RNApolymerasemix',K.CATALYST),('nt','SourceATP/CTP/GTP/UTP,10mM each',K.MOLECULAR),('buffer','Source1xIVTreactionbuffer',K.CONDITION)))
    rna_ivt=g.product('rna_ivt_product','RNAbranch IVTproduct; wholepromoter/template unassembled',K.MOLECULAR)
    g.react('rna_ivt','RNAbranch T7IVT37C15h',S.TRANSCRIPTION,((ds,R.TEMPLATE),(rna_t7,R.CATALYST),(rna_nt,R.STRUCTURAL),(rna_buf,R.CONDITION)),(rna_ivt,))
    rna_clean=_cleanup(g,'rna_ivt_cleanup',rna_ivt,'1.3× AMPure XP')
    fragmented=_one(g,'rna_fragment','Fragment amplified RNA:70°C3min',S.FRAGMENTATION,rna_clean,'RNA fragmentation reagents; STOP mix on ice',R.CONDITION)
    rna_fragment_clean=_cleanup(g,'rna_fragment_cleanup',fragmented,'Source AMPure/STOP mixture')
    random=g.input('pe2_n6','PE2-N6 as printed; effective N count unresolved',K.OLIGO,state=random_primer())
    prime=g.input('primescript_enzyme','Actual separate PrimeScriptRTenzyme supply; allocationviews for two branch RTreactions',K.CATALYST)
    prime_nt=g.input('primescript_nt','Actual separate dNTP supply; source1uL10nM each perRT, printed wording retained',K.MOLECULAR)
    prime_buffer=g.input('primescript_buffer','Actual separate5xPrimeScriptbuffer stock,4uL perRT',K.CONDITION)
    prime_rnase=g.input('primescript_rnaseout','Actual separate RNaseOUTstock,1uL perRT',K.CATALYST)
    rna_cdna=g.product('rna_cdna','Random-primed cDNA from amplified RNA',K.MOLECULAR)
    g.react('rna_rt','INFERRED effective random tract — RNA-branch RT',S.REVERSE_TRANSCRIPTION,
        ((rna_fragment_clean,R.TEMPLATE),(random,R.STRUCTURAL),(prime,R.CATALYST),(prime_nt,R.STRUCTURAL),(prime_buffer,R.CONDITION),(prime_rnase,R.CATALYST)),(rna_cdna,))
    rna_rt_clean=_cleanup(g,'rna_rt_cleanup',rna_cdna,'1.2× AMPure XP')
    rna_contribution=_branch_finish(g,'rna',rna_rt_clean,p1,p2)
    prot_t7,prot_nt,prot_buf=_inventory(g,'protein_ivt_reagent','Actual60uL proteinbranchHiScribeT7IVTmixture',(
        ('enzyme','Source10xdilutedT7RNApolymerasemix',K.CATALYST),('nt','SourceexplicitATP/CTP/GTP/UTP10mM each; prose saysdNTP, sugar-class wording conflict preserved',K.MOLECULAR),('buffer','Source1xIVTreactionbuffer',K.CONDITION)))
    protein_ivt=g.product('protein_ivt_product','Proteinbranch IVTproduct; fullpromoter/template unassembled',K.MOLECULAR)
    g.react('protein_ivt','Proteinbranch T7IVT50C15h as published',S.TRANSCRIPTION,((protein,R.TEMPLATE),(prot_t7,R.CATALYST),(prot_nt,R.STRUCTURAL),(prot_buf,R.CONDITION)),(protein_ivt,))
    protein_clean=_cleanup(g,'protein_ivt_cleanup',protein_ivt,'1.5× AMPure XP; elute20µL')
    backup=g.product('protein_backup','10µL post-IVT protein backup at−80°C',K.MOLECULAR)
    protein_rt_input=g.product('protein_rt_input','10µL post-IVT protein material for RT',K.MOLECULAR)
    g.react('protein_backup_split','Retain post-IVT protein backup',S.SPLIT,
            ((protein_clean,R.TEMPLATE),),(backup,protein_rt_input))
    g.mark_terminal(Terminal(backup.id,Disposition.STORED,'Source explicitly stores half of the20µL eluate at−80°C'))
    specific=g.input('protein_rt_primers','PE2-NNNN-Next1 and PE2-NNNN-BioHash2; full designs unassigned',K.OLIGO)
    protein_cdna=g.product('protein_cdna','Sequence-primed cDNA of amplified protein reporters',K.MOLECULAR)
    g.react('protein_rt','Protein-branch sequence-specific RT',S.REVERSE_TRANSCRIPTION,
            ((protein_rt_input,R.TEMPLATE),(specific,R.STRUCTURAL),(prime,R.CATALYST),(prime_nt,R.STRUCTURAL),(prime_buffer,R.CONDITION),(prime_rnase,R.CATALYST)),(protein_cdna,))
    protein_rt_clean=_cleanup(g,'protein_rt_cleanup',protein_cdna,'1.5× AMPure XP')
    protein_contribution=_branch_finish(g,'protein',protein_rt_clean,p1,p2)
    pooled=g.product('pooled','Distinct selected RNA/protein libraries for source joint sequencing; sample multiplicity unassigned',K.LIBRARY)
    g.react('pool','Pool source RNA/protein library contributions without concatenation',S.POOL,((rna_contribution,R.TEMPLATE),(protein_contribution,R.TEMPLATE)),(pooled,),changes_molecule=False)
    loading=g.product('loading_aliquot','Selected sequencing loading aliquot; amount unassigned',K.LIBRARY)
    remainder=g.product('pool_remainder','Untransferred sequencing pool stock; amount unassigned',K.LIBRARY)
    g.react('select_loading','Select joint sequencing loading aliquot',S.SPLIT,((pooled,R.TEMPLATE),),(loading,remainder),changes_molecule=False)
    g.mark_terminal(Terminal(remainder.id,Disposition.UNRESOLVED,'Sequencing does not imply consumption of the entire library pool.'))
    g.sequence(loading,SequencingPlan(acquisition='unspecified',limitation='NextSeq500 and custom sequencing primers in TableS3. Separate RNA/protein libraries remain distinct in this pool; exact lengths, library strands, run-primer landing and acquired masks unassembled.'))
    g.declare_reviewed_scope((
        'reduce_conditions','fix_conditions','fix_reagent_components','blocking_supply','cell_suspension_supply','optiprep','prepare_cells','droplet_oil','partition','superscript_components','rna_ivt_reagent_components','protein_ivt_reagent_components','primescript_enzyme','primescript_nt','primescript_buffer','primescript_rnaseout','rna_sss_buffer','blocking_supply_components','cell_suspension_supply_components',
        'rna_kapa_components','protein_kapa_components','rna_qc_split','rna_final_qc','rna_select_pool','rna_library_remainder',
        'protein_qc_split','protein_final_qc','protein_select_pool','protein_library_remainder','pool','select_loading','pool_remainder','loading_aliquot',
        'cells','fix_reagent','fix_product','fix','permeabilize_reagent','permeabilize_product','permeabilize',
        'antibody_tags','stained','stain','cell_wash_beads','cell_wash_kept','cell_wash_waste','cell_wash',
        'reduce_product','reduce','bead_primers','uv','released_primers','release_primers',
        'superscript','barcoded_mix','barcode','recover_beads','recover_kept','recover_waste','recover',
        'exo_reagent','exo_product','exo','fractionator','rna_fraction','protein_fraction','fractionate','pe1','pe2',
        'rna_sss_reagent','rna_sss_product','rna_sss','rna_ivt_reagent','rna_ivt_product','rna_ivt',
        'rna_ivt_cleanup_beads','rna_ivt_cleanup_kept','rna_ivt_cleanup_waste','rna_ivt_cleanup',
        'rna_fragment_reagent','rna_fragment_product','rna_fragment','rna_fragment_cleanup_beads',
        'rna_fragment_cleanup_kept','rna_fragment_cleanup_waste','rna_fragment_cleanup','pe2_n6',
        'rna_cdna','rna_rt','rna_rt_cleanup_beads','rna_rt_cleanup_kept','rna_rt_cleanup_waste','rna_rt_cleanup',
        'rna_pcr_input','rna_qpcr','rna_qpcr_split','rna_kapa','rna_amplified','rna_pcr','rna_final_beads',
        'rna_library','rna_final_waste','rna_final_cleanup','protein_ivt_reagent','protein_ivt_product','protein_ivt',
        'protein_ivt_cleanup_beads','protein_ivt_cleanup_kept','protein_ivt_cleanup_waste','protein_ivt_cleanup',
        'protein_backup','protein_rt_input','protein_backup_split','protein_rt_primers','protein_cdna','protein_rt',
        'protein_rt_cleanup_beads','protein_rt_cleanup_kept','protein_rt_cleanup_waste','protein_rt_cleanup',
        'protein_pcr_input','protein_qpcr','protein_qpcr_split','protein_kapa','protein_amplified','protein_pcr',
        'protein_final_beads','protein_library','protein_final_waste','protein_final_cleanup'),SOURCE)
    g.annotate(SourceDiscrepancy(target='pe2_n6',title='PE2-N6 name/length and printed sequence differ',
        claims=(SourceClaim('Main RNA-library Methods and TableS3 name/length','Random hexamer/PE2-N6; Table length38bp.'),
                SourceClaim('Published rendered TableS3,PDF14,printed13','32-base arm followed by five printed Ns:37 characters.')),
        model_choice='Preserve five Ns as written; do not claim effective priming geometry.',
        evidence='Original rendered table inspected; this is not a text-extraction-only discrepancy.',
        uncertainty='Actual manufactured N count/order is unresolved; source follow-up remains open.',slot=AnnotationSlot.BEFORE_WORKFLOW))
    return g.require_publishable()

def render_page():
    from protocol_views import assemble_page,WorkflowSection,OligoSection,OligoEntry
    g=protocol_graph()
    return assemble_page(g,(head('QuRIE-seq — intracellular protein and RNA'),'<div class="wrap" data-protocol-complete="true"><h1>QuRIE-seq</h1>',
        '<p class="research-notes"><a href="01_qurie_seq.html">Research notes</a></p>',info('<a href="https://doi.org/10.1016/j.crmeth.2021.100070">Published QuRIE-seq STARMethods and TableS3</a>'),
        '<h2>Barcode intracellular protein and RNA together</h2>',
        '<p>The selected BJAB preparation fixes and permeabilizes cells for DNA-tagged antibody staining. Cells, loaded beads and RTsupplies are co-encapsulated before DTT reduction and UV primer release. Barcoding copies both RNA and antibody-derived DNA. A0.6x size split retains long RNA-derived material on beads and short protein reporters in supernatant.</p>',
        info('INFERRED — the PE2-N6 random tract is a printed five-N reference whose effective manufactured length is unresolved. Source-qualified preparation with unassembled molecular states: exact antibody-tag/bead/promoter strands, occupancy, complete library products and read-primer geometry remain unknown in this model. No whole inDrop recipe or literal cell/bead array is borrowed. Antibody conjugates enter as prepared source reagents; detailed upstream manufacture is outside the graph.'),
        '<h2>Two physically separate amplification branches</h2>',
        '<p>RNA material undergoes second-strand synthesis, T7IVT at37C, fragmentation and random-primer RT. Protein material undergoes source T7IVT at50C and sequence-specific RT;10uL of the20uL post-IVT eluate is stored as backup. Both branches retain0.5uL qPCR analytical versus9.5uL PCR-template allocations, then separate final QC and library-stock remainders before joint sequencing.</p>',
        '<h2>Source sequence and condition discrepancies</h2>',
        AnnotationSlot.BEFORE_WORKFLOW,
        '<p>The protein-IVT Methods call the nucleotide supply “dNTP” but explicitly list ATP,CTP,GTP,UTP. The model preserves those listed ribonucleotide identities and the literal wording difference, without a deoxynucleotide substitution. The source RNAIVT37C/proteinIVT50C conditions stay separate. Source RT “10nM eachdNTP” wording is retained; no silent dose correction is made.</p>',
        OligoSection('Known source primer orders',(
            OligoEntry('PE1 index pattern',tuple(pe1())),OligoEntry('PE2',(Segment('PE2',PE2),)),
            OligoEntry('PE2-N6 exactly as printed',tuple(random_primer())),
            OligoEntry('CustomRead1',(Segment('CustomRead1',CUSTOM_READ1),)),
            OligoEntry('CustomIndex',(Segment('CustomIndex',CUSTOM_INDEX),)),
            OligoEntry('CustomRead2',(Segment('CustomRead2',CUSTOM_READ2),))),
            note_html='Printed TableS3 orders, not proof of full product or run-primer landing. The PE1 six-base sample index is unassigned; printed PE2-N6 has fiveNs despite its name/length.'),
        '<h2>Sequencing and validation limits</h2>',
        '<p>The source jointly sequences RNA/protein libraries on NextSeq500. Exact run lengths, acquired index masks and cell/UMI feature positions remain unassembled; known custom primer bodies alone do not establish them. Full source-stock and material-fate declarations follow, with unresolved molecular claims visible. No efficiency, purity, specificity or gene/protein measurement is established by the diagrams.</p>',
        WorkflowSection('Full BJAB preparation and molecular limits'),'</div>'))

def selftest():
    from checks import Check
    c=Check();g=protocol_graph()
    c('Both source library contributions reach a joint sequencing pool',len(g.sequencing),1)
    c('Post-IVT protein backup retained',g.terminals['protein_backup'].disposition,Disposition.STORED)
    c('Both qPCR and final characterization fractions accounted',sum(t.disposition is Disposition.QC for t in g.terminals.values()),4)
    c('Source branch-specific IVT temperatures retained',(RNA_IVT_C,PROTEIN_IVT_C),(37,50))
    from storage import sources_root
    f=sources_root()/Path(__file__).resolve().parents[1].name/'j.crmeth.2021.100070_supp_mmc1.pdf.txt'
    if f.exists():
        text=''.join(f.read_text().split())
        for seq in (PE2,CUSTOM_READ1,CUSTOM_INDEX,CUSTOM_READ2):c('Published TableS3 transcription '+seq,seq in text)
        row=f.read_text().split('PE2-N6',1)[1].split('PE1',1)[0]
        c('Source PE2-N6 arm and printed five-N line',il.PE_READ2_ARM in row and 'NNNNN' in row.splitlines())
    else:c.skip('Published primer transcription','External supplement unavailable')
    c('uncertain manufactured PE2-N6 tract defaults to inferred',list(random_primer())[-1].inferred,True)
    c('DTT reduction follows physical encapsulation',g._producer[g.steps['reduce'].inputs[0].material],'partition')
    c('protein backup cannot be pooled',g._consumers['protein_backup'],[])
    c('pool combines actual separate branch contributions',[p.material for p in g.steps['pool'].inputs],['rna_pool_contribution','protein_pool_contribution'])
    c('same actual RTstock supplies components',[p.material for p in g.steps['superscript_components'].inputs],['superscript'])
    c('IVT consumes explicit source ribonucleotide inventories',[any(p.material==b+'_ivt_reagent_nt' and p.role is R.STRUCTURAL for p in g.steps[b+'_ivt'].inputs) for b in ('rna','protein')],[True,True])
    c('untransferred joint pool survives loading',g.terminals['pool_remainder'].disposition,Disposition.UNRESOLVED)
    c.report(exit_on_fail=True)
