"""Published DAb-seq patient branch; source-local capture sites, not a Tapestri recipe."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
import illumina as il
from chemdraw import (Construct, Scene, Segment, feature, complement_segments, revcomp,
                     panel, protocol_graph_panels)
from protocol_graph import ProtocolGraph, MaterialKind as K, InputRole as R, StepKind as S, SequencingPlan, Terminal, Disposition
from pcr import PCRGeometry
from page import head, info

PAGE_FILE='dab.html'
SOURCE='https://doi.org/10.1038/s41467-021-21810-3; published Methods and Supplementary Tables1/2/6'
TAG_PREFIX='GGCTTGTTGTGATTCGACGA'
TAG_SITE='AAGCTTGTTTCTGTGCACTGAG'
RUNX4_WRITTEN='CTCAGTGCACAGAAACAAG*C*T*T'
RUNX4=RUNX4_WRITTEN.replace('*','')
CAPTURE_WRITTEN='/5Biosg/GGCTTGTTGTGATTCGACGA/3C6/'
TAG_WRITTEN='/5AzideN/GGCTTGTTGTGATTCGACGATNNNNN[barcode]NNNNNAAGCTTGTTTCTGTGCACTGAG'
CD3='AACGCTTC'
READ1='GCCTGTCCGCGGAAGCAGTGGTATCAACGCAGAGTAG'
P5_ARM='GCCTGTCCGCGGAAGCAGTGGTATCAACGCAG'
P5_INDEX1=il.P5+'TAGATCGC'+P5_ARM
AB_P7_INDEX1_WRITTEN=il.P7+'TACTACGC'+'GTGACTGGAGTTCCTTGGCACCCGAGAATTCCAGGCTTGTTGTGATTCG*A*C*G'
AB_P7_INDEX1=AB_P7_INDEX1_WRITTEN.replace('*','')

# Source Table2 splits a single molecular identifier around an antibody barcode.
def tag():
    return Construct([
        Segment('5prime prefix',TAG_PREFIX),Segment('prefix boundary T','T'),
        Segment('left UMI half','N'*5,'umi',placeholder=True,
            feature=feature('dab.umi.left','umi','random',group='dab.umi',part='left')),
        Segment('antibody identity','[barcode]','cbc',placeholder=True,length_bp=8,feature=feature('dab.antibody','feature_barcode','whitelist',whitelist=SOURCE+';Table2')),
        Segment('right UMI half','N'*5,'umi',placeholder=True,
            feature=feature('dab.umi.right','umi','random',group='dab.umi',part='right')),
        Segment('RUNX1 capture site',TAG_SITE)],name='Published Table2 antibody-tag pool pattern; UMI halves unspecified')

def runx_binding():
    sc=Scene();sc.strand('tag',list(tag()),label='input antibody tag')
    primer=[Segment('RUNX1_4 primer site',RUNX4,note=RUNX4_WRITTEN)]
    sc.anneal('primer',primer,to='tag',pair=('RUNX1_4 primer site','RUNX1 capture site'),label='bead primer site')
    sc.note('primer','Binding site only: the full5prime barcode/handle extension is not drawn. Three source phosphorothioate bonds remain in the written oligo record')
    return sc

def capture_binding():
    sc=Scene();sc.strand('copied_tag',complement_segments(list(tag())),label='PCR-copied tag strand')
    sc.anneal('probe',[Segment('capture probe',TAG_PREFIX)],to='copied_tag',
        pair=('capture probe',"5prime prefix'"),label='capture oligo',mod5='biotin',mod3='C6 spacer')
    sc.note('probe','Probe anneals to the synthesized complementary tag strand after denaturation; it is not a PCR primer')
    return sc

def _cleanup(g,key,template,reagent,label,kind=K.MOLECULAR):
    beads=g.input(key+'_reagent',reagent,K.SELECTION)
    retained=g.product(key+'_kept',label,kind)
    waste=g.product(key+'_waste','Unretained '+key.replace('_',' ')+' material',K.WASTE)
    g.react(key,key.replace('_',' '),S.CLEANUP,
        ((template,R.TEMPLATE),(beads,R.SELECTION)),(retained,waste))
    g.discard(waste,'Removed at this source-declared recovery/wash step')
    return retained

def _inventory(g,key,label,items):
    stock=g.input(key,label,K.MOLECULAR)
    components=tuple(g.product(key+'_'+suffix,title,kind) for suffix,title,kind in items)
    g.react(key+'_components','Expose source inventory: '+label,S.MAGIC,((stock,R.TEMPLATE),),components,
        changes_molecule=False,reason='Co-present inventory views of the same actual named supply, not fresh duplicate stocks or a physical separation. Aggregate reaction allocations/whole formulations remain unassigned.')
    return components


def protocol_graph():
    g=ProtocolGraph('DAb-seq — patient DNA and antibody libraries',authority=SOURCE,representation='published_instance')
    cells=g.input('cells','Patient mononuclear cells prepared for staining',K.BIOLOGICAL,primary=True)
    tags=g.input('antibody_conjugates','Prepared antibody-DNA panel; source5azide click chemistry precedes this input, residual adduct geometry unassembled',K.OLIGO,state=tag())
    block=g.input('blocking_supply','Actualseparate10uLFcblocking solution supply',K.MOLECULAR)
    dextran=g.input('blocking_dextran','Actualseparate4uL1percentdextransulfate supply',K.MOLECULAR)
    salmon=g.input('blocking_dna','Actualseparate4uL10mg/mlsalmon spermDNA supply',K.MOLECULAR)
    stained=g.product('stained','Source antibody-DNA-stained patient cells',K.MOLECULAR)
    g.react('stain','Block10min and stain30min on ice',S.BINDING,((cells,R.TEMPLATE),(tags,R.STRUCTURAL),(block,R.STRUCTURAL),(dextran,R.STRUCTURAL),(salmon,R.STRUCTURAL)),(stained,))
    wash=g.input('staining_washes_reagent','ActualPBS5percentFBS washsolution; proteincomposition/carryover unassigned',K.MOLECULAR)
    wash_program=g.input('staining_wash_conditions','Sourcefive rounds5mlwash,400gx4min; finalCellBufferresuspension3Mcells/ml',K.CONDITION)
    washed=g.product('staining_washes_kept','Sourcewashed stainedcells; residualFBS/othermaterial unassigned',K.MOLECULAR)
    wash_waste=g.product('staining_washes_waste','Removed excessantibody/washsolution and pelleting losses',K.WASTE)
    g.react('staining_washes','Sourcefivecellwash rounds and resuspension',S.CLEANUP,((stained,R.TEMPLATE),(wash,R.STRUCTURAL),(wash_program,R.CONDITION)),(washed,wash_waste));g.discard(wash_waste,'Sourcepelleting/washsupernatants, not purification by hypotheticalbeads.')
    protease,detergent,lysisbuffer=_inventory(g,'lysis_supply','Actualsource AMLV1 lysisbuffer supply',(
        ('protease','Protease component of same source lysis supply',K.CATALYST),('detergent','Milddetergent component of same source lysis supply',K.MOLECULAR),('buffer','Other lysisbuffer components unassembled',K.CONDITION)))
    oil1=g.input('lysis_oil','Source first-stage droplet continuous phase; formulation/dose unassembled',K.CONDITION)
    droplet1=g.product('lysis_droplets','Aggregate first-stage cell/lysis droplets; occupancy unassigned',K.MOLECULAR)
    phase1=g.product('lysis_phase','First-stage continuous phase; remains with emulsion through transfer',K.MOLECULAR)
    g.react('encapsulate_lysis','Firstencapsulation: cells and source lysis supplies',S.PARTITION,
        ((washed,R.TEMPLATE),(protease,R.CATALYST),(detergent,R.STRUCTURAL),(lysisbuffer,R.CONDITION),(oil1,R.CONDITION)),(droplet1,phase1),changes_molecule=False)
    lysate=g.product('lysate','Source lysed-cell droplet contents; protein digestion products/copresent material unassembled',K.MOLECULAR)
    g.react('lysis','Digest cellular proteins50C1h within first-stage droplets',S.DEPROTEINIZATION,((droplet1,R.TEMPLATE),(protease,R.CATALYST),(detergent,R.STRUCTURAL)),(lysate,))
    heat=g.input('heat','Source80C10min',K.CONDITION)
    inactive=g.product('inactivated','Lysed-cell droplets with protease heat-inactivated',K.MOLECULAR)
    g.react('heat_inactivate','Heat-inactivate protease before barcoding module',S.THERMAL_TREATMENT,((lysate,R.TEMPLATE),(heat,R.CONDITION)),(inactive,))
    beads=g.input('barcoded_beads','Actual AMLV1 loaded hydrogel-bead stock; sourceforward cores includingRUNX1_4, fullbarcodehandle/wholepolymers unassembled',K.MOLECULAR)
    reverse=g.input('reverse_panel','ActualsourceTable1reversepanel withRUNX1_4reverseOMITTED; remaining wholepool unexpanded in capped model',K.OLIGO)
    enzyme,nts,buffer=_inventory(g,'polymerase','Actual MissionBio barcoding polymerasemixture; allocations among source8tubes unassigned',(
        ('enzyme','Polymerase component of same barcoding supply',K.CATALYST),('nt','Nucleotide component of same barcoding supply',K.MOLECULAR),('buffer','Buffer component of same barcoding supply',K.CONDITION)))
    oil2=g.input('barcoding_oil','Source second-stage droplet oil; exact formulation/dose unassigned',K.CONDITION)
    droplet2=g.product('barcoding_droplets','Aggregate source second-stage lysate/bead/PCR droplets',K.MOLECULAR)
    phase2=g.product('barcoding_phase','Co-present second-stage continuous phase',K.MOLECULAR)
    g.react('encapsulate_barcoding','Secondencapsulation: lysed-cell droplets and barcoding supplies',S.PARTITION,
        ((inactive,R.TEMPLATE),(phase1,R.TEMPLATE),(beads,R.STRUCTURAL),(reverse,R.STRUCTURAL),(enzyme,R.CATALYST),(nts,R.STRUCTURAL),(buffer,R.CONDITION),(oil2,R.CONDITION)),(droplet2,phase2),changes_molecule=False,
        reason='Actual second microfluidic stage. Components are co-present views of the same supplies, not fresh additions to sealed tubes. Whole bead/primer designs and occupancy unassembled.')
    # Source has eight technical PCR tubes, not eight invented cells/patients.
    tube_keys=('p1a','p1b','p2a','p2b','p3a','p3b','p4a','p4b')
    tubes={key:g.product(key+'_tube','Anonymous technical PCRtube slot '+key+'; emulsion contents/volume unassigned',K.MOLECULAR) for key in tube_keys}
    g.react('distribute_tubes','Distribute barcoding emulsion into source eight PCRtubes',S.SPLIT,((droplet2,R.TEMPLATE),(phase2,R.TEMPLATE)),tuple(tubes.values()),changes_molecule=False)
    uv=g.input('uv','Source8minUV BlakRayXX15L',K.CONDITION)
    joint_program=g.input('joint_program','Source95C10min;20x95C30s/72C10s/61C4min/72C30s;72C2min;4Chold',K.CONDITION)
    amplified={}
    for key in tube_keys:
        primer=g.product(key+'_released_primers','Co-present cell-barcoded forwardprimer classes in '+key+'; wholebodies unassembled',K.OLIGO)
        material=g.product(key+'_photoreleased','Copresent lysate/matrix material after primer release in '+key,K.MOLECULAR)
        g.react(key+'_photorelease','Photorelease source bead primers in technicaltube '+key,S.CHEMICAL_CONVERSION,((tubes[key],R.TEMPLATE),(uv,R.CONDITION)),(primer,material))
        out=g.product(key+'_amplicons','Joint DNA/tag amplicon classes in technicaltube '+key,K.MOLECULAR)
        g.react(key+'_joint_pcr','Source20cycle joint DNA/antibody barcodingPCR in '+key,S.PCR,
            ((material,R.TEMPLATE),(primer,R.STRUCTURAL),(reverse,R.STRUCTURAL),(enzyme,R.CATALYST),(nts,R.STRUCTURAL),(buffer,R.CONDITION),(joint_program,R.CONDITION)),(out,),
            geometry=PCRGeometry.unavailable(runx_binding(),reason='Source local RUNX1_4 binding only; full cell-barcoded primer/genomic/wholeamplicon geometry unassembled. Reagent ports refer to the same already encapsulated barcoding supplies.'))
        amplified[key]=out
    probe_scene=Scene();probe_scene.strand('capture',[Segment('capture probe',TAG_PREFIX)],mod5='biotin',mod3='C6 spacer')
    probe=g.input('capture_probe','Actualsource /5Biosg/.../3C6/ probe supply;0.6uM in antibody fraction before denaturation',K.OLIGO,state=probe_scene)
    streptavidin=g.input('streptavidin','Actual10uLprewashedstreptavidinbeads perrecoverypool; prepared reagent interface, washdetails outside scope',K.SELECTION)
    capture_wash=g.input('capture_wash','Source15minroomtemperature capture, twoPBSwashes, resuspend30uLwater',K.CONDITION)
    denature=g.input('denature','Source95C5min then snapcool on ice',K.CONDITION)
    p5=g.input('p5','ActualsourcecustomP5primer family; sample memberunassigned',K.OLIGO,state=Construct((Segment('P5',il.P5),Segment('unassignedsourcei5member','N'*8,placeholder=True,feature=feature('dab.sample.i5','sample_index','unknown')),Segment('customP5arm',P5_ARM))))
    dp7=g.input('dna_p7','SourceNexteraP7N7XXprimer family; exact selectedmember unassigned',K.OLIGO)
    ap7=g.input('antibody_p7','ActualsourceantibodyP7family; samplememberunassigned, terminalPSbonds sourcequalified',K.OLIGO,state=Construct((Segment('P7',il.P7),Segment('unassignedsourceantibodyi7member','N'*8,placeholder=True,feature=feature('dab.sample.antibody.i7','sample_index','unknown')),Segment('AbP7arm',AB_P7_INDEX1[len(il.P7)+8:],note=AB_P7_INDEX1_WRITTEN))))
    libenz,libnts,libbuf=_inventory(g,'library_mix','Actual25uL MissionBioLibraryMix allocation per source libraryPCR;50uL totalreaction',(
        ('enzyme','Polymerase component of actual LibraryMix',K.CATALYST),('nt','Nucleotide component of same LibraryMix',K.MOLECULAR),('buffer','Buffer component of same LibraryMix',K.CONDITION)))
    library_program=g.input('library_program','Source95C3min;10DNAor20Abcycles98C20s/62C20s/72C45s;72C2min;4Chold',K.CONDITION)
    phix=g.input('phix','ActualsourcePhiXcontrolDNA, source20percent sequencing spike-in; fraction convention/allocations unassigned',K.LIBRARY)
    dna_plan=SequencingPlan(acquired_roles=('Read 1','Read 2','Index 1 (i7)','Index 2 (i5)'),read_lengths={'Read 1':150,'Read 2':150},limitation='Patient DNA: NovaSeq6000SP2x150, dualindex/customRead1. Terminal streams are library contributions, not a claim of four independent runs; run-grouping, indexmask and fullstrand/landing geometry unassembled.')
    ab_plan=SequencingPlan(acquired_roles=('Read 1','Read 2','Index 1 (i7)','Index 2 (i5)'),read_lengths={'Read 1':38,'Read 2':39},limitation='Patient antibody: NextSeq550Read1=38/Read2=39, dualindex/customRead1; CBCinRead1, Ab8/splitUMI10inRead2. Contributions remain distinct from DNAassay; run-grouping/fullprimer/cyclespans unassembled.')
    for pair in ('p1','p2','p3','p4'):
        joined=g.product(pair+'_paired','Sourcepairedtechnicaltube material '+pair,K.MOLECULAR)
        g.react(pair+'_pool','Pool two technical PCRtubes as one source recoverypair',S.POOL,((amplified[pair+'a'],R.TEMPLATE),(amplified[pair+'b'],R.TEMPLATE)),(joined,),changes_molecule=False,
            reason='Anonymous slots are named by pair membership, not patient or experimental tube IDs. Four independent paired pools; no all-eight pooling followed by splitting.')
        rec=_cleanup(g,pair+'_recover',joined,'ActualsourceExtractionAgent/water; aqueous transfer','Sourcepaired aqueous DNA/tag recovery')
        selection=g.input(pair+'_size_beads','Actual0.75xAMPureXP perrecoveredpair',K.SELECTION)
        dna=g.product(pair+'_dna_beads','SourceDNApanel fraction bound to beads',K.MOLECULAR)
        protein=g.product(pair+'_protein_supernatant','Source short antibody fraction retained in supernatant',K.MOLECULAR)
        g.react(pair+'_size_split','Retain both source sizefractions',S.SPLIT,((rec,R.TEMPLATE),(selection,R.SELECTION)),(dna,protein))
        dna=_cleanup(g,pair+'_dna_wash',dna,'Source two80percentEtOHwashes, elute30uLwater','Source30uL DNApanel eluate')
        premix=g.product(pair+'_probe_mix','Actualcopresent antibodyproduct/probe mixture before heat',K.MOLECULAR)
        g.react(pair+'_add_probe','Add same sourcecaptureprobe before denaturation',S.MAGIC,((protein,R.TEMPLATE),(probe,R.STRUCTURAL)),(premix,),changes_molecule=False,reason='Source0.6uMprobe is added before95C; ensuinghybridization uses the same copresent reagent, not a fresh stock.')
        single=g.product(pair+'_denatured','Denatured antibody/probe mixture',K.MOLECULAR)
        g.react(pair+'_denature','Heat then snapcool actualprobe/product mixture',S.THERMAL_TREATMENT,((premix,R.TEMPLATE),(denature,R.CONDITION)),(single,))
        hybrid=g.product(pair+'_hybrid','Sourcecapture interface on complementary tagcopy; wholebarcodedproduct unassembled',K.MOLECULAR,state=capture_binding())
        g.react(pair+'_hybridize','Hybridize copresent captureprobe after snapcool',S.HYBRIDIZATION,((single,R.TEMPLATE),(probe,R.STRUCTURAL)),(hybrid,))
        ab=g.product(pair+'_captured','Sourcebeadbound antibodyproducts resuspended30uLwater',K.MOLECULAR)
        reject=g.product(pair+'_capture_waste','Unbound antibodyfraction and washmaterial',K.WASTE)
        g.react(pair+'_capture','Capture on streptavidin, isolate and wash',S.AFFINITY_CAPTURE,((hybrid,R.TEMPLATE),(streptavidin,R.SELECTION),(capture_wash,R.CONDITION)),(ab,reject));g.discard(reject,'Actual magnetic separation/twoPBSwashes; efficiency unassigned.')
        for assay,stock,primer,cycles,plan in (('dna',dna,dp7,10,dna_plan),('antibody',ab,ap7,20,ab_plan)):
            key=pair+'_'+assay
            template=g.product(key+'_template','Selected4ngDNA; template+water15uL, sourceeluateallocationvolumeunassigned' if assay=='dna' else 'Actual15uL beadbound antibodytemplate selectedfrom30uL',K.MOLECULAR)
            remainder=g.product(key+'_template_remainder','UntransferredDNAeluate; amount/volumeunassigned' if assay=='dna' else 'Untransferred15uL beadbound antibody stock',K.MOLECULAR)
            g.react(key+'_select_template','Select actualsource libraryPCRtemplate contribution',S.SPLIT,((stock,R.TEMPLATE),),(template,remainder),changes_molecule=False)
            g.mark_terminal(Terminal(remainder.id,Disposition.UNRESOLVED,'DNA4ng in15uL does not mean15uLof30uLeluate. Antibody explicitly uses15uLof30uLbeadsuspension; retainedstockfateunassigned.'))
            product=g.product(key+'_pcr_product','Indexed '+assay+'PCRproduct for source recoverypair '+pair,K.MOLECULAR)
            g.react(key+'_library_pcr','Source'+str(cycles)+'cycle '+assay+'libraryPCR',S.PCR,
                ((template,R.TEMPLATE),(p5,R.STRUCTURAL),(primer,R.STRUCTURAL),(libenz,R.CATALYST),(libnts,R.STRUCTURAL),(libbuf,R.CONDITION),(library_program,R.CONDITION)),(product,),
                geometry=PCRGeometry.unavailable(runx_binding(),reason='Source5uLeachP5/P7at4uMstock and actual25uLLibraryMix; wholeintermediate/barcodeprimer/indexassignments unassembled. No exact final-library molecule is invented.'))
            lib=_cleanup(g,key+'_final_cleanup',product,'Actual0.7xAMPureXP; elute12uLwater','Source12uL '+assay+'library eluate',K.LIBRARY)
            qc=g.product(key+'_qc','SourceBioAnalyzer/Qubit analyticalmaterial; quantitiesunassigned',K.QC)
            retained=g.product(key+'_retained','Library stockretained after analyticalselection',K.LIBRARY)
            g.react(key+'_select_qc','Reserve actualsource final analyticalmaterial',S.SPLIT,((lib,R.TEMPLATE),),(qc,retained),changes_molecule=False)
            g.mark_terminal(Terminal(qc.id,Disposition.QC,'ActualsourceBioAnalyzerHighSensitivityDNA andQubit; detailedassayreactions andaliquots outside this preparation graph.'))
            contribution=g.product(key+'_contribution','Selected library-to-run contribution; volume/run-groupingunassigned',K.LIBRARY)
            stock_rest=g.product(key+'_library_remainder','Untransferred library stock; amountunassigned',K.LIBRARY)
            g.react(key+'_select_contribution','Select library contribution for source assay run',S.SPLIT,((retained,R.TEMPLATE),),(contribution,stock_rest),changes_molecule=False)
            g.mark_terminal(Terminal(stock_rest.id,Disposition.UNRESOLVED,'Run-grouping across technical/patient libraries is unassigned; no whole-library consumption asserted.'))
            control_mix=g.product(key+'_control_mixture','Source20percentPhiX sequencing-mixture interface for selectedlibrary contribution',K.LIBRARY)
            g.react(key+'_control_mix','Include actualsourcePhiX control contribution',S.POOL,((contribution,R.TEMPLATE),(phix,R.TEMPLATE)),(control_mix,),changes_molecule=False,
                reason='Source20percentPhiX; fractionconvention/doses and groupingofsame-assaylibraries unassigned. Each terminal is acontributioninterface, not oneflowcell. PatientDNA/antibodyassays are never combined.')
            loading=g.product(key+'_loading','Selected final sequencing loadingaliquot',K.LIBRARY)
            mix_rest=g.product(key+'_loading_remainder','Untransferred selected sequencingmixture; amountunassigned',K.LIBRARY)
            g.react(key+'_select_loading','Select loadingaliquot without consuming entire mixture',S.SPLIT,((control_mix,R.TEMPLATE),),(loading,mix_rest),changes_molecule=False)
            g.mark_terminal(Terminal(mix_rest.id,Disposition.UNRESOLVED,'Loading/residual volumes unassigned.'))
            g.sequence(loading,plan)
    # Independently curated source inventory: exact technical counts, not graph enumeration.
    g.declare_reviewed_scope((
        'blocking_dextran','blocking_dna','staining_wash_conditions','joint_program','library_program','cells','antibody_conjugates','blocking_supply','stain','staining_washes_reagent','staining_washes',
        'lysis_supply','lysis_supply_components','lysis_oil','encapsulate_lysis','lysis','heat','heat_inactivate',
        'barcoded_beads','reverse_panel','polymerase','polymerase_components','barcoding_oil','encapsulate_barcoding','distribute_tubes','uv',
        'capture_probe','streptavidin','capture_wash','denature','p5','dna_p7','antibody_p7','library_mix','library_mix_components','phix',
    )+tuple(key+suffix for key in ('p1a','p1b','p2a','p2b','p3a','p3b','p4a','p4b') for suffix in ('_photorelease','_joint_pcr'))
     +tuple(pair+suffix for pair in ('p1','p2','p3','p4') for suffix in ('_pool','_recover_reagent','_recover','_size_beads','_size_split','_dna_wash_reagent','_dna_wash','_add_probe','_denature','_hybridize','_capture'))
     +tuple(pair+'_'+assay+suffix for pair in ('p1','p2','p3','p4') for assay in ('dna','antibody') for suffix in ('_select_template','_template_remainder','_library_pcr','_final_cleanup_reagent','_final_cleanup','_select_qc','_qc','_select_contribution','_library_remainder','_control_mix','_select_loading','_loading_remainder','_loading')),SOURCE)
    return g.require_publishable()

def render_page():
    from protocol_views import assemble_page,WorkflowSection,OligoSection,OligoEntry
    return assemble_page(protocol_graph(),(head('DAb-seq — patient DNA and antibody tags'),'<div class="wrap" data-protocol-complete="true"><h1>DAb-seq</h1>',
        '<p class="research-notes"><a href="01_dab_seq.html">Research notes</a></p>',info('<a href="https://doi.org/10.1038/s41467-021-21810-3">Published DAb-seq patient preparation and SupplementaryTables1/2/6</a>'),
        '<h2>Two encapsulations link genotype and protein measurements</h2>',
        '<p>Stained patient cells first enter lysis droplets. After protease treatment and inactivation, lysed-cell droplets enter a second barcoding module with loaded beads, the modified reverseprimer pool and PCRsupplies. UV releases primers before joint DNA/tag amplification.</p>',
        info('Source-qualified patient preparation. Full MissionBioAMLV1 bead handles, remaining whole primer pool, genomic amplicons and final library/primer coordinates are unassembled. Prepared antibody conjugates enter as source reagents; exact residual adduct/protein structures are not drawn. No current vendor recipe, assigned patient genotype or arbitrary barcode is borrowed.'),
        '<h2>Eight technical tubes become four paired recoveries</h2>',
        '<p>The source distributes barcoded droplets into eight PCR tubes and pools them as pairs. The anonymous p1a/p1b through p4a/p4b slots below are named by that pair relation, not actual experimental IDs or eight cells. Four independent aqueous recoveries each retain DNA amplicons on0.75xbeads and short antibody products in supernatant. No all-eight pool followed by splitting is substituted.</p>',
        '<p>The capture probe enters each antibodyfraction before95 °C denaturation/snapcool, then binds the copied complementary tag strand. Streptavidin selection/washes leave30 µLbeadsuspension;15 µL goes to20cycle antibodyPCR. DNA eluate is30 µL, but its4 ng template is prepared in15 µLwater: this does not establish that15 µLofeluate was consumed. DNA-PCR uses10 cycles. Untransferred stocks remain explicit.</p>',
        '<h2>Known source interfaces and primer references</h2>',
        panel(runx_binding().rows(),caption='Source RUNX1_4 core binds tag 3′ site; whole cell-barcoded primer is unassembled. Source PS bonds remain in the recorded oligo; matching core does not establish full amplicon.'),
        panel(capture_binding().rows(),caption='Source5biotin/3C6capture probe binds the complementary tagcopy. It is not a PCR primer and its blocking modification is preserved.'),
        OligoSection('Printed library-primer Index1 references',(
            OligoEntry('P5Index1 reference',(Segment('P5Index1',P5_INDEX1),)),
            OligoEntry('AntibodyP7Index1 reference',(Segment('AbP7Index1',AB_P7_INDEX1,note=AB_P7_INDEX1_WRITTEN),),mods='source last three phosphorothioate bonds; raw ordering in notes'),
            OligoEntry('CustomRead1',(Segment('CustomRead1',READ1),))),
            note_html='These printed Index1 examples are not assigned to any patient or technical pool. Actual selected family members and whole genomic/bead primer designs remain unassigned. The tag contains one split 5+5 UMI around an8 bp antibody barcode.'),
        '<h2>Separate patient sequencing runs and explicit controls</h2>',
        '<p>DNA libraries use NovaSeq6000SP2x150; antibody libraries use NextSeq550Read1=38/Read2=39. Both use dualindex/custom Read1 and source20%PhiX. Cell barcode is parsed fromRead1; antibody8 bp and splitUMI10bp fromRead2. Exact index mask/primer coordinates and grouping of technical/patient libraries into runs remain unassigned. The eight terminal streams represent library contributions, not eight independent flowcells; DNA and antibody assays are never combined.</p>',
        '<p>Final BioAnalyzer/Qubit material, librarycontributions and final loading aliquots leave explicit stock remainders. Sequencing and drawn local pairing do not establish genotype fidelity, capture efficiency, purity or exact final molecules.</p>',
        '<details><summary>Full source preparation, technical arrays and molecular validation</summary>',WorkflowSection('All eight tubes, four recoveries and library contributions'),'</details></div>'))

def selftest():
    from checks import Check
    c=Check();g=protocol_graph()
    c('Four actual paired recoveries retain separate DNA and antibody contributions',len(g.sequencing),8)
    c('Tag capture selects complementary copy, not input polarity',revcomp(RUNX4),TAG_SITE)
    c('C6-blocked capture oligo is not a library-PCR primer',any(p.material=='capture_probe' and p.role is R.STRUCTURAL
        for key in (pair+'_'+assay+'_library_pcr' for pair in ('p1','p2','p3','p4') for assay in ('dna','antibody')) for p in g.steps[key].inputs),False)
    c('Two UMI components share molecular identity',[(s.feature.group,s.feature.part) for s in tag() if s.feature and s.feature.role=='umi'],[('dab.umi','left'),('dab.umi','right')])
    from storage import sources_root
    folder=sources_root()/Path(__file__).resolve().parents[1].name
    f=folder/'s41467-021-21810-3_supp_41467_2021_21810_MOESM1_ESM.pdf.txt'
    if f.exists():
        text=''.join(f.read_text().split())
        for seq in (TAG_WRITTEN,RUNX4_WRITTEN,P5_INDEX1,AB_P7_INDEX1_WRITTEN):c('Source-written sequence/modifications '+seq,seq in text)
    else:c.skip('Published supplement transcription','External supplement unavailable')
    f=folder/'s41467-021-21810-3_PMC7952600.xml.txt'
    if f.exists():
        text=''.join(f.read_text().split())
        for seq in (CAPTURE_WRITTEN,READ1):c('Published Methods transcription '+seq,seq in text)
    else:c.skip('Published Methods transcription','External article unavailable')
    c('first and second encapsulations remain distinct',[(g.steps[k].kind)for k in ('encapsulate_lysis','encapsulate_barcoding')],[S.PARTITION,S.PARTITION])
    c('technical emulsion is allocated into eight real PCRtube slots',len(g.steps['distribute_tubes'].outputs),8)
    c('each source pair pools only its two actual technicalslots',[[p.material for p in g.steps[pair+'_pool'].inputs]for pair in ('p1','p2','p3','p4')],[[pair+'a_amplicons',pair+'b_amplicons']for pair in ('p1','p2','p3','p4')])
    c('probe is introduced before thermaldenaturation',g.steps['p1_denature'].inputs[0].material,'p1_probe_mix')
    c('eight final analytical selections are QC',sum(t.disposition is Disposition.QC for t in g.terminals.values()),8)
    c('source patient assay runplans stay separate',len({id(p)for p in g.sequencing.values()}),2)
    c('whole named LibraryMix supplies only one componentfactory input',[p.material for p in g.steps['library_mix_components'].inputs],['library_mix'])
    c.report(exit_on_fail=True)
