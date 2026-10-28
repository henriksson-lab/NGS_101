"""Published DAb-seq patient branch; source-local capture sites, not a Tapestri recipe."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
import illumina as il
from chemdraw import (Construct, Scene, Segment, feature, complement_segments, revcomp,
                     panel, protocol_graph_panels)
from protocol_graph import ProtocolGraph, MaterialKind as K, InputRole as R, StepKind as S, SequencingPlan
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

def protocol_graph():
    g=ProtocolGraph('DAb-seq — patient DNA and antibody libraries',authority=SOURCE)
    cells=g.input('cells','Patient mononuclear cells prepared for staining',K.BIOLOGICAL,primary=True)
    tags=g.input('antibody_conjugates','Documented antibody-DNA panel; Table2 pattern shown; 5′azide used before click conjugation',K.OLIGO,state=tag())
    stained=g.product('stained','Antibody-DNA-stained cells',K.MOLECULAR)
    g.react('stain','Block and stain cells on ice',S.BINDING,
        ((cells,R.TEMPLATE),(tags,R.STRUCTURAL)),(stained,))
    washed=_cleanup(g,'staining_washes',stained,'Five PBS/5%FBS cell washes','Washed stained cells')
    protease=g.input('protease','Tapestri AMLV1 lysis buffer with protease/mild detergent',K.CATALYST)
    lysate=g.product('lysate','Genomic DNA and tag-containing cell lysate in droplets',K.MOLECULAR)
    g.react('lysis','Digest cellular proteins:50°C1h',S.DEPROTEINIZATION,
        ((washed,R.TEMPLATE),(protease,R.CATALYST)),(lysate,))
    heat=g.input('heat','80°C10min',K.CONDITION)
    inactivated=g.product('inactivated','Lysate with protease heat-inactivated',K.MOLECULAR)
    g.react('heat_inactivate','Heat-inactivate protease',S.THERMAL_TREATMENT,
        ((lysate,R.TEMPLATE),(heat,R.CONDITION)),(inactivated,))
    beads=g.input('barcoded_beads','AMLV1 barcoded hydrogel primers including RUNX1_4 forward site; full bead design unassigned',K.OLIGO)
    uv=g.input('uv','8min UV under Blak-Ray XX15L',K.CONDITION)
    forward=g.product('released_forward','Released cell-barcoded forward primer panel',K.OLIGO)
    g.react('photorelease','Release barcoded PCR primers in droplets',S.CHEMICAL_CONVERSION,
        ((beads,R.TEMPLATE),(uv,R.CONDITION)),(forward,))
    reverse=g.input('reverse_panel','Table1 reverse pool with RUNX1_4 reverse OMITTED',K.OLIGO)
    polymerase=g.input('polymerase','Mission Bio barcoding polymerase mixture',K.CATALYST)
    mixed=g.product('mixed_amplicons','Jointly cell-barcoded DNA amplicons and antibody products',K.MOLECULAR)
    g.react('joint_pcr','20-cycle droplet DNA amplification and antibody capture',S.PCR,
        ((inactivated,R.TEMPLATE),(forward,R.STRUCTURAL),(reverse,R.STRUCTURAL),(polymerase,R.CATALYST)),(mixed,),
        geometry=PCRGeometry.unavailable(runx_binding(),reason='RUNX1_4 tag-site binding is source-derived, but full barcoded primers/genomic targets and whole-PCR geometry remain unassigned.'))
    recovered=_cleanup(g,'emulsion_recovery',mixed,'Mission Bio Extraction Agent; aqueous transfer','Recovered barcoded DNA/tag material')
    selection=g.input('size_beads','0.75× AMPure XP; DNA beads washed twice with80%EtOH',K.SELECTION)
    dna=g.product('dna_fraction','DNA panel amplicons retained on beads',K.MOLECULAR)
    protein=g.product('protein_fraction','Short antibody-tag products retained in supernatant',K.MOLECULAR)
    g.react('size_split','Retain both size fractions',S.SPLIT,
        ((recovered,R.TEMPLATE),(selection,R.SELECTION)),(dna,protein))
    denature=g.input('denature','95°C5min then snap-cool on ice with capture oligo',K.CONDITION)
    single=g.product('denatured_tags','Separated antibody-product strands available for probe annealing',K.MOLECULAR)
    g.react('denature_tags','Denature tag products before probe capture',S.THERMAL_TREATMENT,
        ((protein,R.TEMPLATE),(denature,R.CONDITION)),(single,))
    probe_scene=Scene();probe_scene.strand('capture',[Segment('capture probe',TAG_PREFIX)],label='source capture oligo',mod5='biotin',mod3='C6 spacer')
    probe=g.input('capture_probe','Source /5Biosg/.../3C6/ capture probe',K.OLIGO,state=probe_scene)
    hybrid=g.product('capture_hybrid','Biotin probe hybridized to complementary tag copies',K.MOLECULAR,state=capture_binding())
    g.react('hybridize','Hybridize source capture probe to copied tag strand',S.HYBRIDIZATION,
        ((single,R.TEMPLATE),(probe,R.STRUCTURAL)),(hybrid,))
    streptavidin=g.input('streptavidin','Streptavidin beads;15min room temperature, two PBS washes',K.SELECTION)
    captured=g.product('captured_tags','Bead-bound purified antibody products',K.MOLECULAR)
    rejected=g.product('capture_rejected','Unbound antibody-fraction material and washes',K.WASTE)
    g.react('capture','Streptavidin bead selection of biotin-probe complexes',S.AFFINITY_CAPTURE,
        ((hybrid,R.TEMPLATE),(streptavidin,R.SELECTION)),(captured,rejected))
    g.discard(rejected,'Removed by magnetic separation and PBS washes')
    p5=g.input('p5','Documented custom P5 primer pool; source index1 representative',K.OLIGO,state=Construct([Segment('P5',il.P5),Segment('P5 Index1 sample index','TAGATCGC',feature=feature('dab.sample.i5','sample_index','fixed')),Segment('custom P5 arm',P5_ARM)]))
    dna_p7=g.input('dna_p7','Nextera P7 N7XX primer pool; sample index assignment unassigned',K.OLIGO)
    ab_p7=g.input('antibody_p7','Documented antibody P7 pool; source index1 representative with terminal PS bonds',K.OLIGO,
        state=Construct([Segment('P7',il.P7),Segment('AbP7 Index1 sample index','TACTACGC',feature=feature('dab.sample.antibody.i7','sample_index','fixed')),Segment('AbP7 arm',AB_P7_INDEX1[len(il.P7)+8:],note=AB_P7_INDEX1_WRITTEN)]))
    library_mix=g.input('library_mix','Mission Bio Library Mix',K.CATALYST)
    dna_product=g.product('dna_pcr_product','Amplified indexed DNA panel material',K.MOLECULAR)
    protein_product=g.product('protein_pcr_product','Amplified indexed antibody reporter material',K.MOLECULAR)
    for key,template,p7,out,cycles in (
        ('dna_library_pcr',dna,dna_p7,dna_product,10),
        ('antibody_library_pcr',captured,ab_p7,protein_product,20)):
        g.react(key,'Library PCR:'+str(cycles)+'cycles',S.PCR,
            ((template,R.TEMPLATE),(p5,R.STRUCTURAL),(p7,R.STRUCTURAL),(library_mix,R.CATALYST)),(out,),
            geometry=PCRGeometry.unavailable(runx_binding(),reason='Library primers are published, but complete intermediate strand/barcoded-primer geometry is unassigned in this pass.'))
    dna_lib=_cleanup(g,'dna_final_cleanup',dna_product,'0.7× AMPure XP','Patient genomic DNA library',K.LIBRARY)
    ab_lib=_cleanup(g,'antibody_final_cleanup',protein_product,'0.7× AMPure XP','Patient antibody reporter library',K.LIBRARY)
    roles=('Read 1','Read 2','Index 1 (i7)','Index 2 (i5)')
    g.sequence(dna_lib,SequencingPlan(acquired_roles=roles,read_lengths={'Read 1':150,'Read 2':150},
        limitation='Patient branch: dual-index NovaSeq6000 SP2×150; custom Read1 primer specified. Full strand/primer-coordinate map and index-cycle schedule remain unassigned.'))
    g.sequence(ab_lib,SequencingPlan(acquired_roles=roles,read_lengths={'Read 1':38,'Read 2':39},
        limitation='Patient branch: dual-index NextSeq550 Read1=38/Read2=39; source parses cell barcode inRead1 and antibody8bp plus splitUMI10bp inRead2. Complete terminal primer geometry/cycle spans remain unassigned.'))
    g.declare_reviewed_scope(('cells','antibody_conjugates','stained','stain',
        'staining_washes_reagent','staining_washes_kept','staining_washes_waste','staining_washes',
        'protease','lysate','lysis','heat','inactivated','heat_inactivate','barcoded_beads','uv',
        'released_forward','photorelease','reverse_panel','polymerase','mixed_amplicons','joint_pcr',
        'emulsion_recovery_reagent','emulsion_recovery_kept','emulsion_recovery_waste','emulsion_recovery',
        'size_beads','dna_fraction','protein_fraction','size_split','denature','denatured_tags','denature_tags',
        'capture_probe','capture_hybrid','hybridize','streptavidin','captured_tags','capture_rejected','capture',
        'p5','dna_p7','antibody_p7','library_mix','dna_pcr_product','protein_pcr_product','dna_library_pcr',
        'antibody_library_pcr','dna_final_cleanup_reagent','dna_final_cleanup_kept','dna_final_cleanup_waste',
        'dna_final_cleanup','antibody_final_cleanup_reagent','antibody_final_cleanup_kept',
        'antibody_final_cleanup_waste','antibody_final_cleanup'),SOURCE)
    return g

def render_page():
    return '\n'.join([head('DAb-seq — bounded DNA and protein branches'),'<div class="wrap"><h1>DAb-seq</h1>',
        '<p class="research-notes"><a href="01_dab_seq.html">Research notes</a></p>',
        info('Published patient branch. Complete Tapestri barcoded primers and final library strands remain unassigned. Exact source tag sites are checked below; this does not validate genomic amplification, recovery or genotype fidelity.'),
        protocol_graph_panels(protocol_graph(),cls='long'),
        panel(runx_binding().rows(),caption='Source RUNX1_4 forward binding site on the antibody tag; full5prime barcoded extension is not drawn.'),
        panel(capture_binding().rows(),caption='Biotin/C6 probe binds the copied complementary tag strand after denaturation.'),'</div>'])

def selftest():
    from checks import Check
    c=Check();g=protocol_graph()
    c('Both DNA and antibody libraries retained',len(g.sequencing),2)
    c('Tag capture selects complementary copy, not input polarity',revcomp(RUNX4),TAG_SITE)
    c('C6-blocked capture oligo is not a library-PCR primer',any(p.material=='capture_probe' and p.role is R.STRUCTURAL
        for key in ('dna_library_pcr','antibody_library_pcr') for p in g.steps[key].inputs),False)
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
    c.report(exit_on_fail=True)
