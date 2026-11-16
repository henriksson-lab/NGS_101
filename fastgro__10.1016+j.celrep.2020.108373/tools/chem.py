"""Private source-qualified partial molecular preparation; no full chemistry/absence audit."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Scene,Segment,feature,MolecularState,Workflow,workflow_panel,complement_segments
from protocol_model import role
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,SourceItem,SourceCoverage,SequencingPlan,Terminal,Disposition
class B:
    def __init__(self,title,authority,label):
        self.g=ProtocolGraph(title,authority=authority,representation='published_instance');self.x=self.g.input('source',label,K.BIOLOGICAL,primary=True)
    def stock(self,key,label,kind=K.MOLECULAR,role=R.STRUCTURAL,**kw):return(self.g.input(key,label,kind,**kw),role)
    def cond(self,key,label):return self.stock(key,label,K.CONDITION,R.CONDITION)
    def op(self,key,label,kind=S.MAGIC,inputs=(),product='Source preparation; whole occupied molecules, conversion and side products unassembled',chemical_state=None,**kw):
        y=self.g.product(key+'_product',product,K.MOLECULAR,chemical_state=chemical_state)
        if kind==S.MAGIC:kw.setdefault('reason','Positive source boundary; detailed occupied molecules, molecular products and conversion unassembled. Not source-absence evidence.')
        if kind==S.GENERIC:kw['changes_molecule']=False
        self.g.react(key,label,kind,((self.x,R.TEMPLATE),)+tuple(inputs),(y,),**kw);self.x=y;return y
    def clean(self,key,label):
        y=self.g.product(key+'_kept','Actual retained source workup fraction; carryover and yield unassigned',K.MOLECULAR);w=self.g.product(key+'_waste','Excluded source workup fraction',K.WASTE)
        self.g.react(key,label,S.CLEANUP,((self.x,R.TEMPLATE),self.cond(key+'_condition',label)),(y,w));self.g.discard(w,'Source workup excludes this fraction; no perfect purification claim.');self.x=y;return y
    def aliquot(self,key,label,qc=False):
        y=self.g.product(key+'_selected','Actual selected allocation: '+label,K.QC if qc else K.MOLECULAR);r=self.g.product(key+'_remainder','Untransferred stock; amount and later fate unassigned',K.MOLECULAR)
        self.g.react(key,label,S.SPLIT,((self.x,R.TEMPLATE),),(y,r))
        if qc:self.g.mark_terminal(Terminal(y.id,Disposition.QC,label));self.x=r
        else:self.g.mark_terminal(Terminal(r.id,Disposition.UNRESOLVED,'Untransferred stock retained as unresolved, no invented discard.'));self.x=y
        return self.x
    def finish(self,label,limitation):
        y=self.g.product('library',label,K.LIBRARY);self.g.react('load_library','Selected loading allocation, whole product unassembled',S.GENERIC,((self.x,R.TEMPLATE),),(y,),changes_molecule=False);self.g.sequence(y,SequencingPlan(limitation=limitation));return y

def oligo(seq,name):
    sc=Scene();sc.strand('DNA',[Segment(name,seq)]);return sc

AUTHORITY='Published Cell Reports33,108373(10Nov2020), actual publisher-layout PDF cached as NIHMS1645786-supplement-2, MethodsPDFp16/e3(lines2218–2240). Detailed cited protocol export is separate/version-qualified; no complete equivalence asserted.'
def protocol_graph():
    b=B('fastGRO selected standard HPDP-biotin preparation',AUTHORITY,'Source5–20million cells; upstream biological culture outside selected trace');g=b.g
    g.add_source_item(SourceItem('protocol_differences','Published MethodsPDFp16: freezing pH8.3/2x10^7nuclei per100ul; biotinylation20mMTris/2mMEDTA/25C. Exact cited DOI protocol export ID32136: pH8/2x10^7 per10ul(step8), later100ulnuclei+100ul2xNRO(step13); added recipe implies10mMTris/1mMEDTA and24C.','Selected model uses published Methods values. Actual executed values/final-volume conventions and complete historical revision equivalence unverified. Protocol58 self-refers to57–58/intended200ulDTT eluate; exact second addition not selected or silently reconstructed.'))
    b.clean('pbs','Source two ice-cold PBS washes')
    b.op('swell','Source swelling buffer5min on ice',inputs=(b.cond('swell_buffer','10mMTris7.5/2mMMg/3mMCa/2UmlSUPERasein'),))
    b.clean('glycerol_wash','Source swelling buffer plus10percent glycerol wash')
    b.op('lyse','Source nuclei isolation in1percent Igepal lysis buffer',S.EXTRACTION,inputs=(b.stock('lysis_supplies','Actual source lysis-buffer inventory including protective SUPERasein; packaging unassigned'),b.cond('lysis_conditions','10mMTris7.5/2mMMg/3mMCa/10percentglycerol/1percentIgepal/2UmlSUPERasein')))
    b.clean('nuclei_wash','Source nuclei washed twice in lysis buffer')
    b.op('freeze','Source dry-ice freezing/storage-80C',S.GENERIC,inputs=(b.cond('freeze_buffer','40percentglycerol/5mMMg/printed0.1mM0.5MEDTA/50mMTris8.3; published2x10^7nuclei per100ul; EDTA notation not silently repaired'),),product='Source frozen nuclei with resident engaged RNA-polymerase/RNA/DNA complexes; full molecular state unassembled')
    b.op('thaw','Source thaw nuclei on ice',S.GENERIC)
    # Copresent catalyst is already in the biological nuclei; no fresh purified polymerase addition.
    nuclei=g.product('nuclei_template','Same nuclei DNA/RNA/engaged-complex substrate component',K.MOLECULAR)
    resident=g.product('resident_polymerase','Same nuclei engaged resident transcription polymerase component, not added enzyme',K.CATALYST)
    g.react('resident_components','Inventory decomposition of nuclei, not physical separation or new enzyme addition',S.GENERIC,((b.x,R.TEMPLATE),),(nuclei,resident),changes_molecule=False);b.x=nuclei
    mix=g.input('nro_stock','Actual added pre-warmed NRO reaction buffer; volume EQUAL to nuclei volume, published listed concentrations describe added buffer',K.MOLECULAR)
    nt=g.product('nro_nucleotides','Same NRO buffer500uM ATP/500uM GTP/500uM4-thio-UTP/2uM CTP; final mixed concentrations not called these values',K.MOLECULAR,chemical_state={'supplied_nucleotides':['ATP','GTP','4-thio-UTP','CTP'],'label':'thiolated ribonucleotide, not DNA dNTP/4sU chemical conversion'})
    protect=g.product('nro_superase','Same NRO buffer SUPERasein200Uml protective protein component',K.MOLECULAR)
    cond=g.product('nro_condition','Same buffer10mMTris8/5mMMg/300mMKCl/1mMDTT/1percentSarkosyl; equal-volume addition',K.CONDITION)
    g.react('nro_components','Components of one actual source reaction-buffer addition, not separate fresh doses',S.GENERIC,((mix,R.TEMPLATE),),(nt,protect,cond),changes_molecule=False)
    b.op('runon','Source nuclear run-on7min30C with equal-volume buffer',S.RUN_ON,((resident,R.CATALYST),(nt,R.STRUCTURAL),(protect,R.STRUCTURAL),(cond,R.CONDITION)),product='Actual nuclear run-on preparation with newly transcribed4sUTP-containing RNA classes; exact sequences/positions/label occupancy unassembled')
    b.op('trizol','Source TRIzolLS extraction',S.EXTRACTION,(b.stock('trizol_inventory','Actual TRIzolLS extraction supplies; detailed phases and doses unassembled in selected published model'),b.cond('trizol_conditions','Source TRIzolLS extraction conditions; exact phases/doses unassembled in selected published model')))
    b.clean('nro_precip','Published source ethanol precipitation and water resuspension; exact detailed-protocol isopropanol workup not substituted')
    b.aliquot('rna_qubit','Actual Qubit RNAHS analytical allocation; volume unassigned',qc=True)
    b.aliquot('rna_transfer','Published selected up-to150ug RNA to new tube; untransferred extracted stock retained')
    b.op('spike','Source5–10percent prepared labelled Drosophila RNA spike; selected route has NO nuclei spike',S.GENERIC,(b.stock('spike_rna','Actual prepared labelled DrosophilaS9 RNA spike stock: published50mM4sU5minRT/Directzol/2ug frozen aliquots; amount, percent basis and number of aliquots used unassigned',K.MOLECULAR,R.SPIKE_IN),))
    b.aliquot('unfragmented_qc','Actual unfragmented TapeStation reference allocation; selected volume unassigned',qc=True)
    b.op('sonicate','Source BioruptorUCD200 high1–5cycles30sON/30sOFF',S.FRAGMENTATION,(b.cond('sonic_settings','Published1–5cycles; exact cited protocol one-cycle detail not forced for every sample'),))
    b.aliquot('fragmented_qc','Actual fragmented RNA TapeStation ScreenTape QC allocation; volume unassigned',qc=True)
    b.op('biotinylate','Published HPDP-biotin reaction2h dark25C800rpm',S.CHEMICAL_CONVERSION,(b.stock('hpdp','Actual EZ-link HPDP-biotin200ugml molecular contributor; thiol-linked biotin product class source-positive, full atom map/occupancy unassembled'),b.cond('bio_buffer','Published20mMTris7.5/2mMEDTA8/40percentDMF/25C2h800rpm')),product='Source biotinylated thiol-labelled RNA preparation; conversion specificity/yield and whole structures unassembled',chemical_state={'source_relation':'HPDP/pyridyldithiol-biotin forms reversible disulfide bridge between biotin and thiolated uracil base','linkage':'RNA-S-S-biotin class; full atom map/occupancy unassembled','source':'Published Results500–503; not a biotin-streptavidin covalent bond'})
    b.clean('bio_precip','Published ethanol precipitation and water resuspension after biotinylation')
    beads=g.input('m280_stock','Actual100ul M280 streptavidin Dynabeads/sample',K.MOLECULAR)
    washed=g.product('washed_m280','Same source washed/resuspended100ul M280 capture preparation',K.SELECTION);wash=g.product('m280_wash_removed','Excluded source prewash fraction',K.WASTE)
    g.react('prepare_beads','Source two washes with2volumes buffer, resuspend1volume',S.CLEANUP,((beads,R.TEMPLATE),b.cond('m280_wash_buffer','100mMTris7.5/10mMEDTA8/1MNaCl/.1percentTween20')),(washed,wash));g.discard(wash,'Source prewash excluded.')
    bound=g.product('bound_rna','Actual bead-associated biotin-RNA preparation',K.MOLECULAR);unbound=g.product('capture_unbound','Actual removed capture supernatant',K.WASTE)
    g.react('capture','Source M280 capture15min4Crotation',S.AFFINITY_CAPTURE,((b.x,R.TEMPLATE),(washed,R.SELECTION)),(bound,unbound));g.discard(unbound,'Capture removes this fraction, not a claim of complete depletion.');b.x=bound
    b.clean('hot_wash','Source three65C wash-buffer washes')
    b.clean('rt_wash','Source three room-temperature wash-buffer washes')
    released=g.product('released_rna','Source DTT-released thiol-labelled RNA eluate; biotin adduct reductively removed',K.MOLECULAR,chemical_state={'source_relation':'Reducing DTT cleaves source reversible disulfide biotin-adduct linkage; thiol-labelled RNA released','not_claimed':'No cleavage of biotin-streptavidin bond or conversion to unmodified uridine; full atom map/occupancy unassembled','source':'Published Results500–503/Methods2236–2238'});residual=g.product('residual_beads','Residual bead-associated material after release',K.MOLECULAR)
    g.react('dtt_release','Source100mM DTT release; number/volume of releases not assigned from self-referential protocol58',S.MAGIC,((b.x,R.TEMPLATE),b.stock('dtt','Actual100mM DTT reducing molecular contributor')),(released,residual),reason='Source DTT release chemistry known at class level; whole disulfide products, exact second-elution schedule and quantitative recovery unassembled.');g.mark_terminal(Terminal(residual.id,Disposition.UNRESOLVED,'Residual bead fraction fate/retained RNA unassigned.'));b.x=released
    b.op('column_bind','Source Zymo RNA Clean and Purification column binding',inputs=(b.stock('column_inventory','Actual source column-binding supplies; final buffer formulation unassembled'),))
    b.op('column_dnase','Source in-column DNase removes genomic-DNA traces',S.DIGESTION,(b.stock('dnase','Actual source in-column DNase enzyme',K.CATALYST,R.CATALYST),b.cond('dnase_buffer','Source DNase reaction buffer; published exact volume/conditions unassigned')))
    b.clean('column_wash_elute','Source column purification/water elution after DNase')
    b.aliquot('enriched_qubit','Actual Qubit RNAHS analytical allocation',qc=True)
    b.aliquot('library_rna','Actual enriched RNA allocation for NEBNext UltraII Directional library; amount unassigned')
    kit=g.input('library_supplies','Aggregate source historical NEBNext UltraII Directional RNA-library supplies; packaging/stock count/internal order unassigned',K.MOLECULAR)
    ol=g.product('kit_oligos','Actual kit oligo inventory contributing to source barcoded library; whole bodies/index values and functional partition unassembled',K.OLIGO,chemical_state={'sample_index':feature('fastgro.sample','sample_index','unknown',note='Published barcoded library; whole oligo, index width/value/coordinate and layout unassembled')})
    other=g.product('kit_other','Other actual kit molecular supplies; exact source stock identities and internal order unassembled',K.MOLECULAR)
    g.react('library_components','Inventory bookkeeping, not claim of one prepared kit mix or extra fresh stock',S.GENERIC,((kit,R.TEMPLATE),),(ol,other),changes_molecule=False)
    b.op('library_finish','Source NEBNext UltraII Directional barcoded RNA-library preparation; internal molecular recipe unassembled',inputs=((ol,R.STRUCTURAL),(other,R.STRUCTURAL)))
    b.aliquot('loading','Actual NextSeq500 loading allocation; volume/mass unassigned')
    b.finish('Source selected fastGRO barcoded sequencing library; whole occupied construct unassembled','Published NextSeq500 endpoint; read lengths, primer sites, index values/whole oligos and pooled-sample physical schedule unassembled. No implicit modern kit recipe, RT/PCR cycles or read mask.')
    declare_page_scope(g)
    g.cover_source_item(SourceCoverage('protocol_differences','modelled','freeze'))
    return g

def selftest():
    from checks import Check,Source,have
    from storage import sources_root
    c=Check();g=protocol_graph()
    c('NRO consumes resident polymerase produced from source nuclei',next(i.material for i in g.steps['runon'].inputs if i.role==R.CATALYST),'resident_polymerase')
    c('Resident enzyme not external addition','resident_polymerase' in g.steps['resident_components'].outputs,True)
    c('Labelled nucleotide actual molecular contributor',g.materials['nro_nucleotides'].kind,K.MOLECULAR)
    c('Added NRO buffer not fictitious individual fresh doses',g.steps['nro_components'].inputs[0].material,'nro_stock')
    c('Spike is RNA not nuclei',next(i.material for i in g.steps['spike'].inputs if i.role==R.SPIKE_IN),'spike_rna')
    c('Capture uses physical prepared beads',next(i.material for i in g.steps['capture'].inputs if i.role==R.SELECTION),'washed_m280')
    c('Capture supernatant excluded',g.terminals['capture_unbound'].disposition,Disposition.DISCARD)
    c('DTT released RNA proceeds to column',g.steps['column_bind'].inputs[0].material,'released_rna')
    c('DTT acts on labelledRNA-biotin disulfide class','disulfide' in g.materials['released_rna'].chemical_state['source_relation'],True)
    c('Residual beads not automatic discard',g.terminals['residual_beads'].disposition,Disposition.UNRESOLVED)
    c('DNase after capture/release',g.steps['column_dnase'].inputs[0].material,'column_bind_product')
    c('Library oligos explicitly unassembled',g.materials['kit_oligos'].state,None)
    c('Library kit oligo component consumed',g.steps['library_finish'].inputs[1].material,'kit_oligos')
    c('Both fragment comparison allocations QC',tuple(g.terminals[k].disposition for k in ('unfragmented_qc_selected','fragmented_qc_selected')),(Disposition.QC,Disposition.QC))
    c('Final loading stock remainder retained','loading_remainder' in g.terminals,True)
    c('No automatic LI chemistry stock',any('mts' in m.id for m in g.materials.values()),False)
    p=sources_root()/Path(__file__).resolve().parents[1].name/'j.celrep.2020.108373_supp_NIHMS1645786-supplement-2.pdf.txt'
    if have(c,Source(p,'https://doi.org/10.1016/j.celrep.2020.108373')):
        t=''.join(p.read_text().split());c('Published identity not onlyAM','November10,2020' in t and 'CellReports33,108373' in t,True);c('Actual final standard HPDP route','200mg/mlEZ-linkHPDPBiotin' in t,True)
    c('Typed publication gate',g.publication_problems(),[]);c.report()

# Curator-written scope from published Methods e3/Results and the qualified
# detailed-protocol comparison. Never generated from graph inputs/steps.
def declare_page_scope(g):
    g.declare_reviewed_scope((
        'source','pbs_condition','pbs','swell_buffer','swell',
        'glycerol_wash_condition','glycerol_wash','lysis_supplies',
        'lysis_conditions','lyse','nuclei_wash_condition','nuclei_wash',
        'freeze_buffer','freeze','thaw','resident_components',
        'nro_stock','nro_components','runon','trizol_inventory',
        'trizol_conditions','trizol','nro_precip_condition','nro_precip',
        'rna_qubit','rna_transfer','spike_rna','spike','unfragmented_qc',
        'sonic_settings','sonicate','fragmented_qc','hpdp','bio_buffer',
        'biotinylate','bio_precip_condition','bio_precip','m280_stock',
        'm280_wash_buffer','prepare_beads','capture','hot_wash_condition',
        'hot_wash','rt_wash_condition','rt_wash','dtt','dtt_release',
        'column_inventory','column_bind','dnase','dnase_buffer','column_dnase',
        'column_wash_elute_condition','column_wash_elute','enriched_qubit',
        'library_rna','library_supplies','library_components','library_finish',
        'loading','load_library','library',
    ), AUTHORITY)

def molecular_class(kind):
    sc=Scene()
    if kind in {'engaged','runon'}:
        dna=[role('native transcribed DNA interval',inferred=True)]
        other=complement_segments(dna)
        sc.strand('dna',dna,label='INFERRED genomic DNA')
        sc.anneal('other',other,to='dna',pair=(other[0].name,dna[0].name),unpaired=(other[0].name,),label='DNA complement class')
        rna=[role('pre-existing nascent RNA',inferred=True)]
        if kind=='runon':rna.append(role('run-on RNA extension containing 4-thiouridine',inferred=True))
        sc.strand('rna',rna,label='INFERRED extended RNA'if kind=='runon'else'INFERRED nascent RNA')
        return sc.rows()
    if kind=='library':
        dna=[role('historical directional-kit end / sample index',molecular_feature=feature('fastgro.sample','sample_index','unknown'),inferred=True),
             role('cDNA copied from enriched run-on RNA',inferred=True),role('other historical library end',inferred=True)]
        other=complement_segments(dna);sc.strand('dna',dna,label='INFERRED cDNA library')
        sc.anneal('other',other,to='dna',pair=(other[-1].name,dna[0].name),unpaired=tuple(x.name for x in other),label='DNA complement class')
        return sc.rows()
    rna=[role('run-on RNA fragment containing incorporated 4-thiouridine',inferred=True)]
    state={'isolated':'Extracted labelled RNA; no genomic DNA–RNA covalent join',
           'fragmented':'Sonicated RNA population; positions and end chemistry unassembled',
           'biotin':'RNA–S–S–biotin adduct class; not a phosphodiester backbone modification',
           'bound':'Biotin-associated RNA retained on streptavidin support',
           'released':'DTT-released thiol-labelled RNA; biotin adduct reductively cleaved',
           'clean':'Purified enriched RNA after column DNase; not perfectly pure'}.get(kind,kind)
    sc.strand('rna',rna,label='INFERRED enriched RNA')
    if kind in {'isolated','fragmented','biotin','bound','released','clean'}:
        sc.strand('spike',[role('separate prepared labelled Drosophila RNA spike',inferred=True)],label='Separate RNA spike')
    sc.note('rna',state)
    if kind in {'biotin','bound'}:
        sc.note('rna','HPDP provides the reversible disulfide-linked biotin affinity handle; exact modified sites and atom map remain unassembled.')
    if kind=='bound':sc.note('rna','Streptavidin–biotin is an affinity interaction. It is distinct from the covalent RNA-associated disulfide that DTT reduces.')
    return sc.rows()


def molecular_walkthrough():
    w=Workflow(MolecularState('INFERRED — engaged transcription complexes in source cells',tuple(molecular_class('engaged'))),track_graph=False)
    def step(title,kind,note):w.react(title,molecular_class(kind),name=title,note=note)
    step('1 · Isolate nuclei while preserving engaged transcription complexes','engaged',
         'Start with 5–20 million cells: cold PBS washes, swelling buffer 5 min on ice, glycerol wash, then 1% Igepal lysis and two nuclei washes. Freeze nuclei on dry ice and store at −80°C; thaw on ice for run-on. The source freezing-buffer pH, concentration and EDTA notation differences remain listed below. No purified RNA polymerase is added.')
    step('2 · Nuclear run-on incorporates 4-thio-UTP into new RNA','runon',
         'Mix nuclei with an equal volume of prewarmed run-on buffer and incubate 7 min at 30°C. Resident engaged polymerases extend existing RNA. The added buffer lists ATP, GTP and 4-thio-UTP at 500 µM each and CTP at 2 µM; these are not relabeled as final mixed concentrations. Actual genomic bases, extension length and label occupancy remain unassembled.')
    step('3 · Extract RNA; allocate up to 150 µg and add prepared RNA spike','isolated',
         'TRIzol LS extraction, source ethanol precipitation and water resuspension precede Qubit quantification. Selected route adds 5–10% prepared labelled Drosophila RNA after run-on rather than spiking nuclei before run-on. Percentage basis and exact spike dose remain unassigned. Spike molecules are separate from sample RNA.')
    step('4 · Sonicate RNA and compare fragment QC with an unfragmented aliquot','fragmented',
         'Bioruptor UCD-200 high setting, 1–5 cycles of 30 s on / 30 s off. The run-on label is retained on label-bearing fragments. TapeStation uses separate analytical allocations; its QC material is not returned to the production stock. Fragment boundaries and chemical ends are not invented.')
    step('5 · Couple HPDP-biotin to incorporated thiols','biotin',
         'Published standard route: 200 µg/ml HPDP-biotin, 20 mM Tris pH7.5, 2 mM EDTA and 40% DMF; 2 h at 25°C in the dark, 800 rpm. A reversible disulfide links the RNA-associated thiol and biotin-containing adduct. Ethanol precipitation follows. This is not the separate low-input MTS-biotin recipe.')
    step('6 · Streptavidin capture, hot washes and room-temperature washes','bound',
         'Prepare 100 µl M280 streptavidin beads per sample with two washes; capture 15 min at 4°C with rotation. Wash three times at 65°C and three times at room temperature in the stated Tris/EDTA/NaCl/Tween buffer. Nonselected material leaves with supernatants. Label-specific capture does not prove quantitative recovery or perfect removal of background.')
    step('7 · DTT reduces the disulfide and releases labelled RNA','released',
         'Use source 100 mM DTT. Reduction removes the biotin-associated disulfide adduct and releases RNA from the affinity support; it does not cleave a covalent biotin–streptavidin bond or establish conversion to ordinary unmodified uridine. Residual bead material remains a distinct outcome. The cited export’s self-referential second-release instruction is not silently reconstructed.')
    step('8 · Column-bind, digest DNA on column, wash and elute RNA','clean',
         'The standard route uses Zymo RNA Clean and Purification with in-column DNase before final wash/elution. Quantify recovered RNA by Qubit and allocate it to library preparation. This differs from assigning the low-input route’s earlier TURBO DNase step to the standard workflow.')
    step('9 · Make a barcoded directional cDNA library; sequence on NextSeq500','library',
         'The paper names NEBNext Ultra II Directional RNA preparation. Its historical whole oligos, internal fragmentation/RT/strand-marking schedule, index width/values, primer landing and run lengths remain unassembled. The diagram marks those ends as inferred roles and does not fill them with a modern kit recipe. Loading and untransferred library stock remain separate physical fates.')
    return w

PAGE_FILE = 'fastgro.html'

def render_page():
    from page import head,info,caveat,table
    from protocol_views import assemble_page,WorkflowSection,PlanLimitationSection
    from chemistry_validation import check_graph,Status
    g=protocol_graph();g.require_publishable()
    overview=(
        '<div class="wrap" data-protocol-complete="true"><h1>fastGRO</h1>'
        '<p>Label nascent RNA during nuclear run-on, then enrich the thiolated RNA through a reversible biotin linkage. '
        'The selected standard preparation uses HPDP-biotin and a prepared labelled RNA spike-in.</p>'
        +info('Source: <a href="https://doi.org/10.1016/j.celrep.2020.108373">Barbieri et al., Cell Reports 33, 108373 (2020)</a>, '
              'published Methods page e3 and Results. The cited <a href="https://doi.org/10.17504/protocols.io.bbmgik3w">detailed protocol</a> '
              'was inspected as export ID32136 and remains version-qualified. The low-input MTS-biotin route is separate.')
        +caveat('This page declares the run-on, label/capture/release chemistry and source sequencing endpoint. '
                'The historical NEBNext Ultra II Directional kit continuation is explicit but molecularly unassembled: '
                'whole kit oligos, primer landing, index assignments, internal reaction schedule and exact acquired read lengths are not supplied by this model. '
                'These are review boundaries, not confirmed publication omissions.')
        +'<h2>Label, capture and release</h2>'
        +table(('Stage','Source operation','Chemical meaning'),(
            ('Run-on','Thawed nuclei + an equal volume of prewarmed NRO buffer; 7 min at 30°C.','Resident engaged polymerase incorporates 4-thio-UTP into RNA; no fresh purified polymerase is added.'),
            ('Recover and fragment','TRIzol LS extraction, precipitation, prepared labelled RNA spike and sonication.','RNA spike and sample RNA remain molecular inputs; fragment/Qubit allocations are separate QC fates.'),
            ('Biotinylate','HPDP-biotin, 2 h in the dark at 25°C.','A reversible disulfide links the thiolated uracil class to the biotin-containing reagent. Exact atom map and label occupancy remain unassembled.'),
            ('Capture','Prepared M280 streptavidin beads; cold binding, hot and room-temperature washes.','Affinity capture retains biotin-associated RNA, without asserting complete purity or recovery.'),
            ('Release','100 mM DTT, then column purification and in-column DNase.','Reduction cleaves the disulfide adduct. It does not cleave a covalent biotin–streptavidin bond or necessarily restore unmodified uridine.'),
            ('Library','Source barcoded directional RNA-library preparation and NextSeq500 sequencing.','Historical kit recipe, exact whole oligos and read geometry remain unassembled.'),
        ))
        +'<h2>Known run-on inputs</h2>'
        +table(('Component','Concentration in the added NRO buffer'),(
            ('ATP / GTP / 4-thio-UTP','500 µM each'),('CTP','2 µM'),('SUPERase-in','200 U/ml'),
            ('Other listed components','10 mM Tris pH8, 5 mM MgCl₂, 300 mM KCl, 1 mM DTT and 1% Sarkosyl.'),
        ))
        +info('These describe the added buffer before equal-volume mixing with nuclei. They are not labelled as final mixed reaction concentrations. '
              'Whole source-library oligo sequences are unassembled; no canonical Illumina shell is substituted.')
        +'<h2>Read setup</h2>'
        +table(('Item','Current declaration'),(
            ('Instrument','Published NextSeq500 endpoint.'),('Library','Source barcoded NEBNext Ultra II Directional RNA library.'),
            ('Cycles and indices','Exact read lengths, primer sites, index widths/values and index-read recipe remain unassigned in this model.'),
        ))
        +'<h2>Published Methods and detailed-protocol differences</h2>'
        +caveat('Published Methods use freezing-buffer pH8.3 and 2×10⁷ nuclei per 100 µl. The inspected cited protocol export uses pH8 and '
                '2×10⁷ per 10 µl in step8, then 100 µl nuclei + 100 µl 2×NRO in step13. Published biotinylation specifies 20 mM Tris, '
                '2 mM EDTA and 25°C; the export recipe implies 10 mM Tris, 1 mM EDTA and 24°C. This model selects the published Methods values. '
                'Actual executed values and full revision equivalence remain unresolved. The export’s step58 self-refers to steps57–58; '
                'the model does not invent an exact second DTT-release schedule. The printed freezing-buffer notation “0.1 mM 0.5 M EDTA” is retained without selecting an intended concentration.')
    )
    claims=check_graph(g)
    audit=('<h2>Molecular validation scope</h2>'+table(('Claim status','Computed count'),tuple((s.value,str(sum(c.status==s for c in claims)))for s in Status))
           +'<p>The complete graph below retains unresolved molecular transitions, QC allocations, discarded capture fractions and untransferred stock remainders.</p>')
    return assemble_page(g,[head('fastGRO — nuclear run-on and reversible RNA capture'),overview,
        '<h2>Step by step: extend RNA, label its thiols, capture and release</h2><p>INFERRED — these molecular-class drawings show source-supported changes and retained components. Unknown biological sequences, occupied transcription bubbles and modified-base coordinates remain bracketed; the terminal library is a role-level class, not a validated whole construct.</p>',
        workflow_panel(molecular_walkthrough(),caption='Standard HPDP-biotin route; prepared RNA spike, not nuclei spike'),audit,
        '<details><summary>Complete preparation, allocations and model warnings</summary>',WorkflowSection('Complete preparation and material fates'),'</details>',PlanLimitationSection('library'),'<p><a href="01_fastgro.html">Research note and source comparison</a></p></div>'])

if __name__=='__main__':selftest()
