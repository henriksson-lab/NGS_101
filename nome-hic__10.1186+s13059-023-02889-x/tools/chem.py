"""Protocol-owned selected chemical classes; incomplete source-qualified graph."""
from dataclasses import dataclass
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Scene, Segment, revcomp, feature
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,SourceItem,SequencingPlan,Terminal,Disposition
from pcr import PCRGeometry
from adapter_ligation import AdapterLigationGeometry

def scene(name,text):
    sc=Scene();sc.strand(name,[Segment(name,text,placeholder=True)]);return sc

def clean(g,key,x,label):
    condition=g.input(key+'_condition',label,K.CONDITION)
    kept=g.product(key+'_kept','Source recovered material; no perfect purity or quantitative partition assumed',K.MOLECULAR)
    waste=g.product(key+'_waste','Actually separated cleanup fractions with unquantified molecular carryover',K.WASTE)
    g.react(key,label,S.CLEANUP,((x,R.TEMPLATE),(condition,R.CONDITION)),(kept,waste));g.discard(waste,'Source physical separation, not perfect molecular classification.');return kept

def op(g,key,x,label,kind=S.MAGIC,*,reagent='source mixture',catalyst=False,product=None,state=None,chemical_state=None):
    stock=g.input(key+'_reagent',reagent,K.CATALYST if catalyst else K.CONDITION)
    y=g.product(key+'_product',product or label+'; all physical reaction companions retained until separation',K.MOLECULAR,state=state,chemical_state=chemical_state)
    kwargs={}
    if kind==S.MAGIC:kwargs['reason']='Source operation established; exact internal chemistry/geometry not assembled from unreviewed historical originals.'
    if kind==S.PCR:kwargs['geometry']=PCRGeometry.unavailable(Scene(),reason='Exact source template/primer endpoint geometry unassembled, no ordinary complete library asserted.')
    if kind in (S.LIGATION,S.PROXIMITY_LIGATION):kwargs['geometry']=AdapterLigationGeometry('unavailable','unavailable',(scene('unassembled','[source-known ligation; full substrate/end geometry unassembled]'),),'Unreviewed geometry; not a claim of unavailable disclosure.')
    inputs=[(x,R.TEMPLATE),(stock,R.CATALYST if catalyst else R.CONDITION)]
    if kind in {S.DIGESTION,S.END_REPAIR,S.EXTENSION,S.TAILING,S.PROXIMITY_LIGATION} and not catalyst:
        enzyme=g.input(key+'_enzyme','Explicit catalytic component of source mixture: '+reagent,K.CATALYST);inputs.append((enzyme,R.CATALYST))
    if kind==S.LIGATION:
        adapter=g.input(key+'_adapter','Source ligating molecular adapter; exact historical oligo/end inventory unassembled',K.ADAPTER,state=scene(key+'_adapter','[source adapter geometry unassembled]'));inputs.append((adapter,R.STRUCTURAL))
    if kind==S.PCR:
        for arm in ('f','r'):
            primer=g.input(key+'_primer_'+arm,'Source PCR primer family; exact oligo/endpoint unassembled',K.OLIGO,state=scene(key+'_primer_'+arm,'[source primer body/index geometry unassembled]'));inputs.append((primer,R.STRUCTURAL))
    g.react(key,label,kind,inputs,(y,),**kwargs);return y

def qc_split(g,key,x,label):
    assay=g.product(key+'_assay','Physically separate source analytical aliquot',K.QC)
    kept=g.product(key+'_retained','Unassayed retained preparation for downstream library work',K.MOLECULAR)
    g.react(key,label,S.SPLIT,((x,R.TEMPLATE),),(assay,kept));g.mark_terminal(Terminal(assay.id,Disposition.QC,'Source analytical aliquot never feeds retained library.'));return kept
SOURCE='Fu et al.2023 DOI10.1186/s13059-023-02889-x publishedPMC10018866 mainMethods153/analysis155–163; publishedAdditionalfile1 SupplementaryMethodspp2–6, p3 allocationoriginalrendered.'
ALLOCATION='PublishedSIp3 specifies300ulRNA/1200ulDNA; main153 saysone-thirdofnucleiRNA. Different allocation statements retained, volume homogeneity/nuclei-number relation unverified; selectedgraphusesSI volumes.'
@dataclass(frozen=True)
class FootprintedDNA:
    sequence:str
    endogenous_c:tuple=()
    installed_gpc_c:tuple=()
    def __post_init__(self):
        if not self.sequence or set(self.sequence)-set('ACGT'):raise ValueError('caller-supplied exact selectedDNA required')
        for positions in (self.endogenous_c,self.installed_gpc_c):
            if len(set(positions))!=len(positions) or any(type(i)!=int or not 0<=i<len(self.sequence) or self.sequence[i]!='C' for i in positions):raise ValueError('covalent methylation positions must identifyC')
        if any(i==0 or self.sequence[i-1]!='G' for i in self.installed_gpc_c):raise ValueError('selected M.CviPI target C must followG')
    def footprint(self,positions):
        # Caller supplied reacted sites, never predicted accessibility or complete specificity.
        return FootprintedDNA(self.sequence,self.endogenous_c,tuple(sorted(set(self.installed_gpc_c)|set(positions))))
    def converted(self):
        marks=set(self.endogenous_c)|set(self.installed_gpc_c)
        return ''.join('U' if b=='C' and i not in marks else b for i,b in enumerate(self.sequence))
    def complementary_dna_copy(self):return revcomp(self.converted().replace('U','T'))
    def gch_positions(self):
        # Analysis-only context predicate. It never edits the physical molecule.
        return tuple(i for i in self.installed_gpc_c if i+1<len(self.sequence) and self.sequence[i+1]!='G')

def capture(g,key,x,label,reagent,*,retained_beads=False):
    beads=g.input(key+('_binding_conditions' if retained_beads else '_beads'),reagent,K.SELECTION);bound=g.product(key+'_bound','Source bead-associated selected preparation; no perfect specificity inferred',K.MOLECULAR);other=g.product(key+'_unbound','Physically separate flowthrough with unquantified exact composition/fate',K.MOLECULAR)
    g.react(key,label,S.AFFINITY_CAPTURE,((x,R.TEMPLATE),(beads,R.SELECTION)),(bound,other));g.mark_terminal(Terminal(other.id,Disposition.UNRESOLVED,'Actual separated source companion fraction; no downstream fate or perfect partition inferred.'));return bound

def protocol_graph():
    g=ProtocolGraph('NOMe-HiC — DNA andRNA branches',authority=SOURCE,representation='published_instance')
    x=g.input('cells','SourceIMR90/GM12878 cells withDNA/RNA/protein',K.MOLECULAR,primary=True)
    x=op(g,'crosslink',x,'Source1%formaldehyde10minRT crosslink',S.CROSS_LINKING,reagent='Sourceformaldehyde')
    x=op(g,'quench',x,'Source0.2Mglycine5minRT',S.QUENCHING,reagent='Sourceglycine')
    x=clean(g,'wash_cells',x,'SourceicecoldPBS wash/recovery')
    x=op(g,'isolate_nuclei',x,'SourceHiCnucleiisolation withRNaseinhibitors/DTT ice1h',S.EXTRACTION,reagent='SourceNP40/nucleiisolationbuffer andRNaseinhibitors')
    x=clean(g,'wash_nuclei',x,'Source2500g5min4C pelletandwash')
    x=op(g,'permeabilize',x,'Source0.5%SDS62C10min thenTriton/H2O',S.PERMEABILIZATION,reagent='SourceSDS/Triton')
    x=op(g,'dpnii',x,'Source100UDpnII37Covernight withRNaseinhibitors',S.DIGESTION,reagent='SourceDpnII andbuffer3.1/RNaseinhibitors')
    x=op(g,'heat_inactivate',x,'Source65C10minDpnIIinactivation andcool',S.THERMAL_TREATMENT,reagent='Source65Cincubation thenRTcool10min')
    biotin=g.input('biotin_datp','SourceprintedBiotin-1,4-dATP (notation retained), together with dGTP/dCTP/dTTP; exact adductstructure unassembled',K.MOLECULAR)
    enzyme=g.input('klenow','SourceKlenowM0210L',K.CATALYST);y=g.product('filled','KlenowfilledDNA ends containing sourcebiotin nucleotide, otherDNA andreactioncompanions retained; fullendgeometry unassembled',K.MOLECULAR)
    g.react('fill','Sourcebiotinfill37C90min withRNaseinhibitors',S.EXTENSION,((x,R.TEMPLATE),(biotin,R.STRUCTURAL),(enzyme,R.CATALYST)),(y,))
    x=op(g,'contact_ligation',y,'SourceT4DNA ligasein situ proximalbluntendligationRT4h',S.PROXIMITY_LIGATION,reagent='SourceT4ligaseM0202L/buffer/Triton/BSA',product='Crosslinked nuclei retaining ligatedcontactDNA plusunligatedends andphysicalRNA/protein/reactioncompanions; no perfectligation asserted')
    x=clean(g,'pellet_ligated',x,'Source2500g5min4C nuclei pellet/resuspendGpCbuffer')
    donor=g.input('sam','SourceSAM methyl donor, two1.5ul32mM additions; consumed donor is distinct fromDNA',K.MOLECULAR);enzyme=g.input('mcvi','SourceM.CviPIM0227L',K.CATALYST);y=g.product('footprinted','Nuclei with exogenous selectedGpC5mC andendogenousDNAmarks, contactDNA andRNA/protein retained; SAMreactionproducts retained untilphysicalcleanup',K.MOLECULAR)
    condition=g.input('footprint_buffer','SourceGpCbuffer/sucrose/RNaseinhibitors',K.CONDITION)
    g.react('footprint','SourceM.CviPI200U7.5min37C then100U/SAM7.5min37C',S.CHEMICAL_CONVERSION,((x,R.TEMPLATE),(donor,R.STRUCTURAL),(enzyme,R.CATALYST),(condition,R.CONDITION)),(y,))
    x=op(g,'stop_footprint',y,'Source973.5ulicecoldPBS stop withRNaseinhibitors',S.QUENCHING,reagent='SourcePBS/RNaseinhibitors/DTT')
    rna=g.product('rna_nuclei','SourceSI300ul nuclei-suspension RNA branch; containsDNA/protein too',K.MOLECULAR);dna=g.product('dna_nuclei','SourceSI1200ul nuclei-suspension DNA branch; containsRNA/protein too',K.MOLECULAR)
    g.react('split_branches',ALLOCATION,S.SPLIT,((x,R.TEMPLATE),),(rna,dna))
    rna=op(g,'rna_decrosslink',rna,'SourceMagMAX protease55Covernight90C1h thenRT15min',S.DEPROTEINIZATION,reagent='SourceMagMAXA31881 protease/buffer',product='RNAbranch nucleicacids plusproteolysis/crosslinkreactionproducts retained')
    rna=capture(g,'rna_capture',rna,'SourceMyOneSilane binding/isopropanolRT10min','SourceMyOneSilane37002D/bindingbuffer/isopropanol')
    rna=clean(g,'rna_wash',rna,'SourceRNAwashbuffer/washsolution magneticwash; retainedbead-associatedpreparation')
    rna=op(g,'rna_dnase',rna,'SourceDNaseRT20min onresuspendedbeads',S.DIGESTION,reagent='SourceDNase andbuffer',product='RNAplusdigestedDNAproducts/beads andotherremnants; DNAproducts not removed byreactionalone')
    rna=capture(g,'rna_recapture',rna,'Sourcebindingbuffer/isopropanolrecaptureRT10min','Source200ulbindingbuffer/250ulisopropanol added to retained originalSilane bead suspension; no freshbeads',retained_beads=True)
    rna=clean(g,'rna_elute',rna,'Sourcemagneticwashes/airdry/H2O30ulelution')
    rna=op(g,'rna_library',rna,'Source10ngRNA SMARTerStrandedTotalRNASeqv2PicoMammalian634413',reagent='HistoricalTAKARA634413 kit, internalRT/depletion/PCR/oligogeometryunreviewed',product='SourceRNA-seqlibrarypreparation plus reactioncompanions; exact kitprocess andreadgeometry unassembled')
    rna=qc_split(g,'rna_qc',rna,'SourceRNA-libraryQC aliquot andretainedRNAderivedlibrary')
    g.sequence(rna,SequencingPlan(acquisition='unspecified',limitation='SourceRNAderivedlibrary PE150HiSeqXTen/NovaSeqS6000; actual kitlibrary/readprimers unassembled. RNAfraction distinct fromDNAfraction. '+ALLOCATION))
    dna=op(g,'dna_decrosslink',dna,'SourcePK/SDS55C30min thenNaCl68C4h',S.DEPROTEINIZATION,reagent='SourceproteinaseK/SDS/NaCl',product='DNAbranchDNA plusRNA/protein/crosslinkreactionproducts untilprecipitation')
    dna=clean(g,'dna_precipitate',dna,'SourceEtOH/NaAc minus20Covernight precipitation and80%EtOHwashes')
    dna=op(g,'sonicate',dna,'SourceCovarisM220 mean400bp sonication',S.FRAGMENTATION,reagent='SourceCovaris physicalsonication')
    dna=clean(g,'size_cleanup',dna,'Source1xAMPure sizeselection')
    dna=capture(g,'biotin_capture',dna,'SourceT1streptavidinbiotincaptureRT15min thenwash55C','SourceMyOneT1streptavidin65602')
    dna=op(g,'repair',dna,'SourceonbeadPNK/T4DNApolymerase/KlenowRT30min',S.END_REPAIR,reagent='SourceT4PNK/T4polymerase/Klenow/dNTPs',product='BeadDNA repairpopulation withsourceunligatedfragmentprocessing products retained untilwash; no completebiotinselection claim')
    dna=clean(g,'repair_wash',dna,'Sourcetwowashes55C, retainedbeadDNA')
    dna=op(g,'a_tail',dna,'SourceKlenowexo-minus/dATP37C30min',S.TAILING,reagent='SourceKlenowexo-minusM0212L/dATP')
    dna=clean(g,'tail_wash',dna,'Sourcetwowashes55C retainedbeadDNA')
    dna=op(g,'stubby_ligate',dna,'Source0.9uMcytosinemethylatedtruncatedXgeneUniversalStubbyAdapter QuickLigaseRT15min',S.LIGATION,reagent='SourceIDTmethylatedstubbyadapter/QuickLigase; exact historical oligoendsunassembled')
    dna=clean(g,'adapter_wash',dna,'Source55Cwashbuffer thenTriswashes retainedbeadDNA')
    spike=g.input('lambda','SourceSI0.1%400bpstubby-ligatedunmethylatedlambda control; exact adaptergeometryunassembled',K.MOLECULAR);mix=g.product('spiked','BeadDNA andsourceadaptedlambda mixture',K.MOLECULAR);g.react('spike','Sourcecontrolspikebeforeconversion',S.GENERIC,((dna,R.TEMPLATE),(spike,R.SPIKE_IN)),(mix,))
    dna=op(g,'bisulfite',mix,'SourceEZDNAMethylationGoldD5006 conversion',reagent='HistoricalZymoD5006; exact kitinternalsunreviewed',product='ConvertedDNAcontainingphysicalU atselectedunprotectedC, retainedselected5mC andreactioncompanions; no GCG/CCG physicaldeletion')
    dna=clean(g,'converted_recover',dna,'SourceDNAseparatedfromstreptavidinbeadsafterCTconversion; vendorworkupgeometryunassembled')
    dna=op(g,'index_pcr',dna,'SourceKAPAHiFiUracil+KK2801/truncatedTruSeqindexprimers10–15cyclesMgCl2',S.PCR,reagent='SourceKAPAHiFiUracil+ andsourceindexprimerfamiliesunassembled')
    dna=qc_split(g,'dna_qc',dna,'SourceqPCR/BioAnalyzer analyticalaliquots andretainedDNAderivedlibrary')
    g.sequence(dna,SequencingPlan(acquisition='unspecified',limitation='SourcePE150HiSeqXTen/NovaSeqS6000. Exactstubby/index/primergeometryunassembled; computationalHCG/GCHcontexts,GCG/CCGfilters,SNPimputation do not physicallyalterDNA. SourceRNAseqorPhiX20%balancing spike not part ofsameDNA construct. '+ALLOCATION))
    g.add_source_item(SourceItem('pending','OriginalNOMe-seq/MethylHiC methods andhistoricalMagMAX/SMARTerv2/IDTstubby/TruSeq/Zymoguides','Cited componentpredecessors/vendors not read within cap; fullsource/oligo/end/acquisitiongeometry unassembled. No absenceclaim.'))
    return g

def selftest():
    from checks import Check,Source,have
    from storage import sources_root
    c=Check();d=FootprintedDNA('AGCGCC',(2,)).footprint((2,4))
    c('footprinting retainsoriginalsequence',d.sequence,'AGCGCC')
    c('physicalconversion protectsselectedmarksbutnotunmarkedC',d.converted(),'AGCGCU')
    c('analysisGCHsubsetdoesnotdeleteGCGmarkedC',(d.gch_positions(),d.installed_gpc_c),((4,),(2,4)))
    try:d.footprint((5,));bad=False
    except ValueError:bad=True
    c('GpCconstructorrejectsunqualifiednonGpCsite',bad,True)
    g=protocol_graph();c('methylfootprinting followscontactligation',tuple(g.steps).index('footprint')>tuple(g.steps).index('contact_ligation'),True)
    c('DNA/RNAseparationisphysicalsplit',g.steps['split_branches'].outputs,('rna_nuclei','dna_nuclei'))
    c('RNArecaptureaddsconditionswithoutfreshSilane', 'rna_recapture_beads' not in g.materials and 'rna_recapture_binding_conditions' in g.materials,True)
    c('SAMdonornotinDNAsequenceobject',g.materials['sam'].kind,K.MOLECULAR)
    p=Source(sources_root()/'nome-hic__10.1186+s13059-023-02889-x'/'s13059-023-02889-x_supp_13059_2023_2889_MOESM1_ESM.pdf.txt','get_sources --slugs nome-hic__10.1186+s13059-023-02889-x; originalSIp3rendered')
    if have(c,p,label='publishedallocation'):
        t=p.path.read_text();c('literalSIallocationretained','300µL for RNA-seq library preparation, 1200µL for whole-genome bisulfite' in t,True)
    c.report()
if __name__=='__main__':selftest()
