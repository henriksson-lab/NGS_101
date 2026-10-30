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
SOURCE='Published PMC10307641 overview70–90/Methods188–248/processing256–258; sourceoligooriginalMOESM9 (mainTableS7, workbookheaderTable5); reporting-summaryMOESM2 selectedidentity. DOI10.1038/s41586-023-06228-9.'
MARKERS=(('/5ATTO488N/CTGCTGGAGTTCGTGACC/3ddC/','CTGCTGGAGTTCGTGACC'),('/5ATTO488N/CCTCCTTGTATAAATCCTGGTTGCTGTCTCTTTAT/3ddC/','CCTCCTTGTATAAATCCTGGTTGCTGTCTCTTTAT'))
@dataclass(frozen=True)
class BlockedMarkerDNA:
    body:str
    terminal_base:str='C'
    terminal_sugar:str='2prime3prime-dideoxyribose'
    fluorophore:str='ATTO488N'
    def __post_init__(self):
        if not self.body or set(self.body)-set('ACGT') or self.terminal_base!='C' or self.terminal_sugar!='2prime3prime-dideoxyribose' or self.fluorophore!='ATTO488N':raise ValueError('Explicit sourceATTO-DNAterminalddCnucleotide required')
    @property
    def sequence(self):return self.body+self.terminal_base
    def extend_at_three_oh(self,new_dna):raise ValueError('Source terminalddC lacks3OH required by this selected extension relation; no universal amplification prediction')
    def scene(self,name):
        sc=Scene();sc.strand(name,[Segment(name+'_body',self.body),Segment(name+'_ddc',self.terminal_base,note='Actual terminalddC nucleotide; base retained with2prime3prime-dideoxyribose')],mod5=self.fluorophore,mod3='ddC (no3primeOH)',label='Sourcefluorescent marker, notsampleUMI');return sc
@dataclass(frozen=True)
class SelectedRPF:
    sequence:str
    three_end:str='phosphate'
    def __post_init__(self):
        if not self.sequence or set(self.sequence)-set('ACGU') or self.three_end not in {'phosphate','OH'}:raise ValueError('ExplicitselectedRNA/3end classrequired, populationenddistributionnotpredicted')
    def dephosphorylate(self):return SelectedRPF(self.sequence,'OH')
    def tail(self,adenylates):
        if self.three_end!='OH' or type(adenylates)!=int or adenylates<1:raise ValueError('SelectedpolyAtailingrequires3OH andcallerassignedpositivetail length')
        return SelectedRPF(self.sequence+'A'*adenylates,'OH')

def umi_inventory():
    sc=Scene();sc.strand('umi_adapter',[Segment('umi','N'*12,placeholder=True,feature=feature('riboitp.umi','umi','unknown',note='Source12nt readUMI; physicalkitoligoflanks/orderunassembled')),Segment('unknown_flank','[sourcekitadapterflanks/endgeometryunassembled]',placeholder=True)],label='UMI component inventory only; not an assembledsourceoligo orcomplete read');return sc

def graph_for_mode(mode='single_cell_mnase'):
    if mode not in {'single_cell_mnase','single_cell_rnasei','100_cell_mnase','mouse_mnase_single','mouse_mnase_udi'}:raise ValueError('Selectedsource sample/nuclease/indexbranch required')
    rnasei=mode=='single_cell_rnasei';mouse=mode.startswith('mouse');single_index=mode in {'100_cell_mnase','mouse_mnase_single'}
    g=ProtocolGraph('Ribo-ITP — '+mode,authority=SOURCE,representation='published_instance')
    x=g.input('cells','Source100K562cells' if mode=='100_cell_mnase' else ('Sourcemouseoocyte/embryo collectedindividually' if mouse else 'SourcesingletDAPInegativeK562cellFACScollection'),K.MOLECULAR,primary=True)
    x=op(g,'lyse',x,'Sourcecoldlysisbufferwithcycloheximide, flashfrozenwhere specified',S.EXTRACTION,reagent='SourcelysisbufferTable6deferred; no conventionalmonosome gradientinherited')
    x=op(g,'digest',x,'SourceRNaseI1:300stock37C15min' if rnasei else 'SourceMNase37C30min; dilutionselectedbysamplemode',S.DIGESTION,reagent='SourceRNaseIAmbion' if rnasei else 'SourceMNaseNEB',product='Nuclease-digestedcellularlysate containingRPFs/otherRNA/DNA/protein/debris; exactfragmentends/uniformprotectionnotasserted')
    x=op(g,'quench',x,'Source0.1%SDSstopRNaseI' if rnasei else 'Source10mMEGTAstopMNase',S.QUENCHING,reagent='SourceSDS' if rnasei else 'SourceEGTA')
    inputs=[(x,R.TEMPLATE)]
    for i,(written,body) in enumerate(MARKERS):
        m=BlockedMarkerDNA(body);a=g.input('marker_'+str(i),'SourcefluorescentDNAmarker '+written,K.OLIGO,state=m.scene('marker_'+str(i)),chemical_state=m);inputs.append((a,R.SPIKE_IN))
    xmix=g.product('lysate_markers','Actualdigestedlysate plusblockedfluorescentDNAmarkers',K.MOLECULAR)
    g.react('add_markers','SourcefluorescentmarkeradditionbeforeITP',S.GENERIC,inputs,(xmix,))
    apparatus=g.input('itp_chip','SourcePDMSchip5%/10%polyacrylamide, leading/trailingelectrolytes; fullbuffer/chiporiginalsdeferred',K.SELECTION)
    kept=g.product('itp_eluate','Source10ulbetween-markersRPF-enrichedelution collectedinto3dephosphorylationbuffer; marker/carryoverexclusionnotperfect',K.MOLECULAR);other=g.product('itp_other','Physicalotherchannel/reservoir/washmaterial includingunselectedRNA/DNAmarkers/protein; exactpartitionandlaterfateunresolved',K.MOLECULAR)
    g.react('itp','Sourceelectrophoreticextraction5%gel→10%sizeselection, elutionwellwash/bufferexchange thenbetween-markercollection',S.SELECTION,((xmix,R.TEMPLATE),(apparatus,R.SELECTION)),(kept,other));g.mark_terminal(Terminal(other.id,Disposition.UNRESOLVED,'Actualcompanionfractions; not assumedpurewaste orperfectmarkerremoval.'))
    x=op(g,'dephosphorylate',kept,'SourceD-PlexC05030001 dephosphorylation plus0.5ulT4PNK25min',S.DEPHOSPHORYLATION,catalyst=True,reagent='SourcePNK catalyticcomponent; kitbuffer/othercomponentsunassembled',product='SourceRPF3OHselectedclasses withotherendclasses/markers/companions stillpossible; dephosphorylationdoesnotremoveactualterminalddC')
    x=op(g,'poly_a',x,'SourceFig2cRPFpolyAtailingwithinsingle-tubeD-Plexkit',reagent='HistoricalD-PlexpolyAoperation/internalreagentsunreviewed',product='SourcepolyA-tailedRNAclasses plusuntailedRNA/markers/remnants; taildistributionunassigned')
    umi=g.input('umi_component','SourceUMI-containingadapterinventory; exactD-Plexoligo/polymergeometryunassembled',K.OLIGO,state=umi_inventory())
    rt=g.input('rt_components','SourceD-PlextemplateswitchingRT/primer/TSOmixture; TSOdiluted1:2, actualcomponentsunreviewed',K.CONDITION);cdna=g.product('cdna_mix','SourcecDNApopulationwithUMIadapter plusoriginalRNA/markerDNA/unextendedprimers/remnants; fullTSO/primer/endgeometryunassembled',K.MOLECULAR)
    g.react('rt','SourceD-Plexreverse-transcription/templateswitchinterface',S.MAGIC,((x,R.TEMPLATE),(umi,R.STRUCTURAL),(rt,R.CONDITION)),(cdna,),reason='SourceRT/templateswitchknown, exactkitRTprimer/oligoorder andmolecularendsunreviewed; no importedSMARTrule.')
    amplified_half=g.product('cdna_half','Sourcehalf-cDNAallocation for17cyclePCR',K.MOLECULAR);saved=g.product('cdna_remainder','Sourceotherhalf-cDNA reactionallocation, laterfateunassigned',K.MOLECULAR)
    g.react('halve_cdna','Sourcehalf-cDNAusedforPCR',S.SPLIT,((cdna,R.TEMPLATE),),(amplified_half,saved));g.mark_terminal(Terminal(saved.id,Disposition.UNRESOLVED,'IndependentremainderofphysicalcDNA/RNA/marker/reagentmixture; source laterfate/storageunassigned.'))
    x=op(g,'pcr',amplified_half,'SourceD-Plex17cyclePCR; modulechoice '+('singleindex' if single_index else ('uniquedualindex' if mode=='mouse_mnase_udi' else 'single-cellhumanindexmoduleunassigned')),S.PCR,reagent='SourceD-PlexPCRmix; exactindexmembers/primerendsunassembled',product='PCRDNAcopies plusoriginalRNA/cDNA/blockedmarkers/reagents; no perfectmarkerexclusionasserted')
    assay=g.product('bioanalyzer_aliquot','Source1:4dilutedlibraryanalyticalaliquot; volumeunassigned',K.QC);xkeep=g.product('production_library','RetainedproductionlibraryseparatefromBioanalyzerassayedmaterial',K.MOLECULAR)
    g.react('assay_split','SourceBioanalyzerHSassessment ofdilutedlibrary',S.SPLIT,((x,R.TEMPLATE),),(assay,xkeep));g.mark_terminal(Terminal(assay.id,Disposition.QC,'Analyticallyassayedfractionisnotproductionpool.'))
    other_libraries=g.input('other_libraries','Other source libraries selected for approximately equimolar pooling; actual membership and sequences unassigned',K.MOLECULAR)
    x=g.product('pooled_library','Physical pool of retained production library and other source libraries; exact members/amounts unresolved',K.MOLECULAR)
    g.react('pool','Source approximately equimolar library pooling guided by target-peak concentration',S.POOL,((xkeep,R.TEMPLATE),(other_libraries,R.TEMPLATE)),(x,))
    x=clean(g,'ampure',x,'Source1.8xAMPureXPcleanup')
    if rnasei:selection='SourceRNaseIsingle-cell200bpproduct PAGEpurification'
    elif single_index:selection='SourceBluePippin3%dye-freeBDQ3010tight173–207bp window'
    elif mode=='mouse_mnase_udi':selection='SourceBluePippin3%dye-freeBDQ3010tight183–217bp window'
    else:selection='SourceBluePippinBDQ3010 sizeselection; single-cellhumanindexmodule/windowunassigned'
    x=clean(g,'size_select',x,selection)
    g.sequence(x,SequencingPlan(acquisition='unspecified',limitation='SourceNovaSeq6000; exactrun/readprimer/cyclesunreviewed. Processingextracts12ntUMIanddiscardsfollowing4RT-derivednt computationally; physicalUMI/RTjunctiongeometryunassembled. SourceRNAseqSmartseq3andconventionalSMARTersmRNAlibrariesnotinherited.'))
    g.add_source_item(SourceItem('pending','SourceonlineRiboITPprotocol, Table6buffers, fullSIfigures/vendorD-Plexmoduleguide','Cap3originalsclosed; exacthistoricalkit/revisions, oligos/ends/indices/readgeometry and physicalfractionquantities unreviewed, not absent.'))
    return g

def protocol_graph():return graph_for_mode()
def selftest():
    from checks import Check,Source,have
    from storage import sources_root
    c=Check();marker=BlockedMarkerDNA(MARKERS[0][1]);r=SelectedRPF('ACGU')
    c('ActualddCnucleotideaddsCnotmereannotation',marker.sequence,MARKERS[0][1]+'C')
    try:marker.extend_at_three_oh('A');bad=False
    except ValueError:bad=True
    c('Selected3OHextensioncannoteraseblockedddC',bad,True)
    c('Selected3phosphateremovalretainsRNAsequence',r.dephosphorylate().sequence,r.sequence)
    try:r.tail(2);bad=False
    except ValueError:bad=True
    c('TailconstructorrejectsselectedblockedRPF',bad,True)
    c('TailingpreservesoriginalRPF',r.dephosphorylate().tail(2).sequence,'ACGUAA')
    g=protocol_graph();c('Sourcehalf-cDNAremainderdoesnotfeedPCR',g.steps['pcr'].inputs[0].material,'cdna_half')
    c('ActualpoolcontainsothermolecularDNAinput',g.steps['pool'].kind==S.POOL and tuple(z.material for z in g.steps['pool'].inputs)==('production_library','other_libraries') and g.materials['other_libraries'].kind==K.MOLECULAR,True)
    c('UnusedcDNAhalfhasunresolvedfatenotassumedstorage',g.terminals['cdna_remainder'].disposition,Disposition.UNRESOLVED)
    c('RNaseIusesSDSandPAGEratherthanMNaseEGTA/BluePippin','SDS' in graph_for_mode('single_cell_rnasei').steps['quench'].label and 'PAGE' in graph_for_mode('single_cell_rnasei').steps['size_select'].label,True)
    c('UMIhasstructuredmeaning',umi_inventory().strands['umi_adapter'].segs[0].feature.role,'umi')
    p=Source(sources_root()/'ribo-itp__10.1038+s41586-023-06228-9'/'s41586-023-06228-9_supp_41586_2023_6228_MOESM9_ESM.xlsx.txt','get_sources --slugs ribo-itp__10.1038+s41586-023-06228-9')
    if have(c,p,label='publishedmarkeroriginal'):
        t=p.path.read_text();c('BothwholemodifiedsourceDNAmarkers',all(w in t for w,b in MARKERS),True)
    c.report()
if __name__=='__main__':selftest()
