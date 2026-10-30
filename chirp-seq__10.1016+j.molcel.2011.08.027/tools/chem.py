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
SOURCE='PublishedMolCell44pp667–678 author-labPDF source_pdf_87b24b6fe07e40b97170.pdf, Methods1143–1275; printedDOI/title/dateverifiedoriginalp1. AssociatedMSPMC3249421selectedcomparisononly. DOI10.1016/j.molcel.2011.08.027.'
@dataclass(frozen=True)
class ProbeBackbone:
    nodes:tuple=('target-complementary DNA region (sequence unreviewed)','18-carbon spacer arm (internal composition unreviewed)','biotin')
    def __post_init__(self):
        if self.nodes!=('target-complementary DNA region (sequence unreviewed)','18-carbon spacer arm (internal composition unreviewed)','biotin'):raise ValueError('SourceorderedDNA/spacer/biotin architecture required, notvendorlinkerequivalence')
    def dna_sequence(self):raise ValueError('UnknownDNA plusnon-nucleotidespacer/biotincannotflattenintoDNA')
    def captured_dna(self,dna_id):
        if not dna_id:raise ValueError('IndependentcapturedDNAidentityrequired')
        return {'dna':dna_id,'probe_backbone':self.nodes,'dna_direct_biotin':()}

def graph_for_scope(pool='odd',crosslinker='glutaraldehyde'):
    if pool not in {'odd','even'} or crosslinker not in {'glutaraldehyde','formaldehyde'}:raise ValueError('Sourceprobesplitpool/crosslinkbranchrequired')
    g=ProtocolGraph('ChIRP-seq — '+pool+'/'+crosslinker,authority=SOURCE,representation='published_instance')
    x=g.input('cells','SourceHOTAIR/TERC/roX2cellpreparation; exacttargetcontextunassigned',K.MOLECULAR,primary=True)
    x=op(g,'crosslink',x,'Source1%'+crosslinker+'PBS10minRT',S.CROSS_LINKING,reagent='Source'+crosslinker)
    x=op(g,'quench',x,'Source0.125Mglycine5min',S.QUENCHING,reagent='Sourceglycine')
    x=clean(g,'pellet',x,'SourcePBSwash/scrape/pellet;800gFAor2500gGA')
    x=op(g,'swell',x,'SourceSwellingBuffer10minice/dounce',S.EXTRACTION,reagent='SourceTris/KOAc/MgOAc/NP40/DTT/PMSF/proteaseinhibitor/SuperaseIn')
    x=clean(g,'nuclei',x,'Source2500g5minnucleirecovery')
    x=op(g,'nuclear_lysis',x,'SourceSDSnuclearlysis10minice',S.EXTRACTION,reagent='SourceTris/EDTA/SDS/inhibitors')
    x=op(g,'sonicate',x,'SourceBioruptorsolubilizationandDNA100–500bpfragmentation',S.FRAGMENTATION,reagent='SourceBioruptorphysicalsonication',product='FragmentedcrosslinkedDNA/RNA/proteinpreparation withallreactioncompanions; fullends/occupancyunassembled')
    probe=g.input('probe_pool','Source'+pool+'tilingDNApool with3primebiotin/18carbonspacer, exactTableS4membersunreviewed',K.OLIGO,chemical_state=ProbeBackbone(),state=scene('probe_inventory','[DNA target-complementary pool]—[source18-carbon spacer, chemistry unassembled]—[3primebiotin]'))
    conditions=g.input('hyb_conditions','Sourcetwo-volumeshybridizationbuffer/100pmolprobesper3ml,37C4h',K.CONDITION);y=g.product('hybridized','SourceDNAprobe/RNA/chromatincomplexes plusunboundprobes/othermaterial',K.MOLECULAR)
    g.react('hybridize','Sourcesplit-poolRNAantisensehybridization',S.HYBRIDIZATION,((x,R.TEMPLATE),(probe,R.STRUCTURAL),(conditions,R.CONDITION)),(y,))
    bead=g.input('blocked_beads','SourceC1beads washed/blockedwithyeastRNA+BSAandwashed; exactcarryoverunknown',K.SELECTION);bound=g.product('bound','Sourcebead/probe/RNA/chromatinDNApreparation',K.MOLECULAR);flow=g.product('unbound','Physicalunboundcapturefraction, laterfateunassigned',K.MOLECULAR)
    g.react('capture','Source100ulC1beadsper100pmolprobe,37C30minmagnetcapture',S.AFFINITY_CAPTURE,((y,R.TEMPLATE),(bead,R.SELECTION)),(bound,flow));g.mark_terminal(Terminal(flow.id,Disposition.UNRESOLVED,'Actualremovedcapturefraction; nosourcestatedlaterfateorperfectpurity.'))
    x=clean(g,'wash',bound,'Sourcefive40xbeadvolume2xSSC/0.5%SDSwashes')
    eluents=[]
    for n in (1,2):
        key='rna_elution_'+str(n)
        x=op(g,key,x,'SourceRNaseA100ug/ml/RNaseH0.1U/ul37CelutioninDNAelutionbuffer',S.DIGESTION,reagent='SourceRNaseA/H; sourceDNAelutionNaHCO3/SDS/NaCl',product='ReleasedDNA-containingreactionplusRNAcleavageproducts/beads/probe/remnants; noperfectRNAremovalasserted')
        liquid=g.product(key+'_eluate','SourcecollectedDNA-containingeluent'+str(n),K.MOLECULAR);residue=g.product(key+'_residue','Remainingbead-associatedmaterialafterelution'+str(n),K.MOLECULAR)
        g.react(key+'_collect','Sourcephysicaleluentcollection',S.CLEANUP,((x,R.TEMPLATE),),(liquid,residue));eluents.append(liquid);x=residue
    g.mark_terminal(Terminal(x.id,Disposition.UNRESOLVED,'Sourcepost-elutionbead/probe/residualDNA/RNAcompanion, laterfateunassigned.'))
    x=g.product('pooled_eluate','Actualcombinedfirst+secondDNA-containingeluents, withunremovedreactioncompanions',K.MOLECULAR);g.react('combine_eluents','Sourcecombinesbothelutionsteps',S.POOL,tuple((z,R.TEMPLATE) for z in eluents),(x,))
    if crosslinker=='formaldehyde':x=op(g,'reverse_crosslink',x,'SourceFA65Covernightcrosslinkreversal',S.THERMAL_TREATMENT,reagent='Source65Covernightconditions')
    else:x=op(g,'proteinase',x,'SourceGA/nonreversiblecrosslinkpreparation PK65C45min',S.DEPROTEINIZATION,reagent='SourcePK0.2U/ul',product='ProteolyzedDNA/RNA/proteinremnantmixture; GAdoesnotbecomeFAreversiblecrosslinks')
    x=clean(g,'organic_recover',x,'Sourcephenol/chloroform/isoamylaqueousDNArecoveryandEtOHminus80Covernightprecipitation')
    x=op(g,'dna_library',x,'SourcecitedJohnson2007ChIP-seqDNA-libraryconstructioninterface',reagent='ActualcitedpredecessorDNAendrepair/adapters/PCR/readgeometryunreviewed; no modernChIPkit substituted')
    g.sequence(x,SequencingPlan(acquisition='unspecified',limitation='PublishedGenomeAnalyzerIIx36bpreads; exactlibraryprimers/index/endsunassembled. Odd/evenconcordanceiscomputedcomparison, notphysicalmixingofcapturepools. RNA/proteinvalidationelutionsareseparateassaybranches.'))
    g.add_source_item(SourceItem('pending','ActualpublishedSI/probeTableS4andcitedJohnson2007DNA-libraryoriginal','Mainfinalsourceknownoperationsmodelled; probesequences/fullspacercomposition/end/librarygeometryunreviewed, noabsenceclaim.'))
    return g

def protocol_graph():return graph_for_scope()
def selftest():
    from checks import Check,Source,have
    from storage import sources_root
    c=Check();p=ProbeBackbone();g=protocol_graph();f=graph_for_scope(crosslinker='formaldehyde')
    c('CapturedchromatinDNAdoesnotinheritprobebiotin',p.captured_dna('DNA')['dna_direct_biotin'],())
    try:p.dna_sequence();bad=False
    except ValueError:bad=True
    c('CarbonspacerandbiotincannotflattenintoDNA',bad,True)
    c('SecondRNaseelutionusesremainingbeadpreparation',g.steps['rna_elution_2'].inputs[0].material,'rna_elution_1_residue')
    c('EluentpoolcontainsbothactualDNAfractions',tuple(x.material for x in g.steps['combine_eluents'].inputs),('rna_elution_1_eluate','rna_elution_2_eluate'))
    c('GAbranchdoesnotborrowFAovernightreversal','proteinase' in g.steps and 'reverse_crosslink' not in g.steps and 'reverse_crosslink' in f.steps,True)
    c('Oddandevenareindependentcaptureinputs',graph_for_scope('odd').materials['probe_pool'].label!=graph_for_scope('even').materials['probe_pool'].label,True)
    q=Source(sources_root()/'chirp-seq__10.1016+j.molcel.2011.08.027'/'source_pdf_87b24b6fe07e40b97170.pdf.txt','get_sources --slugs chirp-seq__10.1016+j.molcel.2011.08.027')
    if have(c,q,label='publishedfinal'):
        t=q.path.read_text();c('Finalsourceprints18carbonspacerand36bpacquisition','18-carbon spacer' in t and '36 bp' in t,True)
    c.report()
if __name__=='__main__':selftest()
