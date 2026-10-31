"""Protocol-owned source-qualified incomplete molecular graph."""
from dataclasses import dataclass, replace
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Scene, Segment, revcomp, feature
from protocol_graph import ProtocolGraph, MaterialKind as K, InputRole as R, StepKind as S, SourceItem, SequencingPlan, Terminal, Disposition
from pcr import PCRGeometry

def clean(g,key,x,label):
    condition=g.input(key+'_condition',label,K.CONDITION)
    kept=g.product(key+'_kept','Source recovered preparation; exact carryover unquantified',K.MOLECULAR)
    rest=g.product(key+'_rest','Actual physically separated workup fractions; no perfect molecular partition inferred',K.WASTE)
    ports=[(x,R.TEMPLATE),(condition,R.CONDITION)]
    if any(word in label.lower() for word in ('precip','ethanol','isopropanol')) or key=='cdna_size':
        carrier=g.input(key+'_carrier','Actual source precipitation carrier glycogen(ifRNAquantified) orGlycoBlue; source5–10ug/choice/carryover unassigned',K.MOLECULAR);ports.append((carrier,R.STRUCTURAL))
    g.react(key,label,S.CLEANUP,ports,(kept,rest));g.discard(rest,'Actual source workup separation.');return kept

def op(g,key,x,label,kind=S.MAGIC,*,reagent='source reagents; precise geometry unassembled',catalyst=False,molecular=False,product=None,reason='Source operation known; internal molecular geometry and remnant fates unassembled.'):
    stock=g.input(key+'_reagent',reagent,K.CATALYST if catalyst else K.MOLECULAR if molecular else K.CONDITION)
    out=g.product(key+'_product',product or label+'; accompanying physical reaction remnants retained',K.MOLECULAR)
    ports=[(x,R.TEMPLATE),(stock,R.CATALYST if catalyst else R.STRUCTURAL if molecular else R.CONDITION)]
    if molecular:
        cond=g.input(key+'_conditions',label+'; actualsourceconditions, notadditionalmoleculardose',K.CONDITION);ports.append((cond,R.CONDITION))
    g.react(key,label,kind,ports,(out,),reason=reason);return out

def oligo_scene(name,sequence,*,mod5='',mod3=''):
    sc=Scene();sc.strand(name,[Segment(name,sequence)],mod5=mod5,mod3=mod3);return sc

from adapter_ligation import AdapterLigationGeometry
from circular import CircularizationGeometry
from illumina import P5,P7
SOURCE='Zinshteyn et al. DOI10.1093/nar/gky1172; published PMC6326810 Methods58–80 and original SupplementaryTableS1; final issue2019/online2018.'
LINKER_BODY='CACTCGGGCACCAAGGA'
RT_WRITTEN='/5Phos/AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGC/iSP18/CACTCA/iSp18/TTCAGACGTGTGCTCTTCCGATCTGTCCTTGGTGCCCGAGTG'
@dataclass(frozen=True)
class DNA:
    name:str
    sequence:str
    def __post_init__(self):
        if not self.name or not self.sequence or set(self.sequence)-set('ACGT'):raise ValueError('named exact DNA blocks required')
@dataclass(frozen=True)
class Spacer:
    name:str
    code:str
    def __post_init__(self):
        if not self.name or self.code not in {'/iSP18/','/iSp18/'}:raise ValueError('only actual source-written internal spacer codes supported')
@dataclass(frozen=True)
class Payload:
    name:str
    origin:str
    def __post_init__(self):
        if not self.name or not self.origin:raise ValueError('variable RNA-derived DNA requires molecular provenance')
@dataclass(frozen=True)
class PrimerBackbone:
    nodes:tuple
    circular:bool=False
    def __post_init__(self):
        if not isinstance(self.nodes,tuple) or not self.nodes or not all(isinstance(n,(DNA,Spacer,Payload)) for n in self.nodes) or len({n.name for n in self.nodes})!=len(self.nodes):raise ValueError('ordered uniquely named nucleotide/spacer/payload nodes required')
        if not isinstance(self.nodes[0],DNA) or isinstance(self.nodes[-1],Spacer):raise ValueError('source backbone ends must be nucleotide material')
    def extend_selected(self):
        if self.circular or any(isinstance(n,Payload) for n in self.nodes):raise ValueError('RT extension requires the original linear primer')
        return replace(self,nodes=self.nodes+(Payload('cdna','source complementary RNA copy; selected full-length or stopped endpoint unassembled'),))
    def close(self):
        if self.circular or not isinstance(self.nodes[-1],Payload):raise ValueError('source circle needs extended linear cDNA, not bare primer')
        return replace(self,circular=True)
    def exact_sequence(self):
        if any(not isinstance(n,DNA) for n in self.nodes):raise ValueError('non-nucleotide links and variable copies cannot be erased into an ordinary DNA sequence')
        return ''.join(n.sequence for n in self.nodes)
    def written(self):return ('circle(' if self.circular else '/5Phos/')+''.join(n.sequence if isinstance(n,DNA) else n.code if isinstance(n,Spacer) else '['+n.origin+']' for n in self.nodes)+(')' if self.circular else '')
    def scene(self):
        sc=Scene();segs=[Segment(n.name,n.sequence) if isinstance(n,DNA) else Segment(n.name,'[C18 linker]' if isinstance(n,Spacer) else '['+n.origin+']',placeholder=True,note='Non-nucleotide covalent node, no base length' if isinstance(n,Spacer) else '') for n in self.nodes]
        sc.strand('backbone',segs,label='source circle, opened for display' if self.circular else 'source primer/cDNA',mod5='' if self.circular else 'phosphate');return sc

@dataclass(frozen=True)
class RTMixture:
    selected_cdna:PrimerBackbone
    def scene(self):
        sc=self.selected_cdna.scene()
        sc.strand('original_rna',[Segment('rna','[original modified RNA population retained before NaOH hydrolysis]',placeholder=True)])
        sc.strand('companions',[Segment('carryover','[possible unextended RT primer, enzymes and reagent remnants]',placeholder=True)])
        return sc

def rt_primer():
    return PrimerBackbone((DNA('five','AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGC'),Spacer('spacer1','/iSP18/'),DNA('middle','CACTCA'),Spacer('spacer2','/iSp18/'),DNA('three','TTCAGACGTGTGCTCTTCCGATCTGTCCTTGGTGCCCGAGTG')))

def linker_scene():
    # ddC is a cytidine nucleotide as well as a 3-prime block; it is not deleted from the template.
    sc=Scene();sc.strand('linker',[Segment('body',LINKER_BODY),Segment('terminal_ddc','C',note='source terminal dideoxycytidine nucleotide')],mod5='pre-adenylated',mod3='ddC: no3OH');sc.ligation_end('linker','body',through='terminal_ddc',kind='single-strand',owner='adapter',polarity='5-prime');return sc
@dataclass(frozen=True)
class ProbedRNA:
    sequence:str
    c8_adducts:frozenset[int]
    def __post_init__(self):
        if not self.sequence or set(self.sequence)-set('ACGU') or any(type(i)is not int or not 0<=i<len(self.sequence) or self.sequence[i] not in 'AG' for i in self.c8_adducts):raise ValueError('selected C8 adducts require actual source A/G residues')
    def selected_dna(self,substitutions=None,stop=None):
        # Explicit conditional copying outcomes only; adduct identity is not itself a mutation.
        changes=dict(substitutions or {})
        if set(changes)-self.c8_adducts or any(b not in 'ACGT' for b in changes.values()):raise ValueError('selected mutation outcomes require explicit source-adduct positions and DNA bases')
        if stop is not None and (stop not in self.c8_adducts or changes):raise ValueError('selected stop and mutation classes must be separately supplied')
        positions=range(len(self.sequence)-1,stop if stop is not None else -1,-1)
        return ''.join(changes.get(i,{'A':'T','C':'G','G':'C','U':'A'}[self.sequence[i]]) for i in positions)

def graph_for_branch(branch='map',rt='ssii',sample='ecoli'):
    if branch not in {'map','stop'} or rt not in {'ssii','tgirt'} or sample not in {'ecoli','k562','yeast'} or rt=='tgirt' and sample!='ecoli':raise ValueError('source branch or RT/sample combination unsupported')
    g=ProtocolGraph('LASER-Seq — '+sample+'/'+branch+'/'+rt,authority=SOURCE,representation='published_instance')
    x=g.input('rna_source','Source '+sample+' ribosome/cell preparation; actual adduct sites unassigned',K.MOLECULAR,primary=True)
    x=op(g,'photochemistry',x,'Source NAz/UV: '+('K562300mM,1or3minUV' if sample=='k562' else 'ribosome NAz addition and3minUV after source incubation'),S.CHEMICAL_CONVERSION,reagent='Actual reactive nicotinoyl azide; UV/temperature separateconditions, sourcebranchspecificconcentrations',molecular=True,product='RNA population with selected A/G C8 adducts plus unreacted probe/protein/cell material; no deterministic mutation/stop assignment')
    if sample=='k562':x=op(g,'extract',x,'Source Trizol/DirectZol totalRNA extraction, DNase omitted',reagent='historical ZymoDirectZol extraction procedure; internal mechanism unreviewed')
    else:x=clean(g,'extract',x,'Source isopropanol recovery, PCA extraction until protein interface gone, chloroform and ethanol recovery')
    x=op(g,'fragment',x,'Source10mMZnCl2 fragmentation95°C5min',S.FRAGMENTATION,reagent='ActualsourceZnCl2 molecularcontributor; heat separatecondition',molecular=True,product='Fragmented adduct-bearing RNA plus physical cleavage products; exact termini unassigned')
    x=op(g,'quench',x,'Source ice and immediate EDTA20mM quench',S.QUENCHING,reagent='ActualsourceEDTA molecularquench; ice separatecondition',molecular=True)
    x=clean(g,'rna_size',x,'Source ethanol, denature95°C3min,10%ureaPAGE section100–120bp DNA-ladder equivalent; gel elution/filter/isopropanol')
    x=op(g,'pnk',x,'Source T4PNK37°C1h, no ATP addition in listed mixture',S.DEPHOSPHORYLATION,catalyst=True,reagent='source T4PNK/PNKbuffer; selected ligation-competent3OH inferred from reaction, full terminal distribution unassembled')
    a=g.input('linker','Published oBZ191: /preA/CACTCGGGCACCAAGGA/ddC/; purchased5phosphate then preadenylated in-house',K.ADAPTER,state=linker_scene())
    lig=g.input('rnl2','Source T4RNA ligase2truncated37°C3h',K.CATALYST);joined=g.product('ligated','Source RNA-linker product population plus unligated material; terminal ddC retained as nucleotide/block',K.MOLECULAR)
    g.react('ligate3','Source preadenylated linker ligation',S.LIGATION,((x,R.TEMPLATE),(a,R.STRUCTURAL),(lig,R.CATALYST)),(joined,),geometry=AdapterLigationGeometry.known('single-strand',linker_scene(),statement='Source TableS1 preadenylated5 end and terminalddC; full substrate termini/junction unassembled.'))
    x=clean(g,'lig_cleanup',joined,'Source isopropanol recovery')
    if sample=='ecoli':
        chosen=g.product('rt_selected','Actual selected Ecoli RNA allocation for '+rt+' branch; volumes unassigned',K.MOLECULAR);other=g.product('rt_other','Other actual source EcoliRNA RTbranch allocation; amount/furthertraceunassembled',K.MOLECULAR)
        g.react('rt_allocation','Source EcoliRNA split beforealternativeRTs, notduplicatedwhole stock',S.SPLIT,((x,R.TEMPLATE),),(chosen,other));g.mark_terminal(Terminal(other.id,Disposition.UNRESOLVED,'Other source RTbranch outside selectedtrace.'));x=chosen
    primer=g.input('rt_primer','Published oBZ192 with5phosphate and two non-nucleotide18atom spacers',K.OLIGO,state=rt_primer().scene(),chemical_state=rt_primer())
    if rt=='ssii':
        anneal_buffer=g.input('rt_anneal_buffer','Actual source MaPBasic buffer for6ulRNA/primerannealing preparation',K.CONDITION)
        retained=g.product('retained_rt_primer','Same source oBZ192 carriedthroughannealing; notfreshaddition',K.OLIGO,state=rt_primer().scene(),chemical_state=rt_primer())
        prepared=g.product('annealed_rna','Source RNA/annealingcontext, copresentprimerrepresentedseparately; no physicalfractionation',K.MOLECULAR)
        temp=g.input('rt_anneal_temperature','Source65C4min/55C2min/45C2min/42C2min/ice',K.CONDITION)
        g.react('rt_anneal','SourceRNA+primer annealing BEFOREfresh6ulextensionmix',S.HYBRIDIZATION,((x,R.TEMPLATE),(primer,R.STRUCTURAL),(anneal_buffer,R.CONDITION),(temp,R.CONDITION)),(prepared,retained));x=prepared;primer=retained
    stock=g.input('rt_supplies','Actual source6ulSSIIextensionmix' if rt=='ssii' else 'Actual sourceTGIRT10ulreaction reagentinventory; stockcount/packagingunassigned',K.MOLECULAR)
    cat=g.product('rt_enzyme','Copresentactual'+('SSII' if rt=='ssii' else 'TGIRTIII')+' catalystcomponent',K.CATALYST)
    nts=g.product('rt_nucleotides','Actual source dNTP component ofsameRTsupplies',K.MOLECULAR)
    protective=g.product('rt_protection','Actual copresentSUPERaseIn/protectivereagents, notanotherfreshstock',K.MOLECULAR)
    conditions=g.product('rt_conditions','SameactualRTbuffer/Mn2+forSSII/DTT conditions; source42C1hSSII or60C1hTGIRT',K.CONDITION)
    g.react('rt_components','ActualRTsupplycomponents, noextraindependent doses',S.GENERIC,((stock,R.TEMPLATE),),(cat,nts,protective,conditions),changes_molecule=False)
    cdna=g.product('cdna','Source RT mixture: ordered spacer-containing extended primer, original RNA and reagent remnants retained',K.MOLECULAR,state=RTMixture(rt_primer().extend_selected()).scene(),chemical_state=RTMixture(rt_primer().extend_selected()))
    g.react('rt','Source '+rt+' reverse transcription',S.REVERSE_TRANSCRIPTION,((x,R.TEMPLATE),(primer,R.STRUCTURAL),(cat,R.CATALYST),(nts,R.STRUCTURAL),(protective,R.STRUCTURAL),(conditions,R.CONDITION)),(cdna,))
    x=op(g,'alkaline',cdna,'Source5MNaOH1µl95°C3min RNA hydrolysis',S.CHEMICAL_CONVERSION,reagent='ActualsourceNaOH molecularcontributor; heat separatecondition',molecular=True,product='Spacer-containing cDNA plus hydrolyzed RNA residues and other remnants; hydrolysis alone not physical removal')
    x=clean(g,'cdna_size',x,'Source precipitation,95°C3min denaturation,10%ureaPAGE: '+('full-length products' if branch=='map' else 'products~20nt above primer to~20nt below full-length; untreated controls same class')+'; gel extraction')
    circle=rt_primer().extend_selected().close();circle_stock=g.input('circle_stock','Actual sourceCircLigasepreparation15.5ulTris+2ul10xbuffer+1ul1mMATP+1ul50mMMnCl2+.5ulCircLigaseI',K.MOLECULAR)
    cat=g.product('circligase','Actual copresentCircLigaseI catalyst',K.CATALYST);circle_atp=g.product('circle_atp','ActualcopresentATP substrate fromsameCircLigase mix',K.MOLECULAR);circle_conditions=g.product('circle_conditions','ActualcopresentMn2+/buffer;60C2h/80C20min',K.CONDITION)
    g.react('circle_components','SameactualCircLigase mixture components, notextrafresh doses',S.GENERIC,((circle_stock,R.TEMPLATE),),(cat,circle_atp,circle_conditions),changes_molecule=False)
    out=g.product('circle_material','Source circular cDNA class retaining BOTH source spacers; unclosed/reaction material not assumed absent',K.MOLECULAR,state=circle.scene(),chemical_state=circle)
    g.react('circle','Source cDNA intramolecular circularization',S.CIRCULARIZATION,((x,R.TEMPLATE),(cat,R.CATALYST),(circle_atp,R.STRUCTURAL),(circle_conditions,R.CONDITION)),(out,),geometry=CircularizationGeometry.known('single-strand intramolecular ligation',statement='MainMethods76 and source5phosphate RTprimer; complete circle-to-PCR path unassembled, spacer traversal not assumed.'))
    f=g.input('pcr_f','Source oBZ287 forward primer',K.OLIGO,state=oligo_scene('forward',P5))
    r=g.input('pcr_r','Source barcoded reverse-primer family TableS1; selected member not assigned',K.OLIGO)
    supplies=g.input('pcr_supplies','ActualsourcePCRpreparation supplies; dose/stockcount/freshversusretainedcatalyst provenance unassembled, no borrowedpolymeraseidentity',K.MOLECULAR)
    target=g.product('pcr_template','Same source circularcDNA/primer/reagentcontext; copresentfunctionsrepresentedseparately, notphysicalfractionation',K.MOLECULAR)
    f_retained=g.product('pcr_forward','Same sourceforwardprimer inactualPCRpreparation',K.OLIGO,state=oligo_scene('forward',P5));r_retained=g.product('pcr_reverse','Sameunknownmember ofsourcebarcodedreverseprimerfamily',K.OLIGO)
    poly=g.product('pcr_enzyme','ActualPCRcatalyticfunction; identity/addedversuscarriedprovenanceunassembled',K.CATALYST);pcr_nt=g.product('pcr_nt','ActualPCRnucleotides; dose/sourceprovenanceunassembled',K.MOLECULAR)
    g.react('pcr_components','ActualsourcePCRpreparation inventory, notassertedcommercialpremixorfreshenzymedose',S.GENERIC,((out,R.TEMPLATE),(f,R.STRUCTURAL),(r,R.STRUCTURAL),(supplies,R.STRUCTURAL)),(target,f_retained,r_retained,poly,pcr_nt),changes_molecule=False)
    amp=g.product('amp','Source PCR product mixture plus spacer-bearing circular templates/remnants; no spacer-erased ordinary DNA circle substituted',K.MOLECULAR)
    g.react('pcr','Source8–12cycle PCR',S.PCR,((target,R.TEMPLATE),(f_retained,R.STRUCTURAL),(r_retained,R.STRUCTURAL),(poly,R.CATALYST),(pcr_nt,R.STRUCTURAL)),(amp,),geometry=PCRGeometry.unavailable(Scene(),reason='Exact linker-containing template traversal/polymerase response and primer/end assembly unreviewed.'))
    amp=clean(g,'final_gel',amp,'Source finalPCR gel purification')
    loaded=g.product('loading_selected','Actualsource selectedsequencinglibrary contribution; amountunassigned',K.MOLECULAR);remaining=g.product('loading_remainder','UntransferredfinalDNAstock; amount/fateunassigned',K.MOLECULAR)
    g.react('loading','Source selectedlibrary contribution, noentirestockconsumptionassumption',S.SPLIT,((amp,R.TEMPLATE),),(loaded,remaining));g.mark_terminal(Terminal(remaining.id,Disposition.UNRESOLVED,'UntransferredsourceDNAstock.'))
    g.sequence(loaded,SequencingPlan(acquisition='unspecified',limitation='Published mainMethods76 HiSeq2500 single-end50bp; exact sequencing-primer landing and circle-to-PCR topology remain unassembled.'))
    g.add_source_item(SourceItem('predecessors','Actual directly cited library and synthesis predecessors','MainMethods52 and76 references8,22,23,31; source library/linker traversal detail unreviewed at three-original cap.'))
    return g

def protocol_graph():return graph_for_branch()
def selftest():
    from checks import Check
    c=Check();p=rt_primer();circle=p.extend_selected().close()
    c('RT and closure preserve both ordered non-nucleotide nodes',[(n.name,n.code) for n in circle.nodes if isinstance(n,Spacer)],[('spacer1','/iSP18/'),('spacer2','/iSp18/')])
    try:circle.exact_sequence();bad=False
    except ValueError:bad=True
    c('spacer-containing circle cannot flatten to ordinary DNA',bad,True)
    c('RT mixture retains originalRNA independently before chemical hydrolysis','original_rna' in protocol_graph().materials['cdna'].state.strands,True)
    c('ddC contributes a real terminal cytidine', ''.join(s.top for s in linker_scene().strands['linker'].segs),LINKER_BODY+'C')
    r=ProbedRNA('ACGU',frozenset({2}));c('C8 adduct does not physically replace RNA base',r.sequence,'ACGU')
    c('selected stop gives3prime-side complementary DNA without modified position',r.selected_dna(stop=2),'A')
    c('selected mutation is assigned in copied DNA only',r.selected_dna({2:'T'}),'ATGT')
    try:ProbedRNA('ACGU',frozenset({1}));bad=False
    except ValueError:bad=True
    c('C8 purine operation rejects cytosine as purine substrate',bad,True)
    g=protocol_graph();c('SSII uses primerretainedthroughactualannealing',g.steps['rt'].inputs[1].material,'retained_rt_primer')
    c('RT consumesactualnucleotidecomponent',g.steps['rt'].inputs[3].material,'rt_nucleotides')
    c('Circle consumesactualATP component',g.steps['circle'].inputs[2].material,'circle_atp')
    c('ActualotherEcoliRTbranchretained',g.terminals['rt_other'].disposition,Disposition.UNRESOLVED)
    c('FinalDNAstockremainderretained',g.terminals['loading_remainder'].disposition,Disposition.UNRESOLVED)
    c('both RT readout branches stay source/library incomplete',all(graph_for_branch(b).publication_problems() for b in ('map','stop')),True)
    from checks import Source,have
    from storage import sources_root
    source=Source(sources_root()/'laser-seq__10.1093+nar+gky1172'/'nar_gky1172_supp_gky1172_supplemental_files_180423_LASERseq_SuppInfo_Table_S1_Oligos.xlsx.txt',"python3 tools/get_sources.py --slugs laser-seq__10.1093+nar+gky1172; inspected original selected chemistry sections")
    if have(c,source,label="local original transcription guard"):
        c('original TableS1 retains exact two-spacer RT primer',rt_primer().written() in source.path.read_text(),True)
        c('original TableS1 prints nucleotide ddC with preadenylated linker','/preA/'+LINKER_BODY+'/ddC/' in source.path.read_text(),True)
    c.report()
if __name__=='__main__':selftest()
