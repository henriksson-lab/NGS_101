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
    inputs=[(x,R.TEMPLATE)];kwargs={}
    if kind==S.FRAGMENTATION:
        conditions=g.input(key+'_reagent',reagent,K.CONDITION);inputs.append((conditions,R.CONDITION))
    else:
        group=g.input(key+'_reagent','Actualsource supplies: '+reagent+'; stockcount/packaging unassembled unless explicitlynamedmix',K.MOLECULAR)
        enzyme_kinds={S.DIGESTION,S.END_REPAIR,S.EXTENSION,S.TAILING,S.PROXIMITY_LIGATION,S.LIGATION,S.PCR}
        if catalyst or kind in enzyme_kinds:
            names={'exonucleases':'ActualExoIII+ExoVII catalystcomponents, twoenzymeclasses notonefusedprotein','nick':'ActualUDG+sourceSI EndoIV catalystcomponents; mainEndoVIwordingretained','nick_translate':'ActualfreshDNApolymeraseI','a_tail':'ActualfreshTaqpolymerase','y_adapter':'ActualM0367ligasecomponent','pcr':'ActualQ5UM0597 mastermixpolymerasecomponent','target_pcr':'ActualdistinctQ5M0544 mastermixpolymerasecomponent'}
            enzyme=g.product(key+'_enzyme',names.get(key,'Actualsourcecatalyticcomponents, identitiesunassembled'),K.CATALYST)
            nt_label={'nick_translate':'Actual0.5ul10mMdATP/dGTP/dTTPmix plus0.5ul10mM5-mdCTP','a_tail':'Actual1ul100mMdATP, sourceprintsfinal2mM; totalvolumeconventionunassembled','pcr':'ActualQ5U mastermixnucleotidecomponent','target_pcr':'ActualQ5 mastermixnucleotidecomponent'}.get(key,'Otheractualsourceadditions; no uninspectednucleotideidentity/dose inferred')
            molecules=g.product(key+'_molecules',nt_label,K.MOLECULAR);conditions=g.product(key+'_conditions',label+'; actualsourcebuffer/conditions',K.CONDITION)
            parts=[enzyme,molecules,conditions]
            if kind==S.LIGATION:
                adapter=g.product(key+'_adapter','Actualsource1ul7.5uMmethylatedadapter fromsameaddedreagentgroup; fullhistoricalendsunassembled',K.ADAPTER,state=scene(key+'_adapter','[source adapter full geometry unassembled]'));parts.append(adapter);inputs.append((adapter,R.STRUCTURAL))
            g.react(key+'_components','Actualaddedrecipefunctions; notextra freshdoses orunreportedcommercialpremix',S.GENERIC,((group,R.TEMPLATE),),parts,changes_molecule=False)
            inputs.extend(((enzyme,R.CATALYST),(molecules,R.STRUCTURAL),(conditions,R.CONDITION)))
        else:inputs.append((group,R.STRUCTURAL))
    if kind==S.MAGIC:kwargs['reason']='Sourceoperation established; occupiedmolecules/historicalgeometryunassembled, no sourceabsenceclaim.'
    if kind==S.PCR:
        kwargs['geometry']=PCRGeometry.unavailable(Scene(),reason='Exactsourcefulltemplate/primerendpointgeometryunassembled, noordinarycompletelibraryasserted.')
        primers=[]
        for arm in ('f','r'):
            make=g.product if key=='target_pcr' else g.input
            primer=make(key+'_primer_'+arm,'ActualsourcePCRprimerclass; wholeoligo/indexgeometryunassembled fromhistoricaloriginal',K.OLIGO,state=scene(key+'_primer_'+arm,'[source primer body/index geometry unassembled]'));inputs.append((primer,R.STRUCTURAL));primers.append(primer)
        if key=='target_pcr':
            primerstock=g.input('target_quant_primers','ActualsourceE7603QuantPrimerMix',K.MOLECULAR)
            g.react('target_primer_components','SameactualE7603primerfunctionalclasses, nottwofreshstockadditions',S.GENERIC,((primerstock,R.TEMPLATE),),primers,changes_molecule=False)
    if kind in (S.LIGATION,S.PROXIMITY_LIGATION):kwargs['geometry']=AdapterLigationGeometry('unavailable','unavailable',(scene('unassembled','[source-knownligation; fullsubstrate/endgeometryunassembled]'),),'Available/unassembledgeometry, notclaimofundisclosedsource.')
    y=g.product(key+'_product',product or label+'; sourcecompanions retaineduntilactualseparation',K.MOLECULAR,state=state,chemical_state=chemical_state)
    g.react(key,label,kind,inputs,(y,),**kwargs);return y

def qc_split(g,key,x,label):
    assay=g.product(key+'_assay','Physically separate source analytical aliquot',K.QC)
    kept=g.product(key+'_retained','Unassayed retained preparation for downstream library work',K.MOLECULAR)
    g.react(key,label,S.SPLIT,((x,R.TEMPLATE),),(assay,kept));g.mark_terminal(Terminal(assay.id,Disposition.QC,'Source analytical aliquot never feeds retained library.'));return kept
SOURCE='Yan et al.2022 DOI10.1101/gr.277080.122 published PMC9808626 mainMethods148–164; SupplementalMaterials protocolsteps1–7, originalSIp8 andmainp2088 rendered.'
NAME_BOUNDARY='Mainp2088 printsEndoVI(M0304); SIprotocolp8 printsEndoIV(M0304). SelectedSI name retained; manufacturer original unreviewed. Main adapterE7536 versusSI E7165/E7140 not silently reconciled.'
HAIRPIN='CCACGACGACGACGACGAGCGTTAGGCTCGTCGTCGTCGTCGUGGT'
@dataclass(frozen=True)
class UracilHairpin:
    sequence:str=HAIRPIN
    phosphorylated:bool=True
    def __post_init__(self):
        if self.sequence.count('U')!=1 or set(self.sequence)-set('ACGTU') or not self.phosphorylated:raise ValueError('selected source5phosphate single-uracilDNAhairpin required')
    def scene(self):
        sc=Scene();sc.strand('hairpin',[Segment('hairpin',self.sequence,note='DNA oligo uracil; not evidence ofRNA sugar. Intramolecular stem/loop and ligation junction unassembled.')],mod5='phosphate');return sc
@dataclass(frozen=True)
class OriginalCopyPair:
    original:str
    biological_methyl_c:tuple=()
    def __post_init__(self):
        if not self.original or set(self.original)-set('ACGT'):raise ValueError('explicit selectedoriginal5to3DNA required')
        if len(set(self.biological_methyl_c))!=len(self.biological_methyl_c) or any(type(i)!=int or not 0<=i<len(self.original) or self.original[i]!='C' for i in self.biological_methyl_c):raise ValueError('selected biological5mC positions must beC')
    @property
    def synthetic_copy(self):return revcomp(self.original)
    @property
    def synthetic_methyl_c(self):return tuple(i for i,b in enumerate(self.synthetic_copy) if b=='C')
    def converted_original(self):return ''.join('U' if b=='C' and i not in self.biological_methyl_c else b for i,b in enumerate(self.original))
    def converted_copy(self):return self.synthetic_copy
    def complementary_dna_descendant(self):return revcomp(self.converted_original().replace('U','T'))
    def scene(self):
        sc=Scene();sc.strand('selected_original',[Segment('original',self.original,note='Original biological5mC coordinates:'+str(self.biological_methyl_c))]);sc.strand('selected_protected_copy',[Segment('synthetic',self.synthetic_copy,note='Every copiedC uses5-mdCTP, not original biological methylation')]);sc.strand('companions',[Segment('remnants','[displaced DNA products, uncleaved/other products, enzymes and source reaction companions retained until source cleanup]',placeholder=True)]);return sc

def selected(g,key,x,label):
    take=g.product(key+'_selected',label,x.kind,state=x.state,chemical_state=x.chemical_state);rest=g.product(key+'_remainder','Untransferredsource stock; laterfateunassigned',x.kind,state=x.state,chemical_state=x.chemical_state)
    g.react(key,label,S.SPLIT,((x,R.TEMPLATE),),(take,rest));g.mark_terminal(Terminal(rest.id,Disposition.UNRESOLVED,'Untransferredactualsource stock.'));return take

def graph_for_scope(targeted=False):
    g=ProtocolGraph('Methyl-SNP-seq'+(' — exome branch' if targeted else ' — human WGS'),authority=SOURCE,representation='published_instance')
    x=g.input('genomic_dna','SourceNA12878 human genomicDNA; selectedWGS250bp/exome400–500bp sonication',K.MOLECULAR,primary=True)
    spike=g.input('lambda','Actualtwo1ngunmethylatedlambdaDNA spikecontributions (sourceSI), notassumedonepreviouslymixedtube',K.MOLECULAR)
    lambda1=g.product('lambda1','Actualfirst1nglambdaDNAstockcontribution',K.MOLECULAR);lambda2=g.product('lambda2','Actualsecond1nglambdaDNAstockcontribution',K.MOLECULAR)
    g.react('lambda_inventory','Twoactualsourcecontrolcontributions inventory, notDNAcopying',S.GENERIC,((spike,R.TEMPLATE),),(lambda1,lambda2),changes_molecule=False)
    dna1=g.product('sonication_dna1','Actual2ugNA12878DNA sourceSIsonication1allocation',K.MOLECULAR);dna2=g.product('sonication_dna2','Actual2ugNA12878DNA sourceSIsonication2allocation',K.MOLECULAR)
    g.react('dna_allocate','Actualsource4uginputallocatedtotwoSIsonications; biologicalreplicatemappingunassembled',S.SPLIT,((x,R.TEMPLATE),),(dna1,dna2))
    products=[]
    for n,dna,lam in ((1,dna1,lambda1),(2,dna2,lambda2)):
        mixture=g.product('spiked_dna'+str(n),'Actual2ugDNA+1nglambda in50ulsourceSIsonicationpreparation',K.MOLECULAR)
        g.react('spike'+str(n),'ActualunmethylatedlambdaDNAcontroladdition',S.GENERIC,((dna,R.TEMPLATE),(lam,R.SPIKE_IN)),(mixture,),changes_molecule=False)
        products.append(op(g,'sonicate'+str(n),mixture,'SourceCovarisS2 twoactual50ulSIsonications,250bpWGS or400–500bpexomecontext',S.FRAGMENTATION,reagent='Actualsourcephysicalsonicationconditions; independentproductcontributions'))
    x=g.product('sonicated_dna','Actualsource100ulfragmentedDNAcontributionsforendprep; finalsample/technicalreplicatemappingunassigned',K.MOLECULAR)
    g.react('sonication_collect','Actualtwosonicationproductscontributeto100ulsourceSIendprepinterface',S.POOL,tuple((a,R.TEMPLATE) for a in products),(x,))
    x=op(g,'end_repair',x,'SourceNEBUltraIIE7546 endrepair/dA20C30min65C30min',S.END_REPAIR,reagent='Actual9ulEndPrepEnzymeMix+21ulbuffer+50ulwater to100ulfragmentedDNA; historicalE7546 identities/fullendchemistryunassembled')
    x=clean(g,'end_cleanup',x,'SourceZymoOligoClean80ntprotocol')
    h=UracilHairpin();raw=g.input('hairpin','Sourcewritten5P uracilDNAhairpin, exactlinearbodyretained; physicalfoldunassembled',K.OLIGO,state=h.scene(),chemical_state=h)
    prep=g.product('hairpin_stock','Actualsource20uMannealedhairpinstock; stockvolume/foldoccupancyunassigned',K.ADAPTER,state=h.scene(),chemical_state=h)
    cond=g.input('hairpin_anneal','Source1xannealbuffer,80C2minthenramp25C.1C/s; fullfoldunassembled',K.CONDITION)
    g.react('prepare_hairpin','Sourcehairpinannealingfactory, nofoldbyname assumption',S.MAGIC,((raw,R.TEMPLATE),(cond,R.CONDITION)),(prep,),reason='Knownsource20uMpreparation; intramolecularfold/ends notassembled.')
    hp=selected(g,'hairpin_allocate',prep,'Actual4ul20uMpreparedhairpin contribution; otherstockretained')
    mix=g.input('hairpin_ligation_mix','Actual30ulsource2xNEBM0367ligationmix',K.MOLECULAR)
    cat=g.product('hairpin_ligase','CopresentactualM0367catalystcomponent',K.CATALYST);other=g.product('hairpin_ligation_other','OtheractualM0367molecularcomponents, identitiesunreviewed',K.MOLECULAR);buffer=g.product('hairpin_ligation_buffer','CopresentactualM0367buffer;20C1h',K.CONDITION)
    g.react('hairpin_mix_components','SameactualsourceM0367mixfunctions, notextrafreshdoses',S.GENERIC,((mix,R.TEMPLATE),),(cat,other,buffer),changes_molecule=False)
    y=g.product('dumbbell_mix','Twohairpinclosedselectedclass plusincompleteDNA/freehairpins/remnants; fold/closuregeometryunassembled',K.MOLECULAR)
    g.react('ligate_hairpin','Actual4ulhairpin+26ulDNA+30ul2xmix,20C1h',S.LIGATION,((x,R.TEMPLATE),(hp,R.STRUCTURAL),(cat,R.CATALYST),(other,R.STRUCTURAL),(buffer,R.CONDITION)),(y,),geometry=AdapterLigationGeometry('unavailable','unavailable',(h.scene(),),'Sourcehairpinsequence/5Pknown; preparedfold/fullsubstrateendgeometry available/unassembled.'))
    x=clean(g,'ligation_cleanup',y,'Source0.9xSPRI cleanup')
    x=op(g,'exonucleases',x,'SourceExoIII+ExoVII37C1h digest incompletely protected species',S.DIGESTION,reagent='Actual2.5ulExoIII100Uul+2.5ulExoVII10Uul+5ul10xCutsmart to40ulDNA',product='Closed selectedDNAclass plus digestedDNA products and any incompletely digested material; no perfect enrichment asserted')
    x=clean(g,'exo_cleanup',x,'SourceZymo80ntcleanup afterexonuclease')
    x=selected(g,'nick_allocation',x,'Actual39ulfrom40ulcolumneluate forsource48ulUDG/EndoIVreaction')
    x=op(g,'nick',x,'SourceSI UDG/EndoIV37C15min; printednameboundary retained',S.DIGESTION,reagent='SourceUDGM0280/printedSI EndoIVM0304; '+NAME_BOUNDARY,product='Uracil excision/incision-containing DNA classes, base/cleavageproducts and enzymes retained')
    x=op(g,'nick_translate',x,'SourceDNAPolI/5-mdCTP20C30min added directly without purification',S.EXTENSION,reagent='Actual1ulDNApolI+0.5ul10mMd(A,G,T)TPmix+0.5ul10mM5-mdCTP, freshdirectadditions',product='Translation/replacement reaction mixture; selected original/protectedcopy class shown with displacedDNA products. Opposingnick encounter creates source half-dumbbell classes; fullcut geometry unresolved')
    x=op(g,'a_tail',x,'SourceTaq/dATP65C30min added directly without purification',S.TAILING,reagent='Actual2ulTaqM0267+1ul100mMdATPaddedto50ulPolImix; sourceprintsfinal2mM andnextcleanup50ul, conventionunassembled')
    x=clean(g,'tail_cleanup',x,'SourceZymo80ntcleanup afterTaq; printed50ulcleanup despite3uladdedafter50ulPolImix, volumeconventionunassembled')
    x=op(g,'y_adapter',x,'SourcemethylatedYadapterligation20C1h',S.LIGATION,reagent='SourceSI methylatedadapterE7165/E7140 andM0367; mainE7536 identityunreconciled',product='Y-adapter/one-hairpin linked selected DNA plus remaining source reaction classes; exact adapterjunction unassembled')
    x=clean(g,'y_cleanup',x,'Source0.8xSPRI cleanup')
    x=op(g,'bisulfite',x,'SourceFastBisulfite95C25min, modifiedhistoricalAbcamkit',reagent='Actualsourcebisulfitereagent/kitmolecularsupplies;20ulDNAevaporatedto5ulbeforemixing, notDNA discard; mainab117127/SIab11727 identityunreviewed',product='Selected original unprotectedC→U, biological5mC retained; every syntheticcopyC protectedby5mC. Other chemistry/kit recovery unassembled')
    x=op(g,'pcr',x,'SourceQ5UM0597 PCR; SI8cycles',S.PCR,reagent='SourceQ5U andunassembled indexprimer/template geometry',product='PCR descendants plus original convertedtemplates/primers/remnants; C→U physicaltemplate andTreportedDNA are distinct')
    if targeted:
        x=selected(g,'target_allocation',x,'Actual200–300ngindexedlibrary contribution tosourcecapturepool')
        others=g.input('other_capture_libraries','Actualotherindexedlibrarycontributions tosourcecapturepool; amount/membershipunassigned',K.MOLECULAR)
        pooled=g.product('capture_pool','ActualpooledindexedlibrariesforTwistcapture',K.MOLECULAR)
        g.react('capture_pooling','Sourceactualcapturepool withdistinctlibrarycontributions',S.POOL,((x,R.TEMPLATE),(others,R.TEMPLATE)),(pooled,));x=pooled
        x=op(g,'target_capture',x,'Source200–300ngindexedlibrary pool thenTwistHumanCoreExome102025 enrichment',reagent='HistoricalTwistcapturekit; actual molecularselection/primergeometry unreviewed')
        x=op(g,'target_pcr',x,'SourceNEBNextQ5M0544+E7603 quantprimer mix enrichmentPCR',S.PCR,reagent='SourceQ5M0544; distinct from initialQ5U')
    x=selected(g,'loading',x,'ActualselectedsourceDNA-librarycontribution, volumes/poolingunassigned')
    g.sequence(x,SequencingPlan(acquisition='unspecified',limitation='SourcehumanNovaSeq6000 PE100 (bacterialcomparison NextSeq550PE150); exactlinkedread/hairpin/cyclegeometry unassembled. Cycle100trimming andjointC/Tvariant/methylation calls are computational. '+NAME_BOUNDARY))
    g.add_source_item(SourceItem('pending','HistoricalNEB adapter/endprep/UDG-EndoIV originals, AbcamFastBisulfiteguide andoptionalTwistcaptureguide','Sourceidentified vendor reagents not independentlyread within cap; full hairpin/Yprimer/read geometry unresolved; no absence claim.'))
    return g

def protocol_graph():return graph_for_scope()
def selftest():
    from checks import Check,Source,have
    from storage import sources_root
    c=Check();p=OriginalCopyPair('ACGCCATG',(3,))
    c('originalphysicalCconversion preserves selectedbiological5mC',p.converted_original(),'AUGCUATG')
    c('syntheticcopyCprotected independent oforiginalmethylpattern',p.converted_copy(),revcomp('ACGCCATG'))
    c('changingoriginalbiologicalmark doesnotchangepreservedsequencecopy',OriginalCopyPair(p.original,()).synthetic_copy,p.synthetic_copy)
    c('copiedDNAreadout containsT whilephysicaloriginalcontainsU','U' in p.converted_original() and 'U' not in p.complementary_dna_descendant(),True)
    g=protocol_graph();c('No cleanupbetweenincisionandPolI',g.steps['nick_translate'].inputs[0].material,'nick_product')
    c('No cleanupbetweenPolIandTaq',g.steps['a_tail'].inputs[0].material,'nick_translate_product')
    c('Nickstockusesactual39ulallocation',g.steps['nick'].inputs[0].material,'nick_allocation_selected')
    c('Hairpinligationusesselectedpreparedstock',g.steps['ligate_hairpin'].inputs[1].material,'hairpin_allocate_selected')
    c('PolInucleotidescontributephysically',g.steps['nick_translate'].inputs[2].material,'nick_translate_molecules')
    c('Actualtwosonicationscontributebeforeendprep',len(g.steps['sonication_collect'].inputs),2)
    c('Loadingretainsstockremainder',g.terminals['loading_remainder'].disposition,Disposition.UNRESOLVED)
    c('targetedbranch addsdistinctpostcapturePCR','target_pcr' in graph_for_scope(True).steps and 'target_pcr' not in g.steps,True)
    q=Source(sources_root()/'methyl-snp-seq__10.1101+gr.277080.122'/'gr.277080.122_PMC9808626.xml.txt','get_sources --slugs methyl-snp-seq__10.1101+gr.277080.122')
    if have(c,q,label='publishedhairpin'):
        c('publishedfullhairpinuracilsequence',HAIRPIN in q.path.read_text(),True)
    c.report()
if __name__=='__main__':selftest()
