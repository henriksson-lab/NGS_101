"""Protocol-owned source-qualified incomplete graph; no publication marker."""
from dataclasses import dataclass, replace
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Scene,Segment
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,SourceItem,SequencingPlan,Terminal,Disposition
from adapter_ligation import AdapterLigationGeometry
from pcr import PCRGeometry

def scene(label):
    sc=Scene();sc.strand('preparation',[Segment('population','['+label+']',placeholder=True)]);return sc

def op(g,key,x,label,kind=S.MAGIC,*,reagent='source reagents; geometry unassembled',product=None,state=None,catalyst=False,extra_inputs=()):
    stock=g.input(key+'_stock',reagent,K.CATALYST if catalyst else K.CONDITION)
    out=g.product(key+'_product',product or label+'; molecular mixture/remnants unassembled',K.MOLECULAR,state=state)
    kw={'reason':'Source operation established; exact molecular ends, companion species and geometry unassembled.'} if kind in {S.MAGIC,S.GENERIC} else {}
    g.react(key,label,kind,((x,R.TEMPLATE),(stock,R.CATALYST if catalyst else R.CONDITION))+tuple(extra_inputs),(out,),**kw);return out

def cleanup(g,key,x,label):
    stock=g.input(key+'_condition',label,K.CONDITION)
    out=g.product(key+'_kept','Source recovered preparation; recovery and carryover unquantified',K.MOLECULAR)
    waste=g.product(key+'_waste','Actual source-separated workup fractions',K.WASTE)
    g.react(key,label,S.CLEANUP,((x,R.TEMPLATE),(stock,R.CONDITION)),(out,waste));g.discard(waste,'Source physical cleanup, no perfect partition implied.');return out

def uncertain_ligation():return AdapterLigationGeometry('unavailable','unavailable',(scene('unreviewed exact ends and junction geometry'),),'Known source ligation, complete molecular geometry not assembled.')



def qc(g,key,x,label):
    keep=g.product(key+'_retained','Sourceproductionremainder; analyticalaliquotamountnotinvented',K.MOLECULAR);test=g.product(key+'_assay','Separateanalyticalaliquotforthe sourceassay',K.QC)
    g.react(key,label,S.SPLIT,((x,R.TEMPLATE),),(keep,test));g.mark_terminal(Terminal(test.id,Disposition.QC,'Analyticalmaterialdoesnotbecomethepreparativesequencingtemplate.'));return keep

def unknown_pcr(g,key,x,label,reagent):
    f=g.input(key+'_f','Sourceforwardindexingprimerfamily; memberandfullbodyunreviewed',K.OLIGO);r=g.input(key+'_r','Sourcereverseindexingprimerfamily; memberandfullbodyunreviewed',K.OLIGO);cat=g.input(key+'_enzyme',reagent,K.CATALYST);y=g.product(key+'_product','SourcePCRcopiesplusoriginaltemplates/unextendedprimersandreactionremnants; exactbackbonesunassembled',K.MOLECULAR)
    g.react(key,label,S.PCR,((x,R.TEMPLATE),(f,R.STRUCTURAL),(r,R.STRUCTURAL),(cat,R.CATALYST)),(y,),geometry=PCRGeometry.unavailable(scene('Sourcefulladapter/primer/endgeometryunassembled'),reason='KnownPCRoperationdoesnotestablishcompletefirstcopyorfinalstrandgeometry.'));return y

SOURCE='AssociatedPMC3679896 Methods136–156/Fig1 andpublishedSI MOESM327S14lines481–493; DOI10.1038/ncomms2527. SourcegenomicTAmCcaptureandHiScanSE51route.'
@dataclass(frozen=True)
class SelectedLabelMap:
    dna:str
    marks:tuple
    def __post_init__(self):
        if not self.dna or set(self.dna)-set('ACGT') or len(dict(self.marks))!=len(self.marks) or any(type(i)is not int or not 0<=i<len(self.dna) or self.dna[i]!='C' or m not in {'5mC','5hmC','5gmC','N3-5gmC','SSbiotin-label','thiol-remnant'} for i,m in self.marks):raise ValueError('explicitselectedcytosinesandcovalentclassesrequired')
    def block_existing_hmc(self):return replace(self,marks=tuple((i,'5gmC' if m=='5hmC' else m) for i,m in self.marks))
    def oxidize_and_trap_selected(self):
        if any(m=='5hmC' for _,m in self.marks):raise ValueError('selectedoriginal5hmCmustbeblockedbefore5mCselectivelabeling')
        return replace(self,marks=tuple((i,'N3-5gmC' if m=='5mC' else m) for i,m in self.marks))
    def click_selected(self):return replace(self,marks=tuple((i,'SSbiotin-label' if m=='N3-5gmC' else m) for i,m in self.marks))
    def reduce_selected(self):
        if not any(m=='SSbiotin-label' for _,m in self.marks):raise ValueError('disulfidebiotinsubstraterequired')
        return replace(self,marks=tuple((i,'thiol-remnant' if m=='SSbiotin-label' else m) for i,m in self.marks))

def protocol_graph():
    g=ProtocolGraph('TAmC-Seq',authority=SOURCE,representation='published_instance')
    x=g.input('genomic','SourcesonicatedmESCgenomicDNA; original5mC/5hmCsitesandoccupancyunassigned',K.MOLECULAR,primary=True)
    gt=g.input('block_gt','Actual source betaGT40pmol preblock contribution',K.CATALYST)
    glucose=g.input('block_glucose','Actual source ordinaryUDPglucose300uM preblock contribution',K.MOLECULAR)
    x=op(g,'block',x,'Source3µgDNA/betaGT40pmol/UDPglucose300µM37°C1h',S.CHEMICAL_CONVERSION,reagent='SourceHEPES/MgCl2 and37C1h preblock programme',product='Selectedpreexisting5hmC→ordinaryglucosyl5hmC, original5mC/otherDNA/reagentmaterialretained',extra_inputs=((gt,R.CATALYST),(glucose,R.STRUCTURAL)))
    x=cleanup(g,'block_clean',x,'SourceQiagenDNApurification')
    used=g.product('trap_input','Source2µgtreatedDNAallocatedforcoupledtrap',K.MOLECULAR);rest=g.product('protected_remainder','Unallocatedprotectedpreparation; exactmassafterpurification/recoveryunassigned',K.MOLECULAR);g.react('allocate','Source2µgprotectedDNAreactioninput',S.SPLIT,((x,R.TEMPLATE),),(used,rest));g.mark_terminal(Terminal(rest.id,Disposition.UNRESOLVED,'Notassumedexactly1µg: recoveryisnot100%guaranteed.'))
    tet=g.input('trap_tet','Actual source mTet1 160pmol contribution to coupled reaction',K.CATALYST)
    gttrap=g.input('trap_gt','Actual source betaGT80pmol contribution to same coupled reaction',K.CATALYST)
    sugar=g.input('trap_sugar','Actual source UDP6N3glucose150uM contribution',K.MOLECULAR)
    x=op(g,'oxidize_trap',used,'Source160pmolmTet1/80pmolbetaGT/150µMUDP6N3glucose37°C1h',S.CHEMICAL_CONVERSION,reagent='SourceHEPES/FeII/ascorbate/alphaKG/MgCl2/DTT mixture',product='Selectedoriginal5mC→new5hmCsimultaneouslytrappedasN3-5gmC; preexisting5gmC/otherDNA/offproductsandenzymereagentsremain. NosequentialTET-onlybranchsubstitution',extra_inputs=((tet,R.CATALYST),(gttrap,R.CATALYST),(sugar,R.STRUCTURAL)))
    x=cleanup(g,'trap_clean',x,'SourceDNApurificationbeforeclick')
    linker=g.input('click_linker','Actual source DBCOdisulfidebiotin150uM contribution in40ul click reaction',K.MOLECULAR)
    x=op(g,'click',x,'Source150µMdisulfideDBCO-biotin37°C2h40µl',S.CHEMICAL_CONVERSION,reagent='Source37C2h click programme',product='SelectedSSbiotinN3-5gmCDNAplusunreactedDNA/freeprobeandreagents',extra_inputs=((linker,R.STRUCTURAL),))
    x=cleanup(g,'click_clean',x,'SourceQiagennucleotideremovalkit/H2O50µl')
    beads=g.input('beads','SourceMyOneStreptavidinC1captureperhistoricalmanual',K.SELECTION);bound=g.product('bound','SourcebeadDNAcomplexpluspossiblenonspecificmaterial',K.MOLECULAR);other=g.product('flowthrough','Actualcapture/washfraction; exactcomposition/fateunreviewed',K.MOLECULAR);g.react('capture','SourcebiotinDNApulldown',S.AFFINITY_CAPTURE,((x,R.TEMPLATE),(beads,R.SELECTION)),(bound,other));g.mark_terminal(Terminal(other.id,Disposition.UNRESOLVED,'Captureseparationestablished; completewash/recoveryinventoryunreviewed.'))
    dtt=g.input('dtt_reductant','Actual local source DTT100mM contribution in50ul release solution',K.MOLECULAR)
    x=op(g,'dtt',bound,'Source100mMDTT50µlRT2h',S.CHEMICAL_CONVERSION,reagent='Source roomtemperature2h release programme',product='SelectedreleasedDNAwithcovalentthiollinkerremnantplusbead/biotinDTTproducts; noprematurephysicaldiscard',extra_inputs=((dtt,R.STRUCTURAL),))
    cond=g.input('separate','Sourcesupernatantcollectionafterbeadrelease; exactseparatorunreviewed',K.CONDITION);released=g.product('released','SourcecollectedDNAsupernatant; original5mC-derivedtagremnantretained',K.MOLECULAR);beadrest=g.product('retained_beads','Physicallyseparatedbeads/linkerproducts; retainedDNAfractionunquantified',K.MOLECULAR);g.react('collect','SourcecollectreleasedDNAsupernatant',S.CLEANUP,((x,R.TEMPLATE),(cond,R.CONDITION)),(released,beadrest));g.mark_terminal(Terminal(beadrest.id,Disposition.UNRESOLVED,'Actualsupernatantseparation; nofurtherbeadfateassigned.'));x=released
    x=cleanup(g,'desalt',x,'SourceMicroBioSpin6DTTremoval')
    x=cleanup(g,'minelute',x,'SourceQiagenMinElutePCRpurification/EB')
    selected=g.product('library_input','Source selected25ng enrichedDNA contribution',K.MOLECULAR)
    remainder=g.product('library_input_remainder','Untransferred capturedDNA stock; amount unassigned',K.MOLECULAR)
    g.react('library_allocate','Allocate actual25ng to library, not whole captured stock',S.SPLIT,((x,R.TEMPLATE),),(selected,remainder),changes_molecule=False)
    g.mark_terminal(Terminal(remainder.id,Disposition.UNRESOLVED,'Untransferred captured stock; recovery not normalized.'));x=selected
    x=op(g,'library_endrepair',x,'Source25ngenrichedDNAhistoricalIlluminagenomicChIPlibraryendprep',S.END_REPAIR,reagent='SourcehistoricalIlluminaendprepenzymesunreviewed',catalyst=True)
    x=op(g,'atail',x,'Sourceadenylation',S.TAILING,reagent='SourcehistoricalIlluminaadenylationenzymesunreviewed',catalyst=True)
    adapter=g.input('adapter','SourceIlluminaGenomicDNAadapteroligomix; fullbody/endgeometryunreviewed',K.ADAPTER);cat=g.input('ligase','SourcehistoricalIlluminaligase',K.CATALYST);y=g.product('ligated','Sourceadapter-ligatedenrichedDNAplusunligatedDNA/adaptersandreagents',K.MOLECULAR);g.react('ligate','Sourcegenomicadapterligation',S.LIGATION,((x,R.TEMPLATE),(adapter,R.STRUCTURAL),(cat,R.CATALYST)),(y,),geometry=uncertain_ligation())
    x=cleanup(g,'size',y,'Source2%agarose~200–350bpfragmentgelpurification')
    x=unknown_pcr(g,'pcr',x,'Source18PCRcycles','SourcehistoricalIlluminagenomicChIP-libraryPCRenzymesunreviewed')
    x=qc(g,'bioanalyzer',x,'SourceBioanalyzerDNA1000quality/quantificationaliquot')
    load=g.product('load_selected','Selected production library contribution for actual loading; volume unassigned',K.MOLECULAR)
    rest=g.product('load_remainder','Untransferred production stock; volume unassigned',K.MOLECULAR)
    g.react('load_allocate','Select loading contribution after actualanalyticalQC',S.SPLIT,((x,R.TEMPLATE),),(load,rest),changes_molecule=False)
    g.mark_terminal(Terminal(rest.id,Disposition.UNRESOLVED,'Untransferred production stock; no whole-stock loading assumption.'))
    phix=g.input('phix_spike','Actual source1percent PhiX loadingcontrol library contribution; full control polymer and stock preparation unreviewed',K.LIBRARY)
    pooled=g.product('loading_pool','Actual genomic library/PhiX loading mixture; no covalent fusion or genomic barcode assignment',K.MOLECULAR)
    g.react('phix_pool','Actual1percentPhiX spike in samplelane; dedicatedPhiX lane remains a distinct source control outside selected lane',S.POOL,((load,R.TEMPLATE),(phix,R.TEMPLATE)),(pooled,));x=pooled
    g.sequence(x,SequencingPlan(limitation='ActualHiScan/TruSeqclusterkitv2cBotHS/SBSkitv3HSsingle-end51bp. Separateunenrichedinputand3capturelibrariesreported, notassumedduplicatedDNA. PhiXcontrol lane/1%loading spike separatefromlibraryconstruct. Regionenrichmentisnotbaseconversion; exacthistoricalprimer/endgeometryunreviewed.'))
    g.add_source_item(SourceItem('pending','HistoricalC1capturemanualandIlluminagenomicChIPlibraryprotocol','Mainref26hmC-SeallabelingcomponentexplicitlycitedbutTAmCSIdisulfidelinker/DTT100mMlocalvaluescontrol; noexternalrecipeimportorequalfullprotocolclaim.'))
    return g

def selftest():
    from checks import Check
    c=Check();original=SelectedLabelMap('ACCGC',((1,'5mC'),(2,'5hmC')));blocked=original.block_existing_hmc();tag=blocked.oxidize_and_trap_selected().click_selected();released=tag.reduce_selected()
    c('preexisting5hmCstaysordinaryglucoseblocked',dict(tag.marks)[2],'5gmC')
    c('affinitytagreportsoriginal5mCclass',dict(tag.marks)[1],'SSbiotin-label')
    c('labelingdoesnotchangecytosinenucleotideinventory',released.dna,original.dna)
    c('DTTretainsoriginal5mC-derivedthiolremnant',dict(released.marks)[1],'thiol-remnant')
    try:original.oxidize_and_trap_selected();bad=False
    except ValueError:bad=True
    c('preexistinghmcnotassumedselectivelyexcludedwithoutblock',bad,True)
    g=protocol_graph();c('physicalsupernatantcollectionfollowsDTT',g.steps['collect'].inputs[0].material,'dtt_product')
    c('analyticalBioanalyzermaterialnotsequenced',g.steps['load_allocate'].inputs[0].material,'bioanalyzer_retained')
    c('Library receives actual selected capturedDNA mass',g.steps['library_endrepair'].inputs[0].material,'library_input')
    c('Coupled trap includes both source enzyme contributions',tuple(i.material for i in g.steps['oxidize_trap'].inputs[2:4]),('trap_tet','trap_gt'))
    c('Source PhiX is a physical loading control, not a covalent genomic extension',g.steps['phix_pool'].kind,S.POOL)
    c('Original historical source scope remains uncovered','pending' in g.source_coverage,False)
    c('unassembledhistoricalgeometrypreventspromotion',bool(g.publication_problems()),True);c.report()
if __name__=='__main__':selftest()
