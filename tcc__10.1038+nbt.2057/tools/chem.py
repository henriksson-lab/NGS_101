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

def op(g,key,x,label,kind=S.MAGIC,*,reagent='source reagents; geometry unassembled',product=None,state=None,catalyst=False):
    stock=g.input(key+'_stock',reagent,K.CATALYST if catalyst else K.CONDITION)
    out=g.product(key+'_product',product or label+'; molecular mixture/remnants unassembled',K.MOLECULAR,state=state)
    kw={'reason':'Source operation established; exact molecular ends, companion species and geometry unassembled.'} if kind in {S.MAGIC,S.GENERIC} else {}
    g.react(key,label,kind,((x,R.TEMPLATE),(stock,R.CATALYST if catalyst else R.CONDITION)),(out,),**kw);return out

def cleanup(g,key,x,label,output_kind=K.MOLECULAR):
    stock=g.input(key+'_condition',label,K.CONDITION)
    out=g.product(key+'_kept','Source recovered preparation; recovery and carryover unquantified',output_kind)
    waste=g.product(key+'_waste','Actual source-separated workup fractions',K.WASTE)
    g.react(key,label,S.CLEANUP,((x,R.TEMPLATE),(stock,R.CONDITION)),(out,waste));g.discard(waste,'Source physical cleanup, no perfect partition implied.');return out

def uncertain_ligation():return AdapterLigationGeometry('unavailable','unavailable',(scene('unreviewed exact ends and junction geometry'),),'Known source ligation, complete molecular geometry not assembled.')



def qc(g,key,x,label):
    keep=g.product(key+'_retained','Sourceproductionremainder; analyticalaliquotamountnotinvented',K.MOLECULAR);test=g.product(key+'_assay','Separateanalyticalaliquotforthe sourceassay',K.QC)
    g.react(key,label,S.SPLIT,((x,R.TEMPLATE),),(keep,test));g.mark_terminal(Terminal(test.id,Disposition.QC,'Analyticalmaterialdoesnotbecomethepreparativesequencingtemplate.'));return keep

def unknown_pcr(g,key,x,label,reagent):
    f=g.input(key+'_f','Sourceforwardindexingprimerfamily; memberandfullbodyunreviewed',K.OLIGO);r=g.input(key+'_r','Sourcereverseindexingprimerfamily; memberandfullbodyunreviewed',K.OLIGO);cat=g.input(key+'_enzyme',reagent,K.CATALYST);y=g.product(key+'_product','SourcePCRcopiesplusoriginaltemplates/unextendedprimersandreactionremnants; exactbackbonesunassembled',K.MOLECULAR)
    g.react(key,label,S.PCR,((x,R.TEMPLATE),(f,R.STRUCTURAL),(r,R.STRUCTURAL),(cat,R.CATALYST)),(y,),geometry=PCRGeometry.unavailable(scene('Sourcefulladapter/primer/endgeometryunassembled'),reason='KnownPCRoperationdoesnotestablishcompletefirstcopyorfinalstrandgeometry.'));return y

SOURCE='PublishedMOESM22CORRECTED3July2012 firstpagecorrectionnotice/Methodspp38–45 lines1578–1942;DOI10.1038/nbt.2057. AssociatedAMNIHMS512052overviewonly;uncorrectedhistoricalSIexcluded.'
@dataclass(frozen=True)
class MarkedStrand:
    sequence:str
    biotin_c:tuple=()
    ps_bonds:tuple=()
    inferred:bool=True
    def __post_init__(self):
        if set(self.sequence)-set('ACGT') or len(set(self.biotin_c))!=len(self.biotin_c) or any(type(i)is not int or not 0<=i<len(self.sequence) or self.sequence[i]!='C' for i in self.biotin_c):raise ValueError('Selectedbiotin-dCmustidentifyactualC')
        if len(set(self.ps_bonds))!=len(self.ps_bonds) or any(type(i)is not int or not 0<=i<len(self.sequence)-1 or self.sequence[i+1]!='G' for i in self.ps_bonds):raise ValueError('Source dGTPalphaS installsPSbond5primetoactualG')
    def selected_3prime_resection(self,keep):
        if type(keep)is not int or not 0<=keep<=len(self.sequence) or any(i>=keep-1 for i in self.ps_bonds):raise ValueError('Callerselected3primeresectioncannotcrossorcleaveprotectedPSbond')
        return (MarkedStrand(self.sequence[:keep],tuple(i for i in self.biotin_c if i<keep),self.ps_bonds),MarkedStrand(self.sequence[keep:],tuple(i-keep for i in self.biotin_c if i>=keep),()))
@dataclass(frozen=True)
class ProteinBiotin:
    protein_id:str
    cysteine_sites:tuple
    @property
    def dna_direct_biotin(self):return ()
    def proteolysis_inventory(self):return {'protein_id':self.protein_id,'peptide_biotin_sites':self.cysteine_sites,'dna_biotin_from_protein':()}
def bind(g,key,x,reagent,label):
    out=g.product(key+'_mixture','Sourcebindingmixturewithassociatedandunboundmaterial;notyetphysicallyseparated',K.MOLECULAR)
    g.react(key,label,S.BINDING,((x,R.TEMPLATE),(reagent,R.STRUCTURAL)),(out,));return out
def protocol_graph():
    g=ProtocolGraph('TCC — corrected2012HindIIIprotocol',authority=SOURCE,representation='published_instance')
    x=g.input('cells','Source25millionGM12878cells;noMboI/HiCcontrolrecipeimport',K.MOLECULAR,primary=True)
    x=op(g,'fix',x,'Source1%formaldehydeRT10min',S.CROSS_LINKING,reagent='Sourceformaldehyde/medium')
    x=op(g,'quench',x,'SourceglycineRT15min',S.QUENCHING,reagent='Sourceglycine')
    x=cleanup(g,'cell_pellet',x,'Sourcepellet/supernatantdiscard,flashfreeze/storage')
    x=op(g,'lyse',x,'Sourceonice15minlysis/Dounce20strokes',S.EXTRACTION,reagent='SourceHEPES/NaCl/IGEPAL/proteaseinhibitors')
    x=cleanup(g,'nuclei_wash',x,'Source2500gpellet/twoicecoldwashes')
    x=op(g,'solubilize',x,'SourceSDS65C10minthenRT',S.MAGIC,reagent='SourceSDS/washbuffer')
    ipb=g.input('ipb','SourceIodoacetyl-PEG2-biotin(IPB)protein-cysteinealkylationreagent',K.MOLECULAR);y=g.product('protein_labeled','Sourceprotein-cysteinebiotinylatedchromatinplusunreactedIPB/reagents;notdirectDNA-biotintransfer',K.MOLECULAR)
    label_condition=g.input('label_condition','SourceRT60minIPBalkylationconditions',K.CONDITION)
    g.react('protein_label','SourceIPBRT60minproteinbiotinylation',S.CHEMICAL_CONVERSION,((x,R.TEMPLATE),(ipb,R.STRUCTURAL),(label_condition,R.CONDITION)),(y,));x=op(g,'neutralize_sds',y,'SourceNEBuffer2ice5min/Tritonice10min37C10min',reagent='Sourcebuffer/Triton')
    x=op(g,'hindiii',x,'SourceHindIII/DTT37Covernightdigestion',S.DIGESTION,reagent='SourceHindIII/NEBuffer2/DTT',catalyst=True,product='Crosslinkedprotein-biotinDNAfragmentsplusuncutDNA/freeDNA/reagents;cleavagenotphysicalselection')
    x=cleanup(g,'dialysis',x,'Source20kDcutoffdialysis4hRT/bufferrenewal3h;excessIPBremovedwithoutperfectpurityclaim')
    aliquots=[]
    for i in range(5):aliquots.append(g.product('aliquot_'+str(i),'Sourceoneoffiveequaldialyzedaliquots,approximately500uleach',K.MOLECULAR))
    g.react('five_aliquots','Sourcefiveseparateequalaliquotsforparallelhandling',S.SPLIT,((x,R.TEMPLATE),),tuple(aliquots))
    t1=g.input('t1_stock','Source400ulMyOneT1physicalbeadstock',K.MOLECULAR);t1=cleanup(g,'prepare_t1',t1,'SourceT1stockthreePBSTwashesthenresuspend2mlbeforealiquotallocation')
    bead_parts=[g.product('t1_part_'+str(i),'Source400ulonefifthofwashed2mlT1suspension',K.MOLECULAR) for i in range(5)]
    g.react('allocate_t1','SourcefivealiquotsofpreviouslywashedT1suspension',S.SPLIT,((t1,R.TEMPLATE),),tuple(bead_parts))
    recovered=[]
    for i,a in enumerate(aliquots):
        k='a'+str(i)+'_';beads=bead_parts[i]
        z=bind(g,k+'tether',a,beads,'Sourceproteinbiotin-streptavidintether30minRT;unboundchromatinstillintube')
        blocker=g.input(k+'blocker','SourceIPBneutralizedequimolar2-mercaptoethanol;5ulperaliquot',K.MOLECULAR);z=bind(g,k+'block',z,blocker,'Source15minRTblocksremainingbeadbiotinsitesbeforeDNA-biotinfill')
        z=cleanup(g,k+'wash_tether',z,'SourceTPBS/washbufferphysicalremovalofunboundchromatin/reagents')
        z=op(g,k+'fill',z,'SourceKlenowRT40minbluntfillusingdATP/dTTP/dGTPalphaS/Biotin14dCTP',S.EXTENSION,reagent='SourceKlenow/SpdGTPalphaS/Biotin14dCTP;selectednewDNAmarkpositionsunassembled',catalyst=True,product='SelectedDNAcopieswithPSbond5primetoGandbiotinC3primetoPSplusoriginalDNA/uncopiedends/nucleotides/beads/proteins')
        z=op(g,k+'stop_fill',z,'SourceEDTAquench',S.QUENCHING,reagent='SourceEDTA');z=cleanup(g,k+'wash_fill',z,'Sourcetwowashes/resuspend500ul/transfernew15mltube')
        lig=g.input(k+'ligase','SourceT4DNAligase/ligasebuffer/BSA/Triton',K.CATALYST);out=g.product(k+'ligated','Sourcebead-tetheredligationjunctionclassesplusunligatedDNA/reagents;contact/yieldnotpredicted',K.MOLECULAR)
        g.react(k+'ligate','Sourcebluntligation16C4hwhileprotein-tethered',S.PROXIMITY_LIGATION,((z,R.TEMPLATE),(z,R.STRUCTURAL),(lig,R.CATALYST)),(out,),geometry=AdapterLigationGeometry('unavailable','unavailable',(scene('Sourceknownbluntfilledgenomicends/fullcontactjunctiongeometryunassembled'),),'Genomicbluntendsknown;noactualsequencingadapteratthisproximityligationstep,completegenomicend/productgeometryunassembled.'))
        z=op(g,k+'stop_lig',out,'SourceEDTAstopsbeforemagnetcollection',S.QUENCHING,reagent='SourceEDTA;preventsnewligationduringcollection');z=cleanup(g,k+'aspirate',z,'Sourceaspiratesupernatantafterligationquench')
        z=op(g,k+'pk',z,'SourcePK65CovernightthenadditionalPK65C2h',S.DEPROTEINIZATION,reagent='SourcePK/SDS/extractionbuffer',catalyst=True,product='ReleasedDNAplusprotein-biotinpeptides/beads/reagents;proteinbiotindoesnotbecomeDNAmark')
        z=cleanup(g,k+'collect_dna',z,'Sourcetransfersupernatantfrombeads,PCIextracttwice/chloroformonce')
        z=cleanup(g,k+'ethanol',z,'SourceEtOHminus20Covernight/pellet80%wash/resuspend20ul');recovered.append(z)
    x=g.product('pooled','Sourceallfiverecoveredaliquotspooled100ul',K.MOLECULAR);g.react('pool_five','Sourceactualpooloffiveparallelpreparations',S.POOL,tuple((z,R.TEMPLATE) for z in recovered),(x,))
    x=op(g,'rnase',x,'SourceRNaseA37C30min',S.DIGESTION,reagent='SourceRNaseA',catalyst=True,product='DNAplusdigestedRNAproductsuntilQIAquickcleanup');x=cleanup(g,'qiaquick',x,'SourceQIAquickelute50ul')
    x=qc(g,'quant',x,'Sourceconcentrationmeasurementanalyticalallocation')
    selected=g.product('exo_input','Source5ugpurifiedDNAforExoIII',K.MOLECULAR);rest=g.product('mass_remainder','UnallocatedpurifiedDNAremainder;amount/laterfateunassigned',K.MOLECULAR);g.react('allocate_exo','Source5ugallocationwithoutinventedtotalyield',S.SPLIT,((x,R.TEMPLATE),),(selected,rest));g.mark_terminal(Terminal(rest.id,Disposition.UNRESOLVED,'Sourceunallocatedmassnotimplicitlylost.'))
    x=op(g,'exo',selected,'Source300UExoIII37C1h3primetrimmingofnonligatedends',S.DIGESTION,reagent='SourceExoIII/NEBuffer1;PSprotectionselectednotperfectspecificity',catalyst=True,product='Selectedjunction-biotinDNAplusend-debiotinylatedDNA,biotinylatedmononucleotide/remnantproductsandreagents;noinstantphysicaldiscard')
    x=op(g,'stop_exo',x,'SourceEDTA/NaClthen70C20min',S.QUENCHING,reagent='SourceEDTA/NaCl/heat')
    x=op(g,'shear',x,'SourceCovarisS2100–500bp/180s',S.FRAGMENTATION,reagent='SourceCovarisconditions',product='ShearedDNAclassesplusallExoIIIremnantsandreagentsuntilQIAquick')
    x=cleanup(g,'qiaquick_shear',x,'SourceQIAquickafterCovariselute50ul')
    x=op(g,'end_repair',x,'SourceKlenow/T4DNApolymerase/PNK20C30min',S.END_REPAIR,reagent='Sourceenzymes/dNTPs',catalyst=True);x=cleanup(g,'qiaquick_repair',x,'SourceQIAquickelute40ul')
    x=op(g,'a_tail',x,'SourceKlenowexo-/dATP37C30min',S.TAILING,reagent='SourceKlenowexo-/dATP',catalyst=True);x=op(g,'stop_tail',x,'SourceEDTAquench',S.QUENCHING,reagent='SourceEDTA')
    beads=g.input('c1','Source10ulMyOneC1physicalbeadstock',K.MOLECULAR);beads=cleanup(g,'prepare_c1',beads,'SourceC1twowashesB&W/resuspend50ul2xB&W',output_kind=K.SELECTION)
    bound=g.product('junction_bound','SourceDNA-nucleotide-biotinassociatedC1complexespluspossiblebackground',K.MOLECULAR);no=g.product('junction_unbound','Sourceexcludedunboundcapturefraction',K.WASTE)
    g.react('capture_junction','SourceC1DNA-biotinbindingRT30minandmagnetcollection',S.AFFINITY_CAPTURE,((x,R.TEMPLATE),(beads,R.SELECTION)),(bound,no));g.discard(no,'Sourcephysicalboundfractionselection,notperfectjunctionclassification.')
    x=cleanup(g,'wash_junction',bound,'SourceB&W/TritonthenTriswash')
    adapter=g.input('adapter','SourcehistoricalIlluminapaired-endadapters3uM;exactmember/endgeometryunreviewed',K.OLIGO);lig=g.input('adapter_ligase','SourceT4ligase/rapidbufferRT20min',K.CATALYST);y=g.product('adapter_product','Sourceadapter-bearingbeadboundDNAplusfreeadapters/unligatedDNA/reagents',K.MOLECULAR)
    g.react('adapter_ligate','SourceadapterligationonC1beads',S.LIGATION,((x,R.TEMPLATE),(adapter,R.STRUCTURAL),(lig,R.CATALYST)),(y,),geometry=uncertain_ligation())
    x=op(g,'stop_adapter',y,'SourceEDTAquench',S.QUENCHING,reagent='SourceEDTA');x=cleanup(g,'wash_adapter',x,'SourceB&Wtwice/TEtwice,resuspend30ul')
    chosen=g.product('pcr_input','Source20ulbeadsuspensionforprepPCR',K.MOLECULAR);rest=g.product('bead_remainder','Sourceother10ulbeadsuspension;optionalpilot/laterfateunassigned',K.MOLECULAR);g.react('allocate_pcr','Source30ulpreparation20ulPCRallocation',S.SPLIT,((x,R.TEMPLATE),),(chosen,rest));g.mark_terminal(Terminal(rest.id,Disposition.UNRESOLVED,'Optionalpilot12/15/18cyclecomparisonoutsidechosenprinted15cycleprofile;notinferredconsumed.'))
    x=unknown_pcr(g,'pcr',chosen,'SourcePlatinumPfxprinted15cycleprofile;sourceoptionalpilotoptimizationnotuniversal15assumption','SourcePlatinumPfx/PEprimers/dNTP/MgSO4')
    x=cleanup(g,'gel',x,'Source2%gel350–550bpband/QIAgel/elute50ulEB')
    x=qc(g,'library_qpcr',x,'SourceqPCRquantificationaliquotandretainedlibrary')
    lib=g.product('library','Sourcesize-selectedTCCreadysequencinglibrary;completeadapter/readgeometryunassembled',K.LIBRARY);g.react('library_identity','Sourcefinalphysicalpreparationidentity',S.GENERIC,((x,R.TEMPLATE),),(lib,),changes_molecule=False)
    g.sequence(lib,SequencingPlan(limitation='CorrectedpublishedSIGAIIxpaired-end40bpeachend. ProteinIPBtetheringandnewDNA-biotincaptureareseparatechemicalhandles;computationalcontactclassificationnotphysicalDNAconversion/purity.'))
    g.add_source_item(SourceItem('pending','Finalmainapplicability/fullIlluminaoligos/end/PCR/readgeometryanddirectQuail/predecessororiginals','PublishedSIcorrected3July2012active;earlierincorrectexperimentalprotocolhistoricalexcluded. OtherMboI/non-tetheredHiCsourcebranchesnotimported.'))
    return g
def selftest():
    from checks import Check
    c=Check();a=MarkedStrand('AGCT',(2,),(0,));kept,debris=a.selected_3prime_resection(2)
    c('selectedExoIIIretainsphysicalDNAinventory',kept.sequence+debris.sequence,a.sequence)
    c('biotinreleasedontoDNAremnantnoterased',(kept.biotin_c,debris.biotin_c),((),(0,)))
    c('PSbondretained5primetoG',kept.ps_bonds,(0,))
    try:a.selected_3prime_resection(1);bad=False
    except ValueError:bad=True
    c('resectioncannotcrosssourcePSbarrier',bad,True)
    c('proteinIPBdoesnotbecomeDNA-biotin',ProteinBiotin('protein',('cys1',)).proteolysis_inventory()['dna_biotin_from_protein'],())
    g=protocol_graph();c('allfiveparallelaliquotsphysicallypooled',len(g.steps['pool_five'].inputs),5)
    c('ExoIIIdebrisreachesshearingbeforecleanup',g.steps['shear'].inputs[0].material,'stop_exo_product')
    c('twoindependentbeadpreparationsretainphysicalwashes','prepare_t1' in g.steps and 'prepare_c1' in g.steps,True)
    c('prepPCRdoesnotconsumeentirebeadstock',{i.material for i in g.steps['pcr'].inputs if i.role==R.TEMPLATE},{'pcr_input'})
    c('sourcegeometrylimitsremainprivate',bool(g.publication_problems()),True);c.report()
if __name__=='__main__':selftest()
