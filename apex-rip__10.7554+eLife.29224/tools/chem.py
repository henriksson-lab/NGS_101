"""Protocol-owned source-qualified partial population graph, no finished marker."""
from dataclasses import dataclass,replace
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Scene,Segment
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,SourceItem,SequencingPlan,Terminal,Disposition

def clean(g,key,x,procedure,retained='source-recovered preparation; exact molecular partition/carryover unquantified'):
    condition=g.input(key+'_condition',procedure,K.CONDITION)
    out=g.product(key+'_kept',retained,K.MOLECULAR)
    waste=g.product(key+'_separated','Source physically separated fractions; no perfect molecular recovery/removal claimed',K.WASTE)
    g.react(key,procedure,S.CLEANUP,((x,R.TEMPLATE),(condition,R.CONDITION)),(out,waste));g.discard(waste,'Actual source physical separation only.');return out

def op(g,key,x,label,kind,*,reagent=None,catalyst=False,product=None,reason='',changes_molecule=True,molecular=False):
    inputs=[(x,R.TEMPLATE)]
    if reagent:
        stock=g.input(key+'_stock',reagent,K.CATALYST if catalyst else K.MOLECULAR if molecular else K.CONDITION);inputs.append((stock,R.CATALYST if catalyst else R.STRUCTURAL if molecular else R.CONDITION))
    if molecular and kind in {S.CHEMICAL_CONVERSION,S.QUENCHING,S.CROSS_LINKING}:
        condition=g.input(key+'_conditions',label+'; source conditions, not an extra reagent dose',K.CONDITION);inputs.append((condition,R.CONDITION))
    out=g.product(key+'_product',product or label+'; reaction products/remnants/carryover remain unassigned',K.MOLECULAR)
    g.react(key,label,kind,inputs,(out,),reason=reason,changes_molecule=changes_molecule);return out
SOURCE='Kaewsapsak et al.2017 DOI10.7554/eLife.29224; publisher v1VOR MethodsPartsI–III516–540; SuppFile5construct/qPCRinventory. Selected optimizedBP-firstquench-crosslink route.'
@dataclass(frozen=True)
class CrosslinkedRNP:
    rna_id:str
    protein_ids:frozenset[str]
    biotinylated_proteins:frozenset[str]
    crosslinked_proteins:frozenset[str]
    def __post_init__(self):
        if not self.rna_id or not self.protein_ids or not self.biotinylated_proteins<=self.protein_ids or not self.crosslinked_proteins<=self.protein_ids:raise ValueError('explicit RNA/protein identities and source-qualified selected relationships required')
    def selected_capture_path(self):
        return bool(self.biotinylated_proteins&self.crosslinked_proteins)
    def proteolyzed(self):
        return {'rna_id':self.rna_id,'direct_rna_biotinylation':False,'possible_peptide_remnants':True,'original_protein_ids':self.protein_ids}

def protocol_graph():
    g=ProtocolGraph('APEX-RIP — optimized published route',authority=SOURCE,representation='published_instance')
    current=g.input('cells','Source stable targetedAPEX2expressingHEK293T cells; alternatives/compartments not physically pooled',K.BIOLOGICAL,primary=True)
    current=op(g,'bp',current,'BP500µM30min37°C',S.MAGIC,reagent='source BP/media',reason='Source exposure, distribution unassembled.',molecular=True)
    current=op(g,'protein_label',current,'H2O2final1mM1min initiates source local proteinbiotinylation',S.CHEMICAL_CONVERSION,reagent='Actual source H2O2 oxidant;APEX2 retainedininputcells',molecular=True,product='Selected biotinylated proteins plus RNA/other cell material; not a declaration that RNA is directly biotinylated')
    current=clean(g,'remove_label_medium',current,'Source aspirateBP/H2O2labelingmedium beforequench')
    current=op(g,'radical_quench',current,'Source azide-free Trolox/ascorbatequench1min',S.QUENCHING,reagent='Actual source10mMascorbate/5mMTrolox quenchstock(noazide)',molecular=True)
    current=clean(g,'remove_quench',current,'Source aspirateliquidphase after1minradicalquench')
    current=op(g,'crosslink',current,'Source first0.1percentformaldehydecrosslinkquench1min',S.CROSS_LINKING,reagent='Actual first5mlformaldehyde/ascorbate/Trolox crosslinkquenchsolution',molecular=True,product='Sourcecrosslinkingcellmixture; firstdose retaineduntilstatedaspiration')
    current=clean(g,'remove_first_crosslink',current,'Source aspiratefirstcrosslinkquenchsolution after1min')
    current=op(g,'crosslink_second',current,'Source freshcrosslinkquench9minroomtemperature',S.CROSS_LINKING,reagent='Actual fresh5ml0.1percentformaldehyde/ascorbate/Trolox solution',molecular=True,product='SourceRNA-proteincomplexes, labeled/unlabeledproteins andbackground; nofullcrosslinkyieldassertion')
    current=op(g,'glycine',current,'Glycinefinal125mM5min',S.QUENCHING,reagent='Actual sourceglycine final125mM',molecular=True)
    current=clean(g,'harvest',current,'Source2PBSwashes/scrape/pellet')
    current=op(g,'lyse',current,'RIPAlysis5min4°C and source sonication3x30s',S.MAGIC,reagent='source RIPA/quenchers',reason='Source extraction/sonication envelope, actual break distributions unassembled.')
    current=clean(g,'clarify',current,'Source15,000g5min4°C clarifiedlysate transfer')
    current=op(g,'dilute',current,'Source1mlNLBdilution',S.GENERIC,reagent='source NLB/quenchers',reason='Buffer handling, no new covalent transformation assigned.',changes_molecule=False,molecular=True)
    control=g.product('input_control','Source20%inputaliquot; separately processed control route unassembled',K.MOLECULAR);rest=g.product('ip_input','Remaining sourcepreparation',K.MOLECULAR)
    g.react('split','Allocate20%input separately',S.SPLIT,((current,R.TEMPLATE),),(control,rest));g.mark_terminal(Terminal(control.id,Disposition.UNRESOLVED,'Independent source inputcontrol downstream preparation not fully assembled; not pooled by normalization.'))
    bead=g.input('beads','Source50µlstreptavidinbeadslurry2h4°C;8specifiedwashes',K.SELECTION)
    bound=g.product('bound','Captured protein-biotin/RNAcomplexes plus background',K.MOLECULAR);waste=g.product('ip_waste','Actual source captureflowthrough/wash fractions',K.WASTE)
    g.react('capture','Source streptavidincapture/RIPA-highsalt-urea washes',S.AFFINITY_CAPTURE,((rest,R.TEMPLATE),(bead,R.SELECTION)),(bound,waste));g.discard(waste,'Actual source beadwashing fractions; no perfect specificity.')
    current=op(g,'proteolyze',bound,'ProteinaseK200µg42°C1h/55°C1h in sarcosidebuffer',S.DEPROTEINIZATION,reagent='source ProteinaseK',catalyst=True,product='Released RNA-containing mixture with possible crosslinked peptide remnants/protein fragments/beads; no direct transfer of proteinbiotin to RNA asserted')
    current=clean(g,'rna_clean1',current,'Source RNACleanXP purification')
    current=op(g,'dnase',current,'RQ1RNase-freeDNaseI5U37°C30min',S.DIGESTION,reagent='source DNaseI',catalyst=True,product='RNApreparation with DNA-digestion fragments retained until cleanup')
    current=clean(g,'rna_clean2',current,'Source RNACleanXP cleanup after DNase')
    qa=g.product('rna_qa','ActualsourceRNA Pico/NanoQA contribution; amountunassigned',K.QC)
    stock=g.product('rna_qa_stock','RNAstockafteranalyticalallocation',K.MOLECULAR)
    g.react('rna_qa_allocate','Source RNA concentration/integrity measurement beforeRiboZero',S.SPLIT,((current,R.TEMPLATE),),(qa,stock));g.mark_terminal(Terminal(qa.id,Disposition.QC,'SourceBioanalyzerRNAassay; noheat-coolbeforeloading, notreturnedtopreparation'))
    chosen=g.product('ribo_input','ActualsourceRiboZeroinputRNA, atmost17ul; selectedmass/volumeunassigned',K.MOLECULAR)
    remaining=g.product('rna_remainder','RNAstockoutsideselectedRiboZeroreaction; amount/presence/fateunassigned',K.MOLECULAR)
    g.react('ribo_input_allocate','Select actualRiboZerocontribution',S.SPLIT,((stock,R.TEMPLATE),),(chosen,remaining));g.mark_terminal(Terminal(remaining.id,Disposition.UNRESOLVED,'UntransferredRNAstock'));current=chosen
    current=op(g,'ribo_mix',current,'RiboZeroGoldremoval solution68°C10min/25°C5min',S.MAGIC,reagent='source mass-dependent1or2µlremovalsolution/20µlreaction',reason='Explicit incubation known, internal molecular recognition unassembled.',molecular=True)
    bead=g.input('ribo_beads','Actual32.5ulpreparedRiboZero beads fromsource90ulslurry wash/equilibration; preparationvolumeinterfaceunassembled, notanassumed57.5ulremainder;RT5min/50C5min',K.SELECTION)
    sup=g.product('ribodepleted','Source RNA-containingSUPERNATANT; exact recovery/carryover unquantified',K.MOLECULAR);pellet=g.product('ribo_pellet','Source companion beadfraction; exact RNA/remnant composition unresolved',K.MOLECULAR)
    g.react('ribo_partition','Source recoverSUPERNATANT after magneticbead incubation',S.SPLIT,((current,R.TEMPLATE),(bead,R.SELECTION)),(sup,pellet));g.mark_terminal(Terminal(pellet.id,Disposition.UNRESOLVED,'Source explicit supernatant recovery, companion beadfraction fate not asserted.'))
    current=clean(g,'ribo_cleanup',sup,'Source50ulRNACleanXP cleanup of50ulRNA; beforeelution',retained='SourceRNA-associatedRNACleanbeadmaterial beforeEPFelution, notwater-elutedRNA')
    epf=g.input('epf_mix','Actual19.5ulTruSeqv2 ElutePrimeFragmentmix; innerformulationunassembled',K.MOLECULAR)
    eluted=g.product('epf_rna','SourceEPF-elutedRNA preparation',K.MOLECULAR);residual=g.product('epf_bead_residual','Residualbeadassociatedmaterial; exactcomposition/fateunassigned',K.MOLECULAR)
    g.react('epf_elution','Source elutionintoEPFmix',S.MAGIC,((current,R.TEMPLATE),(epf,R.STRUCTURAL)),(eluted,residual),changes_molecule=False,reason='Actualsourceelutionmix/volume, notanaddedfreshknownRTprimerrecipe.');g.mark_terminal(Terminal(residual.id,Disposition.UNRESOLVED,'Residualcleanupmaterial'));current=eluted
    kit=g.input('truseq','ActualTruSeqRNA sampleprepv2 molecularsupplies fromIncubateRFPstep; uniqueindex/member/oligos/packagingunreviewed',K.MOLECULAR)
    lib=g.product('library','Source uniquely indexedRNA-derived library; exact strand/end/read geometry unassembled',K.LIBRARY)
    g.react('library_prep','Source historicalTruSeqv2kit handoff',S.MAGIC,((current,R.TEMPLATE),(kit,R.STRUCTURAL)),(lib,),reason='Namedkit/version/startpoint known, no internal chemistry or modern substitute inferred.')
    qa=g.product('dna_qa','SourceDNA-HSBioanalyzer libraryQA aliquot; amountunassigned',K.QC)
    stock=g.product('post_qa_library','Sourcepreparedindexedlibrary stockafterQAallocation',K.LIBRARY)
    g.react('dna_qa_allocate','Allocate actualsourceDNA-HS libraryQA',S.SPLIT,((lib,R.TEMPLATE),),(qa,stock));g.mark_terminal(Terminal(qa.id,Disposition.QC,'Actualanalyticalmaterial, notreturnedtopool'))
    selected=g.product('pool_aliquot','Sourcechosenlibrarycontribution forequimolarpool; volume/massunassigned',K.LIBRARY)
    remainder=g.product('library_remainder','Untransferredindexedlibrarystock; amount/presence/fateunassigned',K.LIBRARY)
    g.react('pool_allocate','Select actualequimolarlibrarycontribution',S.SPLIT,((stock,R.TEMPLATE),),(selected,remainder));g.mark_terminal(Terminal(remainder.id,Disposition.UNRESOLVED,'Library stockoutsidepool'))
    others=g.input('other_libraries','Actualselectedcontributionsfromothersourceseparatelypreparedindexedlibraries; equimolar<=10perpool, upstreamstocksunassembled',K.LIBRARY)
    pool=g.product('pool','Source equimolar pooledlibrary preparation',K.LIBRARY);g.react('pool_libraries','Source librarypooling',S.POOL,((selected,R.TEMPLATE),(others,R.TEMPLATE)),(pool,))
    loaded=g.product('loaded_pool','Actualsourcepoolloadingcontribution fortwoHiSeq2500rapidlanes; laneallocationamountsunassigned',K.LIBRARY)
    rem=g.product('pool_remainder','Untransferredpoolstock; amount/presence/fateunassigned',K.LIBRARY)
    g.react('loading_allocate','Select actualsourcepoolloadingcontribution',S.SPLIT,((pool,R.TEMPLATE),),(loaded,rem));g.mark_terminal(Terminal(rem.id,Disposition.UNRESOLVED,'Poolstockoutsideselectedacquisition'))
    g.sequence(loaded,SequencingPlan(acquired_roles=('Read 1','Read 2'),read_lengths={'Read 1':50,'Read 2':50},limitation='Source HiSeq2500rapidPE50thenindex; exact indexcycles/architecture and full primergeometry unassembled.'))
    g.declare_reviewed_scope(('cells','bp_stock','bp','protein_label_stock','protein_label_conditions','protein_label','remove_label_medium_condition','remove_label_medium','radical_quench_stock','radical_quench_conditions','radical_quench','remove_quench_condition','remove_quench','crosslink_stock','crosslink_conditions','crosslink','remove_first_crosslink_condition','remove_first_crosslink','crosslink_second_stock','crosslink_second_conditions','crosslink_second','glycine_stock','glycine_conditions','glycine','harvest_condition','harvest','lyse_stock','lyse','clarify_condition','clarify','dilute_stock','dilute','split','beads','capture','proteolyze_stock','proteolyze','rna_clean1_condition','rna_clean1','dnase_stock','dnase','rna_clean2_condition','rna_clean2','rna_qa_allocate','ribo_input_allocate','ribo_mix_stock','ribo_mix','ribo_beads','ribo_partition','ribo_cleanup_condition','ribo_cleanup','epf_mix','epf_elution','dna_qa_allocate','pool_allocate','other_libraries','pool_libraries','loading_allocate'),SOURCE)
    g.add_source_item(SourceItem('pending_library','Historical TruSeqv2/RiboZerokit originals and full polymer/read assembly','Published Methods536/540'))
    return g

def selftest():
    from checks import Check
    c=Check();x=CrosslinkedRNP('rna1',frozenset({'p1','p2'}),frozenset({'p1'}),frozenset({'p1'}));p=x.proteolyzed()
    c('selected protein-label/crosslink path can recover linked RNA',x.selected_capture_path(),True)
    c('uncrosslinked labeled protein does not supply declared RNAcapture path',replace(x,crosslinked_proteins=frozenset({'p2'})).selected_capture_path(),False)
    c('protein digestion preserves RNAidentity',p['rna_id'],'rna1')
    c('proteinborne biotin is not silently transferred onto RNA',p['direct_rna_biotinylation'],False)
    c('peptide/remnantfates remain possible after proteolysis',p['possible_peptide_remnants'],True)
    c('full sourceknown route retains publicationboundary',bool(protocol_graph().publication_problems()),True);c('Freshcrosslinkdose follows firstsolutionremoval',protocol_graph().steps['crosslink_second'].inputs[0].material,'remove_first_crosslink_kept');c('Pool consumesselectedaliquotnotwholelibrary',protocol_graph().steps['pool_libraries'].inputs[0].material,'pool_aliquot');c.report()
if __name__=='__main__':selftest()
