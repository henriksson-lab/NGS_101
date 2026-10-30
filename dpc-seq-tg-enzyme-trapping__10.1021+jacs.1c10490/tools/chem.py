"""Private Tg DPC-Seq, source-qualified associated-manuscript preparation."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Scene,Segment
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,SourceItem,SequencingPlan,Terminal,Disposition
from adapter_ligation import AdapterLigationGeometry
TITLE='DPC-Seq thymidine-glycol enzyme trapping'
AUTHORITY='AssociatedAM PMC8755629 selected84/98/108/110–116/132–134; onlyonecurrentoriginal. PublishedACSmain403/associatedSIgate; no final-AMequivalenceorSIsequenceclaim.'
def unknown_scene():
    sc=Scene();sc.strand('unassembled',[Segment('unassembled','[source molecular geometry unassembled]',placeholder=True)]);return sc
def protocol_graph():
    g=ProtocolGraph(TITLE,authority=AUTHORITY,representation='published_instance')
    x=g.input('dna','SourceextractedHEK293TgenomicDNA;MonarchT3010/0.1mMBHTextractionoutsidegraph; WT/DKOseparateexperimentalinputs',K.MOLECULAR,primary=True)
    def stock(key,label,kind=K.CATALYST,role=R.CATALYST,**kw):return(g.input(key,label,kind,**kw),role)
    def cond(key,label):return stock(key,label,K.CONDITION,R.CONDITION)
    def op(key,label,kind,inputs=(),product='Source mixture with residualsubstrates/offproducts/freereagents retained untilphysicalworkup',**kw):
        nonlocal x
        if kind==S.GENERIC:kw['changes_molecule']=False
        if kind==S.MAGIC:kw['reason']='Source-specific historicalkitfinish/wholemoleculargeometryunassembled.'
        y=g.product(key+'_product',product,K.MOLECULAR);g.react(key,label,kind,((x,R.TEMPLATE),)+tuple(inputs),(y,),**kw);x=y
    def clean(key,label):
        nonlocal x
        c=g.input(key+'_condition',label,K.CONDITION);a=g.product(key+'_kept','Sourcephysicallyretainedpreparation; recovery/carryoverunassigned',K.MOLECULAR);w=g.product(key+'_waste','Separatedsourceworkupfractions',K.WASTE);g.react(key,label,S.CLEANUP,((x,R.TEMPLATE),(c,R.CONDITION)),(a,w));g.discard(w,'Sourcephysicalworkup,noperfectmolecularspecificity.');x=a
    def aliquot(key,label,rest):
        nonlocal x
        a=g.product(key+'_used',label,K.MOLECULAR);b=g.product(key+'_rest',rest,K.MOLECULAR);g.react(key,label,S.SPLIT,((x,R.TEMPLATE),),(a,b));g.mark_terminal(Terminal(b.id,Disposition.UNRESOLVED,rest));x=a
    op('fragment','SourceNEBNextdsDNAFragmentaseM0348,target100–300bp;sonicationexplicitlyavoided',S.FRAGMENTATION,(stock('fragmentase','SourceFragmentase'),))
    clean('fragment_clean','SourceNEBNextbeadpurificationand sizeselection; detailedfractionsunassembled')
    aliquot('dna_aliquot','Source10ugfragmentedDNAselected','UntransferredfragmentedDNAstock,mass/laterfateunassigned')
    op('pretreat','SourceSMUG1/Fpg/EndoIV25mMTris7.5/50mMNaCl/5mMMgCl2/0.5mMDTT37C60min',S.DIGESTION,(stock('smug1','SourceSMUG1'),stock('fpg','SourceFpg'),stock('endoiv','SourceEndoIV')),product='SourcepretreatedDNAwithprocessedinterferinglesion/APclassesandresidualsubstrates;validationenzymevolumesnotimported, no perfectTgspecificity')
    clean('pretreat_clean','SourceNEBNextSamplePurificationBeadsafterpretreatment')
    nth=stock('nthl1','Source0.8uMHis6-NTHL1in40ul: covalentlyincorporatedproteincontributor,notanewlateenzymestock',K.MOLECULAR,R.STRUCTURAL,chemical_state={'protein':'NTHL1-CterminalHis6','source_fate':'reductively trapped protein-DNA conjugates; exactDNAend/bondgeometryunassembled'})
    op('trap','Source50mMTris7/1mMEDTA/1mMDTT/10mMNaBH3CN37C60min',S.CHEMICAL_CONVERSION,(nth,stock('reductant','Source10mMNaBH3CNmolecularreductant; exactincorporation/bondrelationunassembled',K.MOLECULAR,R.STRUCTURAL),cond('trap_condition','SourceTris/EDTA/DTTreactionbuffer')),product='SourceNTHL1-DNAconjugates, freeNTHL1anduntrappedDNA; originalTgconvertedbytheassay, notnativeendogenousDPCcapture')
    resin=g.input('ni_resin','Source20ulHisPurNi-NTAresin',K.SELECTION);block=g.input('salmon_blocker','Source50ugsssalmon-spermDNA;incomingmolecularblockingstockcancontributebackground',K.MOLECULAR);cc=g.input('capture_condition','Source800ulPBS4Covernight',K.CONDITION);bound=g.product('bound','SourceNi-NTA-boundHisprotein/DPCsplusnonspecific/blockingDNAclasses',K.MOLECULAR);unbound=g.product('capture_supernatant','Source6000gsupernatant',K.WASTE)
    g.react('capture','SourceNi-NTAenrichment',S.AFFINITY_CAPTURE,((x,R.TEMPLATE),(resin,R.SELECTION),(block,R.STRUCTURAL),(cc,R.CONDITION)),(bound,unbound));g.discard(unbound,'Sourceexplicitlydiscardsfirstcapturesupernatant.');x=bound
    clean('capture_wash','Sourcetwice1mlPBSandtwice1ml20mMTris7.5/0.2percentTween;retainresin')
    op('endprep','ON-RESINadd50ulwater/3ulEndPrep/7ulbufferE7103S;20C1h65C30min1200rpmintermixing',S.END_REPAIR,(stock('endprep_enzyme','SourceNEBE7103EndPrep'),),product='Sourceend-preparedDNAonresinwithretainedNTHL1conjugates/freeprotein; wholeendgeometryunassembled')
    adapter=stock('adapter','Source2.5ulNEBNextAdaptorE7103S;historicalwholeordering/stockendsunassembled',K.ADAPTER,R.STRUCTURAL)
    op('ligate','ON-RESINsource30ulLigationMasterMix+2.5uladapter+1ulenhancer20C1h,flip10min',S.LIGATION,(adapter,stock('ligase','SourceNEBligationmastermix')),geometry=AdapterLigationGeometry('unavailable','unavailable',(unknown_scene(),),'Actualhistoricalon-beadligation;adapter/substrateendandsealedproductgeometryunassembled,notconfirmedmissing.'))
    op('user','Source3ulUSER37C30minonresin',S.DIGESTION,(stock('user_enzyme','SourceUSER E7500S'),),product='SourceUSER-processedDNA/adapterclassesandretainedcrosslinkprotein/resin;exactcleavage/endpaths unassembled')
    clean('postligation_wash','SourcetwicePBSandtwicewashingbufferafteron-resinligation/USER')
    op('imidazole','Source20ulPBS500mMimidazole37C30min1200rpmintermixing',S.GENERIC,(cond('imidazole_condition','Sourceimidazoleelution'),))
    eluate=g.product('imidazole_eluate','SourcecollectedDPC-containingeluate',K.MOLECULAR);spent=g.product('spent_resin','Sourceco-separatedresin/residualboundmaterial',K.MOLECULAR);g.react('elute_separate','Physicalsourceeluate/resin separation',S.SPLIT,((x,R.TEMPLATE),),(eluate,spent));g.mark_terminal(Terminal(spent.id,Disposition.UNRESOLVED,'Sourceunreleasedmaterial/resinlaterfateunassigned.'));x=eluate
    op('page_transfer','Source10percentSDS-PAGE90V60min→nitrocelluloseblot;local108procedureexplicitlyinheritedby116',S.GENERIC,(cond('page_condition','SourcePAGE/membranetransfer'),))
    aliquot('band','Sourceexcised50–250kDaDPCbandincut1–2mmnitrocelluloseslices','Othergel/membraneregions anduntransferredmaterial;laterfateunassigned')
    clean('membrane_wash','Sourcewashselectedslices500ul30mMTris8;retainmembraneslices')
    op('proteinase','Source20ulProteinaseK+130ul30mMTris8,55Covernightthen70C20min1200rpmintermixing',S.DEPROTEINIZATION,(stock('proteinase_k','SourceProteinaseK'),),product='SourceDNAreleasedwithproteinfragments/peptide-DNAremnantsandmembraneslices; proteolysisdoesnotchemicallyrestorepristineDNA')
    released=g.product('released_liquid','Source150ulliquidtransferredtoDNALoBindtube',K.MOLECULAR);mem=g.product('residual_membrane','Sourceuntransferredmembrane/residualDNA',K.MOLECULAR);g.react('transfer_150','Physical150ulrelease-fractiontransfer',S.SPLIT,((x,R.TEMPLATE),),(released,mem));g.mark_terminal(Terminal(mem.id,Disposition.UNRESOLVED,'NoassumedcompleteDNAreleasefrommembrane.'));x=released
    clean('release_clean','SourceQIAquickNucleotideRemoval28304')
    q=g.product('bioanalyzer_qc','SourceBioanalyzeranalyticalaliquot;allocationunassigned',K.QC);kept=g.product('qc_kept','ProductionDNAretained',K.MOLECULAR);g.react('qc','SourceBioanalyzerQCseparatefromproduction',S.SPLIT,((x,R.TEMPLATE),),(kept,q));g.mark_terminal(Terminal(q.id,Disposition.QC,'SourceBioanalyzerQC'));x=kept
    op('kit_finish','SourceNEBNextUltraDNAE7103SfinishafterpriorON-BEADadapterligation;remainingoperations/cycles/indexesunassembled',S.MAGIC,(stock('kit_finish_preparation','SourcehistoricalNEBNextUltraDNAfinishpreparation; contributingoligos/nucleotidesunassembled,noinventedsecondligation',K.MOLECULAR,R.STRUCTURAL),),product='Sourcefinishedlibrarieswithunassembledhistoricalkitgeometry;noinventedsecondligationorQ5cyclecount')
    final=g.product('final_qc_kept','Finishedproductionlibraryretained',K.MOLECULAR);fq=g.product('final_bioanalyzer','SourceSECONDpost-finishBioanalyzeranalyticalaliquot;allocationunassigned',K.QC);g.react('final_qc','SourcesecondBioanalyzerassessmentafterkitfinish',S.SPLIT,((x,R.TEMPLATE),),(final,fq));g.mark_terminal(Terminal(fq.id,Disposition.QC,'Sourcepost-kitBioanalyzerQC'));x=final
    lib=g.product('library','SourceHiSeq4000multiplexedsequencingpreparation',K.LIBRARY);g.react('load','Sourcepool/loadingpreparation;exactaliquotsunassigned',S.GENERIC,((x,R.TEMPLATE),),(lib,),changes_molecule=False)
    g.sequence(lib,SequencingPlan(limitation='AssociatedAMsourceHiSeq4000;runlayout/indexgeometryunassembled. Peptide-DNAremnantstallsQ5inmodelassay;complementary-strandamplificationinterpretationnotcompleteproductvalidation. Tgfragmentenrichmentisnotnucleotide-positionstopmapping.'))
    g.add_source_item(SourceItem('pending','PublishedACSmain403andactualassociatedSIPDFgate; final-AMapplicabilityandhistoricalkitgeometryunverified.','Noacceptedabsenceorwhole-sourceaudit.'));return g
def selftest():
    from checks import Check
    c=Check();g=protocol_graph()
    c('GenomicpretreatmentprecedesNTHL1trap',g.steps['trap'].inputs[0].material,'pretreat_clean_kept')
    c('NTHL1isphysicalcontributor',[(p.material,p.role)for p in g.steps['trap'].inputs if p.material=='nthl1'],[('nthl1',R.STRUCTURAL)])
    c('BlockingDNAisphysicalcontributor',[(p.material,p.role)for p in g.steps['capture'].inputs if p.material=='salmon_blocker'],[('salmon_blocker',R.STRUCTURAL)])
    c('Onresinendprepaftercapturewash',g.steps['endprep'].inputs[0].material,'capture_wash_kept')
    c('USERbeforeelution',g.steps['postligation_wash'].inputs[0].material,'user_product')
    c('Noadditionaladapterstockafterprotease',sum(p.material=='adapter'for s in g.steps.values()for p in s.inputs),1)
    c('Selectedbandleavesphysicalremainder',g.terminals['band_rest'].disposition,Disposition.UNRESOLVED)
    c('ProteinaseafterPAGEselection',g.steps['proteinase'].inputs[0].material,'membrane_wash_kept')
    c('150ultransferprecedesQIAquick',g.steps['release_clean'].inputs[0].material,'released_liquid')
    c('Residualmembranefateexplicit',g.terminals['residual_membrane'].disposition,Disposition.UNRESOLVED)
    c('SourceQCisindependent',g.terminals['bioanalyzer_qc'].disposition,Disposition.QC)
    c('OnlysourceFragmentasefragmentation',sum(s.kind==S.FRAGMENTATION for s in g.steps.values()),1)
    c('Molecularreductantpresent',[(p.material,p.role)for p in g.steps['trap'].inputs if p.material=='reductant'],[('reductant',R.STRUCTURAL)])
    c('Kitfinishhasmolecularpreparedcontribution',[(p.material,p.role)for p in g.steps['kit_finish'].inputs if p.material=='kit_finish_preparation'],[('kit_finish_preparation',R.STRUCTURAL)])
    c('SecondQCafteractualfinish',g.steps['final_qc'].inputs[0].material,'kit_finish_product')
    c('Privatepartialscope',bool(g.publication_problems()),True);c.report()
if __name__=='__main__':selftest()
