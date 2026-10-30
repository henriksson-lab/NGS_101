from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Scene,Segment,revcomp,feature
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,SourceItem,SequencingPlan,Terminal,Disposition
from adapter_ligation import AdapterLigationGeometry
from pcr import PCRGeometry
class B:
    def __init__(self,title,authority,start,kind=K.BIOLOGICAL):
        self.g=ProtocolGraph(title,authority=authority,representation='published_instance');self.x=self.g.input('source',start,kind,primary=True)
    def stock(self,key,label,kind=K.CATALYST,role=R.CATALYST,**kw):return(self.g.input(key,label,kind,**kw),role)
    def cond(self,key,label):return self.stock(key,label,K.CONDITION,R.CONDITION)
    def op(self,key,label,kind,inputs=(),geometry=None,product='Source mixture including incomplete/unreacted/off-target classes; occupancy/recovery unassigned',reason=None):
        y=self.g.product(key+'_product',product,K.MOLECULAR);kw={}
        if kind==S.MAGIC:kw['reason']=reason or 'Positive source operation; complete occupied molecular inputs/products and historical formulation unassembled, not a missing-information conclusion.'
        if kind==S.GENERIC:kw['changes_molecule']=False
        if geometry is not None:kw['geometry']=geometry
        self.g.react(key,label,kind,((self.x,R.TEMPLATE),)+tuple(inputs),(y,),**kw);self.x=y;return y
    def clean(self,key,label,extra=()):
        y=self.g.product(key+'_kept','Actual source retained physical fraction; yield/carryover unassigned',K.MOLECULAR);w=self.g.product(key+'_waste','Actual separated source workup fraction',K.WASTE)
        self.g.react(key,label,S.CLEANUP,((self.x,R.TEMPLATE),self.cond(key+'_condition',label))+tuple(extra),(y,w));self.g.discard(w,'Actual workup fraction excluded; no perfect specificity or molecular yield claim.');self.x=y
    def split(self,key,label,selected,remainder):
        y=self.g.product(key+'_selected',selected,K.MOLECULAR);r=self.g.product(key+'_remainder',remainder,K.MOLECULAR)
        self.g.react(key,label,S.SPLIT,((self.x,R.TEMPLATE),),(y,r));self.g.mark_terminal(Terminal(r.id,Disposition.UNRESOLVED,remainder));self.x=y;return y

def unknown_scene():
    s=Scene();s.strand('unassembled',[Segment('source products unassembled','[unassembled]',placeholder=True)]);return s

from nextera import ME,ME_RC,ADAPTOR_S5,S5,S7
from illumina import P5,P7
AUTHORITY='PublishedPMC13645645main243–281/287; same-originalPDFp15; publishedSI1actualTables1–2PDFp11–13 rendered. Two originals; no vendor or predecessor recipe inherited.'
TN5_B_PREFIX='CATCGGCGTACGACT'
A_PREFIX='AGGCCAGAGCATTCG';A_SUFFIX='GTGGCCGATGTTTCG';B_PREFIX='CAAGCGTTGGCTTCTCGCATCT';B_SUFFIX='ATCCACGTGCTTGAG';BARCODE='AACGTGAT'
A1=A_PREFIX+BARCODE+A_SUFFIX;B1=B_PREFIX+BARCODE+B_SUFFIX
N501=P5+'TAGATCGC'+ADAPTOR_S5;N701=P7+'TCGCCTTA'+S7+ME+B_PREFIX
WHOLE_SOURCE={'Tn5ME-A':ADAPTOR_S5,'Tn5MErev':'/5Phos/'+ME_RC,'Tn5ME-B':'/5Phos/'+TN5_B_PREFIX+ME,'N501':N501,'N701':N701,'A1':'/5Phos/'+A1,'B1':B1}
def single(seq,label,mod5=''):
    s=Scene();s.strand('source',[Segment(label,seq)],mod5=mod5);return s

def transposome_duplex(which):
    s=Scene();prefix=S5 if which=='A'else TN5_B_PREFIX
    s.strand('full',[Segment('source prefix',prefix),Segment('source ME',ME)],mod5=''if which=='A'else'P')
    s.anneal('rev',[Segment('source MErev',ME_RC)],to='full',pair=('source MErev','source ME'),mod5='P');return s

def barcode_scene(which):
    s=Scene();prefix,suffix=(A_PREFIX,A_SUFFIX)if which=='A'else(B_PREFIX,B_SUFFIX)
    s.strand('barcode',[Segment('source '+which+' prefix',prefix),Segment('source variable '+which+' spatial axis',BARCODE,feature=feature('spatial_atac_hic.axis_'+which.lower(),'spatial_barcode','fixed',group='spatial_atac_hic.coordinate',part=which.lower(),note='Derived variable window of source50memberfamily; separate orthogonal axis, NOT UMI')),Segment('source '+which+' suffix',suffix)],mod5='P'if which=='A'else'');return s

def pcr_primer(which):
    s=Scene()
    if which=='P5':parts=[Segment('P5',P5),Segment('source N501 sample index','TAGATCGC',feature=feature('spatial_atac_hic.i5','sample_index','fixed',note='Source listed N501 candidate; actual sample assignment unassembled')),Segment('source A extension',ADAPTOR_S5)]
    else:parts=[Segment('P7',P7),Segment('source N701 sample index','TCGCCTTA',feature=feature('spatial_atac_hic.i7','sample_index','fixed',note='Source listed N701 candidate; actual sample assignment unassembled')),Segment('source extension',S7+ME+B_PREFIX)]
    s.strand('primer',parts);return s

def protocol_graph():
    b=B('Spatial-ATAC-Hi-C selected wholeROI contact/accessibility preparation',AUTHORITY,'Actual fresh-frozen7–10um tissue section onPLLslide, upstreamcollection/histology outsidegraph');g=b.g;tissue=b.x
    a=g.input('tn5_a','ActualsourceTn5ME-A oligo',K.OLIGO,state=single(ADAPTOR_S5,'source Tn5ME-A'))
    r=g.input('tn5_rev','ActualsourceTn5MErev',K.OLIGO,state=single(ME_RC,'source Tn5MErev','P'))
    bb=g.input('tn5_b','ActualsourceTn5ME-B',K.OLIGO,state=single(TN5_B_PREFIX+ME,'source Tn5ME-B','P'))
    b.x=a;b.op('transposome_prepare','ActualDiagenodeC01070010 unloadedTn5 assembly perhistoricalvendor instructions',S.MAGIC,((r,R.STRUCTURAL),(bb,R.STRUCTURAL),b.stock('tn5_enzyme','ActualunloadedTn5catalyst'),b.cond('tn5_loading','Completevendorloadingprogramme/stoichiometry/dimerpopulationsunassembled')),reason='Source3oligotypes andunloadedTn5 positive; localMEduplexes pairbycomputation, notassertedfulloccupieddimerdistribution.')
    b.split('transposome_aliquot','Actual7.5ulselectedfrompreparedtransposomestock','Actual7.5ultransposomeinput','Actualuntransferredtransposomestock, quantity/fateunassembled');transposome=b.x;b.x=tissue
    b.op('fix','ThawRT10min;source2percentFA20min',S.CROSS_LINKING,(b.cond('fa','Actual2percentformaldehyde'),))
    b.op('quench','2.5Mglycinestock to.2Mfinal5min',S.QUENCHING,(b.cond('glycine','Actualsource.2Mfinalglycine'),));b.clean('fix_wash','2x1mlDPBS; retainon-slide tissue')
    b.op('permeabilize','30ullysisbuffer4C20min; PDMSreservoir~200ul',S.PERMEABILIZATION,(b.cond('lysis','Source1ul1MTris8/.2ul5MNaCl/2ul10percentNP40/proteaseinhibitor in30ul'),))
    b.op('sds','Source literal36ul10percentSDS62C10min',S.MAGIC,(b.cond('sds_stock','Same-originalpublishedPDFp15 confirms36ul10percentSDS, no silentnormalization'),))
    b.op('triton','Add30ul1percentTriton37C15min',S.GENERIC,(b.cond('triton_stock','ActualsourceTritonaddition'),))
    b.op('restriction','Source20.5ulDpnII/HinfImix37C1h',S.DIGESTION,(b.stock('dpnii','Actual8ulDpnII'),b.stock('hinfi','Actual1.5ulHinfI'),b.cond('restriction_buffer','11ul10xCutSmart/1ul10percentTriton')))
    b.clean('digest_wash','2x200ulwaterRT2–5min')
    b.op('fill_biotin','Source50ulKlenowbiotinfillRT1h',S.MAGIC,(b.stock('klenow','Actual5ul5UulKlenowM0210L'),b.stock('fill_nt','Actual20ul.4mMbiotin14dATP+1uleach10mMdCTP/dGTP/dTTP',K.MOLECULAR,R.STRUCTURAL),b.cond('fill_buffer','5ul10xCutSmart/17ulwater')),reason='Sourcebiotin14dATP incorporation positive; wholeoccupiedrestrictionends/donor/acceptorproducts unassembled, no exactjunction claimed.')
    b.op('contact_ligate','Sourcecontactligation100ul16Covernight BEFORETn5',S.MAGIC,(b.stock('contact_ligase','Actual10ulT4ligaseM0202L'),b.stock('contact_bsa','Actual1ul20mgmlBSA',K.MOLECULAR,R.STRUCTURAL),b.cond('contact_buffer','10ul10xT4buffer/1ul10percentTriton/78ulwater')),reason='Physical on-slidecontactligation; complete occupied DNAends/junction classes unassembled. No HiCARfront-end recipe inherited.')
    b.clean('pre_tn5_wash','200ulpre-tagmentwashRT5min; source1MTris7.4/5MNaCl/1MMgCl2/10percentTween20/5percentBSA/water descriptors; finalformulationunassembled')
    b.op('tagment','Source100ultransposition37C1h aftercontactligation',S.MAGIC,((transposome,R.STRUCTURAL),b.cond('tagment_buffer','50ul2xtagmentbuffer/33ulDPBS/1ul10percentTween20/3.5ulwater')),reason='Sourceactual7.5ulpreparedtransposome; completevendorbuffer/enddistribution/loaded-dimerproductsunassembled, notborrowedcanonicalTn5products.')
    b.op('tn5_stop','Source100ul40mMEDTART5min',S.QUENCHING,(b.cond('edta','Actualsource40mMEDTAaddedsolution'),));b.clean('tn5_wash','Removeallreservoirliquid;2x500ulNEBuffer3.1/5min')
    for axis in ('A','B'):
        letter=axis.lower();stock=g.input('barcode_'+letter+'_plate','Actual50-channel preparedbarcode'+axis+'stockarray; exactA1/B1representativecomponentdrawn, not50solutionsmixedtogether',K.MOLECULAR,state=barcode_scene(axis),chemical_state={'whole_source_representative':WHOLE_SOURCE[axis+'1'],'source_preanneal':'10ul100uMbarcode+10ul100uMlinker1or2+20ul2xannealbuffer→25uMbarcode','whole_linker':'unassembled','spatial_assignment':'Orthogonal50channels; wholeROI remains physicallyononeslide'})
        b.op('barcode_'+letter,'Actual'+axis+'50channelin-situligation37C30min, alignedbrightfieldimage',S.MAGIC,((stock,R.STRUCTURAL),b.stock('barcode_'+letter+'_ligase','ActualsourceT4ligaseinpreparedchannelmix'),b.cond('barcode_'+letter+'_conditions','Perchannel5ul=2ulligationmix(27ulT4buffer/72.4ulwater/5.4ul5percentTriton/11ulligase)+2ulNEBuffer3.1+1ul25uMpreannealedbarcode')),reason='Source spatial-channel ligation, separate orthogonalA/Bidentities; wholelinkerstrands/occupiedsubstrategeometry/productsunassembled. Notphysicalpixelaliquots or pooled50barcodes inonevessel.')
        b.clean('barcode_'+letter+'_wash','RemovePDMSchip/clamp; waterwash/airdry/imaging; retainwholeROI')
    b.op('roi_lyse','100ul2mgmlPK/TE8/1percentSDS65C1honslide',S.MAGIC,(b.stock('pk','Actual2mgmlProteinaseK'),b.cond('roi_buffer','Actual100ulreversecrosslinkmix')),reason='ActualwholeROIlysis; spatialpixels distinguished bybarcodes later, no physical2500-tubesplit.')
    b.op('reverse_crosslink','Collectlysateintube65Covernight',S.THERMAL_TREATMENT,(b.cond('reverse_program','Actualsource65Covernight'),));b.clean('dna_purify','ZymoChIPDNA/elute50ulEB')
    bound=g.product('t1_bound','ActualT1retainedbiotinDNApreparation, occupancy/purityunassigned',K.MOLECULAR);flow=g.product('t1_flow','Actualnonretainedcapturefraction',K.WASTE)
    g.react('biotin_capture','50ulDNA+50ulT1beads/2xB&W1hrotatorUU40',S.AFFINITY_CAPTURE,((b.x,R.TEMPLATE),b.stock('t1','Actual50ul10mgmlT1stock washed2x200ul1xB&W then50ul2xB&W',K.SELECTION,R.SELECTION),b.cond('capture_buffer','2xB&W10mMTris7.5/1mMEDTA/2MNaCl; actualvolumes50+50')),(bound,flow));g.discard(flow,'Actualnonretainedfraction');b.x=bound;b.clean('t1_wash','2x200ulTWB55C2minmix/magneticdiscard; resuspend30ulEB')
    b.split('pcr_beads','Actual10ulaliquotof30ulbeadsuspension','Actual10ulPCRbeadinput','Actual20uluntransferredbeadstock, fateunassigned')
    forward=g.input('p5_primer','Source listedN501 representativecustomP5candidate; actualsamplechoiceunassembled',K.OLIGO,state=pcr_primer('P5'),chemical_state={'source_ordering':N501,'actual_sample_assignment':'unassembled'})
    reverse=g.input('i7_primer','Source listedN701 representativecustomi7candidate; actualsamplechoiceunassembled',K.OLIGO,state=pcr_primer('i7'),chemical_state={'source_ordering':N701,'actual_sample_assignment':'unassembled'})
    q5stock=g.input('q5_stock','Actual one25ul2xQ5 High-Fidelity master-mix aliquot; whole historical formulation/doses unassembled',K.MOLECULAR)
    q5enzyme=g.product('q5_enzyme','Copresent Q5 DNApolymerase component of SAME25ulmastermix, moleculardoseunassembled',K.CATALYST)
    q5nt=g.product('q5_nucleotides','Copresent dNTP component of SAME25ulmastermix, individualconcentrationsunassembled',K.MOLECULAR)
    q5buffer=g.product('q5_buffer','Copresent reactionbuffer component of SAME25ulmastermix, completeformulationunassembled',K.CONDITION)
    g.react('q5_components','Expose copresent enzyme/nucleotide/buffer in ONE25ulstock, not physical fractions or fresh additions',S.GENERIC,((q5stock,R.TEMPLATE),),(q5enzyme,q5nt,q5buffer),changes_molecule=False)
    b.op('pcr','25ul2xQ5+2.5ul25uMP5+2.5ul25uMi7+10ulbeads+10ulwater;72C5min/98C30s/10–13cycles98C10s63C30s72C1min',S.MAGIC,((forward,R.STRUCTURAL),(reverse,R.STRUCTURAL),(q5enzyme,R.CATALYST),(q5nt,R.STRUCTURAL),(q5buffer,R.CONDITION)),reason='Sourcewholeprimercandidatesretained, actualsampleindexchoice/occupiedwholeprimerlanding/PCRproductsunassembled; noprimerproductproofinvented.')
    b.clean('pcr_clean','Source1xAMPure/elute50ulwater')
    b.op('double_size','Sourceprinted0.55xthen0.15xAMPure double-size-selection',S.MAGIC,(b.stock('size_beads','ActualsourceAMPurestockallocations',K.SELECTION,R.SELECTION),),reason='Positive source double-size-selection; exactintermediatefraction/transfer schedule unassembled, notsilently interpreted as0.70xfinal.')
    lib=g.product('joint_library','Actualsinglecontact/accessibilityspatialDNA library; no physicalATACversusHiClibrarysplit',K.LIBRARY);g.react('load','SourceIlluminaXplus2x150paired-end',S.GENERIC,((b.x,R.TEMPLATE),),(lib,),changes_molecule=False);g.sequence(lib,SequencingPlan(read_lengths={'R1':150,'R2':150},limitation='SourceR2containsDNA+spatialbarcodes; pipelinecomputationallysplitstoR1DNA/R2barcodes andoriginalR1→R3DNA. Exactwholelinkers/readlanding/indexmaskunassembled;2500pixels notphysicalaliquots.'))
    g.add_source_item(SourceItem('printed_reagent_descriptors','ActualmainPDFp15confirms36ul10percentSDS; prewashlists1MTris/5MNaCl/1MMgCl2/10percentTween/5percentBSAwithoutvolumes. Modelretainsprinteddescriptors, finalwashformulationunassembled.','No silentunitfix or provenpublicationerrorclaim; sourceonlyprimer/template/componentgeometrypartial.'));return g

def selftest():
    from checks import Check,Source,have
    from storage import sources_root
    c=Check();g=protocol_graph();transposome_duplex('A');transposome_duplex('B')
    c('Contactligationbeforetransposition',g.steps['pre_tn5_wash'].inputs[0].material,'contact_ligate_product');c('Transpositionusesactualselectedstock',g.steps['tagment'].inputs[1].material,'transposome_aliquot_selected');c('BarcodeBfollowsphysicalAworkup',g.steps['barcode_b'].inputs[0].material,'barcode_a_wash_kept');c('WholeROIlysisfollowsbothaxes',g.steps['roi_lyse'].inputs[0].material,'barcode_b_wash_kept');c('PCRusesselected10ul',g.steps['pcr'].inputs[0].material,'pcr_beads_selected');c('Retains20ulstock',g.terminals['pcr_beads_remainder'].disposition,Disposition.UNRESOLVED);c('OneactualQ5physicalstock',g.steps['q5_components'].inputs[0].material,'q5_stock');c('PCRconsumessameQ5components',{x.material for x in g.steps['pcr'].inputs if x.material.startswith('q5_')},{'q5_enzyme','q5_nucleotides','q5_buffer'});c('Onephysicaljointlibrary',set(g.sequencing),{'joint_library'})
    root=sources_root()/Path(__file__).resolve().parents[1].name;p=root/'s41592-026-03217-4_supp_41592_2026_3217_MOESM1_ESM.pdf.txt'
    if have(c,Source(p,'https://doi.org/10.1038/s41592-026-03217-4')):
        t=''.join(p.read_text().split())
        for k,seq in WHOLE_SOURCE.items():c('Originalwholemodifiedordering '+k,seq in t,True)
    p=root/'s41592-026-03217-4_PMC13645645.xml.txt'
    if have(c,Source(p,'https://doi.org/10.1038/s41592-026-03217-4')):
        t=''.join(p.read_text().split());c('Actual30ulbeadstock','resuspendedin30μlofEBbuffer'in t or'resuspendedin30µlofEBbuffer'in t,True)
    c('Partialnotpublicready',bool(g.publication_problems()),True);c.report()
if __name__=='__main__':selftest()
