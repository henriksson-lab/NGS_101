"""Private source-owned PolII/CPD PADD-seq preparation, selected detailed SI route."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Scene,Segment,revcomp
from illumina import NEBNEXT_ARM_READ2,TRUSEQ_READ1
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,SourceItem,SequencingPlan,Terminal,Disposition
from adapter_ligation import AdapterLigationGeometry
from pcr import PCRGeometry
TITLE='PADD-seq selected PolII/CPD route with NVP-2'
AUTHORITY='Published PMC10085696 main91–115/129; published SI TableS1 originalp2 and SupplementaryProtocol160–413, PCRtableoriginalp16 viewed; two originals. H3K9ac/cisplatin and Input/ChIP-seq controls separate.'
AD1T=NEBNEXT_ARM_READ2+'A'
AD1B='NNNNNGACTGGTTCCAATTGAAAGTGCTCTTCCGATCT'
O3P=revcomp(NEBNEXT_ARM_READ2)
AD2T=revcomp(TRUSEQ_READ1)
AD2B=TRUSEQ_READ1+'NNNNN'
RAW={'ad1t':'/5Phos/'+AD1T+'/3SpC3/','ad1b':AD1B[:-1]+'*'+AD1B[-1],'o3p':O3P,'ad2t':'/5Phos/'+AD2T+'/3SpC3/','ad2b':AD2B+'/3SpC3/'}
def unknown_scene():
    sc=Scene();sc.strand('unassembled',[Segment('unassembled','[source whole product geometry unassembled]',placeholder=True)]);return sc
def inventory(name):
    seq={'ad1t':AD1T,'ad1b':AD1B,'o3p':O3P,'ad2t':AD2T,'ad2b':AD2B}[name]
    sc=Scene();parts=[Segment('random splint window','N'*5,placeholder=True),Segment('fixed adapter',seq[5:])]if name=='ad1b'else([Segment('fixed adapter',seq[:-5]),Segment('random splint window','N'*5,placeholder=True)]if name=='ad2b'else[Segment('source ordering',seq)])
    sc.strand('oligo',parts,mod5='P'if name in{'ad1t','ad2t'}else'',mod3='C3 spacer'if name in{'ad1t','ad2t','ad2b'}else'');return sc
def adapter_scene(number):
    sc=Scene()
    if number==1:
        n=14
        if AD1T[:n]!=revcomp(AD1B[23:37]):raise ValueError('Source Ad1 terminal stem mismatch')
        sc.strand('top',[Segment('Ad1 stem top',AD1T[:n]),Segment('Ad1 top fork',AD1T[n:])],mod5='P',mod3='C3 spacer')
        sc.anneal('bottom',[Segment('Ad1 random5','N'*5,placeholder=True),Segment('Ad1 bottom fork',AD1B[5:23]),Segment('Ad1 stem bottom',AD1B[23:37]),Segment('Ad1 terminal T',AD1B[-1:])],to='top',pair=('Ad1 stem bottom','Ad1 stem top'),unpaired=('Ad1 random5','Ad1 bottom fork'))
        sc.overhang('bottom','Ad1 terminal T',polarity='3-prime',owner='adapter',text='Computed source Ad1 3primeT; source DNA dA-tailing known')
    else:
        sc.strand('top',[Segment('Ad2 paired donor',AD2T)],mod5='P',mod3='C3 spacer')
        sc.anneal('bottom',[Segment('Ad2 paired splint',TRUSEQ_READ1),Segment('Ad2 random splint','N'*5,placeholder=True)],to='top',pair=('Ad2 paired splint','Ad2 paired donor'),mod3='C3 spacer')
        sc.ligation_end('top','Ad2 paired donor',kind='splinted',polarity='5-prime',owner='adapter')
        sc.note('top','Source random5nt splint pool contacts an unassigned extension-product terminus; no exact occupied splint/substrate junction or UMI meaning asserted.')
    return sc
def geometry(n):
    sc=adapter_scene(n)
    if n==1:
        sc.anneal('source_dna_3a',[Segment('sourceDNAunknownbody','[source fragment body unassigned]',placeholder=True),Segment('sourceDNA3A','A')],to='bottom',pair=('sourceDNA3A','Ad1 terminal T'),unpaired=('sourceDNAunknownbody',))
        sc.ligation_end('source_dna_3a','sourceDNA3A',kind='cohesive',polarity='3-prime',owner='substrate')
        sc.note('source_dna_3a','Only source dA-tail/adapterT contact is asserted; other genomicstrand/5primeend and full sealed product remain unassembled. Not added to independentadapterstockinventory.')
    return AdapterLigationGeometry.known('cohesive'if n==1 else'splinted',sc,unknown_scene(),statement='Computed exact local source adapter pairing/end and,forAd1,source3Aprime-tail contact only; complete genomic sealed products, splint occupancy and yields unassembled.')

def protocol_graph():
    g=ProtocolGraph(TITLE,authority=AUTHORITY,representation='published_instance')
    x=g.input('cells','Source three15cmdishes~3e7humanfibroblasts; selected PolII/CPD+NVP-2 route',K.BIOLOGICAL,primary=True)
    def inp(key,text,kind=K.CATALYST,role=R.CATALYST,**kw):return(g.input(key,text,kind,**kw),role)
    def condition(key,text):return inp(key,text,K.CONDITION,R.CONDITION)
    def op(key,label,kind,ins=(),product='Source reaction mixture with residual substrates/free reagents retained until workup',geometry=None):
        nonlocal x
        kw={}
        if kind==S.GENERIC:kw['changes_molecule']=False
        if geometry is not None:kw['geometry']=geometry
        y=g.product(key+'_product',product,K.MOLECULAR);g.react(key,label,kind,((x,R.TEMPLATE),)+tuple(ins),(y,),**kw);x=y
    def clean(key,label):
        nonlocal x
        c=g.input(key+'_condition',label,K.CONDITION);a=g.product(key+'_kept','Physically retained source fraction; recovery/carryover unassigned',K.MOLECULAR);w=g.product(key+'_waste','Source discarded workup fraction',K.WASTE);g.react(key,label,S.CLEANUP,((x,R.TEMPLATE),(c,R.CONDITION)),(a,w));g.discard(w,'Source physically removes this fraction; no ideal molecular specificity.');x=a
    def split(key,label,rest,disp=Disposition.UNRESOLVED):
        nonlocal x
        a=g.product(key+'_used',label,K.MOLECULAR);r=g.product(key+'_rest',rest,K.QC if disp==Disposition.QC else K.MOLECULAR);g.react(key,label,S.SPLIT,((x,R.TEMPLATE),),(a,r));g.mark_terminal(Terminal(r.id,disp,rest));x=a
    def ready_adapter(n):
        names=('ad1t','ad1b')if n==1 else('ad2t','ad2b');a,b=[g.input(name,'Source '+RAW[name],K.OLIGO,state=inventory(name),chemical_state={'source_ordering':RAW[name],'polymer':'DNA','phosphorothioate_after':(len(AD1B)-2,)if name=='ad1b'else()})for name in names]
        c=g.input(f'ad{n}_anneal_condition','Source100uM each oligo in50ul0.5mMTris8/100mMNaCl/0.05mMEDTA,boil2min/naturalcool/store-20C; selectedaliquot stock concentration bridge unassigned',K.CONDITION)
        bulk=g.product(f'ad{n}_annealed_stock','Source complete prepared adapter stock before aliquot; dilution to40uM statedligationstock unassembled',K.MOLECULAR,state=adapter_scene(n));g.react(f'ad{n}_anneal','Source independent adapter annealing',S.HYBRIDIZATION,((a,R.TEMPLATE),(b,R.STRUCTURAL),(c,R.CONDITION)),(bulk,))
        y=g.product(f'ad{n}_ready','Source selected '+('50pmol/1.25ul40uM'if n==1 else'40pmol/1ul40uM')+' adapter aliquot; preparedstock dilution bridge unassembled',K.ADAPTER,state=adapter_scene(n));rest=g.product(f'ad{n}_remainder','Untransferred prepared adapter stock; quantity and laterfate unassigned',K.MOLECULAR);g.react(f'ad{n}_aliquot','Source actual adapter aliquot, not entire50ulannealingstock',S.SPLIT,((bulk,R.TEMPLATE),),(y,rest));g.mark_terminal(Terminal(rest.id,Disposition.UNRESOLVED,'Source stores prepared stock at-20C; remaining stock quantity/lateruseunassigned.'));return y
    ad1=ready_adapter(1);ad2=ready_adapter(2)
    op('nvp','Source250nMNVP-2 5minbeforeUV; recovery time afterUV selected byexperiment, unassigned here',S.GENERIC,(condition('nvp_condition','Source NVP-2 medium'),))
    op('uv','Source remove medium,20Jm2 254nmUV-C,restoreNVP-medium andsourceincubation',S.CHEMICAL_CONVERSION,(condition('uv_condition','Source UV irradiation'),),product='Source UV-damaged cells; lesion positions/occupancy unassigned')
    op('crosslink','Source1percentformaldehydeRT10min',S.CROSS_LINKING,(condition('formaldehyde','Source formaldehyde'),))
    op('quench','Source150mMglycineRT5min',S.QUENCHING,(condition('glycine','Source glycine'),))
    clean('cell_wash','Source twiceicePBS,collectcells andpellet; source supernatantsdiscarded')
    clean('lysis1','Source2.1mllysisbuffer1/ice10min/850g10min4C; retainpellet anddiscard supernatant')
    clean('lysis2','Source2.1mllysisbuffer2/ice10min/850g10min4C; retainpellet anddiscard supernatant')
    op('lysis3','Source630ul1.5percentSDS-RIPA/proteaseinhibitors,ice30min',S.PERMEABILIZATION,(condition('lysis3_condition','Source statedSDS-RIPA'),),product='Source mixed lysate with soluble/insoluble material retained together throughsonication; no physical separation untilclarification')
    op('sonicate','Source~140ulaliquots Q80030percentamplitude4.5min30son/off4C,combine;target300–600bp',S.FRAGMENTATION,(condition('sonicate_condition','Source Q800 programme'),))
    clean('clarify','Source20000g10min4C; transferchromatin supernatant, separatedpellet')
    split('qubit_chromatin','Productionchromatinretained','Source1ulQubitBR analyticalaliquot',Disposition.QC)
    split('chip_aliquot','Source167ugchromatinselected;14volSDS-freeRIPA gives0.1percentSDS','Remainingchromatin; allocation toindependentInput/Damagecontrols unassigned')
    ab=inp('rpb1','SourceSI25ulantiRPB1A304-405A; main97gives25ug, stock concentration relation unassigned',K.SELECTION,R.SELECTION)
    op('rpb_bind','Source0.1percentBSA/0.1ugultRNA/proteaseinhibitors+antiRPB1,4Crotation1–2h',S.BINDING,(ab,condition('rpb_condition','Source blocking/incubationconditions')))
    # Prepared bead stock is made independently, not added as fresh stock twice.
    pa=g.input('protein_a','Source60ulProteinAagarose slurry',K.MOLECULAR);pc=g.input('protein_a_condition','Sourcepellet/two1mlRIPAwashes;470ulRIPA+20.8ulinhibitor+17.2ul3percentBSA+1.7ul20ugultRNA,4C~1h',K.CONDITION);pready=g.product('protein_a_ready','Prepared sourceProteinAbeads inblockingmixture',K.SELECTION);pw=g.product('protein_a_waste','Discardedsourcebeadstock/washsupernatants',K.WASTE)
    g.react('protein_a_prepare','Source independent washed/blockedProteinAstock',S.CLEANUP,((pa,R.TEMPLATE),(pc,R.CONDITION)),(pready,pw));g.discard(pw,'Source explicitlyremoves wash supernatants.')
    op('chip_capture','Source transferblockedProteinAmixture toantibody/chromatin,4Covernight',S.BINDING,((pready,R.SELECTION),),product='Source ProteinA/antiRPB1/chromatin bindingmixture withunboundchromatin beforepelleting')
    clean('chip_wash','Source750g5min/discardsupernatant;3x1mlRIPA,3xRIPA500,3xLiCl,1xTE; retainedProteinAcomplexes')
    op('chip_elute','Source240uldirectelutionbuffer10mMTris8.8/5mMEDTA/300mMNaCl/1percentSDS65C15min1500rpm',S.GENERIC,(condition('chip_elute_condition','Source directelutionbuffer'),))
    clean('chip_separate','Source20000g1minRT; transferDNA-containingeluateawayfrombeads')
    op('rnase','Source9.6ulRNaseA37C30min',S.DIGESTION,(inp('rnase_a','SourceRNaseA'),))
    op('proteinase','Source6ulProteinaseK55C2h',S.DEPROTEINIZATION,(inp('pk','SourceProteinaseK'),))
    op('reverse_crosslink','Source65C~10hovernightreversecrosslink',S.THERMAL_TREATMENT,(condition('reverse_condition','Source heatprogramme'),))
    clean('chip_clean','Sourcephenol/chloroform+ethanolprecipitation afterreversecrosslink')
    split('chip_dna_aliquot','Source50–100ngChIPDNAselected (SIminimum25ng)','OtherChIPDNA inclQubitHS analyticalallocation unassigned')
    op('endprep','SourceNEBNextUltraIIendrepair/dAtailing; historicalprogrammeunassembled',S.END_REPAIR,(inp('endprep_enz','SourceNEBNextUltraIIendprepmix'),))
    op('ligate1','Source50pmolAd1 at4Covernight withNEBligationmodule',S.LIGATION,((ad1,R.STRUCTURAL),inp('ligase1','SourceNEBfirstligationmix')),geometry=geometry(1))
    clean('ad1_clean','Source0.8xFragSelectcleanup/elute12.2ul0.1xTE')
    op('denature1','Sourceadd5ul8Murea,98C2min/ice1min',S.THERMAL_TREATMENT,(condition('denature1_condition','Sourceurea/heat'),))
    beads=g.input('cpd_bead_stock','Source1.25ulProteinGDynabeads+1.25ulanti-rabbitDynabeads',K.MOLECULAR);anti=g.input('cpd_antibody_stock','Source0.5ulantiCPD+0.5ulrabbitanti-mouseIgG,source50percentglycerolstocks',K.SELECTION);bc=g.input('cpd_prepare_condition','Source2x190ulIPbufferwashes,5ulbuffer+0.25ulboiledsalmonDNA,4Covernight,1wash,resuspend5ulbuffer+0.5ulboiledsalmonDNA',K.CONDITION);br=g.product('cpd_beads_ready','Source5.5ulblockedantiCPDbeadpreparation',K.SELECTION);bw=g.product('cpd_prepare_waste','Source removedbeadloadingwashfractions',K.WASTE)
    g.react('cpd_prepare','Source independentblockedantiCPDbeadpreparation',S.CLEANUP,((beads,R.TEMPLATE),(anti,R.SELECTION),(bc,R.CONDITION)),(br,bw));g.discard(bw,'Source physicalbeadstock/loadingwash fractions.')
    op('cpd_bind','Sourceadd2.5ul8xIPbuffer+5.5ulpreparedantiCPDbeads,4Covernightrotation',S.BINDING,((br,R.SELECTION),),product='Source damage-antibodybindingmixture, notperfectlesionselection')
    clean('cpd_wash','Source magneticcapture;sequentialwashU/I/II/III/IV/TEatRT, physicallyretainCPDbeads')
    op('cpd_elute','Source100ul1percentSDS/TE65C10min1100rpm',S.GENERIC,(condition('cpd_elution_condition','SourceSDSelution'),))
    clean('cpd_clean','Sourceisolateeluate/phenolchloroform/ethanol;dissolve6ul0.1xTE')
    op('extend','Source6ulDNA+7.5ulQ5mix+1.5ul30pmolO3P;98C50s65C5min',S.EXTENSION,(inp('o3p','SourceO3P '+O3P,K.OLIGO,R.STRUCTURAL,state=inventory('o3p')),inp('q5_extension','SourcefreshNEBNextUltraIIQ5')),product='SourceO3Pextensionstrands, lesion-stalledandfulllength/unextendedclasses plusoriginaldamagedtemplates; no perfectstopfrequency or fragment disappearance')
    op('exoi','Source1ulExoI37C15min; nucleaseactionnotphysicalcleanup',S.DIGESTION,(inp('exoi_enzyme','SourceExoI;SIprintsExo1,mainExoI'),))
    clean('extension_clean','Sourceadd25ulTE,1.1xFragSelect45.6ul/elute100ul;ethanolconcentrate6.5ul')
    op('denature2','Source98C3min/ice1min',S.THERMAL_TREATMENT,(condition('denature2_condition','Sourceheatprogramme'),))
    op('ligate2','Source1ul40uMAd2+7.5ul2xInstantStickyEndsLigase,4Covernight',S.LIGATION,((ad2,R.STRUCTURAL),inp('ligase2','SourcefreshInstantStickyEndsLigase')),geometry=geometry(2))
    clean('ad2_clean','Sourceadd35ulTE,1.1xFragSelect55ul/elute15.5ul')
    production=g.product('qc_used','Source15ul productionligationproduct',K.MOLECULAR);analytical=g.product('qc_rest','Source0.3ul analytical22cycleQ5PCR/5percentPAGE',K.QC);unused=g.product('qc_unused','Nominal0.2uluntransferredligationeluate;laterfateunassigned',K.MOLECULAR)
    g.react('qc','Physical15.5ul eluate allocation:15ulproduction,0.3ulQC,nominal0.2ulunused',S.SPLIT,((x,R.TEMPLATE),),(production,analytical,unused));g.mark_terminal(Terminal(analytical.id,Disposition.QC,'Source0.3ulQCworkflow'));g.mark_terminal(Terminal(unused.id,Disposition.UNRESOLVED,'Untransferrednominal0.2ul, not analyticalQC'));x=production
    op('pcr','Source15ulDNA+5uleachi5/i7atLITERAL10mM+25ulQ5mix;98C30s,15x98C10s65C75s,65C5min',S.PCR,(inp('i5','SourceNEBMultiplexi5primer stock,exactorder/indexunassembled',K.OLIGO,R.STRUCTURAL),inp('i7','SourceNEBMultiplexi7primerstock,exactorder/indexunassembled',K.OLIGO,R.STRUCTURAL),inp('q5_pcr','SourcefreshNEBNextUltraIIQ5')),geometry=PCRGeometry.unavailable(unknown_scene(),reason='Sourcehistoricalmultiplexprimersandwholeproductsunassembled;notconfirmedmissing.'))
    c=g.input('size1_condition','Sourceadd50ulTEthen40ulbeads,RT>=5min/magnet',K.CONDITION);sup=g.product('size1_supernatant','Source transferredsupernatant afterfirstsizebinding',K.MOLECULAR);big=g.product('size1_beads','Separatedfirstbead-boundfraction;laterfateunassigned',K.MOLECULAR)
    g.react('size1','Source firstsizecut retainsSUPERNATANT',S.CLEANUP,((x,R.TEMPLATE),(c,R.CONDITION)),(sup,big));g.mark_terminal(Terminal(big.id,Disposition.UNRESOLVED,'Source transferssupernatant; firstboundfractionlaterfateunassigned.'));x=sup
    clean('size2','Sourceadd50ulfreshbeadstotransferredsupernatant,RT>=5min;discardsecondsupernatant,two200ul70percentEtOHwashes,elute16ulTris8;target~300–1000bp')
    lib=g.product('library','Source physicallysizeselectedPADDlibrary',K.LIBRARY);g.react('load','SourceQubitHS/loadingpreparation',S.GENERIC,((x,R.TEMPLATE),),(lib,),changes_molecule=False)
    g.sequence(lib,SequencingPlan(limitation='SourceSI PE150IlluminaNovaSeq,~20Mrawpairedreads. Damageinference2ntupstreamofmappedreads is computational; no read-endlanding/wholelibrary claim. Input/ChIPcomparatorsphotorepairnotimported intoPADD.'))
    g.add_source_item(SourceItem('pending','Full historicalkit/UDI/readgeometry, antibodymass/volume stockbridge and suspiciousliteral10mMprimerunitsunassembled.','Source-positive selectedscope; no confirmedabsence/fullreactionreview.'));return g
def selftest():
    from checks import Check,Source,have
    from storage import sources_root
    c=Check();g=protocol_graph()
    c('ChIP beforedamageIP',g.steps['cpd_bind'].inputs[0].material,'denature1_product')
    c('Reversecrosslinkbeforelibraryprep',g.steps['chip_clean'].inputs[0].material,'reverse_crosslink_product')
    c('PreparedProteinAonlyonce',sum(p.material=='protein_a_ready'for s in g.steps.values()for p in s.inputs),1)
    c('CPDcaptureusespreparedstock',[p.material for p in g.steps['cpd_bind'].inputs if p.role==R.SELECTION],['cpd_beads_ready'])
    c('O3Paddedonce',sum(p.material=='o3p'for s in g.steps.values()for p in s.inputs),1)
    c('ExoIfollowsextension',g.steps['exoi'].inputs[0].material,'extend_product')
    c('Heatbeforesecondligation',g.steps['ligate2'].inputs[0].material,'denature2_product')
    c('Ad2sourceblockingretained',inventory('ad2b').strands['oligo'].mod3,'C3 spacer')
    c('FirstsizefractionisdiscardedNOTsilently',g.terminals['size1_beads'].disposition,Disposition.UNRESOLVED)
    c('Secondsizecutreceivesfirstsupernatant',g.steps['size2'].inputs[0].material,'size1_supernatant')
    c('QCseparatefromproductionPCR',g.terminals['qc_rest'].disposition,Disposition.QC)
    c('UnusedeluateisnotQC',g.terminals['qc_unused'].disposition,Disposition.UNRESOLVED)
    for n in(1,2):c('Adapterselectedfrompreparedstock '+str(n),g.steps[f'ad{n}_aliquot'].inputs[0].material,f'ad{n}_annealed_stock')
    src=sources_root()/Path(__file__).resolve().parents[1].name/'nar_gkad008_supp_gkad008_supplemental_file.pdf.txt'
    if have(c,Source(src,'https://doi.org/10.1093/nar/gkad008')):
        text=''.join(src.read_text().split())
        for name,s in RAW.items():c('Wholemodifiedsourceordering '+name,s in text,True)
    c('Privateincompletescope',bool(g.publication_problems()),True);c.report()
if __name__=='__main__':selftest()
