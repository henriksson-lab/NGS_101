"""Private source-qualified partial preparation; no complete molecular/absence claim."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Scene,Segment,feature,revcomp
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,SourceItem,SourceCoverage,SequencingPlan,Terminal,Disposition
class B:
    def __init__(self,title,authority,label,kind=K.MOLECULAR):
        self.g=ProtocolGraph(title,authority=authority,representation='published_instance');self.x=self.g.input('source',label,kind,primary=True)
    def stock(self,key,label,kind=K.MOLECULAR,role=R.STRUCTURAL,**kw):return(self.g.input(key,label,kind,**kw),role)
    def cat(self,key,label):return self.stock(key,label,K.CATALYST,R.CATALYST)
    def cond(self,key,label):return self.stock(key,label,K.CONDITION,R.CONDITION)
    def op(self,key,label,kind=S.MAGIC,inputs=(),product='Actual source preparation; whole occupied substrates/products/side products and conversion unassembled',**kw):
        y=self.g.product(key+'_product',product,K.MOLECULAR)
        if kind==S.MAGIC:kw.setdefault('reason','Positive source operation; whole occupied molecules, end distributions, detailed molecular products and conversion unassembled, not a source-absence claim.')
        if kind==S.GENERIC:kw['changes_molecule']=False
        self.g.react(key,label,kind,((self.x,R.TEMPLATE),)+tuple(inputs),(y,),**kw);self.x=y;return y
    def clean(self,key,label,inputs=()):
        y=self.g.product(key+'_kept','Actual retained physical source fraction, carryover/recovery unassigned',K.MOLECULAR);w=self.g.product(key+'_waste','Actual excluded workup fraction',K.WASTE)
        self.g.react(key,label,S.CLEANUP,((self.x,R.TEMPLATE),self.cond(key+'_condition',label))+tuple(inputs),(y,w));self.g.discard(w,'Actual source workup fraction excluded; no quantitative purification claim.');self.x=y;return y
    def aliquot(self,key,label,kept='Selected source aliquot',rest='Untransferred stock, amount/fate unassigned',qc=False):
        y=self.g.product(key+'_selected',kept,K.MOLECULAR);r=self.g.product(key+'_remainder',rest,K.MOLECULAR)
        self.g.react(key,label,S.SPLIT,((self.x,R.TEMPLATE),),(y,r))
        if qc:self.g.mark_terminal(Terminal(y.id,Disposition.QC,kept));self.x=r;return r
        self.g.mark_terminal(Terminal(r.id,Disposition.UNRESOLVED,rest));self.x=y;return y




from adapter_ligation import AdapterLigationGeometry
from illumina import P5,P7,TRUSEQ_READ1
AUTHORITY='PublishedPMC6821191 FICCMethods/introductionFig1/computation plusoriginalSupplementaryMethodsDOCXfullselectedrecipe. EndogenousHEK293TRMT2A/5FU only; noHAP1/iCLIP/mutantmiCLIP/LCMSrecipeinheritance. Twooriginals; cappedpartialnotabsenceaudit.'
SOURCE_ORDERS={'l3':'rAppAGATCGGAAGAGCGGTTCAG/ddC/','rt':'5Phos/NNXXXXNNNAGATCGGAAGAGCGTCGTGGATCCTGAACCGC','bamhi':'GTTCAGGATCCACGACGCTCTTCAAAA','p5':'AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT','p3':'CAAGCAGAAGACGGCATACGAGATCGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT'}
L3='AGATCGGAAGAGCGGTTCAG';RT_FIXED=revcomp(TRUSEQ_READ1[-18:])+'GGATCC'+'TGAACCGC'
# Literal curator-reviewed inventories: adding graph targets does not extend review.
REVIEWED_TRAP = '''source fu_reagent fu_conditions lysis_buffer turbo_dnase rnase_i
clarify_condition protein_g_beads anti_trmt2a antibody_buffer ip_wash_condition
fu_exposure lysis nuclease clarify prepare_antibody_beads ip ip_wash'''.split()
REVIEWED_RNA = '''pnk3_enzyme pnk3_conditions pnk3_wash_condition l3_adapter l3_ligase
l3_conditions l3_wash_condition radio_pnk_enzyme gamma_atp radio_conditions
radio_supernatant_condition heat_conditions gel_transfer_conditions pk pk_urea_conditions
phenol_condition rna_glycoblue rna_precip_condition pnk3 pnk3_wash l3_ligation l3_wash
radio_pnk radio_supernatant heat_elute gel_transfer membrane_select pk_urea phenol rna_precip'''.split()
REVIEWED_DNA = '''rt_primer rt_dntp ssii rt_conditions cdna_glycoblue cdna_precip_condition
cdna_gel_conditions cdna_extract_condition size_glycogen size_precip_condition circligase
circle_conditions bamhi_oligo bamhi_anneal_conditions bamhi linear_glycoblue
linear_precip_condition p5 p3 accuprime_master pcr_conditions sybr library_gel_conditions
library_extract_condition library_glycogen library_precip_condition rt_anneal rt_components
rt cdna_precip cdna_gel cdna_band cdna_extract size_precip circle bamhi_anneal bamhi_digest
linear_precip pcr_input accuprime_components pcr library_gel library_band library_extract
library_precip loading load library'''.split()
def rt_parts():
    return [Segment('umi_first','NN',feature=feature('ficc.molecule','umi','random',group='ficc.molecule',note='Sourcefive randomdeduplicationbases,firsttwohere')),Segment('sample_index','XXXX',placeholder=True,feature=feature('ficc.sample','sample_index','fixed',note='Sourcefixed4ntsampledemultiplexingcode; actualsequenceunassigned')),Segment('umi_last','NNN',feature=feature('ficc.molecule','umi','random',group='ficc.molecule',note='Remainingthree randomdeduplicationbases')),Segment('fixed_adapter_prefix',RT_FIXED[:6]),Segment('bamhi_binding_window',RT_FIXED[6:29]),Segment('fixed_rt_suffix',RT_FIXED[29:])]
def oligo_scene(name):
    s=Scene()
    if name=='l3':
        s.strand('oligo',[Segment('l3_body',L3),Segment('additional_ddC','C')],mod5='App',mod3='ddC');s.ligation_end('oligo','l3_body',kind='single-strand',polarity='5-prime',owner='adapter',text='Sourcepreadenylated5primeL3donor; terminaladditionalddCnucleotideat3prime. WholeoccupiedRNA/sealedjunctionunassembled.')
    elif name=='rt':s.strand('oligo',rt_parts(),mod5='P')
    elif name=='p5':s.strand('oligo',[Segment('canonicalP5',P5),Segment('canonicalRead1body',TRUSEQ_READ1[4:])])
    elif name=='p3':s.strand('oligo',[Segment('canonicalP7',P7),Segment('sourcecustombody',SOURCE_ORDERS['p3'][len(P7):])])
    else:s.strand('oligo',[Segment('bamhi_complement',SOURCE_ORDERS['bamhi'][:23]),Segment('sourceunpairedtail',SOURCE_ORDERS['bamhi'][23:])])
    return s

def local_bamhi_interface():
    s=Scene();s.strand('primer_derived_local_cdna',[Segment('bamhi_binding_window',RT_FIXED[6:29])]);s.note('primer_derived_local_cdna','Local23ntRT-derivedtargetwindowonly; omittedknownRTflanksremaininthewholeRTinventory. SourceBamHIoligoAAAAisnotforcedtopairwiththeRTprefix; fullcircleoccupiedgeometryunassembled.')
    s.anneal('bamhi_oligo',[Segment('bamhi_complement',SOURCE_ORDERS['bamhi'][:23]),Segment('sourceunpairedtail',SOURCE_ORDERS['bamhi'][23:])],to='primer_derived_local_cdna',pair=('bamhi_complement','bamhi_binding_window'))
    return s

def oligo_input(b,key,name):
    chem={'source_ordering':SOURCE_ORDERS[name]}
    if name=='l3':chem.update({'five_prime':'adenylated phosphate','three_prime':'ddCblock','terminal_nucleotide':'additionalphysicalC','source_polymer':'DNA20nt;5adenylationnotanextraRNApolymerresidue'})
    if name=='rt':chem.update({'five_prime':'P','identifier_scope':'N2+N3UMI; X4fixedsamplecode, actualsampleunassigned'})
    return b.g.input(key,'SourceFICC '+name,K.OLIGO,state=oligo_scene(name),chemical_state=chem)

def protocol_graph():
    b=B('FICC-Seq selectedendogenousHEK293TRMT2A catalytictrap preparation',AUTHORITY,'SourceexponentiallygrowingHEK293cellswithendogenousTRMT2A; noengineeredmutantinput',K.BIOLOGICAL);g=b.g
    g.add_source_item(SourceItem('catalytic_trap','MainFig1/introduction5FU→FUTPincorporation/normalmethyltransfer→failureofreleasefrombaseC6.','Intracellularpathways/adductatommap/occupiedsites/convertedfractionunassembled. NoUV/AzaC/mutantmiCLIPorliteralall-RNAmodificationclaim.'))
    g.add_source_item(SourceItem('source_units','Supplementtwiceprints400mlTEforgelcrushing/extractionfollowedbycollectioninto1.5mltube; otherworkupvolumesinul.','Preserveliteral400ml; executedvolumeunassigned, no silent400ulrepair. Listedligationmix12ulnotnormalizedbyinventedwater.'))
    adduct=g.product('fu_trapped_cells','Source5FU-exposedcell/RNA/TRMT2Aadductclasses; fullsites/conversion/offtargetsunassembled',K.MOLECULAR,chemical_state={'reaction_class':'SourceTRMT2A-RNAcovalentcatalyticintermediateatincorporatedfluorouridineresidue','enzyme_provenance':'Endogenousproteininsideactualsourcecells, notfreshpurifiedenzymeaddition','metabolic_substrate':'5FU-derivedFUTPincorporationintosubsetofnewRNA; allproductsunassembled'})
    g.react('fu_exposure','Source100uM5FU24hthenharvest',S.MAGIC,((b.x,R.TEMPLATE),b.stock('fu_reagent','Actual100uM5FUmediumcontribution'),b.cond('fu_conditions','Source24hcelltreatment')),(adduct,),reason='Positiveintracellularcatalytictrap; metabolicconversion/fullatommap/occupancy/yieldunassembled. NotordinaryUVcrosslinking.');b.x=adduct
    b.op('lysis','Source0.8mllysisbuffer4C30minintermittentpipetting',inputs=(b.cond('lysis_buffer','SourceTris7.4/100mMNaCl/NP40/SDS/deoxycholate'),))
    b.op('nuclease','Source2ulTurboDNase+2ulRNaseI1to200/37C3min/ice5min',S.DIGESTION,(b.cat('turbo_dnase','Actual2ulTurboDNaseAM2239'),b.cat('rnase_i','Actual2ulRNaseI1to200AM2295')));b.clean('clarify','Source13000rpm15min4Cinsolublematerialexcluded')
    raw_beads=g.input('protein_g_beads','Actual100ulProteinGDynabeadstock; physicalsupport/proteinpreparationsubstrate',K.MOLECULAR);ab=g.input('anti_trmt2a','Actual5ugantiTRMT2AOTI1C3antibody',K.MOLECULAR);loaded=g.product('prepared_antibody_beads','ActualwashedProteinG/antiTRMT2Apreparedbeads',K.SELECTION);wash=g.product('antibody_wash','Actualbeadpreparationwashfractions',K.WASTE)
    g.react('prepare_antibody_beads','Sourcebeadwash3x/antibody0.5mlbuffer30min4C/wash3x',S.MAGIC,((raw_beads,R.TEMPLATE),(ab,R.STRUCTURAL),b.cond('antibody_buffer','Actualsource0.5mllysisbuffer')),(loaded,wash),reason='Sourcephysicalaffinitybeadpreparation; stoichiometry/occupiedconjugatesunassembled.');g.discard(wash,'Actualantibodybeadpreparationwashfractions')
    bound=g.product('ip_bound','ActualselectedTRMT2A-RNA/beadpreparation',K.MOLECULAR);sup=g.product('ip_supernatant','ActualexcludedIPsupernatant',K.WASTE);g.react('ip','SourceendogenousTRMT2AIP2h4Cgentlerotation',S.AFFINITY_CAPTURE,((b.x,R.TEMPLATE),(loaded,R.SELECTION)),(bound,sup));g.discard(sup,'ActualunboundIPfractionexcluded, nopurecapture/yieldclaim');b.x=bound;b.clean('ip_wash','Source5xlysisbuffer/2xPNKwash')
    b.op('pnk3','Source20ulPNKpH6.5noaddedATP37C20min',S.DEPHOSPHORYLATION,(b.cat('pnk3_enzyme','Actual.5ulT4PNK'),b.cond('pnk3_conditions','Actual4ul5xPNK6.5/.5ulRNaseinhibitor/15ulwater')));b.clean('pnk3_wash','Source3xPNKwash')
    l3=oligo_input(b,'l3_adapter','l3')
    b.op('l3_ligation','Source5ul4xbuffer/1ulligase/.5ulinhibitor/1.5ul20uML3/4ulPEG400/16Covernight1100rpm',S.LIGATION,((l3,R.STRUCTURAL),b.cat('l3_ligase','ActualT4RNAligaseNEB, exactsubtypeunassigned'),b.cond('l3_conditions','Source4xbufferTris/Mg/DTT/PEG400; listed12ul, no furtherwateraddedinmodel')),geometry=AdapterLigationGeometry.known('single-strand',oligo_scene('l3'),statement='Sourceknown5App/3ddCL3localend; healedRNA3OHconditional, wholeRNA/adduct/sealedjunctionunassembled.'));b.clean('l3_wash','Source3xPNKwash')
    b.op('radio_pnk','SourcefreshPNK/gamma32PATP37C5min5endlabel',S.PHOSPHORYLATION,(b.cat('radio_pnk_enzyme','Actual1ulPNK'),b.stock('gamma_atp','Actual2ulgamma32P-ATPphosphatedonor'),b.cond('radio_conditions','Actual2ul10xPNKbuffer/15ulwater')));b.clean('radio_supernatant','Actualsupernatantremovedbeforeloadingbuffer')
    eluate=g.product('gel_input','Actualheat-elutedRNPpreparation',K.MOLECULAR);beads=g.product('heat_residual_beads','Actualresidualaffinitybeads/material',K.MOLECULAR);g.react('heat_elute','Source20ul1xNuPAGEloadingbuffer70C5min',S.MAGIC,((b.x,R.TEMPLATE),b.cond('heat_conditions','Actualsourceheatloadingbuffer')),(eluate,beads),reason='Sourcegel-loadedeluate; exactbeadretention/recoveryandwholeRNPunassembled.');g.mark_terminal(Terminal(beads.id,Disposition.UNRESOLVED,'Residualbeadfateoutsidepartialtrace'));b.x=eluate
    b.op('gel_transfer','Source4to12percentBisTris/MOPS→nitrocellulosetransfer/filmexposure',inputs=(b.cond('gel_transfer_conditions','Actualgel/membrane/autoradiographyprogramme'),));selected=g.product('membrane_slice','ActualexcisedradioactiveTRMT2A-RNPnitrocellulosepiece',K.MOLECULAR);excluded=g.product('membrane_rest','Unselectedmembrane/othermaterial',K.MOLECULAR);g.react('membrane_select','Sourceradioactivity-guidedRNPmembraneexcision',S.SPLIT,((b.x,R.TEMPLATE),),(selected,excluded));g.mark_terminal(Terminal(excluded.id,Disposition.UNRESOLVED,'Unselectedmembranefateunassigned'));b.x=selected
    b.op('pk_urea','Source200ugPK200ul37C20minthen200ul7Mureamix/further20min',S.DEPROTEINIZATION,(b.cat('pk','Actual200ugPK'),b.cond('pk_urea_conditions','SourceTris/NaCl/EDTA; secondadditiongives3.5Murea')),product='ReleasedRNA/adapter/peptideremnantclasses; completeadductproteolysis/yieldunassembled')
    b.clean('phenol','Source400ulsolution+400ulphenol/chloroform30C5min/spin; retainaqueousphase')
    b.clean('rna_precip','Source40ul3MNaOAc5.5/1mlEtOH/overnightminus20/pellet/wash/resuspend5ul',inputs=(b.stock('rna_glycoblue','Actual1ulGlycoBluecarrier'),))
    primer=oligo_input(b,'rt_primer','rt');nt=g.input('rt_dntp','Actual1uldNTPmixNEB',K.MOLECULAR)
    mix=g.product('rt_annealed_mixture','ActualRNA+sourceRTprimer+dNTPmixtureafter65C5min/ice5min',K.MOLECULAR);g.react('rt_anneal','SourceactualRNA/primer/dNTPpreannealing',S.MAGIC,((b.x,R.TEMPLATE),(primer,R.STRUCTURAL),(nt,R.STRUCTURAL)),(mix,),reason='SourcewholeRNA/adapterprimercomplexunassembled; allactualcontributingmoleculesretainedbeforeheating.')
    template=g.product('rt_template','Templatecomponentofsameannealedmixture',K.MOLECULAR);retained=g.product('rt_retained_primer','SameactualexternalRTprimercomponentalreadyannealed',K.OLIGO,state=oligo_scene('rt'),chemical_state={'source_ordering':SOURCE_ORDERS['rt'],'provenance':'sameexternalprimeralreadyadded, notnewstock'});retained_nt=g.product('rt_retained_nt','SameactualdNTPcomponentalreadyadded',K.MOLECULAR);g.react('rt_components','Same-tubecomponentbookkeeping,nophysicalsplitorfreshprimer/NTaddition',S.GENERIC,((mix,R.TEMPLATE),),(template,retained,retained_nt),changes_molecule=False);b.x=template
    b.op('rt','SourceSSII25C5min/42C20min/50C40min/80C5min',inputs=((retained,R.STRUCTURAL),(retained_nt,R.STRUCTURAL),b.cat('ssii','Actual.5ulSuperscriptII'),b.cond('rt_conditions','Actual7ulwater/.5ulinhibitor/1ul.1MDTT/4ul5xbuffer')),product='SourceRT-truncatedcDNA/RNA/peptideremnantclasses; fulloccupiedprimingandproductjunctionunassembled')
    b.clean('cdna_precip','Source380ulwater/40ulNaOAc/1mlEtOH/overnightminus20/pellet/wash/resuspend6ul',inputs=(b.stock('cdna_glycoblue','Actual1ulGlycoBluecarrier'),))
    b.op('cdna_gel','Source2ul2xTBEureabuffer/70C3min/6percentTBEurea40min180V/select70to150nt',inputs=(b.cond('cdna_gel_conditions','Actualsourceprogramme; fullRTproductendsunassembled'),))
    b.aliquot('cdna_band','Actual70to150ntcDNAgelbandexcision','Actualselected70to150ntcDNAgelpiece','UnselectedcDNAgelbands/material; fateunassigned')
    b.clean('cdna_extract','Sourcegelcrush2h37C1100rpm/SpinXglasswoolflowthrough; literal400mlTE/executedvolumeunassigned')
    b.clean('size_precip','Source40ulNaOAc/1mlEtOH/overnightminus20/pellet/wash/resuspend6.5ul',inputs=(b.stock('size_glycogen','Actual1ulglycogencarrier'),))
    b.op('circle','Source6.5ulcDNA+.8ulbuffer/.4ul50mMMn/.3ulCircLigaseII→8ul60C1h',inputs=(b.cat('circligase','Actual.3ulCircLigaseII'),b.cond('circle_conditions','Actualsource8ulprogramme')),product='SourcecircularcDNA/unreacted/offproductclasses; wholecirclejunctionunassembled')
    bam=oligo_input(b,'bamhi_oligo','bamhi');local_bamhi_interface()
    b.op('bamhi_anneal','Source26ulwater/3ulFastDigestbuffer/1ul10uMBamHIoligo;95C2min/70stepsminus1C',inputs=((bam,R.STRUCTURAL),b.cond('bamhi_anneal_conditions','Actualsourceannealingprogramme; local23bpduplexconstructorverified')),product='Sourceannealedcircle/BamHIoligoclasses; localRT-derivedsiteknown, wholecircleunassembled')
    b.op('bamhi_digest','Source2ulFastDigestBamHI37C30min',S.DIGESTION,(b.cat('bamhi','Actual2ulFastDigestBamHI'),),product='SourcelinearizedcDNA/unreacted/offproductclasses; fullcirclecut/retainedoligodistributionunassembled')
    b.clean('linear_precip','Source350ulTE/40ulNaOAc/1mlEtOH/overnightminus20/pellet/wash/resuspend20ul',inputs=(b.stock('linear_glycoblue','Actual1ulGlycoBluecarrier'),))
    b.aliquot('pcr_input','Actual10of20ullinearizedcDNAforPCR','Actual10ulPCRtemplate','Actual10uluntransferredlinearizedcDNAstock')
    pf=oligo_input(b,'p5','p5');pr=oligo_input(b,'p3','p3');master=g.input('accuprime_master','Actual20ulAccuprimeSupermix',K.MOLECULAR);enz=g.product('accuprime_enzyme','SameactualSupermixenzymecomponent',K.CATALYST);dn=g.product('accuprime_nt','SameactualSupermixnucleotidecomponent',K.MOLECULAR);buf=g.product('accuprime_buffer','SameactualSupermixbuffercomponent',K.CONDITION);g.react('accuprime_components','Actualmastercomponentbookkeeping, nofreshseparateadditions',S.GENERIC,((master,R.TEMPLATE),),(enz,dn,buf),changes_molecule=False)
    b.op('pcr','Source40ulPCR94C2min/25cycles94C15s65C30s68C30s/68C3min',inputs=((pf,R.STRUCTURAL),(pr,R.STRUCTURAL),(enz,R.CATALYST),(dn,R.STRUCTURAL),(buf,R.CONDITION),b.cond('pcr_conditions','Actual1uleach10uMprimer/8ulwater')),product='SourcePCRlibrary/unreacted/offproducts; fullduplex/readlandingunassembled')
    b.op('library_gel','Source6percentTBEgel/SYBRGreenIstain/select150to200nt',inputs=(b.stock('sybr','ActualSYBRGreenIstainmolecularcontribution'),b.cond('library_gel_conditions','Actualsourcebandselection')),product='SourceSYBRstainedselectedlibraryband; dye-carryoverunassigned')
    b.aliquot('library_band','Actual150to200ntlibrarybandexcision','Actualselected150to200ntlibrarygelpiece','Unselectedlibrarygelbands/material; fateunassigned')
    b.clean('library_extract','Sourcegelcrush/SpinXglasswoolflowthrough; literal400mlTE2h37C/executedvolumeunassigned')
    b.clean('library_precip','Source40ulNaOAc/1mlEtOH/overnightminus20/pellet/wash/resuspend20ul',inputs=(b.stock('library_glycogen','Actual1ulglycogencarrier'),))
    b.aliquot('loading','ActualHiSeq4000loadingallocation; amountunassigned')
    lib=g.product('library','SourceFICCSeqDNAsequencinglibraryclass; wholemolecules/readlandingunassembled',K.LIBRARY);g.react('load','SourceHiSeq4000singleend50cycles',S.GENERIC,((b.x,R.TEMPLATE),),(lib,),changes_molecule=False);g.sequence(lib,SequencingPlan(read_lengths={'R1':50},limitation='PublishedSupplementaryMethodssingleend50cycles. N2+N3randomUMI/X4samplecode sourcequalified, actualfixedcode/readprimerbindingandwholeduplexunassembled. Truncation/readassignmentcomputationalnotphysicalbaseconversion.'))
    g.cover_source_item(SourceCoverage('catalytic_trap','modelled','fu_exposure'))
    g.cover_source_item(SourceCoverage('source_units','modelled','cdna_extract'))
    g.declare_reviewed_scope(REVIEWED_TRAP,'Published PMC6821191 main52/56/78 and SupplementaryMethods FICC opening: endogenous HEK293/5FU, actual prepared antibody bead capture, not mutant/UV/iCLIP comparison.')
    g.declare_reviewed_scope(REVIEWED_RNA,'Published SupplementaryMethods FICC: source PNK/ligation/radioactive donor, RNP membrane selection, proteolysis and recovered RNA; whole adducts/recovery unassembled.')
    g.declare_reviewed_scope(REVIEWED_DNA,'Published SupplementaryMethods FICC: five whole source oligos, RNA/primer/dNTP annealing, SSII, size selection/circle/BamHI, 10of20ul PCR, gel selection and actual HiSeq4000SE50 contribution. Literal400ml TE retained; actual extraction volume unassigned.')
    return g

def selftest():
    from checks import Check,Source,have
    from storage import sources_root
    import zipfile,xml.etree.ElementTree as ET
    c=Check();g=protocol_graph()
    c('AppadaptercontainsadditionalphysicalddC', ''.join(x.top for x in oligo_scene('l3').strands['oligo'].segs),L3+'C')
    c('RTpolymerretainsallsourceplaceholderandfixedbases',''.join(x.top for x in oligo_scene('rt').strands['oligo'].segs),SOURCE_ORDERS['rt'].split('/')[1])
    c('CanonicalP5stillwholepublishedbody',''.join(x.top for x in oligo_scene('p5').strands['oligo'].segs),SOURCE_ORDERS['p5'])
    c('CanonicalP7prefixretainswholecustomP3',''.join(x.top for x in oligo_scene('p3').strands['oligo'].segs),SOURCE_ORDERS['p3'])
    c('KnownBamHIlocalduplexactuallyconstructed',len(local_bamhi_interface().strands),2)
    c('ActualRNAenters65Cannealing',g.steps['rt_anneal'].inputs[0].material,'rna_precip_kept')
    c('RTusesproducedprimercomponentnotnewstock',g.steps['rt'].inputs[1].material,'rt_retained_primer')
    c('RTusesproducedNTcomponentnotnewstock',g.steps['rt'].inputs[2].material,'rt_retained_nt')
    c('RadioPNKphosphatedonorisstructural',g.steps['radio_pnk'].inputs[2].role,R.STRUCTURAL)
    c('PKprecedesureaextractionclass',g.steps['pk_urea'].kind,S.DEPROTEINIZATION)
    c('OnlyselectedcDNAbandentersextraction',g.steps['cdna_extract'].inputs[0].material,'cdna_band_selected')
    c('Onlyselectedlibrarybandentersextraction',g.steps['library_extract'].inputs[0].material,'library_band_selected')
    c('PCRtakesselected10ulof20stock',g.steps['pcr'].inputs[0].material,'pcr_input_selected')
    c('Pooledloadingwhole-stockconsumptionavoided',g.steps['load'].inputs[0].material,'loading_selected')
    p=sources_root()/Path(__file__).resolve().parents[1].name/'nar_gkz658_supp_gkz658_supplemental_files_Supplementary_Methods.docx'
    if have(c,Source(p,'https://doi.org/10.1093/nar/gkz658')):
        with zipfile.ZipFile(p)as z:text=''.join(ET.fromstring(z.read('word/document.xml')).itertext())
        for name,raw in SOURCE_ORDERS.items():c('OriginalDOCXwholemodifiedordering '+name,raw in text,True)
        c('Literal400mlTEsourcewordingretained','400ml TE buffer' in text,True)
        c('SourceSE50','single end using 50 cycles' in text,True)
    c('Source-reviewed preparation satisfies normal public gate',g.publication_problems(),[])
    g.input('unreviewed_addition','Regression-only stock not in curated source inventory',K.CONDITION)
    c('New source input cannot auto-publish',any('unreviewed_addition' in x and 'source-inventory' in x for x in g.publication_problems()),True)
    c.report()

PAGE_FILE = 'ficc_seq.html'
def render_page():
    from page import head,info
    from chemdraw import panel
    from protocol_views import assemble_page,WorkflowSection,PlanLimitationSection
    g=protocol_graph().require_publishable()
    parts=[head('FICC-Seq'),'<div class="wrap" data-protocol-complete="true"><h1>FICC-Seq</h1>',
      '<p>Fluorouracil-induced catalytic crosslinking captures RNA substrates of endogenous TRMT2A. Metabolically incorporated fluorouridine traps the enzyme–RNA catalytic intermediate; reverse transcription can stop at the recovered peptide–RNA crosslink.</p>',
      info('Source: <a href="https://doi.org/10.1093/nar/gkz658">FICC-Seq: a method for enzyme-specified profiling of methyl-5-uridine in cellular RNA</a>, published article and Supplementary Methods. This schematic selects endogenous TRMT2A in HEK293 cells treated with 100 µM 5-fluorouracil for 24 hours. HAP1, mutant miCLIP, UV-iCLIP and LC–MS comparisons are separate branches; their recipes are not substituted here.'),
      '<h2>Source volume ambiguity</h2>',info('The Supplementary Methods print “400ml TE buffer” twice: after the 70–150-nt cDNA gel and after the final 150–200-nt library gel. Each extraction subsequently collects flow-through into a 1.5-ml tube. The model preserves 400 ml as printed and leaves the executed volume unresolved; it does not silently replace it with 400 µl. The listed ligation additions total 12 µl, with no extra water invented.'),
      '<h2>Known oligos and identifiers</h2><p>The preadenylated L3 adapter ends in an additional physical ddC nucleotide. The 5′-phosphorylated RT primer carries a split five-base UMI (N2 + N3) and a fixed four-base sample code (X4). The actual sample-code value is unassigned. The BamHI oligo pairs with a known local RT-derived window; a complete occupied circle is not inferred.</p>']
    for name,caption in [('l3','Source L3 adapter: 5′ adenylation and terminal physical ddC.'),('rt','Source RT primer: N2–X4–N3 identifier arrangement and 5′ phosphate.'),('bamhi','Source BamHI annealing oligo; its terminal AAAA is retained.'),('p5','Source P5 amplification primer.'),('p3','Source custom P3 amplification primer.')]:
        parts.append(panel(oligo_scene(name).rows(),caption=caption))
    parts += [panel(local_bamhi_interface().rows(),caption='Known local 23-bp BamHI-oligo pairing; whole circle, remote flanks and occupied products remain unassembled.'),
      '<h2>Preparation and reads</h2><p>Trim lysate RNA and capture endogenous TRMT2A complexes on antibody beads. Dephosphorylate RNA 3′ ends, ligate L3 and radiolabel 5′ ends with fresh PNK and gamma-32P-ATP. Resolve the complexes, transfer to nitrocellulose and excise the radioactive region. Proteolyze, recover RNA and reverse transcribe with the source primer. Select 70–150-nt cDNA, circularize, anneal the BamHI oligo and linearize. Use 10 µl of the 20-µl stock for PCR, select the 150–200-nt library band and allocate material for sequencing.</p>',
      info('Intracellular 5FU metabolism, adduct occupancy, peptide remnants, complete RNA/cDNA/circle structures, ligation junctions, yields and primer landing remain unresolved. The local duplex checks complementarity only. The source library is sequenced single-end for 50 cycles on HiSeq 4000; exact identifier cycle spans and read orientation are not computed from a complete construct. Computational duplicate removal does not physically purify the preparation.'),PlanLimitationSection('library'),
      '<details><summary>Full preparation and chemistry warnings</summary>',WorkflowSection('Published endogenous TRMT2A preparation'),'</details>',
      '<p><a href="01_ficc-seq.html">Source notes and model limits</a></p></div>']
    return assemble_page(g,parts)
if __name__=='__main__':selftest()
