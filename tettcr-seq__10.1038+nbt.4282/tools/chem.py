from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Construct, Scene, Segment, feature, workflow_from_protocol_graph, workflow_panel
from protocol_graph import ProtocolGraph, MaterialKind as K, InputRole as R, StepKind as S, SourceItem, SequencingPlan
from page import head,info
PAGE_FILE='chem.html'

def operations(g):
    def op(key,x,label,reagent=None,state=None,kind=S.MAGIC,reagent_kind=K.MOLECULAR,reagent_role=R.STRUCTURAL):
        inputs=[(x,R.TEMPLATE)]
        if reagent:
            z=g.input(key+'_reagent',reagent,reagent_kind);inputs.append((z,reagent_role))
        y=g.product(key+'_product',label+'; preparation includes unremoved reaction components',K.MOLECULAR,state=state)
        kw={'reason':'Actualsourceoperation; occupiedmolecules/internalrecipe unassembled, no genericrecipe inheritance.'} if kind==S.MAGIC else {}
        if kind==S.GENERIC:kw['changes_molecule']=False
        g.react(key,label,kind,inputs,(y,),**kw);return y
    def clean(key,x,label):
        z=g.input(key+'_selection',label,K.SELECTION)
        y=g.product(key+'_kept','Physically retained source fraction: '+label,K.MOLECULAR)
        w=g.product(key+'_removed','Physically removed fractions: '+label,K.WASTE)
        g.react(key,label,S.CLEANUP,((x,R.TEMPLATE),(z,R.SELECTION)),(y,w));g.discard(w,'Source physical workup only; no perfect specificity or purity assertion.');return y
    return op,clean

def end(g,x,label):g.sequence(x,SequencingPlan(acquisition='unspecified',limitation=label))
def incomplete(g):
    g.declare_reviewed_scope(('linker','fluor_streptavidin','conjugation_supplies','conjugate','conj_cleanup_selection','conj_cleanup','reporter_template','anneal_condition','encoding_nt','anneal_reporter','encoding_mix','encoding_components','encoding_enzyme','encoding_buffer','encode','stop_extension_reagent','stop_extension','ivtt_template','ivtt_f','ivtt_r','ivtt_pcr_supplies','ivtt_pcr_components','ivtt_polymerase','ivtt_nt','ivtt_buffer','ivtt_pcr','ivtt_water','ivtt_dilute','ivtt_allocation','ivtt_stock','ivtt_mix_allocation','ivtt','pmhc','uv_condition','uv_exchange','tetramerize','panel_allocate','other_panel','panel_pool','stain_buffer','panel_concentrate_selection','panel_concentrate'),SOURCE+'; fixed selectedfactory/panel targets only; unreadlinkedTCR/Solulink/pMHC/fullgeometry scope notcovered.')
    g.add_source_item(SourceItem('unassembled','Unreviewed linked originals and full molecular geometry',LIMIT));return g


from dataclasses import dataclass
import illumina as il
from protocol_graph import Terminal,Disposition
TITLE='TetTCR-seq — DNA-encoded pMHC binding and paired TCR recovery'
SOURCE='Published university mirror Zhang_Nat_Biotech_2018.pdf printedDOI/date/title andOnlineMethods973–1064;publishedTable10MOESM15byte-identicalaliasofNIHMSsupp11;publishedSI MOESM4Fig3/controlthreshold'
LIMIT='Private partial source model. TemplateDNA feeds physically separate IVTT and reporter-encoding branches;DNA barcodes are not made by reverse-transcribing antibodies. Linker5amineC12/iSp18 andsynthetic12NreporterMIDremain distinctfromwelltags and5baseNHNHNwindows. ExactSolulink/fulltetramer/duplex/TCR-parent/readgeometryunassembled.'
LINKER_HEAD='TAGTACTCAGAGGTTGATCTACATTG'
LINKER_SEED='GACGATGACGACAAG'
TEMPLATE='CGAGGTGCTTCGTTAAACAGCGTTGATACCCAGAGCAACCAGTTTCTTGTCGTCATCGTC'
P1='TCAGAGGTTGATCTACATTG'
P2='AGCGAGGTGCTTCGTTA'
@dataclass(frozen=True)
class LinkerDNA:
    attached:bool=False
    encoded:bool=False
    def scene(self):
        sc=Scene();sc.strand('linker',[Segment('source5DNAflank',LINKER_HEAD),Segment('synthetic reporterMID','N'*12,placeholder=True,feature=feature('tettcr.reporter_mid','umi','random',note='Synthetic linkerMID counts reporter species;notTCR-RNAUMI orcellidentity.')),Segment('source3annealingseed',LINKER_SEED)],label='Original sourceDNA linker strand',mod5=('protein–linker attachment unassembled; iSp18 spacer' if self.attached else 'amine-C12; iSp18 spacer'))
        sc.note('linker','iSp18 is a nonnucleotide spacer,not18DNA bases. Original5amine isattachmentreactant,not assertedfreeafterconjugation. Encodedreporter='+str(self.encoded)+'.')
        if self.encoded:sc.strand('encoded',[Segment('encodedHCVpeptideDNA/copyregion','[peptide-associated DNA andcopiedduplexextentunassembled]',placeholder=True,feature=feature('tettcr.hcv_klv','feature_barcode','unknown'))],label='SourceKlenowencodesreporter;wholecopygeometrynotreconstructed')
        sc.linker_dna=self;return sc

def template():return Construct([Segment('sourceHCV-KLV(WT)DNAtemplate',TEMPLATE,feature=feature('tettcr.hcv_template','feature_barcode','fixed',note='Source peptide-encoding DNA template example;not whole protein sequence.'))])
def well_primer(which):
    if which=='f':head='GACGTGTGCTCTTCCGATCT';bc='ATCACG';arm='TACTCAGAGGTTGATCTACATTG'
    elif which=='r':head=il.TRUSEQ_READ1[4:];bc='TCAAG';arm=P2
    else:raise ValueError(which)
    return Construct([Segment('sourcepartialIlluminahandle',head),Segment('sourceNHNHNclass','NHNHN',placeholder=True,note='Printed source pattern retained;function unassigned,not labelled another UMI.'),Segment('sourcewellbarcodeexample',bc,feature=feature('tettcr.well_'+which,'cell_barcode','fixed',note='Source primer member1;experimentalwellassignmentunassigned.')),Segment('sourcebindingarm',arm)])
def illumina_primer(which):
    if which=='f':return Construct([Segment('sourceILLU_f',il.P5+'TCTTTCCCTACACGAC')])
    if which=='r':return Construct([Segment('sourceP7andfixedAA',il.P7+'AA'),Segment('sixXsourcewindow','N'*6,placeholder=True,feature=feature('tettcr.library_index','sample_index','unknown',note='Printed6Xwindowonly;nofullacquiredindexboundaryassigned.')),Segment('sourceIlluminabindingarm','GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT')])
    raise ValueError(which)

def protocol_graph():
    g=ProtocolGraph(TITLE,authority=SOURCE,representation='published_instance');op,clean=operations(g)
    def select(k,x,label):
        take=g.product(k+'_selected',label,K.MOLECULAR);rest=g.product(k+'_remainder','Untransferred source stock; amount/fateunassigned',K.MOLECULAR)
        g.react(k,label,S.SPLIT,((x,R.TEMPLATE),),(take,rest));g.mark_terminal(Terminal(rest.id,Disposition.UNRESOLVED,'Actualuntransferredstock.'));return take
    linker=g.input('linker','Source5amineC12/iSp18DNA linker',K.OLIGO,state=LinkerDNA().scene())
    protein=g.input('fluor_streptavidin','SourcePE-orAPC-streptavidin;conjugationkit-specificgeometryunreviewed',K.MOLECULAR)
    conj=g.product('linker_conjugate','SourceDNA–fluorescentstreptavidinconjugatesplusunreactedcomponents',K.MOLECULAR,state=LinkerDNA(True).scene())
    conjugation=g.input('conjugation_supplies','ActualsourceSolulinkconjugation supplies; identities/doses/revisionunreviewed, no reactionrecipeborrowed',K.MOLECULAR)
    g.react('conjugate','SourceSolulinkcoupling;fullsite/chemistryunreviewed',S.MAGIC,((protein,R.TEMPLATE),(linker,R.STRUCTURAL),(conjugation,R.STRUCTURAL)),(conj,),reason='Sourceconjugationreported;manufacturerrecipeandexactsiteunknown,notinvented.')
    conj=clean('conj_cleanup',conj,'SixVivaspin100kDawashes,concentrate,0.2umfilter;retainconjugatefraction')
    reporter_template=g.input('reporter_template','SourceHCV-KLV(WT)DNAoligoaliquotforreporterencoding',K.OLIGO,state=template())
    annealed=g.product('annealed_reporter','Source linker/template annealingpreparation;exactcompletegeometryunassembled',K.MOLECULAR)
    conditions=g.input('anneal_condition','55C5min→25Catminus0.1C/s;actualCutsmartbuffer',K.CONDITION)
    nt=g.input('encoding_nt','Actual250uMdNTPs addedBEFORElinker/templateannealing',K.MOLECULAR)
    retained_nt=g.product('encoding_retained_nt','SameactualdNTPcomponent carriedthroughannealing, notnewaddition',K.MOLECULAR)
    g.react('anneal_reporter','AnnealpeptideDNAtoattachedlinker withcopresentdNTPs',S.HYBRIDIZATION,((conj,R.TEMPLATE),(reporter_template,R.STRUCTURAL),(nt,R.STRUCTURAL),(conditions,R.CONDITION)),(annealed,retained_nt))
    mix=g.input('encoding_mix','Actualsource1ulextensionmixture:.1ul10xCutsmart+.125ulKlenowexo-minus;remainingvolumeunassigned',K.MOLECULAR)
    enzyme=g.product('encoding_enzyme','CopresentactualKlenowexo-minus catalyst',K.CATALYST);buffer=g.product('encoding_buffer','Copresentactualsourcebuffer;37C1h',K.CONDITION)
    g.react('encoding_components','Sameactualextensionmixcomponents, noextra doses',S.GENERIC,((mix,R.TEMPLATE),),(enzyme,buffer),changes_molecule=False)
    encoded=g.product('encode_product','Actualsource encodedreporterconjugatepreparation; fulloccupiedduplex/junctionunassembled',K.MOLECULAR,state=LinkerDNA(True,True).scene())
    g.react('encode','Sourceexo-minusKlenow37C1h extendsattachedlinker',S.MAGIC,((annealed,R.TEMPLATE),(retained_nt,R.STRUCTURAL),(enzyme,R.CATALYST),(buffer,R.CONDITION)),(encoded,),reason='Actualsourceextension; wholeprimer/template/productgeometryunassembled, nofreshdNTPaddition.')
    encoded=op('stop_extension',encoded,'StopsourceextensionwithEDTA;no physicalpurification','ActualEDTAmolecularquench',state=LinkerDNA(True,True).scene(),reagent_kind=K.MOLECULAR,reagent_role=R.STRUCTURAL)
    ivtt_dna=g.input('ivtt_template','PhysicallyseparatealiquotofsamepeptideDNAidentityforIVTT',K.OLIGO,state=template())
    f=g.input('ivtt_f','Actual1uMIVTT_f sourcePCRprimer; wholebodyavailableinTable10/unexpandedhere',K.OLIGO);r=g.input('ivtt_r','Actual1.05uMIVTT_r sourcePCRprimer; wholebodyavailable/unexpanded',K.OLIGO)
    supplies=g.input('ivtt_pcr_supplies','Aggregateactualsource ExTaqHS/dNTP/buffer supplies; stockcount/packagingunassigned',K.MOLECULAR)
    pol=g.product('ivtt_polymerase','ActualsourceExTaqHS.0375Uulcomponent',K.CATALYST);nt=g.product('ivtt_nt','Actualsource400uMdNTPcomponent',K.MOLECULAR);buf=g.product('ivtt_buffer','ActualsourcePCRconditions95C3min/30×95C20s52C40s72C45s/72C5min',K.CONDITION)
    g.react('ivtt_pcr_components','ActualPCRsupplycomponents, notassumedpremix ornewdoses',S.GENERIC,((supplies,R.TEMPLATE),),(pol,nt,buf),changes_molecule=False)
    ivtt=g.product('ivtt_pcr_product','ActualsourceIVTTtemplatePCRDNAplusremainingreagents; wholeampliconunassembled',K.MOLECULAR)
    g.react('ivtt_pcr','Source30cycleIVTTtemplatePCR',S.MAGIC,((ivtt_dna,R.TEMPLATE),(f,R.STRUCTURAL),(r,R.STRUCTURAL),(pol,R.CATALYST),(nt,R.STRUCTURAL),(buf,R.CONDITION)),(ivtt,),reason='ActualsourcePCR; wholeinputprimer/amplicongeometryunassembled.')
    water=g.input('ivtt_water','Actual73.3ulwater dilutionaddition; finaldilutionfactorunassigned',K.CONDITION)
    diluted=g.product('ivtt_diluted','Actual dilutedIVTTtemplatePCRstock',K.MOLECULAR);g.react('ivtt_dilute','SourcedilutePCRproductwith73.3ulwater',S.GENERIC,((ivtt,R.TEMPLATE),(water,R.CONDITION)),(diluted,),changes_molecule=False)
    ivtt=select('ivtt_allocation',diluted,'Actual1uldilutedPCRDNA contribution toIVTT')
    stock=g.input('ivtt_stock','Actualsource20ul1.5xPURExpresspreparedmastermix:10ulA/7.5ulB/.8ulRF1+2+3/.25ulenterokinase/.25ulRNaseinhibitor/1.2ulwater',K.MOLECULAR)
    selected=select('ivtt_mix_allocation',stock,'Actual2ulpreparedPURExpressmastermix; nominal18ulstockrestnotdiscarded')
    peptide=g.product('ivtt_product','ActualcompletedIVTTmixture includesDNAtemplate/RNA/peptide/enzymes/nucleotide/aminoacid/reagentremainders; no perfectconversion',K.MOLECULAR)
    g.react('ivtt','Source1ultemplate+2ulactualIVTTmix,30C4h',S.MAGIC,((ivtt,R.TEMPLATE),(selected,R.STRUCTURAL)),(peptide,),reason='ActualsourceRNA/proteinproducer; PURExpressfunctionaltranslation/transcriptionstock known, detailedkitmolecularcomposition/occupiedproducts unassembled. NoT7/ribosomerecipeinherited.')
    mhc=g.input('pmhc','SourcebiotinylatedpMHCwithUV-labilepeptide',K.MOLECULAR)
    exchanged=g.product('exchanged','UV-exchangedpMHCpreparation;IVTTtemplates/RNA/peptide/reagentcarryoverunquantified',K.MOLECULAR)
    uv=g.input('uv_condition','SourceUVexchange60minonice,then4Catleast12h',K.CONDITION)
    g.react('uv_exchange','AddpMHCdirectlytocompletedIVTTreaction;sourceUVexchange,notDNAphotocleavage',S.GENERIC,((peptide,R.TEMPLATE),(mhc,R.STRUCTURAL),(uv,R.CONDITION)),(exchanged,))
    tetramer=g.product('tetramers','SourcematchedDNA-barcodedpMHCtetramerpreparationplusunremovedreactioncomponents',K.MOLECULAR,state=LinkerDNA(True,True).scene())
    g.react('tetramerize','AddmatchedDNA-streptavidintopMHC1:6.7,4C1h',S.GENERIC,((exchanged,R.TEMPLATE),(encoded,R.STRUCTURAL)),(tetramer,))
    own=select('panel_allocate',tetramer,'Actual.65ulown tetramerstockcontribution topanel')
    other_panel=g.input('other_panel','Actual.65uleachotherpaneltetramercontributions; membership/sourcefactorydetailsunexpanded, notcopiesofownmember',K.MOLECULAR)
    panel=g.product('panel_combined','ActualcombinedDNA-pMHCtetramerpanelpreparation',K.MOLECULAR)
    g.react('panel_pool','Actualphysicalpoolofdistincttetramerstockcontributions',S.POOL,((own,R.TEMPLATE),(other_panel,R.TEMPLATE)),(panel,))
    stainbuffer=g.input('stain_buffer','Actual400ulsourcePBS/5mMEDTA/2percentFBS/100ugmlsalmonDNA/100uMD-biotin/.05percentazidebuffer; concentratorpreblockhistoryunassembled',K.MOLECULAR)
    buffered=g.product('buffered_panel','Actualsourcepanelplusstainingbuffermaterial',K.MOLECULAR)
    g.react('panel_buffer','Source400ulstainingbufferaddition',S.GENERIC,((panel,R.TEMPLATE),(stainbuffer,R.STRUCTURAL)),(buffered,),changes_molecule=False)
    tetramer=clean('panel_concentrate',buffered,'Actual100kDaVivaspinconcentration7000g10minorlonger to~50ul; preblockprotocolnotexpanded')
    cells=g.input('cells','SourceCD8Tcells;RosetteSep/Ficollpreparation',K.BIOLOGICAL,primary=True)
    stained=g.product('stained','Tetramer-stainedcellsplusunboundreporterpanel/cellreagents',K.MOLECULAR)
    g.react('stain','Source~60nMeachtetramer+CD8Ab,4C1h',S.GENERIC,((cells,R.TEMPLATE),(tetramer,R.STRUCTURAL)),(stained,))
    x=clean('stain_wash',stained,'SourcecellwashbeforeantiPE/antiAPCmicrobeadbinding')
    x=op('microbead_bind',x,'SourceantiPE/antiAPCmicrobeadbinding thenwash','Sourcemicrobeads;exactbindingchemistryunassembled')
    beads=g.input('ls_column','SourceLSmagneticcolumnselection',K.SELECTION)
    enriched=g.product('enriched','Elutedtetramer-enrichedcellfraction',K.MOLECULAR)
    depleted=g.product('depleted','Retainedtetramer-depletedflowthroughcontrolfraction',K.MOLECULAR)
    g.react('enrich','PhysicalLScolumnfractionation;bothenrichedandflowthroughfractionsretained',S.AFFINITY_CAPTURE,((x,R.TEMPLATE),(beads,R.SELECTION)),(enriched,depleted))
    depleted=op('flow_control_stain',depleted,'Sourceflowthroughphenotypingwithsamepanel;notassignedaDNA/RNAlibrary','SourceFACSAbpanel')
    g.mark_terminal(Terminal(depleted.id,Disposition.UNRESOLVED,'Sourceflowthroughphenotypingprepared; no inventedstorage/discard/NGSendpoint.'))
    x=clean('facs_wash',enriched,'WashintoFACSbuffer;sourcephenotypicpanelstaining;sourceindividualwellmaybetetramerpositiveornegative')
    chosen=g.product('single_well','OnesourceFACS-sortedcellin4uLlysisbuffer;representativewell',K.MOLECULAR)
    other=g.product('not_chosen','Otherenrichedcells/unselectedsortfractions',K.MOLECULAR)
    g.react('single_sort','Physicalsinglecellallocation;controlidentitynotassignedbybindingprediction',S.SPLIT,((x,R.TEMPLATE),),(chosen,other));g.mark_terminal(Terminal(other.id,Disposition.UNRESOLVED,'Sourceunselectedcellsnotautomaticallydiscarded.'))
    x=op('tcr_rt',chosen,'SourceTCRRTinsamewellwhileDNA-BCreporterspresent;DNAbarcodenotcreatedbyRT','CitedTCRprotocolRTreagents/primersunreviewed')
    p1=g.input('p1','SourceP1primerforDNA-BCfirstPCR',K.OLIGO,state=Construct([Segment('P1',P1)]))
    p2=g.input('p2','SourceP2primerforDNA-BCfirstPCR',K.OLIGO,state=Construct([Segment('P2',P2)]))
    first=g.product('first_pcr','JointfirstPCRpopulation:TCRcDNAandDNA-BCampliconsplusoriginalreporter/cellreagents',K.MOLECULAR)
    joint=g.input('joint_supplies','ActualsourcecitedfirstPCRsupplies; parentrecipe/nucleotide/catalystprovenanceunassembled',K.MOLECULAR)
    g.react('joint_pcr','SourcefirstTCRPCRwithP1/P2each100nMforconcurrentDNA-BCamplification',S.MAGIC,((x,R.TEMPLATE),(p1,R.STRUCTURAL),(p2,R.STRUCTURAL),(joint,R.STRUCTURAL)),(first,),reason='ActualsourcePCRproducer; parentrecipeunreviewed, no genericchemistry-freePCRorborrowedpolymeraseidentity.')
    dna_part=g.product('dna_aliquot','Source1uLfirstPCRallocatedtoDNA-BCsecondPCR',K.MOLECULAR)
    tcr_part=g.product('tcr_aliquot','SourceparallelTCRPCRaliquot;volumeunassembled',K.MOLECULAR)
    rest=g.product('first_remaining','UnallocatedfirstPCRmaterial',K.MOLECULAR)
    g.react('split_libraries','SourceparallelwellsforDNA-BCandTCRpreparation;notmodalitypurification',S.SPLIT,((first,R.TEMPLATE),),(dna_part,tcr_part,rest));g.mark_terminal(Terminal(rest.id,Disposition.UNRESOLVED,'UntransferredfirstPCRremainderunassigned.'))
    tcr=op('tcr_later',tcr_part,'CitedTCRprotocolsubsequentPCRs;exactparentrecipeunreviewed','SourceTCRnestedPCR/indexingcomponents')
    tcr=clean('tcr_pool_gel',tcr,'Sourcepoolmultiplecells,gelelectrophoresisandgelextraction')
    end(g,tcr,'SourceIlluminaMiSeqV2;fullTCRlibrary/runlength/indexgeometryunassembled. PairedTCRchainsassociatedbysinglecellidentity,notcovalentlyfused.')
    f=g.input('well_f','SourceV1F_rxn2_1example;notassignedtoexperimentalwell',K.OLIGO,state=well_primer('f'))
    r=g.input('well_r','SourceV1R_rxn2_1example;notassignedtoexperimentalwell',K.OLIGO,state=well_primer('r'))
    second=g.product('second_pcr','DNA-BCampliconswithsourcewellbarcodesandpartialadapters;remnantsretained',K.MOLECULAR)
    pol=g.input('ex_taq','SourceExTaqHS;5uLreaction10cycles95/55/72',K.CATALYST)
    g.react('dna_pcr2','Source10cycleDNA-BCsecondPCRwithwellprimers',S.MAGIC,((dna_part,R.TEMPLATE),(f,R.STRUCTURAL),(r,R.STRUCTURAL),(pol,R.CATALYST)),(second,),reason='Actualsource10cyclePCR; fullnucleotide/sourcebuffer/occupiedampliconunassembled, notgenericno-changeoperation.')
    f3=g.input('illu_f','SourceILLU_f',K.OLIGO,state=illumina_primer('f'))
    r3=g.input('illu_r','SourceILLU_r6Xwindow;acquiredindexboundaryunassembled',K.OLIGO,state=illumina_primer('r'))
    third=g.product('third_pcr','DNA-BCIlluminaadapterPCRpopulation',K.MOLECULAR)
    g.react('dna_pcr3','Source5cyclethirdPCRsameprogramwithILLU_f/ILLU_r',S.MAGIC,((second,R.TEMPLATE),(f3,R.STRUCTURAL),(r3,R.STRUCTURAL)),(third,),reason='Actualsource5cyclePCR; additionalreactioninventory/occupiedampliconunassembled, no namednewenzymeinvented.')
    lib=clean('dna_pool_gel',third,'Poolmultiplewells;gelelectrophoresisandgelextraction')
    end(g,lib,'SourceDNA-BClibrarysequencedatleast6000reads/cell;exactrun/read/indexgeometryunassembled. Synthetic12NMIDisnotlaterNHNHNpatternorwellbarcode.')
    return incomplete(g)
def render_page():
    g=protocol_graph();return head(TITLE)+'<main class="wrap"><h1>'+TITLE+'</h1>'+info(LIMIT)+'<p><a href="01_tettcr.html">Research notes</a></p>'+''.join(workflow_panel(workflow_from_protocol_graph(g,x)) for x in ('dna_pool_gel_kept','tcr_pool_gel_kept'))+'</main>'
def selftest():
    from checks import Check
    from storage import sources_root
    c=Check();g=protocol_graph()
    c('ReporterMIDexistsbeforePCR',g.materials['linker'].state.strands['linker'].segs[1].feature.role,'umi')
    c('LinkeriSp18notcountedasDNA',len(g.materials['linker'].state.strands['linker'].segs),3)
    c('IVTTandreporterencodingusedistinctDNAaliquots',g.steps['ivtt_pcr'].inputs[0].material!=g.steps['anneal_reporter'].inputs[1].material,True)
    c('EDTAquenchdoesnotphysicallypurifyencodedDNA',g.steps['tetramerize'].inputs[1].material,'stop_extension_product')
    c('Flowthroughcontrolretained',g.terminals['flow_control_stain_product'].disposition,Disposition.UNRESOLVED)
    c('JointfirstPCRbeforelibrarysplit',g.steps['split_libraries'].inputs[0].material,'first_pcr')
    c('DNAsecondPCRusesallocated1uL',g.steps['dna_pcr2'].inputs[0].material,'dna_aliquot')
    c('Reporterextensionusespreannealednucleotides',g.steps['encode'].inputs[1].material,'encoding_retained_nt')
    c('IVTTusesselected1uldilutedDNA',g.steps['ivtt'].inputs[0].material,'ivtt_allocation_selected')
    c('IVTTpreparedstockremainderretained',g.terminals['ivtt_mix_allocation_remainder'].disposition,Disposition.UNRESOLVED)
    c('Actualpanelpoolhasdistinctsourcecontributions',len(g.steps['panel_pool'].inputs),2)
    c('ChemicalPCRnotgenericno-change',g.steps['dna_pcr2'].kind,S.MAGIC)
    c('Modelprivate',g.publishable,False)
    p=sources_root()/Path(__file__).resolve().parents[1].name/'nbt.4282_supp_NIHMS1508189-supplement-11.xlsx.txt'
    if p.exists():
        import re
        text=p.read_text();line=next(line for line in text.splitlines() if line.startswith('DNA Linker\t'))
        c('PublishedsyntheticMIDhas12Ns',len(re.findall(r'\(N(?::25252525)?\)',line)),12)
        c('PublishedamineandnonDNAspacer',('/5AmMC12/' in line,'/iSp18/' in line),(True,True))
        c('PublishedtemplateDNA',TEMPLATE in text,True)
        c('PublishedwellprimerprintedNHNHN',all(''.join(s.top for s in well_primer(k).segments).replace(' ','') in text.replace(' ','') for k in ('f','r')),True)
    else:c.skip('PublishedTable10transcription','Externaloriginalabsent')
    c.report()
