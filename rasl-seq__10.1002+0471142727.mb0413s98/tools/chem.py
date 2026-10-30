"""RASL-seq associated-author-manuscript basic purified-RNA route, private partial graph."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Scene,Segment,feature
from illumina import P5,P7
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,SourceItem,SequencingPlan,Terminal,Disposition
from adapter_ligation import AdapterLigationGeometry
from pcr import PCRGeometry
TITLE='RASL-seq: selected purified-RNA route'
AUTHORITY='AssociatedauthorMS PMC3325489.1/NIHMS368182, sourceplanning45–47/materials63–85/basicprotocol95–197. Oneoriginal; source-qualified historicalprotocol, notfinalpublished-body equivalence orabsenceaudit.'
INDEX_ARM='TAGCTGATGCTACGACCACAGG'

def seq(name):return {'p5':P5[:23],'barcode':P7[:-2]+'N'*7+INDEX_ARM,'oligodt':'T'*25,'index_read':INDEX_ARM}[name]

def inventory(name):
    sc=Scene()
    parts=[Segment(name,seq(name))]if name!='barcode'else[Segment('source_P7_prefix',P7[:-2]),Segment('sample_index','N'*7,placeholder=True,feature=feature('rasl_sample','sample_index','unknown'),note='Source7ntfixedsamplecode;actualmemberunassigned,noUMI'),Segment('source_index_arm',INDEX_ARM)]
    sc.strand(name,parts,label='SourcebiotinylatedoligodT25;biotinpositionunassigned'if name=='oligodt'else'Sourcehistoricalprimer');return sc

def unknown():
    sc=Scene();sc.strand('unknown',[Segment('unassembled','[unassembled]',placeholder=True)]);return sc

def probe_end_inventory():
    sc=Scene();sc.strand('donor',[Segment('source_3side_target','N'*20,placeholder=True,note='Source20nttargetarm;actualbasesunassigned'),Segment('source_index_arm','[index primer arm unassembled]',placeholder=True)],mod5='P')
    sc.strand('acceptor',[Segment('source_P5_arm','[P5 primer arm unassembled]',placeholder=True),Segment('source_5side_target','N'*20,placeholder=True,note='Source20nttargetarm;actualbasesunassigned')])
    sc.strand('rna_template',[Segment('source_RNA','[target RNA unassembled]',placeholder=True)],label='Independent RNA template inventory; no pairing assigned')
    return sc

def protocol_graph():
    g=ProtocolGraph(TITLE,authority=AUTHORITY,representation='published_instance');x=g.input('rna','Sourceapproximately1ugisolatedtotalRNAin20ulwater;chosenbasicroute,cellMELTalternativenotmerged',K.MOLECULAR,primary=True)
    def term(m,why,disp=Disposition.STORED):g.mark_terminal(Terminal(m.id,disp,why))
    def wash(key,label):
        nonlocal x
        b=g.input(key+'_condition',label,K.CONDITION);y=g.product(key+'_kept','Sourcephysicallyretainedbead-associatedpreparation;recovery/carryoverunassigned',K.MOLECULAR);w=g.product(key+'_waste','Actualaspiratedwash/supernatantmaterial',K.WASTE);g.react(key,label,S.CLEANUP,((x,R.TEMPLATE),(b,R.CONDITION)),(y,w));g.discard(w,'Actualaspirationnotperfectmolecularselectivity');x=y
    cocktail=g.input('binding_cocktail','Source20ul2xbindingcocktailwithprobes5nMeach/biotin(dT)25 100nM/5ulwashedSTVbeads/TrisNaClEDTASDS;literalstocknotassumedfinalconcentrations',K.MOLECULAR)
    probes=g.product('binding_probes','DNAprobe-poolcomponentofonebindingcocktail;two20nttargetarmsperspecificpair,wholeorderingsunassigned;3sideprobe5Pperplanning45',K.OLIGO,chemical_state={'ligation_donor':'5primephosphateonsource3sidetargetprobe;noorientationinferredfromambiguousupstreamterminology'})
    dt=g.product('binding_dt','BiotinylatedoligodT25componentofsamebindingcocktail;biotinattachmentpositionunassigned',K.OLIGO,state=inventory('oligodt'),chemical_state={'biotin_position':'unassigned'})
    beads=g.product('binding_beads','WashedSTVbeadcomponentofsamecocktail;notnewsecondaddition',K.SELECTION);b=g.product('binding_buffer','Buffercomponentofsamecocktail',K.CONDITION)
    g.react('binding_components','Representonecocktailcomponentinventory;notphysicalfractionation',S.GENERIC,((cocktail,R.TEMPLATE),),(probes,dt,beads,b),changes_molecule=False)
    y=g.product('hybridized','SourceRNA/probepairs/oligodT/beadmixture;targetstoichiometryandalloff-targetclassesunassembled',K.MOLECULAR);g.react('hybridize','Source65C10minthen45C60minannealing',S.HYBRIDIZATION,((x,R.TEMPLATE),(probes,R.STRUCTURAL),(dt,R.STRUCTURAL),(beads,R.SELECTION),(b,R.CONDITION)),(y,));x=y
    wash('capture_washes','Sourcemagneticstandaspirationsleaveapproximately10ul;5x50ulwashingbufferthen1x50ulligationbuffer')
    template=g.product('retained_rna','Same-tubeRNA/oligodT/beadtemplatecomponent;notphysicalseparation',K.MOLECULAR);ret=g.product('retained_probes','Same-tubeco-retainedDNAprobecomponent;notnewfreeprobeaddition',K.OLIGO)
    g.react('retained_components','Representretainedligationtemplatesandprobes',S.GENERIC,((x,R.TEMPLATE),),(template,ret),changes_molecule=False)
    lig=g.input('ligase','Sourceadd20ul1xligationbuffercontaining5WeissUT4DNAligase;residualwashvolumeunassigned',K.CATALYST);y=g.product('ligated','SourceRNA-templatedligatedDNAprobeproductsplusRNA/beads/unligated/randomligationcompanions',K.MOLECULAR)
    g.react('ligate','SourceT4DNAligation37C60min',S.LIGATION,((template,R.TEMPLATE),(ret,R.STRUCTURAL),(lig,R.CATALYST)),(y,),geometry=AdapterLigationGeometry('unavailable','unavailable',(probe_end_inventory(),),'Sourceadjacent20ntDNAprobesonRNAanddonor5Pknown;actualtarget/wholeprobeorderingsandfullRNA-splintedendstateunassembled.Notasourceabsenceclaim.'));x=y
    wash('ligation_washes','Source3washingbufferwashes')
    water=g.input('resuspend_water','Source30ulwaterresuspension',K.CONDITION);hot=g.product('released','Source65C5minreleaseofDNAprobesfromRNA;completeeluatecompositionunassigned',K.MOLECULAR);g.react('release','Sourcewaterresuspension/65C5min',S.THERMAL_TREATMENT,((x,R.TEMPLATE),(water,R.CONDITION)),(hot,))
    selected=g.product('five_ul','Source5ulsupernatanttransferredaftermagneticstand',K.MOLECULAR);backup=g.product('backup_plate','Originalplatewithremainingmaterialincludingbeads;approximately25uluntransferredliquidif30ulretained',K.MOLECULAR);g.react('aliquot','Source5ulreleased-probeallocation',S.SPLIT,((hot,R.TEMPLATE),),(selected,backup));term(backup,'Sourceoriginalplatestoredminus20Casbackup');x=selected
    pcrmix=g.input('pcr_mix','Source10ulPCRcocktailAmpliTaqGold.45U/buffer/MgCl2/dNTP/P5.6uM/barcodedprimer.6uM;actualstock/finalconcentrationsnotreinterpreted',K.MOLECULAR)
    p5=g.product('pcr_p5','P5componentofsamePCRcocktail',K.OLIGO,state=inventory('p5'));idx=g.product('pcr_barcode','BarcodeprimercomponentofsamePCRcocktail',K.OLIGO,state=inventory('barcode'));e=g.product('pcr_enzyme','AmpliTaqGoldcomponentofsamePCRcocktail',K.CATALYST);b=g.product('pcr_buffer','Buffer/dNTPcomponentofsamePCRcocktail',K.CONDITION);g.react('pcr_components','Representoneactual10ulPCRcocktailaddition',S.GENERIC,((pcrmix,R.TEMPLATE),),(p5,idx,e,b),changes_molecule=False)
    y=g.product('pcr_product','SourceindexedDNAampliconsplusprimer/backgroundclasses;completeproductgeometryunassembled',K.MOLECULAR);g.react('pcr','Source94C10minthen20to25cycles94C30s58C30s72C30s',S.PCR,((x,R.TEMPLATE),(p5,R.STRUCTURAL),(idx,R.STRUCTURAL),(e,R.CATALYST),(b,R.CONDITION)),(y,),geometry=PCRGeometry.unavailable(unknown(),reason='SourcePCRprimerorderingsknownbuttarget-specificDNAprobes/fulllandingandsealedproductunassembled.'));x=y
    gel=g.product('gel_sample','Sourceanalyticalgelallocationforselectedcheckedwell;volumeunassigned',K.QC);production=g.product('production','SourcePCRproductionmaterialafteranalyticalallocation',K.MOLECULAR);g.react('gel','Sourcegelcheckforselectedsubsetofwells,notallwells',S.SPLIT,((x,R.TEMPLATE),),(production,gel));term(gel,'Analyticalmaterialnotreturned',Disposition.QC);x=production
    other=g.input('other_wells','ActualotherPCRwellproductsfromsource1to4plates;samplecomposition/amountunassembled',K.MOLECULAR);pooled=g.product('pooled','ActualsourcepooledPCRproducts;selectedwellplusotherwells',K.MOLECULAR);g.react('pool','SourcePCRsamplepool',S.POOL,((x,R.TEMPLATE),(other,R.TEMPLATE)),(pooled,));x=pooled;wash('column','SourceNucleoSpinExtractIIpurificationawayfromfreeprimers')
    q=g.product('qubit','SourceQubitanalyticalallocation;volumeunassigned',K.QC);prod=g.product('postqc','SourceproductionpoolafterQubit',K.MOLECULAR);g.react('qubit_split','SourceQubitallocation',S.SPLIT,((x,R.TEMPLATE),),(prod,q));term(q,'Actualanalyticalallocation',Disposition.QC)
    lib=g.product('library','SourcehistoricalRASL-seqpool;completeprimerlanding/flowcellgeometryunassembled',K.LIBRARY);g.react('identity','Physicallibraryidentity',S.GENERIC,((prod,R.TEMPLATE),),(lib,),changes_molecule=False);g.sequence(lib,SequencingPlan(limitation='SourceHiSeq2000/Solexa:40ntligated-proberegionreadthenstrip-and-reprimeindexarmtoread7ntsamplecode;notmodernpairedenddefault. HistoricalP7prefix22ntandP5prefix23ntpreserved. NoUMIorRTinthischosenroute.'))
    g.add_source_item(SourceItem('selected',AUTHORITY,'Actualtarget-specificprobeDNAandRNA-splintedligationgeometryunassembled;noabsenceclaim.'));return g

def selftest():
    from checks import Check
    c=Check();g=protocol_graph()
    c('WholehistoricalP5',seq('p5'),'AATGATACGGCGACCACCGAGAT')
    c('Wholehistoricalbarcodeprimer',seq('barcode'),'CAAGCAGAAGACGGCATACGAGNNNNNNNTAGCTGATGCTACGACCACAGG')
    c('Wholeindexarm',seq('index_read'),'TAGCTGATGCTACGACCACAGG')
    c('Biotinpositionnotinvented',g.materials['binding_dt'].chemical_state['biotin_position'],'unassigned')
    c('Ligateusesretainednotnewprobe',[(i.material,i.role)for i in g.steps['ligate'].inputs if i.role==R.STRUCTURAL],[('retained_probes',R.STRUCTURAL)])
    c('KnownDNAprobe5Pdrawn',probe_end_inventory().strands['donor'].mod5,'P')
    c('Originalplatestockkept',g.terminals['backup_plate'].disposition,Disposition.STORED)
    c('SamebulkprimersenterPCR',[(i.material,i.role)for i in g.steps['pcr'].inputs if i.role==R.STRUCTURAL],[('pcr_p5',R.STRUCTURAL),('pcr_barcode',R.STRUCTURAL)])
    c('NoRToperationinvented',any(s.kind==S.REVERSE_TRANSCRIPTION for s in g.steps.values()),False)
    c('PoolincludesotherphysicalDNA',len([i for i in g.steps['pool'].inputs if i.role==R.TEMPLATE]),2)
    c('Qubitmaterialnotreturned',g.terminals['qubit'].disposition,Disposition.QC)
    from checks import Source,have
    from storage import sources_root
    source=Source(sources_root()/'rasl-seq__10.1002+0471142727.mb0413s98'/'0471142727.mb0413s98_PMC3325489.xml.txt','https://pmc-oa-opendata.s3.amazonaws.com/PMC3325489.1/PMC3325489.1.xml via tools/get_sources.py; materials85')
    if have(c,source):
        import re
        normalized=re.sub('[^ACGTNX]','',source.path.read_text().upper()).replace('X','N')
        for name in ('p5','barcode','index_read'):c('Whole'+name+'againstselectedoriginal',seq(name)in normalized,True)
    c('Privatepartialgraph',bool(g.publication_problems()),True);c.report()
if __name__=='__main__':selftest()
