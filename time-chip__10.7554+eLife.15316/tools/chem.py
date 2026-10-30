"""Published time-ChIP0hr arm; local partial preparation graph, not a full library model."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Scene,Segment
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,Terminal,Disposition,SequencingPlan,SourceItem

CONTROLS={'GC30': 'CAATTTGTGTTAGAATCGTTTGATTATACATAGTTTTCTGCTAATGACTCAACATATGCTATCTTTCTATATTCAGAGTCACCTTCATAATTGTCTATCATAGTTCGTATTATACAGTATATTAGTCGTGACTATTAATGTGGAATACCA', 'GC40': 'CCTATTCTCCTACGATACTATATATGATCCGCGTTATGGAACGTCCATATAAGTAGAAGTAAGTAGTTGATTGGTACTACGTAGAGCCAAGTGTCCGTTATACCCGAACGTGCTACACTTGCACATCAGGTAATTATAGAACGTCATTAA', 'GC50': 'CCTCATAGATACTTGTCCTTGCGAGGAAGTTTCGTTCGCTTGGCTTCTGCAAGAAGGCCGAAGCGCGCCTTAGCTAAGTTGTATCGGGGTCTTATACCTTGTAGCGACGGATAATGTGCGCACTAAGCAATAACGTTGCAACGACCCTAG', 'GC60': 'CGGCTACCGCTGCCATGCCGTATTATCGATAGGAACACGGGTCAGGGAGGACCTGCTCCACTGCCACGTTTGCAGCGACATCTCGGTCGTCGTCGCATAGGGTTGCTTCGGTTAGCGCCTACTGTGCGCTATTCCAAGTTCGCCGCCGAG', 'GC70': 'CGGCGGAATTGCCAGCACCGCGTGGGCTCTGGGAGCTTCACCCGCCCAACACGGTCACGCCGCCGGATACTGCGTGGATTGGCTCGGCCAGCGCTCCACCGGCTCCCTGCTGCCCGCGGTTGGGTTTTGTGGCTGCGCGGTGAGGCCAGA'}
SOURCE='PublishedDeaton2016 eLife15316v2 mainMethods349–355/385, originalPDFp18 units; publishedS5wholeworkbook. Only source0hr H3.3ESC arm here; othertimes/variantsnotanotherprotocol.'

def control_state(name):
 sc=Scene();sc.strand(name+' listed strand; partner and end chemistry unassigned',[Segment(name,CONTROLS[name])]);return sc

def cleanup(g,key,x,label):
 c=g.input(key+'_buffer',label,K.CONDITION);y=g.product(key+'_kept','Source retained preparation; yield/carryover unassigned',K.MOLECULAR);w=g.product(key+'_waste','Physical removed fraction; complete molecular inventory unassembled',K.WASTE)
 g.react(key,label,S.CLEANUP,((x,R.TEMPLATE),(c,R.CONDITION)),(y,w));g.discard(w,'Physical source wash/separation, not computational filtering.');return y

def treatment(g,key,x,label,kind=S.MAGIC,reagent='Source condition',role=R.CONDITION):
 c=g.input(key+'_stock',reagent,K.CATALYST if role==R.CATALYST else K.MOLECULAR if role==R.STRUCTURAL else K.CONDITION);y=g.product(key+'_product',label+'; exact states and residual components unassembled',K.MOLECULAR)
 kw={'reason':'Source operation known, complete source molecular states and transfer geometry not assembled.'}if kind==S.MAGIC else {}
 inputs=[(x,R.TEMPLATE),(c,role)]
 if role==R.STRUCTURAL:
  medium=g.input(key+'_medium','Sourcecomplete medium andincubationconditions; distinctsolventcomponent ofCP-biotinpulse',K.CONDITION);inputs.append((medium,R.CONDITION))
 g.react(key,label,kind,inputs,(y,),**kw);return y

def protocol_graph():
 g=ProtocolGraph('time-ChIP — source H3.3ESC0hr arm',authority=SOURCE,representation='published_instance')
 x=g.input('induced_cells','Source independently cultured H3.3_SNAP_HA ESC cells;500ng/ml doxycycline48h induction upstream. Genetic integration/selection/culture not molecularlyassembled.',K.BIOLOGICAL,primary=True)
 x=treatment(g,'pulse',x,'CP-biotin pulse10uM40min37C in50ulcomplete medium; source0.5–2x10^8cells/50ul literal scale retained, no silent unit correction',S.CHEMICAL_CONVERSION,'Source CP-biotin labeling reagent in medium; tag covalently contributes to protein, complete protein-DNA complex states unassembled',R.STRUCTURAL)
 for n in (1,2):x=cleanup(g,'pulse_wash'+str(n),x,'Source10mlmedium wash'+str(n))
 x=treatment(g,'free_label_release',x,'Resuspend10mlcomplete medium20min37C toremoveunboundCP-biotin',S.THERMAL_TREATMENT,'Sourcecomplete medium and37C incubation')
 x=cleanup(g,'final_medium_wash',x,'Source final mediumwash beforetimeallocation')
 zero=g.product('zero_hour','Source one-fifth labeled cells at0hr',K.BIOLOGICAL);chase=g.product('chase_cohort','Remainingfour-fifthslabeledcells platedfor3/6/12hrharvests; allocation/proliferation notassembled',K.BIOLOGICAL)
 g.react('time_allocation','Source one-fifth0hr harvest; remainingcellsplated',S.SPLIT,((x,R.TEMPLATE),),(zero,chase));g.mark_terminal(Terminal(chase.id,Disposition.UNRESOLVED,'Other three source temporal arms and their separate libraries outside selected0hrgraph; not discarded.'))
 x=treatment(g,'nuclei',zero,'Extract/count nuclei andselectequalcellnumberper timepoint; exacttransfer/recoveryunassembled',S.EXTRACTION,'Source500ul20mMHepespH7.7/20mMKCl/0.5mMEDTA/300mMNaCl/1mMDTT/proteaseinhibitors after nuclei extraction')
 x=treatment(g,'mnase',x,'Source120–140UMNase,3mMCaCl2,20min25C; nativechromatin, notcrosslinktitrationbranch',S.DIGESTION,'SourceMNaseWorthington',R.CATALYST)
 x=treatment(g,'quench',x,'StopMNasewithEDTA/EGTA',S.QUENCHING,'SourceEDTA/EGTA; fullfinalconcentrationsunassigned')
 detergent=g.input('detergent','Source0.05%NP40',K.CONDITION)
 soluble=g.product('soluble','Source12000g5min supernatant',K.MOLECULAR);pellet=g.product('pellet','Physical centrifugation pellet',K.WASTE)
 g.react('solubilize_spin','Source0.05%NP40 solubilization followedby12000g5min',S.SOLUBILIZATION,((x,R.TEMPLATE),(detergent,R.CONDITION)),(soluble,pellet));g.discard(pellet,'Source downstream uses resulting supernatant; completepelletinventory unassembled.')
 used=g.product('capture_aliquot','Source400ul supernatant to capture',K.MOLECULAR);rest=g.product('supernatant_rest','Untransferredsupernatant; recoveredvolumeunassigned',K.MOLECULAR)
 g.react('capture_allocation','Source400ulactualsupernatanttransfer',S.SPLIT,((soluble,R.TEMPLATE),),(used,rest));g.mark_terminal(Terminal(rest.id,Disposition.UNRESOLVED,'No perfect500ulrecovery/100ulrestassumption.'))
 beads=g.input('beads','Source15ulstreptavidinmagneticbeadsprebocked50mgBSA/200ugmlyeasttRNA; literal50mgretained',K.SELECTION)
 bound=g.product('bound','Sourcebead-associatednativehistone/DNAcomplexes; capture/carryover notperfect',K.MOLECULAR);unbound=g.product('unbound','Sourcecaptureflowthrough',K.WASTE)
 g.react('capture','Source1h4C streptavidincapture',S.AFFINITY_CAPTURE,((used,R.TEMPLATE),(beads,R.SELECTION)),(bound,unbound));g.discard(unbound,'Physicalsupernatant removal forsource beadwashes; notassumedallunlabeledDNA.')
 x=bound
 for n in (1,2,3):x=cleanup(g,'bead_wash'+str(n),x,'Source5minwash:20mMHepespH7.7/20mMKCl/500 M NaCl LITERALPUBLISHEDVALUE/3mMCaCl2/5mMEGTA/0.5%NP40/1mMDTT/proteaseinhibitors; NaCl formulationunresolved, notnormalizedto500mM')
 rnase=g.input('rnase','SourceDNase-freeRNaseA',K.CATALYST);buffer=g.input('elution_buffer','Source250ul10mMTrispH8/100mMNaCl/10mMEDTA/0.5%SDS',K.CONDITION)
 eluted=g.product('eluate','SourceDNA-bearingeluate includingreleasedmaterial; completeprotein/RNA/beadfates notassembled',K.MOLECULAR);beadrest=g.product('spent_beads','Physicaluntransferredbead-associatedfraction; capture/releaseyieldunknown',K.WASTE)
 g.react('elute','SourceDNAelution10min30C withRNaseA; noproteinaseKorspecificdeproteolysisinvented',S.MAGIC,((x,R.TEMPLATE),(rnase,R.CATALYST),(buffer,R.CONDITION)),(eluted,beadrest),reason='Sourceelution/RNaseAknown;wholecomplexreleaseandfractionstates unassembled.');g.discard(beadrest,'Separatephysicalspentcapturematerial, noquantitative releaseassertion.')
 x=cleanup(g,'phenol_chloroform',eluted,'Sourcephenol/chloroformextraction')
 x=cleanup(g,'chloroform',x,'Sourcechloroformextraction')
 x=cleanup(g,'column',x,'SourceQiagenPCRpurificationcolumn')
 retained=g.product('quantified','RetainedpurifiedDNA',K.MOLECULAR);qc=g.product('picogreen_sample','SourcePicoGreenDNAquantificationmaterial; actualamountunassigned',K.QC)
 g.react('picogreen_allocation','SourcePicoGreenquantification; materialallocationnotassembled',S.SPLIT,((x,R.TEMPLATE),),(retained,qc));g.mark_terminal(Terminal(qc.id,Disposition.QC,'Separateassaymaterial; nosilentreturntoreaction.'))
 stocks=[g.input('control_'+name.lower(),'SourceS5syntheticDNAcontrol'+name+'; listedstrand, partnerandendsunassigned',K.MOLECULAR,state=control_state(name))for name in CONTROLS]
 pool=g.product('control_pool','Source5controls equimolarstock',K.MOLECULAR);g.react('controls_pool','SourceS5equimolarcontrolmixture',S.POOL,tuple((v,R.TEMPLATE)for v in stocks),(pool,))
 spike=g.product('spike_aliquot','SourcecontrolDNAamount=0hrrecoveredDNA/10000, sameabsoluteamountforeachtimepoint',K.MOLECULAR);spikerest=g.product('spike_remainder','Untransferredcontrolstock',K.MOLECULAR)
 g.react('spike_allocation','SourcefixedcontrolDNAallocationperlibrary',S.SPLIT,((pool,R.TEMPLATE),),(spike,spikerest));g.mark_terminal(Terminal(spikerest.id,Disposition.UNRESOLVED,'Unusedstockamount/fateunassigned.'))
 library=g.product('library','Sourcebarcodedtime-ChIPlibrary; exactBowman2013adapter/ends/stepsunassembled',K.LIBRARY)
 g.react('library_preparation','Bowman2013delegatedDNA-librarypreparationwithphysicalDNAspike; stageofspikeadditionunassembled',S.MAGIC,((retained,R.TEMPLATE),(spike,R.SPIKE_IN)),(library,),reason='ActualsourceexplicitlydelegateslibrarypreparationtoBowman2013; directlylinkedpredecessorunreadwithinthispass. Cannotinvent oligos/ends/cyclesfromsharedIlluminaname.')
 g.sequence(library,SequencingPlan(limitation='Publishedpaired-end50HiSeq2500,4barcodedsamples/lane; exactindexcycles/primers/wholelibrarygeometryunassembled. NoTruSeqRNAkitborrowedfromseparateRNA-seqbranch.'))
 g.add_source_item(SourceItem('pending','WholeBowman2013librarychemistry/source5controlsphysicalends/other3timepointlibraryarms','Selectedsourcepreparativefront-end only; notfullsourceorabsenceaudit.'))
 return g

def selftest():
 from checks import Check,Source,have
 from storage import sources_root
 import zipfile,xml.etree.ElementTree as ET
 c=Check();g=protocol_graph()
 c('KnownNaClsourceambiguitynotpromotedtofinishedgraph',bool(g.publication_problems()),True)
 c('Actual400ultransferhasseparatestockfate','supernatant_rest'in g.terminals,True)
 c('Otherphysicaltimecohortisnotdiscarded',g.terminals['chase_cohort'].disposition,Disposition.UNRESOLVED)
 src=Source(sources_root()/'time-chip__10.7554+eLife.15316'/'eLife.15316_elife-15316-supp5-v2.xlsx','OriginalS5Sheet1controlrows')
 if have(c,src,label='PublishedS5original'):
  ns={'s':'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
  with zipfile.ZipFile(src.path)as z:
   strings=[''.join(el.itertext())for el in ET.fromstring(z.read('xl/sharedStrings.xml')).findall('s:si',ns)];rows=[]
   for row in ET.fromstring(z.read('xl/worksheets/sheet1.xml')).findall('s:sheetData/s:row',ns):
    vals=[]
    for cell in row.findall('s:c',ns):
     v=cell.find('s:v',ns);vals.append(strings[int(v.text)]if cell.get('t')=='s'and v is not None else v.text if v is not None else '')
    rows.append(vals)
  actual={r[0]:r[2]for r in rows if r and r[0].startswith('GC')}
  c('Allfivewholecontrolsagreewithindependentoriginalworkbook',CONTROLS,actual)
 src=Source(sources_root()/'time-chip__10.7554+eLife.15316'/'eLife.15316_v2.xml.txt','PublishedMethods353/385')
 if have(c,src,label='Publishedmain'):
  t=''.join(src.path.read_text().split());c('Literal500MNaClnotextractionrepair','500MNaCl'in t,True);c('SourcePE50readout','Paired-end50'in t,True)
 c.report()
if __name__=='__main__':selftest()
