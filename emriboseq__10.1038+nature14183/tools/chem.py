"""Original IonTorrent emRiboSeq; RNaseH2 capture, no borrowed Illumina endpoint."""
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent
sys.path[:0]=[str(HERE),str(HERE.parents[1]/'lib')]
from chemdraw import Construct,Scene,Segment,oligo,revcomp
from adapter_ligation import AdapterLigationGeometry
from pcr import PCRGeometry
from protocol_graph import ProtocolGraph,InputRole as R,MaterialKind as K,StepKind as S,SequencingPlan,Terminal,Disposition
from protocol_model import molecule,role
import seqprimers as sp
PAGE_FILE='emriboseq.html';TITLE='emRiboSeq — original IonTorrent RNaseH2 end capture';NOTES='01_emriboseq.html'
SOURCE='<a href="https://doi.org/10.1038/nature14183">Reijns et al.(2015)</a>,associated accepted manuscriptPMC4374164,library preparation/design/analysisMethods.'
SUMMARY='Attach trP1, block old 3′ ends with ddATP, then nick with RNase H2. A phosphorylated adapter ligates new 3′-OH ends through a biotinylated random splint. Release the nonbiotinylated strand, copy with primer A, size-select and perform Ion Torrent emulsion PCR.'
CAVEAT='INFERRED — strand assemblies derive from printed original Ion adapters; genomic positions and run primers are unassigned. This is PGM/Proton with OneTouch emPCR,not a later Illumina EndoSeq recipe. No pre-sequencing library amplification is inserted before source second-strand synthesis and gel selection. Random6 is splint degeneracy,notUMI. The printed trP1 bottom has a blocked terminal base outside its maximal paired core; that base is retained rather than silently corrected. Full VOR and source-linked originals remain unreviewed.'
TRP1_CORE='CCTCTCTATGGGCAGTCGGTGAT';TRP1_TOP=TRP1_CORE+'T';TRP1_BOTTOM='ATCACCGACTGCCCATAGAGAGGC'
A_TOP='CTGAGTCGGAGACACGCAGGGATGAGATGG';A_BOTTOM=revcomp(A_TOP)+'N'*6;PRIMER_A='CCATCTCATCCCTGCGTGTCTCCGAC'
CAPTURED=Construct((Segment('trP1top',TRP1_TOP),role('original DNA upstream of RNaseH2nick;genomicsequenceunassigned'),Segment('covalentA-top',A_TOP)),name='INFERRED — nonbiotinylated captured strand')
FINAL_LIBRARY=Construct((Segment('A-derivednewstrandend',revcomp(A_TOP),'r1',inferred=True),role('reversecomplementofcapturedgenomicDNA'),Segment('trP1reversecomplement',revcomp(TRP1_TOP),'r2',inferred=True)),name='INFERRED — original Ion library reference')
FINAL_CAPTION='INFERRED — newly synthesized complement of captured strand; no biotin or random-splint6 is covalently assigned to it. Source physical terminal chemistry and exact Ion run primers remain separately qualified.'
SEQ_PRIMERS=();SEQUENCING_UNAVAILABLE='Original Ion Torrent PGM or Proton with OneTouch emulsion PCR. Exact flow and read-length configuration, loaded run primers and emPCR primer polymers remain unassembled.'
SEQUENCING_ENDING='The source maps the rNMP one base upstream and on the opposite strand from the read 5′ end. No library amplification preceded second-strand synthesis and selection; platform emulsion PCR remains separate. Reads were not deduplicated.'
READ_SPECS=(sp.read_spec('Read 1','unassigned','original Ion readout; source opposite-strand upstream coordinate assignment','include'),)
def state(label,*,debris=''):
 sc=Scene();sc.strand('dna',[role(label)],label='DNA target class')
 if debris:sc.strand('other',[role(debris)],label='separate material until cleanup')
 return sc

def trp1_scene():
 sc=Scene();sc.strand('top',[Segment('trP1core',TRP1_CORE),Segment('trP1T','T')],label='trP1-top;finalbondphosphorothioate')
 sc.anneal('bottom',[Segment('bottomcore',TRP1_BOTTOM[:-1]),Segment('blockedterminalC',TRP1_BOTTOM[-1])],to='top',pair=('bottomcore','trP1core'),label='trP1-bottom',mod5='5Phos',mod3='dideoxy',unpaired=('blockedterminalC',))
 sc.ligation_end('top','trP1T',kind='cohesive',polarity='3-prime',owner='adapter');return sc

def a_scene():
 sc=Scene();sc.strand('top',[Segment('Atop',A_TOP)],label='A-top',mod5='5Phos',mod3='dideoxy')
 sc.anneal('bottom',[Segment('Abottomcore',revcomp(A_TOP)),Segment('random6','N'*6)],to='top',pair=('Abottomcore','Atop'),label='A-bottom',mod5='biotin',mod3='C3')
 sc.ligation_end('top','Atop',kind='splinted',polarity='5-prime',owner='adapter');return sc

def oligos():return tuple(oligo(n,[Segment(n,x)],mods=m) for n,x,m in [('trP1top',TRP1_TOP,'finalbondphosphorothioate'),('trP1bottom',TRP1_BOTTOM,'5Phos;3-prime dideoxy'),('Atop',A_TOP,'5Phos;3-prime dideoxy'),('Abottom',A_BOTTOM,'5-primebiotin;3-primeC3'),('primerA',PRIMER_A,'')])
SOURCE_OLIGO_EXCEPTIONS={'trp1_top':'Source spells internal bond as -phosphorothioate-T; shared parser may not rejoin this English-word modification. Local source test checks full printed23base-phosphorothioate-T,retainingterminalT.'}
def protocol_graph():
 g=ProtocolGraph(TITLE,authority='Reijns2015 associated manuscript original Ionlibrary Methods; all printed sourceoligos')
 dna=g.input('dna','purified yeast genomic DNA with embedded rNMP',K.MOLECULAR,state=state('genomic DNA;positionsunassigned'),primary=True)
 def op(i,label,kind,x,condition,outstate,*,rk=K.CONDITION,rr=R.CONDITION):
  r=g.input(i+'_reagent',condition,rk);y=g.product(i+'_dna',label,K.MOLECULAR,state=outstate);g.react(i,label,kind,((x,R.TEMPLATE),(r,rr)),(y,));return y
 def clean(i,x,label,*,outstate=None,kind=K.MOLECULAR):
  r=g.input(i+'_condition',label,K.CONDITION);y=g.product(i+'_retained','recoveredDNA',kind,state=x.state if outstate is None else outstate);w=g.product(i+'_waste','purificationremovedfractions',K.WASTE)
  g.react(i,label,S.CLEANUP,((x,R.TEMPLATE),(r,R.CONDITION)),(y,w));g.discard(w,'Actualsourcepurification;noquantitativerecoveryassumption.');return y
 current=op('sonicate','Sonicate to about 400 bp',S.FRAGMENTATION,dna,'BioruptorPlus',state('sonicated genomic DNA'))
 current=clean('fragment_cleanup',current,'ethanolprecipitation;1.2xAMPure')
 current=op('endprep','NEBNext end repair and dA-tailing',S.END_REPAIR,current,'sourceNEBNextEndRepair/dA-Tailing',state('sourceA-tailedgenomicfragments'),rk=K.CATALYST,rr=R.CATALYST)
 current=clean('endprep_cleanup',current,'1.2xAMPure;subsequentbeadscaptureNaCl/PEG')
 top=g.input('trp1_top','trP1top;terminalbondphosphorothioate',K.OLIGO,state=Construct([Segment('trP1top',TRP1_TOP)]))
 bottom=g.input('trp1_bottom','trP1bottom;5Phos,3dideoxy',K.OLIGO,state=Construct([Segment('trP1bottom',TRP1_BOTTOM)]))
 trp1=g.product('trp1','sourceannealedtrP1adapterpair',K.ADAPTER,state=trp1_scene())
 g.react('anneal_trp1','AnnealprintedtrP1pair;95C5min,gradualcooling',S.HYBRIDIZATION,((top,R.TEMPLATE),(bottom,R.STRUCTURAL)),(trp1,))
 ligase=g.input('ligation_mix','NEBNextQuickLigation,14–18h16C',K.CATALYST)
 joined=g.product('trp1_ligated','trP1-ligatedDNA',K.MOLECULAR,state=state('trP1-bearing genomic DNA;endassemblyunassigned'))
 g.react('trp1_ligate','LigatetrP1beforeendblocking',S.LIGATION,((current,R.TEMPLATE),(trp1,R.STRUCTURAL),(ligase,R.CATALYST)),(joined,),geometry=AdapterLigationGeometry.known('cohesive',trp1_scene(),statement='Source trP1top3-T;blockedbottomterminalCretainedasprinted;genomicsourceA-tailed.'))
 current=op('block','Block old free 3′ ends',S.TAILING,joined,'TdT/ddATP,20Uperµg,37C2h',state('oldfree3-primeendsdd-blocked;notnewRNaseH2nicks',debris='excessddATPandreactionproducts'),rk=K.CATALYST,rr=R.CATALYST)
 current=clean('block_cleanup',current,'AMPurepurification;removebeads',outstate=state('old3-primeendclassdd-blocked'))
 current=op('nick','RNase H2 nicking',S.DIGESTION,current,'RNaseH2,10pmolperµg,37C2h;80C20min',state('new3-OHupstreamofrNMP;newdownstream5-P;oldendsblocked'),rk=K.CATALYST,rr=R.CATALYST)
 current=clean('nick_cleanup',current,'1.8xAMPure')
 current=op('sap','Remove 5′ phosphates',S.DEPHOSPHORYLATION,current,'SAP5U,37C1hperµg;65C15min',state('new3-OHretained;5-phosphatesremoved',debris='phosphataseproducts'),rk=K.CATALYST,rr=R.CATALYST)
 current=clean('sap_cleanup',current,'AMPurepurification',outstate=state('new3-OHvsoldblockedends;5-OHclass'))
 current=op('denature','Denature and snap-cool',S.THERMAL_TREATMENT,current,'95C5min;snapcool',current.state)
 adapter=g.input('a_adapter','sourceAadapter/biotinsplint',K.ADAPTER,state=a_scene())
 complex_state=Scene();complex_state.strand('captured',list(CAPTURED),label='INFERRED — covalentnonbiotinylatedstrand',mod3='dideoxy');complex_state.strand('splint',[Segment('biotinsplint',A_BOTTOM)],label='separate noncovalent biotin splint',mod5='biotin',mod3='C3')
 captured=g.product('captured','INFERRED — A-top-ligatedDNAwithbiotinylatedsplint',K.MOLECULAR,state=complex_state,inferred=True)
 g.react('a_ligate','INFERRED — A-top5-PjoinsnewDNA3-OH',S.LIGATION,((current,R.TEMPLATE),(adapter,R.STRUCTURAL),(ligase,R.CATALYST)),(captured,),geometry=AdapterLigationGeometry.known('splinted',a_scene(),statement='Source6randomAbottomsplintsnewDNA3-OHtoAtop5-P;Atop3-ddandAbottom3-C3blocked.'))
 beads=g.input('streptavidin','M280streptavidinbeads;sourcevendorguidelines',K.SELECTION)
 bound=g.product('bound','biotin-splint-associatedDNAcapturedonbeads',K.MOLECULAR,state=complex_state)
 flow=g.product('flowthrough','nonretainedcapturematerial',K.MOLECULAR,state=state('unretained material;exactcompositionunassigned'))
 g.react('capture_beads','Capture biotinylated adapter complexes',S.AFFINITY_CAPTURE,((captured,R.TEMPLATE),(beads,R.SELECTION)),(bound,flow));g.mark_terminal(Terminal(flow.id,Disposition.UNRESOLVED,'Nonretainedfractionnotsequenced;detailedphysicaldispositionnotreconstructed.'))
 naoh=g.input('naoh','0.15MNaOHstrandrelease',K.CONDITION)
 released=g.product('released','INFERRED — releasednonbiotinylatedcapturedDNA',K.MOLECULAR,state=CAPTURED,inferred=True)
 splint_state=Scene();splint_state.strand('splint',[Segment('Abottom',A_BOTTOM)],label='biotinylatedsplintretainedwithbeads',mod5='biotin',mod3='C3');splint_state.note('splint','Biotin-streptavidinassociation;beadattachmentnotacovalentDNAnucleotide.')
 splint=g.product('bead_splint','biotinylatedsplintretainedwithbeads',K.MOLECULAR,state=splint_state)
 g.react('release','Release nonbiotinylated strand; retain biotin splint',S.STRAND_SEPARATION,((bound,R.TEMPLATE),(naoh,R.CONDITION)),(released,splint));g.mark_terminal(Terminal(splint.id,Disposition.UNRESOLVED,'SourceNaOHreleasesnonbiotinylatedstrand;retainedbiotin-beadfractionnotsequenced,finaldisposalunassigned.'))
 current=clean('release_cleanup',released,'ethanolprecipitationofreleasedssDNA')
 primer=g.input('primer_a','sourceprimerA',K.OLIGO,state=Construct([Segment('primerA',PRIMER_A)]));pol=g.input('phusion','PhusionFLASHforsecondstrandsynthesis',K.CATALYST)
 copied=g.product('copied','INFERRED — second-strandlibraryreference',K.MOLECULAR,state=FINAL_LIBRARY,inferred=True)
 g.react('second_strand','INFERRED — primerAsecondstrandsynthesis,nopre-sequencingPCRseries',S.EXTENSION,((current,R.TEMPLATE),(primer,R.STRUCTURAL),(pol,R.CATALYST)),(copied,))
 current=clean('gel',copied,'2%E-GelEX200–300bpselection;BioanalyzerQC')
 f=g.input('ion_bead_primer','OneTouchemPCRbeadprimer;exactpolymerunassigned',K.OLIGO,state=molecule('beadprimer',(role('exactIonbeadprimerunassigned'),)));r=g.input('ion_other_primer','OneTouchemPCRsecondprimer;exactpolymerunassigned',K.OLIGO,state=molecule('emPCRprimer',(role('exactsecondprimerunassigned'),)));mix=g.input('ion_mix','sourceIonTorrentOneTouchemPCRmixture',K.CATALYST)
 library=g.product('ion_library','sourceIonemPCRsequencinglibrary',K.LIBRARY,state=FINAL_LIBRARY,inferred=True)
 g.react('emulsion_pcr','Ion Torrent OneTouch emulsion PCR after library QC',S.PCR,((current,R.TEMPLATE),(f,R.STRUCTURAL),(r,R.STRUCTURAL),(mix,R.CATALYST)),(library,),geometry=PCRGeometry.unavailable(FINAL_LIBRARY,reason='SourceOneTouchplatformspecified;emPCRoligos/conditionsandbead-strandgeometryunassembled.'))
 g.sequence(library,SequencingPlan(library=FINAL_LIBRARY,limitation=SEQUENCING_UNAVAILABLE,acquired_roles=('Read 1',)))
 g.declare_reviewed_scope(('dna','sonicate_reagent','sonicate_dna','sonicate','fragment_cleanup_condition','fragment_cleanup_retained','fragment_cleanup_waste','fragment_cleanup','endprep_reagent','endprep_dna','endprep','endprep_cleanup_condition','endprep_cleanup_retained','endprep_cleanup_waste','endprep_cleanup','trp1_top','trp1_bottom','trp1','anneal_trp1','ligation_mix','trp1_ligated','trp1_ligate','block_reagent','block_dna','block','block_cleanup_condition','block_cleanup_retained','block_cleanup_waste','block_cleanup','nick_reagent','nick_dna','nick','nick_cleanup_condition','nick_cleanup_retained','nick_cleanup_waste','nick_cleanup','sap_reagent','sap_dna','sap','sap_cleanup_condition','sap_cleanup_retained','sap_cleanup_waste','sap_cleanup','denature_reagent','denature_dna','denature','a_adapter','captured','a_ligate','streptavidin','bound','flowthrough','capture_beads','naoh','released','bead_splint','release','release_cleanup_condition','release_cleanup_retained','release_cleanup_waste','release_cleanup','primer_a','phusion','copied','second_strand','gel_condition','gel_retained','gel_waste','gel','ion_bead_primer','ion_other_primer','ion_mix','ion_library','emulsion_pcr'),g.authority)
 return g.require_publishable()
def render_page():
 from batch_page import render
 return render(sys.modules[__name__])
def selftest():
 import re
 from checks import Check,Source,have
 from storage import sources_root
 c=Check();source=Source(sources_root()/HERE.parent.name/'nature14183_PMC4374164.xml.txt','originalIonOligonucleotidesandadapterdesign')
 if have(c,source,label='originalIonoligos'):
  t=re.sub(r'\s+','',source.path.read_text())
  for n,x in [('trP1top',TRP1_CORE+'-phosphorothioate-T'),('trP1bottom',TRP1_BOTTOM+'-dideoxy'),('Atop',A_TOP+'-dideoxy'),('Abottom',A_BOTTOM+'-C3'),('primerA',PRIMER_A)]:c('source '+n,x in t,True)
 protocol_graph();render_page();c.report()
