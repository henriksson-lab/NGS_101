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
 sc.anneal('bottom',[Segment('Abottomcore',revcomp(A_TOP)),Segment('random6','N'*6,placeholder=True,note='Source random splint degeneracy; not a UMI')],to='top',pair=('Abottomcore','Atop'),label='A-bottom',mod5='biotin',mod3='C3')
 sc.ligation_end('top','Atop',kind='splinted',polarity='5-prime',owner='adapter');return sc

def oligos():return tuple(oligo(n,[Segment(n,x)],mods=m) for n,x,m in [('trP1top',TRP1_TOP,'finalbondphosphorothioate'),('trP1bottom',TRP1_BOTTOM,'5Phos;3-prime dideoxy'),('Atop',A_TOP,'5Phos;3-prime dideoxy'),('Abottom',A_BOTTOM,'5-primebiotin;3-primeC3'),('primerA',PRIMER_A,'')])
SOURCE_OLIGO_EXCEPTIONS={'trp1_top':'Source spells internal bond as -phosphorothioate-T; shared parser may not rejoin this English-word modification. Local source test checks full printed23base-phosphorothioate-T,retainingterminalT.'}

def inventory(g,key,label,nt_label):
 stock=g.input(key,label,K.MOLECULAR)
 enzyme=g.product(key+'_enzyme','Enzyme class from '+label,K.CATALYST)
 nts=g.product(key+'_nts',nt_label,K.MOLECULAR)
 buffer=g.product(key+'_buffer','Actual buffer/condition classes; full formulation unassembled',K.CONDITION)
 g.react(key+'_components','Copresent source supply classes; not physical separation or fresh duplicate additions',S.GENERIC,((stock,R.TEMPLATE),),(enzyme,nts,buffer),changes_molecule=False)
 return enzyme,nts,buffer

def allocation(g,key,x,label,qc=False):
 selected=g.product(key+'_selected',label,K.QC if qc else K.MOLECULAR,state=x.state,chemical_state=x.chemical_state)
 rest=g.product(key+'_remainder','Untransferred source DNA stock; amount/fate unassigned',K.MOLECULAR,state=x.state,chemical_state=x.chemical_state)
 g.react(key,label,S.SPLIT,((x,R.TEMPLATE),),(selected,rest))
 g.mark_terminal(Terminal(selected.id if qc else rest.id,Disposition.QC if qc else Disposition.UNRESOLVED,'Analytical allocation excluded from production'if qc else'Actual untransferred stock remains outside selected processing'))
 return rest if qc else selected

def protocol_graph():
 g=ProtocolGraph(TITLE,authority='Reijns2015 associated manuscript original Ionlibrary Methods; all printed sourceoligos')
 dna=g.input('dna','purified yeast genomic DNA with embedded rNMP',K.MOLECULAR,state=state('genomic DNA;positionsunassigned'),primary=True)
 def op(i,label,kind,x,condition,outstate,*,rk=K.CONDITION,rr=R.CONDITION):
  if i in ('endprep','block'):
   enzyme,nts,buffer=inventory(g,i+'_reagent','Actual '+condition+' aggregate supplies; original stock packaging unassigned','Actual ddATP contribution'if i=='block'else'Actual end-repair/dA-tailing nucleotide classes; exact mixture unassembled')
   ports=((enzyme,R.CATALYST),(nts,R.STRUCTURAL),(buffer,R.CONDITION))
  else:ports=((g.input(i+'_reagent',condition,rk),rr),)
  y=g.product(i+'_dna',label,K.MOLECULAR,state=outstate);g.react(i,label,kind,((x,R.TEMPLATE),)+ports,(y,));return y
 def clean(i,x,label,*,outstate=None,kind=K.MOLECULAR):
  r=g.input(i+'_condition',label,K.CONDITION);y=g.product(i+'_retained','recoveredDNA',kind,state=x.state if outstate is None else outstate,chemical_state=x.chemical_state);w=g.product(i+'_waste','purificationremovedfractions',K.WASTE)
  g.react(i,label,S.CLEANUP,((x,R.TEMPLATE),(r,R.CONDITION)),(y,w));g.discard(w,'Actualsourcepurification;noquantitativerecoveryassumption.');return y
 current=op('sonicate','Sonicate to about 400 bp',S.FRAGMENTATION,dna,'BioruptorPlus',state('sonicated genomic DNA'))
 current=clean('fragment_cleanup',current,'ethanolprecipitation;1.2xAMPure')
 current=allocation(g,'nanodrop_qc',current,'Actual NanoDrop analytical allocation before end preparation',qc=True)
 current=allocation(g,'endprep_input',current,'Actual up to 5 µg DNA selected for end preparation')
 current=op('endprep','NEBNext end repair and dA-tailing',S.END_REPAIR,current,'sourceNEBNextEndRepair/dA-Tailing',state('sourceA-tailedgenomicfragments'),rk=K.CATALYST,rr=R.CATALYST)
 current=clean('endprep_cleanup',current,'1.2xAMPure;subsequentbeadscaptureNaCl/PEG')
 top=g.input('trp1_top','trP1top;terminalbondphosphorothioate',K.OLIGO,state=Construct([Segment('trP1top',TRP1_TOP)]))
 bottom=g.input('trp1_bottom','trP1bottom;5Phos,3dideoxy',K.OLIGO,state=Construct([Segment('trP1bottom',TRP1_BOTTOM)]))
 trp1=g.product('trp1','sourceannealedtrP1adapterpair',K.ADAPTER,state=trp1_scene())
 g.react('anneal_trp1','AnnealprintedtrP1pair;95C5min,gradualcooling',S.HYBRIDIZATION,((top,R.TEMPLATE),(bottom,R.STRUCTURAL)),(trp1,))
 ligase,ligation_other,ligation_buffer=inventory(g,'ligation_mix','Actual NEBNext Quick Ligation aggregate supplies; stock packaging unassigned','Other actual ligation molecular supplies; identities/doses unassembled')
 joined=g.product('trp1_ligated','trP1-ligatedDNA',K.MOLECULAR,state=state('trP1-bearing genomic DNA;endassemblyunassigned'))
 g.react('trp1_ligate','LigatetrP1beforeendblocking',S.LIGATION,((current,R.TEMPLATE),(trp1,R.STRUCTURAL),(ligase,R.CATALYST),(ligation_other,R.STRUCTURAL),(ligation_buffer,R.CONDITION)),(joined,),geometry=AdapterLigationGeometry.known('cohesive',trp1_scene(),statement='Source trP1top3-T;blockedbottomterminalCretainedasprinted;genomicsourceA-tailed.'))
 current=op('block','Block old free 3′ ends',S.TAILING,joined,'TdT/ddATP,20Uperµg,37C2h',state('oldfree3-primeendsdd-blocked;notnewRNaseH2nicks',debris='excessddATPandreactionproducts'),rk=K.CATALYST,rr=R.CATALYST)
 current=clean('block_cleanup',current,'AMPurepurification;removebeads',outstate=state('old3-primeendclassdd-blocked'))
 current=op('nick','RNase H2 nicking',S.DIGESTION,current,'RNaseH2,10pmolperµg,37C2h;80C20min',state('new3-OHupstreamofrNMP;newdownstream5-P;oldendsblocked'),rk=K.CATALYST,rr=R.CATALYST)
 current=clean('nick_cleanup',current,'1.8xAMPure')
 current=op('sap','Remove 5′ phosphates',S.DEPHOSPHORYLATION,current,'SAP5U,37C1hperµg;65C15min',state('new3-OHretained;5-phosphatesremoved',debris='phosphataseproducts'),rk=K.CATALYST,rr=R.CATALYST)
 current=clean('sap_cleanup',current,'AMPurepurification',outstate=state('new3-OHvsoldblockedends;5-OHclass'))
 current=op('denature','Denature and snap-cool',S.THERMAL_TREATMENT,current,'95C5min;snapcool',current.state)
 adapter=g.input('a_adapter','sourceAadapter/biotinsplint',K.ADAPTER,state=a_scene())
 complex_state=Scene();complex_state.strand('captured',list(CAPTURED),label='INFERRED — covalentnonbiotinylatedstrand',mod3='dideoxy');complex_state.strand('splint',[Segment('biotinsplint',A_BOTTOM)],label='separate noncovalent biotin splint',mod5='biotin',mod3='C3')
 captured=g.product('captured','INFERRED — A-top-ligatedDNAwithbiotinylatedsplint',K.MOLECULAR,state=complex_state,inferred=True,chemical_state={'source_parent_end':'Original captured strand retains source A-top3prime dideoxy block; whole covalent junction context inferred'})
 g.react('a_ligate','INFERRED — A-top5-PjoinsnewDNA3-OH',S.LIGATION,((current,R.TEMPLATE),(adapter,R.STRUCTURAL),(ligase,R.CATALYST),(ligation_other,R.STRUCTURAL),(ligation_buffer,R.CONDITION)),(captured,),geometry=AdapterLigationGeometry.known('splinted',a_scene(),statement='Source6randomAbottomsplintsnewDNA3-OHtoAtop5-P;Atop3-ddandAbottom3-C3blocked.'))
 beads=g.input('streptavidin','M280streptavidinbeads;sourcevendorguidelines',K.SELECTION)
 bound=g.product('bound','biotin-splint-associatedDNAcapturedonbeads',K.MOLECULAR,state=complex_state,chemical_state=captured.chemical_state)
 flow=g.product('flowthrough','nonretainedcapturematerial',K.MOLECULAR,state=state('unretained material;exactcompositionunassigned'))
 g.react('capture_beads','Capture biotinylated adapter complexes',S.AFFINITY_CAPTURE,((captured,R.TEMPLATE),(beads,R.SELECTION)),(bound,flow));g.mark_terminal(Terminal(flow.id,Disposition.UNRESOLVED,'Nonretainedfractionnotsequenced;detailedphysicaldispositionnotreconstructed.'))
 naoh=g.input('naoh','Actual 0.15 M NaOH strand-release reagent',K.MOLECULAR)
 released=g.product('released','INFERRED — releasednonbiotinylatedcapturedDNA',K.MOLECULAR,state=CAPTURED,inferred=True,chemical_state=captured.chemical_state)
 splint_state=Scene();splint_state.strand('splint',[Segment('Abottom',A_BOTTOM)],label='biotinylatedsplintretainedwithbeads',mod5='biotin',mod3='C3');splint_state.note('splint','Biotin-streptavidinassociation;beadattachmentnotacovalentDNAnucleotide.')
 splint=g.product('bead_splint','biotinylatedsplintretainedwithbeads',K.MOLECULAR,state=splint_state)
 g.react('release','Release nonbiotinylated strand; retain biotin splint',S.STRAND_SEPARATION,((bound,R.TEMPLATE),(naoh,R.STRUCTURAL),(g.input('release_conditions','Source 0.15 M NaOH denaturation; actual additional temperature/time unassigned',K.CONDITION),R.CONDITION)),(released,splint));g.mark_terminal(Terminal(splint.id,Disposition.UNRESOLVED,'SourceNaOHreleasesnonbiotinylatedstrand;retainedbiotin-beadfractionnotsequenced,finaldisposalunassigned.'))
 current=clean('release_cleanup',released,'ethanolprecipitationofreleasedssDNA')
 primer=g.input('primer_a','sourceprimerA',K.OLIGO,state=Construct([Segment('primerA',PRIMER_A)]))
 pol,flash_nt,flash_buffer=inventory(g,'phusion','Actual Phusion FLASH High-Fidelity PCR Master Mix for second-strand synthesis; amount unassigned','Nucleotides copresent in the same actual prepared Phusion FLASH stock')
 copied=g.product('copied','INFERRED — second-strandlibraryreference',K.MOLECULAR,state=FINAL_LIBRARY,inferred=True)
 g.react('second_strand','INFERRED — primerAsecondstrandsynthesis,nopre-sequencingPCRseries',S.EXTENSION,((current,R.TEMPLATE),(primer,R.STRUCTURAL),(pol,R.CATALYST),(flash_nt,R.STRUCTURAL),(flash_buffer,R.CONDITION)),(copied,))
 current=clean('gel',copied,'2% E-Gel EX 200–300bp selection')
 current=allocation(g,'bioanalyzer_qc',current,'Actual Bioanalyzer analytical allocation before OneTouch emPCR',qc=True)
 current=allocation(g,'em_pcr_input',current,'Actual selected source library contribution for OneTouch emPCR')
 f=g.input('ion_bead_primer','OneTouchemPCRbeadprimer;exactpolymerunassigned',K.OLIGO,state=molecule('beadprimer',(role('exactIonbeadprimerunassigned'),)));r=g.input('ion_other_primer','OneTouchemPCRsecondprimer;exactpolymerunassigned',K.OLIGO,state=molecule('emPCRprimer',(role('exactsecondprimerunassigned'),)));mix,ion_nt,ion_buffer=inventory(g,'ion_mix','Actual OneTouch emPCR aggregate supplies; original stock packaging unassigned','Actual emPCR nucleotide class; whole formulation/doses unassembled')
 library=g.product('ion_library','sourceIonemPCRsequencinglibrary',K.LIBRARY,state=FINAL_LIBRARY,inferred=True)
 g.react('emulsion_pcr','Ion Torrent OneTouch emulsion PCR after library QC',S.PCR,((current,R.TEMPLATE),(f,R.STRUCTURAL),(r,R.STRUCTURAL),(mix,R.CATALYST),(ion_nt,R.STRUCTURAL),(ion_buffer,R.CONDITION)),(library,),geometry=PCRGeometry.unavailable(FINAL_LIBRARY,reason='SourceOneTouchplatformspecified;emPCRoligos/conditionsandbead-strandgeometryunassembled.'))
 g.sequence(library,SequencingPlan(library=FINAL_LIBRARY,limitation=SEQUENCING_UNAVAILABLE,acquired_roles=('Read 1',)))
 g.declare_reviewed_scope(('dna','sonicate_reagent','sonicate_dna','sonicate','fragment_cleanup_condition','fragment_cleanup_retained','fragment_cleanup_waste','fragment_cleanup','endprep_reagent','endprep_dna','endprep','endprep_cleanup_condition','endprep_cleanup_retained','endprep_cleanup_waste','endprep_cleanup','trp1_top','trp1_bottom','trp1','anneal_trp1','ligation_mix','trp1_ligated','trp1_ligate','block_reagent','block_dna','block','block_cleanup_condition','block_cleanup_retained','block_cleanup_waste','block_cleanup','nick_reagent','nick_dna','nick','nick_cleanup_condition','nick_cleanup_retained','nick_cleanup_waste','nick_cleanup','sap_reagent','sap_dna','sap','sap_cleanup_condition','sap_cleanup_retained','sap_cleanup_waste','sap_cleanup','denature_reagent','denature_dna','denature','a_adapter','captured','a_ligate','streptavidin','bound','flowthrough','capture_beads','naoh','released','bead_splint','release','release_cleanup_condition','release_cleanup_retained','release_cleanup_waste','release_cleanup','primer_a','phusion','copied','second_strand','gel_condition','gel_retained','gel_waste','gel','ion_bead_primer','ion_other_primer','ion_mix','ion_library','emulsion_pcr'),g.authority)
 g.declare_reviewed_scope(('release_conditions','endprep_reagent_components', 'block_reagent_components', 'ligation_mix_components', 'phusion_components', 'ion_mix_components', 'nanodrop_qc', 'endprep_input', 'bioanalyzer_qc', 'em_pcr_input'),'Same accepted PMC4374164 Methods149–159: actual NanoDrop/up-to5ug contribution, TdT/ddATP, Quick Ligation supplies, Phusion FLASH prepared master mix, Bioanalyzer and OneTouch emPCR supply classes. Bookkeeping views do not invent stock packaging or extra doses.')
 return g.require_publishable()
def render_page():
 from page import head,caveat,table
 from protocol_views import assemble_page,WorkflowSection
 from chemdraw import panel
 g=protocol_graph()
 return assemble_page(g,(head(TITLE),'<main class="wrap" data-protocol-complete="true"><h1>'+TITLE+'</h1><p>'+SUMMARY+'</p><p>'+SOURCE+' · <a href="01_emriboseq.html">Research note and source boundaries</a></p>',
 caveat('<strong>INFERRED — strand assembly and genomic context.</strong> This page follows the associated accepted manuscript, not a fully verified final-version recipe. The printed original Ion adapters are retained; captured genomic sequence, complete strand assembly and run primers remain unassigned. Random N6 is splint degeneracy, not a UMI. No library-PCR series is inserted before second-strand synthesis; platform OneTouch emPCR remains separate.'),
 '<h2>Source preparation</h2><ol><li>Sonicate DNA to about 400 bp, purify and quantify. Select up to 5 µg for end repair/dA-tailing, then ligate trP1 before blocking old free 3′OH ends with TdT/ddATP.</li><li>RNaseH2 creates new 3′OH ends upstream of embedded ribonucleotides. Remove 5′ phosphates with SAP, then denature and snap-cool.</li><li>Ligate A-top 5′ phosphate to a new DNA 3′OH through the random splint. A-top is 3′ dideoxy-blocked; A-bottom has 5′ biotin and a 3′ C3 block.</li><li>Capture on M-280 beads. NaOH releases nonbiotinylated DNA while the biotinylated splint remains with beads.</li><li>PrimerA and the actual Phusion FLASH master mix synthesize the complement. Select 200–300 bp DNA, allocate Bioanalyzer QC material and carry a selected library contribution into OneTouch emPCR.</li></ol>',
 '<h2>Known adapter geometry</h2>',panel(trp1_scene().rows(),caption='Source printed trP1 pair: terminal T and blocked bottom C retained; the terminal C outside the maximal paired core is not silently removed'),panel(a_scene().rows(),caption='Source A adapter and noncovalent biotinylated random splint; whole occupied genomic binding context remains inferred'),
 '<h2>Read setup and interpretation</h2><p>Ion Torrent PGM or Proton after OneTouch emulsion PCR. Exact flow/read-length configuration, emPCR primer bodies and loaded run primers remain unassembled. The source maps a ribonucleotide one base upstream and on the opposite strand from the read 5′ end; that coordinate is distinct from the captured upstream DNA 3′OH. The downstream ribonucleotide is not assigned to the captured fragment’s terminus.</p>',
 '<h2>Source scope</h2><p>The original trP1 bottom strand contains a blocked terminal base outside its maximal exact duplex. This page keeps that printed base and its separate unpaired position. Modifications on the original adapter/template are not copied onto newly synthesized DNA. Full VOR, RNaseH2 preparation reference and later detailed emRiboSeq/EndoSeq methods remain outside this source pass. No Illumina recipe is substituted.</p>',
 '<details><summary>Full molecular workflow, physical fractions and unresolved chemistry</summary>',WorkflowSection('INFERRED strand contexts — accepted-source Ion preparation'),'</details></main>'))
def selftest():
 import re
 from checks import Check,Source,have
 from storage import sources_root
 c=Check();source=Source(sources_root()/HERE.parent.name/'nature14183_PMC4374164.xml.txt','originalIonOligonucleotidesandadapterdesign')
 if have(c,source,label='originalIonoligos'):
  t=re.sub(r'\s+','',source.path.read_text())
  for n,x in [('trP1top',TRP1_CORE+'-phosphorothioate-T'),('trP1bottom',TRP1_BOTTOM+'-dideoxy'),('Atop',A_TOP+'-dideoxy'),('Abottom',A_BOTTOM+'-C3'),('primerA',PRIMER_A)]:c('source '+n,x in t,True)
 g=protocol_graph();render_page()
 c('ddATP is a molecular contributor, not hidden in TdT',g.materials['block_reagent_nts'].kind,K.MOLECULAR)
 c('Phusion NT comes from actual master-mix stock',g.steps['phusion_components'].inputs[0].material,'phusion')
 c('Second strand uses same-stock nucleotides',any(i.material=='phusion_nts'for i in g.steps['second_strand'].inputs),True)
 c('Original captured3prime block retained through cleanup',g.materials['release_cleanup_retained'].chemical_state,g.materials['captured'].chemical_state)
 c('Template3primeblock not propagated to copied DNA',g.materials['copied'].chemical_state,None)
 c('Bioanalyzer allocation excluded from production',g.terminals['bioanalyzer_qc_selected'].disposition,Disposition.QC)
 c('PlatformPCR gets actual selected library portion',g.steps['emulsion_pcr'].inputs[0].material,'em_pcr_input_selected')
 c('OneTouch NT class is explicit',any(i.material=='ion_mix_nts'for i in g.steps['emulsion_pcr'].inputs),True)
 c('Publication gate',g.publication_problems(),[])
 c.report()
