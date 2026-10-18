"""Source-owned native 2023 spatial epigenome/transcriptome graph.

All molecular choices below belong to this declared assay. The 2023 paper explicitly
assigns RT/ligations/cDNA processing to both assays; no later protocol is substituted.
"""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Construct,Scene,Segment,feature,revcomp,oligo,complement_segments
from protocol_model import role
from protocol_graph import ProtocolGraph,MaterialKind as M,InputRole as I,StepKind as K,SequencingPlan,Terminal,Disposition
from adapter_ligation import AdapterLigationGeometry
from multimodal_spatial import HybridizationLigationRound
from pcr import PCRGeometry
import illumina as il
import nextera as nx
import seqprimers as sp

CUTTAG = True
TITLE = 'spatial-CUT&amp;Tag–RNA-seq'
NOTES = '01_spatial_cuttag_rna_seq.html'
SOURCE = 'Zhang et al., <a href="https://doi.org/10.1038/s41586-023-05795-1">Nature2023, published Methods244–292 and SupplementaryTables4–6</a>.'
SUMMARY = 'Independent tissue DNA and RNA-derived cDNA receive the same two channel coordinates; streptavidin separates biotinylated originals before distinct native library preparations.'
CAVEAT = ('The graph uses the2023 publication. A1/B1 and N501/N701 are source-printed reference members; actual study-specific primer/index assignments remain unverified. The DNA reference readout is INFERRED, including compatible standardNextera sequencing primers. Native transposition gaps, first-cycle completion and exact RNA-library primer ends remain unassigned. Biotin stays on the original RT-primer cDNA; PCR copying does not reproduce that modification. Residual RNA/DNA and other molecular classes have unresolved recovery unless a physical separation is specified.')
FINAL_LIBRARY = None
SEQ_PRIMERS = ()
SEQUENCING_ENDING = 'Source native acquisition: NovaSeq6000 paired-end150 for both libraries. DNARead1 maps genomic DNA; RNARead1 maps transcript-derived DNA; Read2 supplies spatial identifiers and RNAUMI where present. Native index-read acquisition is unassigned; the DNA compatible-primer panel is a conditional reference.'
SEQUENCING_UNAVAILABLE = ''
A_FEATURE=feature('spatial_a','spatial_barcode','whitelist',whitelist='Zhang2023 Table5 A1–A100')
B_FEATURE=feature('spatial_b','spatial_barcode','whitelist',whitelist='Zhang2023 Table6 B1–B100')
UMI=feature('spatial_rna_umi','umi','random')
I5=feature('reference_i5','sample_index','fixed',note='Table4 N501 reference, not study assignment')
I7=feature('reference_i7','sample_index','fixed',note='Table4 N701 reference, not study assignment')
RT_HANDLE='CATCGGCGTACGACT'
RT_UMI='N'*10
RT_POLYT='T'*15+'VN' # /iBiodT/ is an additional physical T before these15T.
TN5_A=nx.S5+nx.ME
TN5_B=RT_HANDLE+nx.ME
TN5_REV=nx.ME_RC
A_PREFIX='AGGCCAGAGCATTCG';A1='AACGTGAT';A_SUFFIX='GTGGCCGATGTTTCG'
B_PREFIX='CAAGCGTTGGCTTCTCGCATCT';B1='AACGTGAT';B_SUFFIX='ATCCACGTGCTTGAG'
LINKER1='AGTCGTACGCCGATGCGAAACATCGGCCAC'
LINKER2='CGAATGCTCTGGCCTCTCAAGCACGTGGAT'
TSO_WRITTEN='AAGCAGTGGTATCAACGCAGAGTGAATrGrG+G'
TSO_DNA='AAGCAGTGGTATCAACGCAGAGTGAATGGG'
PCR2='AAGCAGTGGTATCAACGCAGAGT'
N501=il.P5+'TAGATCGC'+nx.S5+nx.ME
N701=il.P7+'TCGCCTTA'+nx.S7+nx.ME+B_PREFIX
# The obligatory substrate5′P in this neutral nick constructor is local source
# evidence: RT/Tn5ME-B carry5P, and Table5 barcodeA carries5P for the second nick.
A_ROUND=HybridizationLigationRound('barcode A',RT_HANDLE,A_PREFIX,A_SUFFIX,LINKER1,A_FEATURE,phosphorylated_5=True,barcode_length=8)
B_ROUND=HybridizationLigationRound('barcode B',A_PREFIX,B_PREFIX,B_SUFFIX,LINKER2,B_FEATURE,phosphorylated_5=False,barcode_length=8)

def ss(label,parts,**kw):
 sc=Scene();sc.strand('strand',list(parts),label=label,**kw);return sc

def barcode_parts(stage):
 a=[Segment('A outgoing arm',A_PREFIX),Segment('A barcode','N'*8,placeholder=True,feature=A_FEATURE),Segment('A ligating arm',A_SUFFIX)]
 b=[Segment('B PCR arm',B_PREFIX),Segment('B barcode','N'*8,placeholder=True,feature=B_FEATURE),Segment('B ligating arm',B_SUFFIX)]
 return (b if stage>=2 else [])+(a if stage>=1 else [])

def rt_parts(*,copied=False):
 return [Segment('RT handle',RT_HANDLE),Segment('RNA UMI',RT_UMI,'umi',placeholder=True,feature=UMI),Segment('unmodified copied T' if copied else 'internal biotin-dT','T',note='ordinary DNA copy; biotin not copied' if copied else 'physical /iBiodT/'),Segment('dT15','T'*15),Segment('VN anchor','VN',placeholder=True)]

def native_state(stage=0,*,tagged=False,rt=False,switched=False,rna_only=False,dna_only=False):
 sc=Scene()
 if not rna_only:
  dna=([Segment('transferred B5P handle',RT_HANDLE),Segment('transferred mosaic end',nx.ME),role('genomic DNA; opposite transferred/nontransfer/gap geometry unassigned')] if tagged else [role('original tissue genomic DNA')])
  sc.strand('DNA',barcode_parts(stage)+dna,label='native DNA population; disclosed5′ arm' if tagged else 'native genomic DNA',mod5='p' if tagged and stage<2 else '')
  if tagged:sc.strand('opposite DNA',[role('opposite genomic strand and nontransfer/gap geometry unassigned')],label='independent complementary DNA component')
 if not dna_only:
  sc.strand('RNA',[role('original RNA and possible RNA-derived fragments; recovery unassigned')],label='independent RNA population')
  if rt:
   parts=barcode_parts(stage)+rt_parts()+[role('first-strand transcript cDNA; native hybrid alignment unassigned')]
   if switched:parts.append(Segment('DNA copy of TSO',revcomp(TSO_DNA),inferred=True,note='normalDNA copied from modified TSO; RNA/LNA chemistry not copied'))
   sc.strand('cDNA',parts,label='original biotinylated ss cDNA'+(' on streptavidin' if rna_only else ''),mod5='p' if stage<2 else '')
   sc.note('cDNA','Physical biotin remains at RT-primer dT. Source ss/hybrid second-strand inventory is not replaced by a completed amplicon.')
 if tagged or rt:
  sc.strand('carried classes',[role('possible unconverted DNA/RNA/nonproductive and reaction constituents; partition unassigned')],label='independent co-carried populations')
 for key in ('DNA','cDNA'):
  if key not in sc.strands:continue
  if stage>=1:sc.junction(key,'A ligating arm','transferred B5P handle' if key=='DNA' else 'RT handle','sealed ligation A')
  if stage>=2:sc.junction(key,'B ligating arm','A outgoing arm','sealed ligation B')
 return sc

def transposome_scene():
 sc=Scene()
 for key,seq,phos in [('A',TN5_A,''),('B',TN5_B,'p')]:
  arm=seq[:-len(TN5_REV)]
  sc.strand('transfer'+key,[Segment(key+' outer arm',arm),Segment(key+' mosaic end',nx.ME)],label='Tn5ME-'+key,mod5=phos)
  sc.anneal('nontransfer'+key,[Segment(key+' nontransfer',TN5_REV)],to='transfer'+key,pair=(key+' nontransfer',key+' mosaic end'),label='Tn5MErev',mod5='p')
 return sc

def dna_reference():
 # Mathematical N501/N701 product reference only. Native PCR geometry below
 # receives actual source-bounded gapped intermediates, never this finalduplex.
 body=[Segment('Tn5 A outer',nx.S5),Segment('Tn5 A ME',nx.ME),role('genomic insert; length unassigned')]
 right=[Segment('Tn5 B ME rc',nx.ME_RC),Segment('shared handle rc',revcomp(RT_HANDLE)),Segment('A ligating arm rc',revcomp(A_SUFFIX)),Segment('A barcode rc','N'*8,placeholder=True,feature=A_FEATURE),Segment('A outgoing arm rc',revcomp(A_PREFIX)),Segment('B ligating arm rc',revcomp(B_SUFFIX)),Segment('B barcode rc','N'*8,placeholder=True,feature=B_FEATURE),Segment('B PCR arm rc',revcomp(B_PREFIX))]
 parts=[Segment('P5',il.P5),Segment('reference i5','TAGATCGC',feature=I5),*body,*right,Segment('R2 ME rc',nx.ME_RC),Segment('S7 rc',nx.S7_RC),Segment('reference i7 rc',revcomp('TCGCCTTA'),feature=I7),Segment('P7 rc',il.P7_RC)]
 for s in parts:s.inferred=True
 return Construct(parts,name='INFERRED N501/N701 representative DNA product')
DNA_REFERENCE=dna_reference()
DNA_READ_PRIMERS=(sp.NEXTERA['R1'],sp.NEXTERA['R2'])
FINAL_CAPTIONS={'dna_library':'INFERRED — N501/N701 representative DNA product and compatible standardNextera readout. The source prints these ordered oligos; actual sample assignment, native first-cycle completion and loaded sequencing primers remain unverified. Genomic length is unassigned.'}

def dna_copies():
 sc=Scene.duplex(list(DNA_REFERENCE),label='INFERRED representative DNA PCR daughter',unpaired=tuple(s.name+"'" for s in DNA_REFERENCE if s.is_role_token()))
 return sc

def rna_copies():
 parts=barcode_parts(2)+rt_parts(copied=True)+[role('transcript cDNA copy'),Segment('normalDNA TSO-derived end',revcomp(TSO_DNA),inferred=True)]
 sc=Scene.duplex(parts,label='RNA-derived PCR daughters; original biotin not copied',unpaired=('transcript cDNA copy\'',))
 sc.strand('possible free classes',[role('possible freeRNA/fragments/primers/non-targetDNA; component-specific partition unassigned')],label='independent possible freeco-carriedclasses')
 sc.note('top','These are copies initiated by PCR Primer1/2; neither primer introduces RT-primer biotin or TSO RNA/LNA modifications.')
 return sc

def rna_tagged(stage):
 sc=ss('native RNA-derived '+stage,[role('native NexteraXT partial ends/gaps; ordered primer geometry unassigned'),role('fragmented RNA-derived DNA; tag-bearing subset source-dependent')])
 sc.strand('unplaced identifiers',[role('retained B feature where present',molecular_feature=B_FEATURE),role('retained A feature where present',molecular_feature=A_FEATURE),role('retained RNA UMI where present',molecular_feature=UMI)],label='identifier inventory; native strand placement unassigned')
 sc.strand('possible free classes',[role('possible freeRNA/fragments/primers/non-targetDNA; component-specific recovery unassigned')],label='independent possible freeco-carriedclasses')
 sc.note('strand','Other fragments/originals are not physically removed by tagmentation or PCR; recovery remains unassigned.')
 return sc

def oligos():
 for name,seq,mods in [('Tn5ME-A',TN5_A,''),('Tn5ME-B',TN5_B,'/5Phos/'),('Tn5MErev',TN5_REV,'/5Phos/'),('A1',A_PREFIX+A1+A_SUFFIX,'/5Phos/'),('B1',B_PREFIX+B1+B_SUFFIX,''),('linker1',LINKER1,''),('linker2',LINKER2,''),('PCR Primer1',B_PREFIX,''),('PCR Primer2',PCR2,''),('N501 reference',N501,''),('N701 reference',N701,'')]:yield f'{name}:5′-{mods}{seq}-3′'
 yield 'RT primer:5′-/5Phos/'+RT_HANDLE+RT_UMI+'/iBiodT/'+RT_POLYT+'-3′'
 yield 'Template switch primer:5′-'+TSO_WRITTEN+'-3′;terminal rG/rG/LNA-G'

def protocol_graph():
 g=ProtocolGraph(TITLE,authority=SOURCE)
 def inp(key,label,kind=M.CONDITION,state=None,inferred=False):return g.input(key,label,kind,state=state,inferred=inferred,source=SOURCE)
 def op(key,label,kind,x,state,reagents=(),**kw):
  out=g.product(key+'_out',label,M.MOLECULAR,state=state)
  g.react(key,label,kind,((x,I.TEMPLATE),*reagents),(out,),**kw);return out
 def unresolved(x,why):g.mark_terminal(Terminal(x.id,Disposition.UNRESOLVED,why))
 def recover(key,label,x):
  beads=inp(key+'_selection',label,M.SELECTION);kept=g.product(key+'_out','retained '+x.label,M.MOLECULAR,state=x.state)
  waste=g.product(key+'_removed','physically removed supernatant/washes',M.WASTE)
  g.react(key,label,K.CLEANUP,((x,I.TEMPLATE),(beads,I.SELECTION)),(kept,waste),changes_molecule=False);g.discard(waste,'Physical non-retained fraction; complete species-specific depletion is not inferred.');return kept
 def amp(key,label,x,state,fwd,rev,pol,conditions=()):
  out=g.product(key+'_out',label,M.MOLECULAR,state=state,inferred=True)
  rest=g.product(key+'_originals','original templates/nonselected molecular classes',M.MOLECULAR,state=x.state)
  g.react(key,label,K.PCR,((x,I.TEMPLATE),(fwd,I.STRUCTURAL),(rev,I.STRUCTURAL),(pol,I.CATALYST),*tuple((c,I.CONDITION) for c in conditions)),(out,rest),geometry=PCRGeometry.unavailable(x.state,reason='Actual incoming native ss/gapped source material shown; complete first-cycle/native end geometry unassigned. A conditional completed-product reference is a separate output.'))
  unresolved(rest,'PCR does not erase original templates or physically remove co-carried molecular classes; later recovery unassigned.');return out
 tissue=g.input('tissue','fresh-frozen tissue section',M.BIOLOGICAL,state=native_state(),primary=True)
 fix=inp('fixative','0.2% formaldehyde5min;glycine quench')
 x=op('fix','Fix and quench tissue',K.CROSS_LINKING,tissue,tissue.state,((fix,I.CONDITION),))
 x=recover('fix_wash','Two source tissue washes after fixation',x)
 perm=inp('perm','source assay-specific digitonin/NP40 permeabilization')
 x=op('permeabilize','Permeabilize fixed tissue',K.PERMEABILIZATION,x,x.state,((perm,I.CONDITION),))
 if CUTTAG:
  antibody=inp('primary_antibody','selected histone-mark primary antibody;4°C overnight',M.MOLECULAR)
  x=op('bind_primary','Bind primary antibody',K.BINDING,x,x.state,((antibody,I.STRUCTURAL),))
  secondary=inp('secondary','guinea-pig anti-rabbit secondary;30minRT',M.MOLECULAR)
  x=op('bind_secondary','Bind secondary antibody',K.BINDING,x,x.state,((secondary,I.STRUCTURAL),));x=recover('antibody_wash','Source wash after secondary antibody',x)
 adapters=inp('transposome','source Tn5ME-A and5P-Tn5ME-B each paired with5P-MErev',M.ADAPTER,transposome_scene())
 enzyme=inp('tn5','pA-Tn5' if CUTTAG else 'Tn5',M.CATALYST)
 if CUTTAG:
  x=op('tether_tn5','Bind pA-Tn5 adaptor complex1hRT in300mMNaCl',K.BINDING,x,x.state,((adapters,I.STRUCTURAL),(enzyme,I.CATALYST)));x=recover('tn5_wash','Source300-wash beforeMg activation',x)
 activation=inp('activation','10mMMgCl2;37°C1h' if CUTTAG else 'source transposition mixture;37°C30min')
 tagged=g.product('tagged','tagmented DNA and independent RNA',M.MOLECULAR,state=native_state(tagged=True))
 other=g.product('tag_other','other transposition/unconverted DNA classes',M.MOLECULAR,state=x.state)
 g.react('tagment','Tag native DNA',K.TAGMENTATION,((x,I.TEMPLATE),(adapters,I.STRUCTURAL),(enzyme,I.CATALYST),(activation,I.CONDITION)),(tagged,other));unresolved(other,'Off-products/unconverted DNA are not chemically discarded; later recovery unassigned.')
 edta=inp('edta','40mMEDTA5minRT withRNase inhibitor')
 x=op('stop_tn5','Stop transposition',K.QUENCHING,tagged,tagged.state,((edta,I.CONDITION),));x=recover('tag_wash','Two0.5×DPBS-RI washes',x)
 rtprimer=inp('rt_primer','5P RT primer withN10 andphysical internalbiotin-dT',M.OLIGO,ss('native RT primer',rt_parts(),mod5='p'))
 rt=inp('rt_enzyme','MaximaHMinus+dNTPs;RT30min then42°C90min',M.CATALYST)
 x=op('rt','Copy RNA in situ',K.REVERSE_TRANSCRIPTION,x,native_state(tagged=True,rt=True),((rtprimer,I.STRUCTURAL),(rt,I.CATALYST)));x=recover('rt_wash','NEBuffer3.1/RNaseinhibitor tissue wash',x)
 ligase=inp('ligase','T4DNA ligase;37°C30min each channel set',M.CATALYST)
 for key,round_,stage in [('a',A_ROUND,1),('b',B_ROUND,2)]:
  adapter=inp('barcode_'+key,'source barcode/linker preannealed plate;representative1',M.ADAPTER,round_.adapter_scene())
  geom=AdapterLigationGeometry.known('splinted',round_.adapter_scene(),x.state,statement='Source Tables4–6 plus computed antiparallel nick; preceding5P is source RT/Tn5ME-B forA andsource5P-A forB.')
  y=g.product('lig_'+key+'_out','source channel-'+key+' ligated populations',M.MOLECULAR,state=native_state(stage,tagged=True,rt=True))
  residual=g.product('lig_'+key+'_residual','unligated/partially ligated molecular classes',M.MOLECULAR,state=x.state)
  g.react('lig_'+key,'Ligate spatial barcode'+key.upper(),K.LIGATION,((x,I.TEMPLATE),(adapter,I.STRUCTURAL),(ligase,I.CATALYST)),(y,residual),geometry=geom);unresolved(residual,'Ligation is not physical purification; subsequent recovery of unligated classes unassigned.')
  x=recover('lig_'+key+'_wash','Source post-ligation tissue wash',y)
 prot=inp('proteinase','proteinaseK/SDS58°C2h then65°Covernight',M.CATALYST)
 x=op('reverse_crosslink','Lyse and reverse-crosslink',K.DEPROTEINIZATION,x,x.state,((prot,I.CATALYST),));x=recover('lysate_cleanup','ZymoDCC-5;elute100µl',x)
 strept=inp('streptavidin','DynabeadsMyOneStreptavidinC1;1hRT',M.SELECTION)
 dna=g.product('dna_supernatant','DNA-enriched supernatant; possible carriedRNA/classes unassigned',M.MOLECULAR,state=native_state(2,tagged=True,rt=False))
 cdna=g.product('rna_bound','biotinylated original cDNA-enriched bead fraction',M.MOLECULAR,state=native_state(2,tagged=False,rt=True,rna_only=True))
 g.react('separate','Affinity enrich original biotinylated cDNA; retainDNA supernatant',K.AFFINITY_CAPTURE,((x,I.TEMPLATE),(strept,I.SELECTION)),(dna,cdna))
 dna=recover('dna_cleanup','DNA supernatantZymoDCC-5;elute20µl',dna)
 p5=inp('dna_p5','source P5PCR primer;N501 conditional reference',M.OLIGO,ss('Table4N501 reference',[Segment('N501',N501)]),inferred=True)
 p7=inp('dna_p7','source indexed i7PCR primer;N701 conditional reference',M.OLIGO,ss('Table4N701 reference',[Segment('N701',N701)]),inferred=True)
 pol=inp('dna_pol','NEBNextMasterMix;source assay-specific PCR programme',M.CATALYST)
 if CUTTAG:
  hold=inp('cut_hold','58°C5min sourceprehold')
  dna=op('dna_prehold','Source58°Cprehold withprimers/polymerase',K.THERMAL_TREATMENT,dna,dna.state,((hold,I.CONDITION),(p5,I.STRUCTURAL),(p7,I.STRUCTURAL),(pol,I.CATALYST)))
 initial=inp('dna_initial','72°C5min beforecycling')
 dna=op('dna_initial_extension','Source72°C initial extension; native gap completion unassigned',K.EXTENSION,dna,dna.state,((initial,I.CONDITION),(p5,I.STRUCTURAL),(p7,I.STRUCTURAL),(pol,I.CATALYST)))
 dna=amp('dna_pcr','INFERRED sourceDNA PCR;12cycles' if CUTTAG else 'INFERRED sourceDNA PCR;initial5cycles',dna,dna_copies(),p5,p7,pol)
 if not CUTTAG:
  allocation=inp('qpcr_rule','5µl qPCR aliquot;45µl remaining preamp')
  qc=g.product('qpcr_input','5µl forqPCR',M.MOLECULAR,state=dna.state);rest=g.product('dna_main','45µl main preamp',M.MOLECULAR,state=dna.state)
  g.react('qpcr_sample','Sample actualqPCR aliquot',K.SAMPLING,((dna,I.TEMPLATE),(allocation,I.CONDITION)),(qc,rest),changes_molecule=False)
  newp5=inp('qpcr_p5','source newP5PCRprimer;ordered identity unassigned',M.OLIGO)
  qc=amp('qpcr','20cycle sourceqPCR withSYBRGreen;determineone-thirdplateau',qc,dna_copies(),newp5,p7,pol);unresolved(qc,'Measurement products not physically recombined with45µl mainlibrary; postmeasurementfate unassigned.')
  threshold=inp('additional_cycles','sourceqPCR-selected additionalcyclecount;sample-specific')
  dna=amp('dna_additional','INFERRED amplify45µl mainfraction tosourceqPCR-selected cyclecount',rest,dna_copies(),p5,p7,pol,conditions=(threshold,))
 dna=recover('dna_final_cleanup','1.3×AMPureXP' if CUTTAG else '1×AMPureXP',dna)
 library=g.product('dna_library','INFERRED representative nativeDNA library',M.LIBRARY,state=dna.state,inferred=True)
 g.react('retain_dna','RetainDNA library',K.GENERIC,((dna,I.TEMPLATE),),(library,),changes_molecule=False)
 g.sequence(library,SequencingPlan(library=DNA_REFERENCE,primers=DNA_READ_PRIMERS,acquired_roles=('Read 1','Read 2'),read_lengths={'Read 1':150,'Read 2':150},limitation='INFERRED N501/N701 compatible-standard-primer reference. Native study primer/index member, loaded sequencing-primer identities and first-cycle gap completion remain unverified; sourcePE150known, native indexacquisition unassigned.'))
 cdna=recover('bound_wash','TwoB&W/Tween washes thenTris/Tween',cdna)
 tso=inp('tso','physical sourceTSO withterminalrGrG+G',M.OLIGO,ss('modifiedTSO',[Segment('TSO',TSO_DNA,note='terminalRNA/RNA/LNA-G; modificationsnotcopied')]))
 cdna=op('switch','On-beadtemplate switch;MaximaHMinus/RT30min/42°C90min',K.TEMPLATE_SWITCHING,cdna,native_state(2,tagged=False,rt=True,switched=True,rna_only=True),((tso,I.STRUCTURAL),(rt,I.CATALYST)))
 cdna=recover('switch_wash','Tris/Tween thenwater bead washes',cdna)
 f=inp('rna_f','sourcePCRPrimer1',M.OLIGO,ss('PCRPrimer1',[Segment('P1',B_PREFIX)]));r=inp('rna_r','sourcePCRPrimer2',M.OLIGO,ss('PCRPrimer2',[Segment('P2',PCR2)]));kapa=inp('kapa','KAPAHiFiHotStart;source5initialcycles',M.CATALYST)
 original_bound=cdna
 first_mixture=rna_copies()
 first_mixture.strand('bead originals',[role('original biotinylated cDNA/RNA and possible co-carried classes; remain on beads')],label='original bead-associated population retained through first5cycles')
 first=g.product('rna_first_mixture','first5cycle PCR mixture with original beads and free daughterDNA',M.MOLECULAR,state=first_mixture)
 g.react('rna_first_pcr','Amplify actual bead-boundss/hybrid cDNA;5cycles',K.PCR,((cdna,I.TEMPLATE),(f,I.STRUCTURAL),(r,I.STRUCTURAL),(kapa,I.CATALYST)),(first,),geometry=PCRGeometry.unavailable(cdna.state,reason='Actual incoming on-bead ssDNA/RNA-hybrid original is explicit; first-cycle completion unassigned. Original bead-associated material remains in the output until physical removal.'))
 cdna=first
 bead_rule=inp('bead_removal_selection','Source removebeads after5PCRcycles',M.SELECTION)
 free=g.product('free_cdna','free daughterDNA inPCRsolution',M.MOLECULAR,state=rna_copies());beads=g.product('removed_beads','removed originalbiotinDNA/matrix and possibleco-carriedclasses',M.MOLECULAR,state=original_bound.state)
 g.react('remove_beads','Physically remove originalcDNA beads',K.SELECTION,((cdna,I.TEMPLATE),(bead_rule,I.SELECTION)),(free,beads));unresolved(beads,'Source removes beads; subsequentstorage/disposal unassigned. Originalbiotin/RNA classes are not chemically erased by PCR.')
 qpol=inp('rna_qpol','KAPA/SYBRGreen;up to15cycles stopatplateau',M.CATALYST)
 cdna=amp('rna_more_pcr','Amplify freecDNA daughters tosourceplateau',free,rna_copies(),f,r,qpol);cdna=recover('rna_preamp_cleanup','0.8×AMPureXP;elute20µl',cdna)
 mass=inp('rna_1ng_rule','1ng selected cDNA input intoXT')
 aliquot=g.product('rna_1ng','1ng forXT',M.MOLECULAR,state=cdna.state);excess=g.product('rna_preamp_rest','remainingpreampnotallocated',M.MOLECULAR,state=cdna.state)
 g.react('sample_rna','Take1ng fornativeXTpreparation',K.SAMPLING,((cdna,I.TEMPLATE),(mass,I.CONDITION)),(aliquot,excess),changes_molecule=False);unresolved(excess,'Unsampledsourcepreamp remainder; no physical discard prescribed.')
 xt=inp('xt','nativeATM-loaded adapters;ordered strands unassigned',M.ADAPTER);xtpol=inp('xt_tn5','NexteraXTATM/TagmentDNAbuffer55°C5min',M.CATALYST)
 tagged=op('rna_tagment','Tagment amplifiedcDNA withsourceXT',K.TAGMENTATION,aliquot,rna_tagged('gappedtagmentationproducts'),((xt,I.STRUCTURAL),(xtpol,I.CATALYST)))
 nt=inp('nt','NTbuffer5minRT')
 tagged=op('rna_stop','Neutralize nativeXT',K.QUENCHING,tagged,tagged.state,((nt,I.CONDITION),))
 fp=inp('rna_index_p5','sourceP5PCRprimer orderedstrandunassigned',M.OLIGO);rp=inp('rna_index_p7','sourceindexedP7PCRprimer orderedstrandunassigned',M.OLIGO);npm=inp('npm','sourcePCRmix;12cycle programme begins95°C30s; no extra72°Cprehold imposed',M.CATALYST)
 tagged=amp('rna_index_pcr','Native12cycleRNA-library PCR',tagged,rna_tagged('indexedRNA-librarydaughterDNA'),fp,rp,npm);tagged=recover('rna_final_cleanup','0.7×AMPureXP',tagged)
 lib=g.product('rna_library','native RNA-derived library;exactorderedendsunassigned',M.LIBRARY,state=tagged.state)
 g.react('retain_rna','RetainRNA library',K.GENERIC,((tagged,I.TEMPLATE),),(lib,),changes_molecule=False)
 g.sequence(lib,SequencingPlan(acquired_roles=('Read 1','Read 2'),read_lengths={'Read 1':150,'Read 2':150},limitation='SourcePE150;RNARead1maps transcript-derivedDNA,Read2parses spatialA/BandUMI wherepresent. Native finalorderedP5/P7, strand-specificprimerlanding and indexacquisition remainunassigned.'))
 reviewed=list(REVIEWED_SHARED)
 reviewed+=list(REVIEWED_CUT if CUTTAG else REVIEWED_ATAC)
 g.declare_reviewed_scope(reviewed,SOURCE)
 return g.require_publishable()

# Explicit source inventory, not generated from graph nodes. Each named recovery
# contributes its known selection reagent and physical operation.
RECOVERY_SHARED=('fix_wash','tag_wash','rt_wash','lig_a_wash','lig_b_wash','lysate_cleanup','dna_cleanup','dna_final_cleanup','bound_wash','switch_wash','rna_preamp_cleanup','rna_final_cleanup')
REVIEWED_SHARED=('tissue','fixative','fix','perm','permeabilize','transposome','tn5','activation','tagment','edta','stop_tn5','rt_primer','rt_enzyme','rt','ligase','barcode_a','lig_a','barcode_b','lig_b','proteinase','reverse_crosslink','streptavidin','separate','dna_p5','dna_p7','dna_pol','dna_initial','dna_initial_extension','dna_pcr','dna_library','retain_dna','tso','switch','rna_f','rna_r','kapa','rna_first_pcr','bead_removal_selection','remove_beads','rna_qpol','rna_more_pcr','rna_1ng_rule','sample_rna','xt','xt_tn5','rna_tagment','nt','rna_stop','rna_index_p5','rna_index_p7','npm','rna_index_pcr','rna_library','retain_rna')+tuple(x for key in RECOVERY_SHARED for x in (key+'_selection',key))
REVIEWED_CUT=('primary_antibody','bind_primary','secondary','bind_secondary','antibody_wash_selection','antibody_wash','tether_tn5','tn5_wash_selection','tn5_wash','cut_hold','dna_prehold')
REVIEWED_ATAC=('qpcr_rule','qpcr_sample','qpcr_p5','qpcr','additional_cycles','dna_additional')
