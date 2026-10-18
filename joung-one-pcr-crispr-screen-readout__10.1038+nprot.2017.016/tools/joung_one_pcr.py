from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Construct,Scene,Segment,feature,revcomp,complement_segments
from protocol_model import role
from protocol_graph import ProtocolGraph,MaterialKind as M,InputRole as I,StepKind as K,SequencingPlan,Terminal,Disposition
from pcr import PCRGeometry,PCRPrimer
import illumina as il
import seqprimers as sp
import crispr as cr
TITLE='Joung one-PCR CRISPR screen readout'
NOTES='01_joung-one-pcr.html'
SOURCE='Joung et al., <a href="https://doi.org/10.1038/nprot.2017.016">NatureProtocols2017, published manuscriptPMC5526071 Steps32–33/37/58–64 andTable3</a>.'
SUMMARY='Ten forward-primer species occupy separate PCR reactions. Each reaction simultaneously enriches the integrated GeCKO guide and adds sequencing arms; purified reaction populations are then pooled.'
CAVEAT='Selected genomic-DNA GeCKO readout; SAM and the plasmid gel-excision branch are distinct. A non-G-starting20nt guide with a transcription-startG and deposited v2 scaffold is the representative model. Source-loaded sequencing-primer identity is unverified; the standardTruSeq sites are a compatible reference. Twelve cycles are the printed default; source discussion recommends sample-specific titration. Actual genomic flanks, copy numbers and off-product recovery remain unassigned.'
SEQUENCING_ENDING='Source MiSeq/NextSeq acquisition: Read1=80, Index1=8; 5%/20% independentPhiX control respectively. No Read2 is assigned. Ten different forward staggers change guide start cycles; the representative geometry below usesFwd1.'
SEQ_PRIMERS=(sp.TRUSEQ['R1'],sp.TRUSEQ['I1'])
GUIDE=feature('guide','guide_barcode','unknown')
I7=feature('sample_i7','sample_index','whitelist',whitelist='Joung2017Table3 KO reverse1–8')
FORWARD_PRIMERS=('AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCTTAAGTAGAGGCTTTATATATCTTGTGGAAAGGACGAAACACC', 'AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCTATCATGCTTAGCTTTATATATCTTGTGGAAAGGACGAAACACC', 'AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCTGATGCACATCTGCTTTATATATCTTGTGGAAAGGACGAAACACC', 'AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCTCGATTGCTCGACGCTTTATATATCTTGTGGAAAGGACGAAACACC', 'AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCTTCGATAGCAATTCGCTTTATATATCTTGTGGAAAGGACGAAACACC', 'AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCTATCGATAGTTGCTTGCTTTATATATCTTGTGGAAAGGACGAAACACC', 'AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCTGATCGATCCAGTTAGGCTTTATATATCTTGTGGAAAGGACGAAACACC', 'AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCTCGATCGATTTGAGCCTGCTTTATATATCTTGTGGAAAGGACGAAACACC', 'AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCTACGATCGATACACGATCGCTTTATATATCTTGTGGAAAGGACGAAACACC', 'AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCTTACGATCGATGGTCCAGAGCTTTATATATCTTGTGGAAAGGACGAAACACC')
REVERSE_PRIMERS=('CAAGCAGAAGACGGCATACGAGATTCGCCTTGGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTCCGACTCGGTGCCACTTTTTCAA', 'CAAGCAGAAGACGGCATACGAGATATAGCGTCGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTCCGACTCGGTGCCACTTTTTCAA', 'CAAGCAGAAGACGGCATACGAGATGAAGAAGTGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTCCGACTCGGTGCCACTTTTTCAA', 'CAAGCAGAAGACGGCATACGAGATATTCTAGGGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTCCGACTCGGTGCCACTTTTTCAA', 'CAAGCAGAAGACGGCATACGAGATCGTTACCAGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTCCGACTCGGTGCCACTTTTTCAA', 'CAAGCAGAAGACGGCATACGAGATGTCTGATGGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTCCGACTCGGTGCCACTTTTTCAA', 'CAAGCAGAAGACGGCATACGAGATTTACGCACGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTCCGACTCGGTGCCACTTTTTCAA', 'CAAGCAGAAGACGGCATACGAGATTTGAATAGGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTCCGACTCGGTGCCACTTTTTCAA')
FWD_SITE='GCTTTATATATCTTGTGGAAAGGACGAAACACC'
REVERSE_SITE='CCGACTCGGTGCCACTTTTTCAA'
# Deposited Addgene52961 map confirms this specific invariant scaffold component.
# Its empty-vector filler is not a synthesized GeCKO member or a screen genome.
SCAFFOLD=cr.SCAFFOLD_V1

def final_library(member=0):
 fwd=FORWARD_PRIMERS[member];r=REVERSE_PRIMERS[0]
 prefix=il.TRUSEQ_P5_FULL
 if not fwd.startswith(prefix) or not fwd.endswith(FWD_SITE):raise ValueError('source forward decomposition')
 tail=r[len(il.P7)+8:]
 if tail!=il.TRUSEQ_READ2+REVERSE_SITE:raise ValueError('source reverse decomposition')
 end=SCAFFOLD.index(revcomp(REVERSE_SITE))+len(REVERSE_SITE)
 # The source reverse primer binds scaffold positions50–72, leaving3 terminal scaffold bases outside PCR.
 return Construct([Segment('P5',il.P5),Segment('Read1 arm',prefix[len(il.P5):]),Segment('source stagger',fwd[len(prefix):-len(FWD_SITE)]),Segment('hU6 site',FWD_SITE),Segment('representative transcriptionG','G',inferred=True),Segment('guide','N'*20,placeholder=True,feature=GUIDE),Segment('scaffold through reverse site',SCAFFOLD[:end]),Segment('Read2 arm rc',revcomp(il.TRUSEQ_READ2)),Segment('i7 rc',revcomp(r[len(il.P7):len(il.P7)+8]),feature=I7),Segment('P7 rc',il.P7_RC)],name='representative source Fwd'+str(member+1)+'/KORev1 daughter')
FINAL_LIBRARY=final_library()

def genomic_state():
 sc=Scene();sc.strand('native',[role('native genomic DNA containing integrated GeCKO guide cassettes; no flowcell arms; actual genomic/template coordinates unassigned')]);sc.strand('other',[role('possible co-extracted RNA and other DNA; no purity assumed')]);return sc

def oligos():
 for n,x in enumerate(FORWARD_PRIMERS,1):yield f'NGS-Lib-Fwd-{n}:5′-{x}-3′'
 yield 'Representative NGS-Lib-KO-Rev-1:5′-'+REVERSE_PRIMERS[0]+'-3′; source eight-index family'

def protocol_graph():
 g=ProtocolGraph(TITLE,authority=SOURCE)
 gdna=g.input('genomic_dna','screen genomicDNA; all sourcegDNA scaled across reactions',M.BIOLOGICAL,state=genomic_state(),primary=True)
 alloc=g.input('allocation_rule','one reaction per each of10 Fwd species; ≤2.5µg gDNA/50µl; scale reaction count to totalgDNA',M.CONDITION)
 aliquots=[g.product('aliquot_'+str(n+1),'genomicDNA allocated toFwd'+str(n+1)+' reaction class',M.MOLECULAR,state=gdna.state) for n in range(10)]
 g.react('allocate','Allocate sourcegDNA across all10 forward-primer reaction classes',K.SPLIT,((gdna,I.TEMPLATE),(alloc,I.CONDITION)),tuple(aliquots),changes_molecule=False)
 rev=g.input('reverse','selected KORev1 reference member; actual sample member unassigned',M.OLIGO,state=Construct([Segment('KORev1',REVERSE_PRIMERS[0])]))
 pol=g.input('polymerase','NEBNextHighFidelity2×PCRMasterMix',M.CATALYST)
 column=g.input('columns','ZymoSpinV withreservoir/DNABindingBuffer andwashbuffer',M.SELECTION)
 pcrs=[]
 for n,aliquot in enumerate(aliquots):
  tag=str(n+1);f=g.input('forward_'+tag,'NGS-Lib-Fwd-'+tag,M.OLIGO,state=Construct([Segment('Fwd'+tag,FORWARD_PRIMERS[n])]))
  out=g.product('pcr_'+tag+'_out','Fwd'+tag+' PCR daughters and possible co-carried originals',M.MOLECULAR,state=final_library(n))
  original=g.product('pcr_'+tag+'_originals','native templates and unselected molecular classes',M.MOLECULAR,state=aliquot.state)
  g.react('pcr_'+tag,'INFERRED — sourcePCR:95°C5min;12cycles98°C20s/60°C15s/72°C15s;72°C1min',K.PCR,((aliquot,I.TEMPLATE),(f,I.STRUCTURAL),(rev,I.STRUCTURAL),(pol,I.CATALYST)),(out,original),geometry=PCRGeometry.unavailable(aliquot.state,reason='Actual native genomicDNA is the input; a completed adapter-bearing daughter is not its template. Exact genomic boundaries and initial binding are unassigned.'))
  g.mark_terminal(Terminal(original.id,Disposition.UNRESOLVED,'PCR does not remove original genomes/RNA/other DNA; their subsequent class-specific recovery is unassigned.'))
  clean=g.product('clean_'+tag,'Fwd'+tag+' sourcecolumn retained population',M.MOLECULAR,state=out.state)
  waste=g.product('column_waste_'+tag,'physically removed flowthrough/washes',M.WASTE)
  g.react('purify_'+tag,'Source5volbindingbuffer;two2ml washes;dryspin;150µl elution',K.CLEANUP,((out,I.TEMPLATE),(column,I.SELECTION)),(clean,waste),changes_molecule=False)
  g.discard(waste,'Sourcecolumn flowthrough andwashes; doesnot prove complete class-specific depletion.')
  pcrs.append(clean)
 mixed=Scene()
 for n in range(10):mixed.strand('Fwd'+str(n+1),list(final_library(n)),label='independent daughterpopulation Fwd'+str(n+1))
 mixed.strand('possible_carried',[role('possible residual originalDNA/RNA/oligos; column recovery unassigned')])
 pool=g.product('library','pooled purified source forward-stagger populations',M.LIBRARY,state=mixed)
 g.react('pool','Pool purified PCRs ofone sample (Step64)',K.POOL,tuple((x,I.TEMPLATE) for x in pcrs),(pool,),changes_molecule=False)
 g.sequence(pool,SequencingPlan(library=FINAL_LIBRARY,primers=SEQ_PRIMERS,acquired_roles=('Read 1','Index 1 (i7)'),read_lengths={'Read 1':80,'Index 1 (i7)':8},limitation='Fwd1/KORev1 source representative; allten independentFwd populations remain in the pool. Compatible standardTruSeq primers are a reference, not independently verified loaded oligos. SourcePhiX control remains an independent sequencing mixture.'))
 reviewed=['genomic_dna','allocation_rule','allocate','reverse','polymerase','columns','pool','library']
 reviewed += [key+str(n) for n in range(1,11) for key in ('forward_','pcr_','purify_')]
 g.declare_reviewed_scope(reviewed,SOURCE)
 return g.require_publishable()

FINAL_CAPTION='INFERRED — Fwd1/KORev1 representative daughter with non-G-starting20nt guide and compatible standardTruSeq readout; other nine staggers remain independent populations.'

SEQUENCING_INTRO='Compatible standardTruSeq primer geometry for the representative Fwd1/KORev1 daughter; actual source acquisition80+8.'
