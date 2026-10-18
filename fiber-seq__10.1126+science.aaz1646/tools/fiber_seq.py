"""Protocol-owned native molecules and source-bounded graph."""
from chemdraw import Scene
from protocol_model import role
from protocol_graph import ProtocolGraph,MaterialKind as M,InputRole as I,StepKind as K,SequencingPlan,Terminal,Disposition
FINAL_LIBRARY=None
SEQ_PRIMERS=()
SEQUENCING_UNAVAILABLE=''
SHOW_READ_CONTENT=False
import seqprimers as sp
def population(*rows):
 sc=Scene()
 for key,parts,label in rows:sc.strand(key,tuple(parts),label=label)
 return sc
def recover(g,key,x,label):
 reagent=g.input(key+'_selection',label,M.SELECTION)
 kept=g.product(key+'_retained','recovered '+x.label,M.MOLECULAR,state=x.state)
 waste=g.product(key+'_waste','physical non-retained supernatant/washes',M.WASTE)
 g.react(key,label,K.CLEANUP,((x,I.TEMPLATE),(reagent,I.SELECTION)),(kept,waste),changes_molecule=False)
 g.discard(waste,'Actual non-retained fraction; complete species-specific depletion is not asserted.')
 return kept
def unresolved(g,x,reason):g.mark_terminal(Terminal(x.id,Disposition.UNRESOLVED,reason))
def op(g,key,x,label,kind,state,reagent,port=I.CATALYST):
 out=g.product(key+'_out',label,M.MOLECULAR,state=state)
 g.react(key,label,kind,((x,I.TEMPLATE),(reagent,port)),(out,),**({'changes_molecule':False} if kind in (K.GENERIC,K.QUENCHING) else {}))
 return out
TITLE='Fiber-seq — original2020 native chromatin marking'
NOTES='01_fiber-seq.html'
SOURCE='Stergachis et al. Science2020 doi:10.1126/science.aaz1646; published coauthor-hosted article and ScienceSupplementaryMaterials with printedDOI/title, AdMTase-seq/Fiber-seq methods.'
SUMMARY='Original Hia5 modifies accessible adenines in native chromatin. Two source extraction/precipitation rounds and AMPure cleanup precede mechanical shearing; marked original genomic strands remain distinct from any sequencing daughter.'
CAVEAT='Selected original200U/10min/37C route, not later livingbench25C/commercialHia5/SMRTbellkit3.0. Original2020 supplementaryMethods names the sequencing system but no preparationkit; complete kit-specific repair, hairpin/end geometry and native closure remain a bounded unread-source assignment. No fixed methylation pattern, copying step,18base hairpin or A/T junction is invented.'
LIMIT='Native2020 kit-specific hairpin/primer/closure geometry remains unassigned. PacBio reads original modified genomicDNA; methylation calls and computational read filters are distinct from physical purification. No fixed acquired readlength is assigned.'
SEQUENCING_ENDING=LIMIT
READ_SPECS=(sp.read_spec('Sequencing output',None,'Original nativeDNA and kinetic methylation information; SequelII System2.0 selected. SourceSystem1.0/Sequel/control alternatives are separate experiments.','include'),)
def fiber(marked=False):
 label='original native genomicDNA; Hia5-associatedm6dA positions/stochastic efficiency unassigned' if marked else 'original native genomicDNA; native marks unassigned'
 return population(('strandA',(role(label+' strandA'),),'original genomic strandA'),('strandB',(role(label+' strandB'),),'opposite original genomic strandB; independent methylation pattern'),('other',(role('chromatin proteins/RNA or their processing products; physical partition unassigned'),),'independent possible co-carried constituents'))
def protocol_graph():
 g=ProtocolGraph(TITLE,authority=SOURCE)
 nuclei=g.input('nuclei','source2.5millionS2 or300000K562 nuclei; selectedS2 route',M.BIOLOGICAL,state=fiber(),primary=True)
 hia5=g.input('hia5','original in-houseHia5,200U;37C10min',M.CATALYST)
 sam=g.input('sam','source0.8mMSAM/BufferA',M.CONDITION)
 marked=g.product('marked','original Hia5-marked chromatinDNA',M.MOLECULAR,state=fiber(True))
 g.react('mark','Deposit m6dA at accessible native adenines; positions unassigned',K.CHEMICAL_CONVERSION,((nuclei,I.TEMPLATE),(hia5,I.CATALYST),(sam,I.CONDITION)),(marked,))
 sds=g.input('sds','sourceSDSstop1% thenadditionalSDS/BufferA',M.CONDITION)
 x=op(g,'stop',marked,'SDS quench original marking reaction',K.QUENCHING,marked.state,sds,I.CONDITION)
 enzymes=g.input('extract_enzymes','RNaseA+ProteinaseK;50C1h',M.CATALYST)
 x=op(g,'deproteinize',x,'Digest proteins/RNA; hydrolysis is not physical removal',K.DEPROTEINIZATION,fiber(True),enzymes)
 x=recover(g,'phenol1',x,'First1:1phenol/chloroform/isoamyl extraction; physically transferupper aqueousphase')
 x=recover(g,'precip1',x,'FirstNaOAc/2.5volEtOH overnight−20C;pellet/70%EtOHwash;resuspend200µlBufferA')
 rnase=g.input('rnase2','secondRNaseA37C1h',M.CATALYST)
 x=op(g,'rna_digest2',x,'Second RNaseA digestion; hydrolysates remain until separation',K.DIGESTION,fiber(True),rnase)
 x=recover(g,'phenol2',x,'Second1:1phenol/chloroform/isoamyl extraction;upperaqueous transfer')
 x=recover(g,'chloroform',x,'Separatechloroform/isoamyl extraction;upperaqueous transfer removesresidualphenol')
 x=recover(g,'precip2',x,'SecondNaOAc/2.5volEtOH−20C1h;pellet/70%EtOHwash;resuspend100µlwater')
 x=recover(g,'ampure',x,'SourceAMPure repurification; exact ratio/species partition unassigned')
 shear=g.input('shear_conditions','Megaruptor10–15kb target;FragmentAnalyzer is analyticalQC',M.CONDITION)
 x=op(g,'shear',x,'Mechanically shear original markedDNA to source10–15kb target',K.FRAGMENTATION,fiber(True),shear,I.CONDITION)
 kit=g.input('native_kit','original2020 PacBio preparation referenced withoutkitname',M.CONDITION)
 lib=g.product('library','native original markedDNA PacBiolibrary; exact topology/ends unassigned',M.LIBRARY,state=fiber(True),inferred=True)
 g.react('native_library_boundary','Original2020 PacBiolibrary construction; exactkit unavailable in inspectedscope',K.MAGIC,((x,I.TEMPLATE),(kit,I.CONDITION)),(lib,),reason='ScienceSM names SequelII System2.0, not a preparationkit. Laterkit3.0 cannot establish original2020 repair/hairpins/nuclease order. Original modifiedDNA population retained without fabricated canonical kit steps.')
 g.sequence(lib,SequencingPlan(limitation=LIMIT))
 g.declare_reviewed_scope(REVIEWED,SOURCE)
 return g.require_publishable()
REVIEWED=('nuclei','hia5','sam','mark','sds','stop','extract_enzymes','deproteinize','phenol1_selection','phenol1','precip1_selection','precip1','rnase2','rna_digest2','phenol2_selection','phenol2','chloroform_selection','chloroform','precip2_selection','precip2','ampure_selection','ampure','shear_conditions','shear','native_kit','native_library_boundary','library')
