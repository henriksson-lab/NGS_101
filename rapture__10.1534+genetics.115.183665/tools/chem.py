"""Protocol-owned, source-qualified incomplete preparation graph."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Scene,Segment,revcomp,feature
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,Terminal,Disposition,SequencingPlan,SourceItem
from pcr import PCRGeometry
from adapter_ligation import AdapterLigationGeometry,reported_adapter_end

def unknown(label):
 s=Scene();s.strand('preparation',[Segment('unassembled','['+label+']',placeholder=True)]);return s

def op(g,key,x,label,kind=S.MAGIC,reagent='Source reagents; complete formulation/geometry not assembled'):
 role=R.CATALYST if kind in {S.DIGESTION,S.EXTENSION,S.END_REPAIR,S.TAILING,S.NICK_SEALING} else R.STRUCTURAL if kind==S.HYBRIDIZATION else R.CONDITION
 st=g.input(key+'_stock',reagent,K.CATALYST if role==R.CATALYST else K.OLIGO if role==R.STRUCTURAL else K.CONDITION);y=g.product(key+'_product',label+'; all substrates/products/companion species and reaction remnants retained until stated workup',K.MOLECULAR)
 kw={'reason':'Known source operation; exact molecular states, ends and geometry not yet assembled.'} if kind in {S.MAGIC,S.GENERIC} else {}
 g.react(key,label,kind,((x,R.TEMPLATE),(st,role)),(y,),**kw);return y

def cleanup(g,key,x,label):
 st=g.input(key+'_condition',label,K.CONDITION);y=g.product(key+'_kept','Source recovered preparation; recovery/carryover not presumed perfect',K.MOLECULAR);w=g.product(key+'_waste','Actual separated physical workup fractions',K.WASTE);g.react(key,label,S.CLEANUP,((x,R.TEMPLATE),(st,R.CONDITION)),(y,w));g.discard(w,'Source physical separation, not computational filtering.');return y

def split(g,key,x,labels):
 ys=[g.product(key+'_'+str(i),lab,K.MOLECULAR) for i,lab in enumerate(labels)];g.react(key,key,S.SPLIT,((x,R.TEMPLATE),),ys);return ys

def assay(g,key,x,label):
 y,a=split(g,key,x,('Retained preparation after source assay allocation; assay amount unassigned',label));g.mark_terminal(Terminal(a.id,Disposition.QC,'Separate analytical material; no implicit return to preparation.'));return y

def pcr(g,key,x,label,primers):
 ps=[g.input(key+'_p'+str(i),lab,K.OLIGO,state=state) for i,(lab,state) in enumerate(primers)];e=g.input(key+'_enzyme',label+' source catalyst',K.CATALYST);y=g.product(key+'_product','New DNA copies plus original DNA/free primers/reagents',K.MOLECULAR);geom=PCRGeometry.unavailable(unknown('Whole operational PCR geometry not assembled'),reason='Source primer identities do not by themselves establish all inward sites/products.',symmetric_species=len(ps)==1);g.react(key,label,S.PCR,((x,R.TEMPLATE),*((p,R.STRUCTURAL) for p in ps),(e,R.CATALYST)),(y,),geometry=geom);return y

def oligo_scene(name,seq,mod5='',mod3=''):
 s=Scene();s.strand(name,[Segment(name,seq)],mod5=mod5,mod3=mod3);return s

def finish(g,x,limitation):
 lib=g.product('library','Source recovered sequencing preparation; complete molecular state unassembled',K.LIBRARY);g.react('library_identity','Source recovered sequencing library',S.GENERIC,((x,R.TEMPLATE),),(lib,),changes_molecule=False);g.sequence(lib,SequencingPlan(limitation=limitation));g.add_source_item(SourceItem('pending','Full source-specific substrate/adapter/end/read geometry and historical kit originals','Selected original scope only; ordinary unresolved work is not confirmed source absence.'));return g

from illumina import P5,P7,TRUSEQ_READ1
SOURCE='Published Genetics2016 main PMC4788223, NewRADMethods270–278/readlayout280–281/Fig1 caption151–152 (primary browser original; standard local acquisition blocked). TraditionalRAD comparator not imported.'
POSTCAP_F=P5+TRUSEQ_READ1[4:21]
POSTCAP_R=P7[:21]
POSTCAP_F_WRITTEN=POSTCAP_F[:-1]+'*'+POSTCAP_F[-1]
POSTCAP_R_WRITTEN=POSTCAP_R[:-1]+'*'+POSTCAP_R[-1]
def ps_oligo(name,written):
 seq=written.replace('*','');s=Scene();s.strand(name,[Segment(name,seq,note='Source terminal phosphorothioate bond marked * in '+written)]);return s

def protocol_graph():
 g=ProtocolGraph('Rapture — source newRAD plus capture route',authority=SOURCE,representation='published_instance')
 x=g.input('gdna','Source50ng gDNA of one indexed well; tissue extraction and actual DNA coordinates unassembled',K.MOLECULAR,primary=True)
 x=op(g,'sbfi_digest',x,'Source SbfI-HF2.4U/12ul1xNEBuffer4,37C1h;80C20min heat inactivation',S.DIGESTION,'Source SbfI-HF and NEBuffer4; no source cut-coordinate prediction')
 adapter=g.input('biotin_well_adapter','Source2ul50nM indexed SbfI/PstI BIOTINYLATED RAD adapter;8bp fixed well barcode, full FileS1 strands/modifications unreviewed',K.ADAPTER,state=unknown('Full source biotinylated adapter/end geometry requires actual FileS1'))
 ligase=g.input('t4_ligase','Source2ul ligation mix1.28ulH2O/0.4ulNEBuffer4/0.16ul100mMrATP/0.16ulT4DNAligase',K.CATALYST);y=g.product('ligation_mixture','Source tagged/untagged genomic fragments plus free adapters/reagents; exact junctions/unconverted species unassembled',K.MOLECULAR);g.react('ligation','Source20C1h;65C15min ligase inactivation',S.LIGATION,((x,R.TEMPLATE),(adapter,R.STRUCTURAL),(ligase,R.CATALYST)),(y,),geometry=AdapterLigationGeometry.unavailable('ActualsourceFileS1 unread; no conventionalRAD or modernvendor adapter stock substituted.'))
 used,rest=split(g,'well_aliquot',y,('Source5ul selected well ligation to plate pool','Untransferred indexed well stock reserved for additional library preparation'));g.mark_terminal(Terminal(rest.id,Disposition.UNRESOLVED,'Explicit source reserved well aliquot; not discarded by read filtering.'))
 other=g.input('other_wells','Other95 source indexed well ligation aliquots5ul each; actual sample/barcode assignments unreviewed',K.MOLECULAR,primary=True);pool=g.product('plate_pool','Source96-well physical pool; independent sample DNA/index classes retained',K.MOLECULAR);g.react('plate_pooling','Source equal well aliquots pooled platewise',S.POOL,((used,R.TEMPLATE),(other,R.TEMPLATE)),(pool,))
 x=cleanup(g,'ampure_pool',pool,'Source1xAMPureXP purification;resuspend210ul lowTE')
 x=op(g,'shear',x,'Source BioruptorNGS9cycles30s on/90s off; additional source-dependent cycles conditional',S.FRAGMENTATION)
 x=assay(g,'fragment_analyzer',x,'Separate source fragment-analyzer QC aliquot; additional shearing decision not automatically predicted')
 x=cleanup(g,'rad_streptavidin',x,'Source30ul M280streptavidin beads washed2x100ul2xBW thenmixed equal2xBW/samplevol;RT20min;2x150ul1xBW,2x56CBW,2xNEBuffer4 washes; retain bound RAD-tag fraction')
 enz=g.input('release_sbfi','Source2ulSbfI-HF in40ul1xNEBuffer4,37C1h',K.CATALYST);released=g.product('released_rad','Source liberated barcoded DNA in recovered supernatant; exact remaining ends/adapter fractions unassembled',K.MOLECULAR);beads=g.product('release_beads','Physically retained beads, biotinylated adapter pieces and unreleased/carryover fractions',K.WASTE);g.react('release','Source SbfI processing and physical supernatant recovery',S.MAGIC,((x,R.TEMPLATE),(enz,R.CATALYST)),(released,beads),reason='Main explicitly source restriction release from beads; exact FileS1 base/end inventory and physical recovery imperfect/unassembled.');g.discard(beads,'Physical separation after restriction release, not presumed DNA hydrolysis of all unwanted fragments.')
 x=cleanup(g,'ampure_release',released,'Source1xAMPureXP purification;elute55.5ul lowTE')
 x=op(g,'historical_neb_library',x,'Source NEBNextUltra DNA LibraryPrep forIllumina with NOmodifications; plate index installation and PCR source historical kit recipe unassembled')
 other=g.input('other_plates','Independent source plate-indexed libraries; actual selected pool and pre-capture aliquots unassigned',K.MOLECULAR,primary=True);pool=g.product('plate_library_pool','Source physically combined plate libraries for capture',K.MOLECULAR);g.react('pool_libraries','Source pooled indexed plate libraries before one capture',S.POOL,((x,R.TEMPLATE),(other,R.TEMPLATE)),(pool,))
 pre,capture=split(g,'pre_post_capture',pool,('Independent source pre-capture comparison aliquot','Source library aliquot entering MYbaits capture'));g.mark_terminal(Terminal(pre.id,Disposition.UNRESOLVED,'Actual separate pre-capture sequencing arm; not returned to captured library.'))
 x=op(g,'mybaits_hybridization',capture,'Source MYcroarray500RAD-locus bait pool, supplied MYbaits protocol; selected recommended bait concentration case,one-fifth comparator separate',S.HYBRIDIZATION,'Source biotinylated oligonucleotide bait collection; actual500 FileS2 sequences unreviewed')
 x=cleanup(g,'mybaits_capture',x,'Source MYbaits target capture and physical wash/recovery; exact historical kit workup unassembled')
 x=pcr(g,'postcapture_pcr',x,'Source final capture-library amplification with universal source primers; program/kit stock unassigned in inspected main',(('Source universal forward PS primer',ps_oligo('postcapF',POSTCAP_F_WRITTEN)),('Source universal reverse PS primer',ps_oligo('postcapR',POSTCAP_R_WRITTEN))))
 return finish(g,x,'Source HiSeq2500PE100. Source8bp hamming well sample barcode plus separate plate index; exact full adapter/primer binding/indexcycle geometry unassembled. NewRAD wellbarcode/partialsite may occur on either read; do not impose traditionalRAD R1-only demultiplex. Random shear endpoints are not syntheticUMIs; read duplicate removal computational. Historical NEB andMYbaits source recipes deferred, not borrowed from current vendor kits.')

def selftest():
 from checks import Check,Source,have
 from storage import sources_root
 c=Check();g=protocol_graph();c('Source reserved well stock remains physical','well_aliquot_1' in g.terminals,True);c('Separate pre-capture sequencing arm retained','pre_post_capture_0' in g.terminals,True);c('Release is recovered from source bound fraction',g.steps['release'].inputs[0].material,'rad_streptavidin_kept');c('Historical library kit remains explicit gap',g.steps['historical_neb_library'].kind,S.MAGIC);c('Full read/end geometry not invented',bool(g.publication_problems()),True)
 src=Source(sources_root()/'rapture__10.1534+genetics.115.183665'/'genetics.115.183665_PMC4788223.xml.txt','PublishedmaincaptureMethods278:wholeuniversalprimers')
 if have(c,src,label='Locallycachedpublishedmain'):
  t=''.join(src.path.read_text().split());c('Whole final-capture primers and PS marks sourceexact',all(x in t for x in (POSTCAP_F_WRITTEN,POSTCAP_R_WRITTEN)),True)
 c.report()
if __name__=='__main__':selftest()
