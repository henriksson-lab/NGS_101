"""Protocol-owned, source-qualified incomplete preparation graph."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Scene,Segment,feature
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,Terminal,Disposition,SequencingPlan,SourceItem
from pcr import PCRGeometry

def unknown(label):
 s=Scene();s.strand('preparation',[Segment('unassembled','['+label+']',placeholder=True)]);return s

def op(g,key,x,label,kind=S.MAGIC,reagent='Source reagents; complete formulation/geometry not assembled'):
 role=R.CATALYST if kind in {S.DIGESTION,S.DEPROTEINIZATION,S.EXTENSION,S.END_REPAIR,S.TAILING,S.NICK_SEALING} else R.STRUCTURAL if kind==S.HYBRIDIZATION else R.CONDITION
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

def finish(g,x,limitation):
 lib=g.product('library','Source recovered sequencing preparation; complete molecular state unassembled',K.LIBRARY);g.react('library_identity','Source recovered sequencing library',S.GENERIC,((x,R.TEMPLATE),),(lib,),changes_molecule=False);g.sequence(lib,SequencingPlan(limitation=limitation));g.add_source_item(SourceItem('pending','Full source-specific substrate/adapter/end/read geometry and historical kit originals','Selected original scope only; ordinary unresolved work is not confirmed source absence.'));return g

from illumina import P5,TRUSEQ_READ1
SOURCE='Published XYZeq main PMC8059935, Results136/Fig1A-B and Methods210–222/238; rendered original Fig1B p3. Selected tissue route; complete supplementary cocktail/template-switch/tagmentation geometry remains unassembled.'
def rt_primer():
 s=Scene();s.strand('well_rt',[Segment('R1 3prime partial body',TRUSEQ_READ1[-22:]),Segment('UMI','N'*10,placeholder=True,feature=feature('xyz_umi','umi','random')),Segment('spatial microwell barcode','N'*16,placeholder=True,feature=feature('xyz_spatial','spatial_barcode','unknown')),Segment('oligo dT','T'*18)]);return s

def pcr_index_primer():
 s=Scene();s.strand('indexing_primer',[Segment('P5',P5),Segment('PCR well barcode i5; exact value and molecular length unassembled','[i5]',placeholder=True,feature=feature('xyz_pcr_well','cell_barcode','unknown')),Segment('R1',TRUSEQ_READ1)]);return s

def protocol_graph():
 g=ProtocolGraph('XYZeq — source-qualified tissue route',authority=SOURCE,representation='published_instance')
 x=g.input('organ','Source harvested tissue organ with resident RNA; selected tissue experiment, not 10x comparator',K.BIOLOGICAL,primary=True)
 x=op(g,'prefreeze',x,'Source30min ice-cold DMSO-free freezing medium; then30min ice-cold DSP plus10%FCS',S.CROSS_LINKING,'Source DSP crosslinking reagent and FCS; complete molecular adduct states unassembled')
 x=op(g,'neutralize',x,'Source neutralize in ice-cold20mMTrisHClpH7.5; airtightcryomold,slowfreeze overnight−80C',S.QUENCHING,'Source Tris quench/freezing conditions; actual fixative disposition unassembled')
 x=op(g,'section',x,'Source25um tissue section placed on preloaded microwell array at−80C; orientation recorded and chamber gasket sealed',S.GENERIC,'Source preprintedspatial oligo array; distinct RT wells are physical partitions')
 rp=g.input('spatial_rt_primer','Source main214 oligo dT; UMI/spatial meaning verified Fig1B, no concrete spatial whitelist assigned',K.OLIGO,state=rt_primer())
 control=g.input('hek293t_control','Source DSP-fixed HEK293T cells,5ul at10x10^6cells/ml added to RT digestion mix for tissue experiments; actual per-well distribution variable',K.BIOLOGICAL)
 enz=g.input('rt1_cocktail','Source RT/digestion/permeabilization cocktail; complete enzyme/buffer identity not assembled from main',K.CATALYST)
 y=g.product('spatial_rt_product','Retained source cells/tissue/cDNA hybrids/free primers/cocktail and control cells; chemical products unassembled',K.MOLECULAR)
 g.react('spatial_rt','Gradual stepwise warming to42C then20min in chamber; ramps not invented',S.REVERSE_TRANSCRIPTION,((x,R.TEMPLATE),(rp,R.STRUCTURAL),(control,R.SPIKE_IN),(enz,R.CATALYST)),(y,))
 x=cleanup(g,'release',y,'Source chip in50ml1xSSC plus25%FCS,vortex,1000rcf10min,remove excess volume; collected cells retained, carryover not assumed zero')
 x=op(g,'filter_dapi',x,'Source recovered cells filtered and stainedDAPI before BDAria sorting',S.GENERIC,'Physical filtering/DAPI staining; DAPI is not a DNA-sequence barcode')
 sorted_cells,unsorted=split(g,'sort_to_plate',x,('Source sorted cells partitioned to96wellplates preloaded5ul secondRT mix','Uncollected sorting/filtering material; exact source recovery/fate unassigned'))
 g.mark_terminal(Terminal(unsorted.id,Disposition.UNRESOLVED,'Source does not establish exact sort recovery here. This is a physical fraction, not computational count filtering.'))
 rp2=g.input('second_rt_primer','Source second RT primer(s), exact identities/geometry unassembled',K.OLIGO);e2=g.input('second_rt_catalyst','Source second RT enzyme/mix; identity unassembled',K.CATALYST);x=g.product('second_rt_product','Second RT preparation retaining original sorted material/reagents; complete product unassembled',K.MOLECULAR);g.react('rt2','Source second RT42C1.5h in plate wells; full mix/lysis/template-switch chemistry unassembled',S.REVERSE_TRANSCRIPTION,((sorted_cells,R.TEMPLATE),(rp2,R.STRUCTURAL),(e2,R.CATALYST)),(x,))
 x=pcr(g,'well_index_pcr',x,'Source2xKapaHotStartReadyMix PCR with P5-i5-R1 indexing primer; companion primer/program unassembled',(('Main222 indexing primer; i5 identifies PCR well',pcr_index_primer()),('Source companion PCR primer not assembled',None)))
 pooled=g.product('pooled_pcr','Physical pool of source PCR well contents, free primers/reagents retained',K.MOLECULAR);g.react('pool','Source PCR plate contents physically pooled into2mlEppendorf tubes',S.GENERIC,((x,R.TEMPLATE),),(pooled,),reason='Known physical multiwell pool; input node currently aggregates well classes, so exact separate-input pooling relation is unassembled.',changes_molecule=False);x=pooled
 x=cleanup(g,'cdna_spri',x,'Source AMPureXP SPRI cDNA purification; exact ratio/elution not specified in selected main scope')
 enzyme=g.input('nextera_catalyst','Source tagmentation transposase, exact formulation/version unassembled',K.CATALYST);adapter=g.input('nextera_adapter','Source tagmentation adapter stock, exact loaded oligo/strand ends unassembled',K.ADAPTER)
 y=g.product('tag_product','Source tagmented cDNA/reagents; exact strand ends and retained molecules unassembled',K.MOLECULAR)
 g.react('tagmentation','Source cDNA tagmented before final Nextera-p7 index PCR',S.TAGMENTATION,((x,R.TEMPLATE),(enzyme,R.CATALYST),(adapter,R.STRUCTURAL)),(y,))
 x=pcr(g,'final_pcr',y,'Source amplified with Illumina Nextera library p7 index; full partner primer/cycles unknown in selected main scope',(('Source Nextera p7 index primer',None),('Source paired finalPCR primer',None)))
 x=assay(g,'final_qc',x,'Source BioAnalyzer and Qubit analytical aliquots; amounts unassigned')
 return finish(g,x,'Source NovaSeq/HiSeq4000 R1=26,R2=98,I1=8,I2=8. Fig1B cellular identifier combines spatial16 and PCRwelli5; UMI10 separate. No complete primer landing/run-cycle span validation or invented template-switch route. DecontX/UMI deduplication is computational, not physical removal from this preparation.')

def selftest():
 from checks import Check,Source,have
 from storage import sources_root
 c=Check();g=protocol_graph()
 c('IndexPCR precedes physical pool',list(g.steps).index('well_index_pcr')<list(g.steps).index('pool'),True)
 c('SecondRT remains separate from arrayRT',all(k in g.steps for k in ('spatial_rt','rt2')),True)
 c('SourcePCRwelli5 has cell meaning',pcr_index_primer().strands['indexing_primer'].segs[1].feature.role,'cell_barcode')
 c('Complete unknown chemistry prevents public promotion',bool(g.publication_problems()),True)
 src=Source(sources_root()/'xyzeq__10.1126+sciadv.abg4755'/'sciadv.abg4755_PMC8059935.xml.txt','Publishedmain214/222 andFig1B renderedp3')
 if have(c,src,label='Publishedmain'):
  t=''.join(src.path.read_text().split());c('WholeRT fixed pieces and source16bp placeholder retained',TRUSEQ_READ1[-22:]+'N'*10+'[16–basepairuniquespatialbarcode]'+'T'*18 in t,True)
  c('Wholeindexprimer source body matches importedcanonicalparts',P5+'[i5]'+TRUSEQ_READ1 in t,True)
  c('Selected read cycles source transcription','read1:26cycles;read2:98cycles;index1:8cycles;index2:8cycles'in t,True)
 c.report()
if __name__=='__main__':selftest()
