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
 role=R.CATALYST if kind in {S.DIGESTION,S.DEPROTEINIZATION,S.TRANSCRIPTION,S.EXTENSION,S.END_REPAIR,S.TAILING,S.NICK_SEALING} else R.STRUCTURAL if kind==S.HYBRIDIZATION else R.CONDITION
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

from dataclasses import dataclass
from illumina import P5,P7,TRUSEQ_SMALL_RNA_RA5,TRUSEQ_SMALL_RNA_RA3
SOURCE='Published PMC8292322.1 main56–66/136–176 and published MOESM1 SI p7 whole oligo table (rendered). PubMed identity counts as original1; three-original cap closed. Selected eight-NPOM-cage RT primer, not all optimization controls.'
CAGED_WRITTEN='GCCGGTAATACGACTCACTATAGGGtttGAGttCtACAGTCCGACGATCNNNNNNTCGAAGttTTTTTTTTTTTTTTTTTTTTTTV'
RANDOM_HEX='GCCTTGGCACCCGAGAATTCCA'+'N'*6
PCR1=P5+'GTTCAGAGTTCTACAGTCCGA'
PCR2=P7+'CGTGAT'+'GTGACTGGAGTTCCTTGGCACCCGAGAATTCCA'
@dataclass(frozen=True)
class SourceCagedODN:
 written:str=CAGED_WRITTEN
 @property
 def cage_positions(self):return tuple(i for i,b in enumerate(self.written) if b=='t')
 def scene(self):
  # Lowercase t is source NPOM-dT, never silently erased from the reagent identity.
  s=Scene();prefix,tail=self.written.split('NNNNNN');segs=[]
  for i,b in enumerate(prefix):segs.append(Segment('NPOM-caged dT' if b=='t' else 'source backbone',b.upper(),note='Source NPOM-caged dT at this position; pairing/read-through not predicted' if b=='t' else ''))
  segs.append(Segment('source random molecular identifier','N'*6,placeholder=True,feature=feature('pic_umi','umi','random')))
  segs.append(Segment('source selected inline sample barcode',tail[:6],feature=feature('pic_rt_sample','sample_index','fixed')))
  for b in tail[6:]:segs.append(Segment('NPOM-caged dT' if b=='t' else 'source oligo-dT / anchor',b.upper(),placeholder=b=='V',note='Source NPOM-caged dT at this position' if b=='t' else ''))
  s.strand('selected eight-cage RT ODN',segs);return s

def pcr_reverse():
 s=Scene();s.strand('RNA PCR primer2',[Segment('source P7',P7),Segment('source sequencer-index window','CGTGAT',feature=feature('pic_pcr_sample','sample_index','fixed')),Segment('source remainder','GTGACTGGAGTTCCTTGGCACCCGAGAATTCCA')]);return s

def rt(g,key,x,label,primer):
 pr=g.input(key+'_primer',label+' source RT primer',K.OLIGO,state=primer);e=g.input(key+'_enzyme',label+' SuperscriptII/dNTP/DTT/RNaseOut/buffer stock',K.CATALYST);y=g.product(key+'_product','Source RNA/cDNA hybrids plus original RNA, unextended primers, enzyme/reagents and biological companions',K.MOLECULAR);g.react(key,label,S.REVERSE_TRANSCRIPTION,((x,R.TEMPLATE),(pr,R.STRUCTURAL),(e,R.CATALYST)),(y,));return y

def allocated(g,key,x,label):
 used,rest=split(g,key,x,(label,'Untransferred source stock; actual remainder/fate unassigned'));g.mark_terminal(Terminal(rest.id,Disposition.UNRESOLVED,'Exact source aliquot is not assumed to consume entire stock.'));return used

def protocol_graph():
 g=ProtocolGraph('PIC — selected eight-cage transcript preparation',authority=SOURCE,representation='published_instance')
 x=g.input('specimen','Source fixed-position cultured cells or tissue specimen containing RNA; actual chosen sample/ROI unassigned',K.BIOLOGICAL,primary=True)
 x=op(g,'fix_permeabilize',x,'Source3.7%formaldehyde10minRT;5%Triton3min;0.1NHCl5min;1MTrispH8 neutralization10min')
 x=cleanup(g,'rt_rinse',x,'Source PBS-DEPC washes/65C5min specimen treatment and cold rinse before RT')
 x=rt(g,'in_situ_rt',x,'Source primer mix65C then4C before SSII; cagedODN250ng,dNTP;1xFirstStrand/DTT/RNaseOut/SSII100U;42C60min,70CPBS10min',SourceCagedODN().scene())
 x=op(g,'image_and_uncage',x,'Source optional immunostaining/imaging guides ROI;340–380nm15min UV uncages NPOM-dT; nonirradiated molecules/background remain physically in whole specimen')
 x=op(g,'lysis',x,'Source0.1%Tween/400ugml ProteinaseK PBS55C30min',S.DEPROTEINIZATION,'Source ProteinaseK and lysis buffer; DNA/RNA are not declared digested by proteolysis')
 other=g.input('other_barcode_lysates','Actual independently prepared source replicates carrying distinct RT sample barcodes; identities unassigned',K.MOLECULAR,primary=True);pooled=g.product('pooled_lysates','Source physically combined distinct-barcode lysates and companions',K.MOLECULAR);g.react('pool_replicates','Source pooling after lysis, before hybrid purification',S.POOL,((x,R.TEMPLATE),(other,R.TEMPLATE)),(pooled,));x=cleanup(g,'minelute',pooled,'Source QiagenMinElute PCR purification of cDNA:mRNA hybrids; H2Oelution')
 x=allocated(g,'second_strand_allocation',x,'Source15ul recovered hybrids for second-strand reaction')
 x=op(g,'second_strand',x,'Source E.coli DNApolymerase,DNAligase,RNaseH,dNTP/FirstStrand+SecondStrandbuffers;16C2h; uncaged-dependent promoter copying, not perfect gate prediction',S.EXTENSION,'Source E.coli DNApolymerase3U,DNAligase0.8U,RNaseH0.16U and buffers/dNTPs')
 x=cleanup(g,'ampure_ds',x,'Source AMPureXP purification of double-stranded cDNA; H2Oelution')
 x=allocated(g,'ivt_allocation',x,'Source6.4ul ds-cDNA eluate into T7IVT')
 x=op(g,'ivt',x,'Source MEGAscriptT7:1.6ul each A/G/C/UTP plus1.6ul10xbuffer+1.6ulenzyme,37C17h; aRNA copies plus DNA/reagents retained',S.TRANSCRIPTION,'Source T7 RNA polymerase and NTPs; exact transcribed molecular endpoints unassembled')
 x=op(g,'turbo_dnase',x,'Source1ulTURBODNase2U/ul,37C15min; DNA hydrolysis products remain until physical workup',S.DIGESTION,'Source TURBO DNase')
 x=op(g,'exosap',x,'Source6ulExoSAP-IT37C15min; exact historical reagent/end products unassembled',S.DIGESTION,'Source ExoSAP-IT PCR Product Cleanup Reagent')
 x=op(g,'fragment',x,'Source5.5ul200mMTrisacetatepH8.1/500mMKOAc/150mMMgOAc stock,94C3min;2.75ul0.5MEDTAquench',S.FRAGMENTATION)
 x=cleanup(g,'rnaclean',x,'Source RNACleanXP purification of fragmented aRNAs; H2Oelution')
 x=allocated(g,'second_rt_allocation',x,'Source5ul aRNA eluate for randomhex RT')
 x=rt(g,'randomhex_rt',x,'Source randomhexRT250ng/dNTP65C5min then4C before SSII;25C10min/42C1h',oligo_scene('source randomhexRT',RANDOM_HEX))
 x=allocated(g,'pcr_allocation',x,'Source9ul RT product enters final PCR')
 x=pcr(g,'library_pcr',x,'Source RNA PCR1/2 primers1.8ul each10uM/22.5ulPhusionHFmastermix/9.9ulH2O;98C30s;11x98C10s/60C30s/72C30s;72C10min',(('Source RNA PCR1',oligo_scene('RNA PCR1',PCR1)),('Source RNA PCR2',pcr_reverse())))
 x=cleanup(g,'library_size',x,'Source AMPureXP purification/size selection250–500bp; no perfect size threshold asserted')
 x=assay(g,'bioanalyzer_qc',x,'Separate source Bioanalyzer2100 high-sensitivityDNA QC aliquot')
 x=assay(g,'quant_qc',x,'Separate source Clontech library-quantification aliquot')
 return finish(g,x,'Source IlluminaHiSeq1500; main data analysis extracts R1 pattern NNNNNNCCCCCC (6nt molecular UMI then6nt RT sample barcode), R2 gene mapping. PCR2 has sourceCGTGAT sequencer-index window; exact run/index-cycle settings and molecular landing geometry unassembled. RNA PCR/IVT products not native-RNA sequences. ROI coordinates are image selections, not assumed encoded spatial grid. Duplicate/UMI filtering computational, not physical DNA removal.')

def selftest():
 from checks import Check,Source,have
 from storage import sources_root
 c=Check();g=protocol_graph();c('NPOM source identity retains eight positions',len(SourceCagedODN().cage_positions),8);c('All four source aliquots retain physical stock',all(k+'_1' in g.terminals for k in ('second_strand_allocation','ivt_allocation','second_rt_allocation','pcr_allocation')),True);c('Proteolysis is deproteinization',g.steps['lysis'].kind,S.DEPROTEINIZATION);c('No complete source-molecular/read claim',bool(g.publication_problems()),True)
 src=Source(sources_root()/'pic__10.1038+s41467-021-24691-8'/'s41467-021-24691-8_supp_41467_2021_24691_MOESM1_ESM.pdf.txt','PublishedSIp7 whole oligo table rendered')
 if have(c,src,label='PublishedSI'):
  t=''.join(src.path.read_text().split());c('Whole selected eight-cage modification pattern exact',CAGED_WRITTEN in t,True);c('Whole randomhex and both final PCR primers exact',all(x in t for x in (RANDOM_HEX,PCR1,PCR2)),True)
 c.report()
if __name__=='__main__':selftest()
