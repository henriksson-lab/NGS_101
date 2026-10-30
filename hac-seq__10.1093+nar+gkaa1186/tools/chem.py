"""Protocol-owned, source-qualified incomplete preparation graph."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Scene,Segment
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,Terminal,Disposition,SequencingPlan,SourceItem
from pcr import PCRGeometry

def unknown(label):
 s=Scene();s.strand('preparation',[Segment('unassembled','['+label+']',placeholder=True)]);return s

def op(g,key,x,label,kind=S.MAGIC,reagent='Source reagents; complete formulation/geometry not assembled'):
 role=R.CATALYST if kind in {S.DIGESTION,S.DEPROTEINIZATION,S.DEPHOSPHORYLATION,S.PHOSPHORYLATION,S.EXTENSION,S.END_REPAIR,S.TAILING,S.NICK_SEALING} else R.STRUCTURAL if kind==S.HYBRIDIZATION else R.CONDITION
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

SOURCE='Published HAC-seq main PMC7969016 Methods65/77/81/85,Results95/103/Fig2caption107 and originalFig2p6. Source three groups not collapsed. WorkbookTableS1 opening is cleavage-ratio data, not an operational primer table; SI DOCXcaption only. Whole kit/end/read geometry unassembled.'
GROUPS=('Ctrl','HAC','DM-HAC')
def demethylate(g,key,x,label):
 e=g.input(key+'_alkb','Source160pmolAlkB-WT+200pmolAlkB-D135S andRNasin; enzymeproduction/whole molecularcofactorstates unassembled',K.CATALYST)
 conditions=g.input(key+'_buffer','Source100ulreaction,50ul2xbuffer stock:0.6MKCl/4mMMgCl2/0.1mM(NH4)2Fe(SO4)2.6H2O/0.6mM2ketoglutarate/4mMascorbicacid/0.1mgmlBSA/0.1MMESpH5.0;stocknotfinal concentrations',K.CONDITION)
 y=g.product(key+'_product','Source treated RNA with variable demethylation/unconverted modifications,originalsubstrates/enzyme/buffer/byproducts; no perfect m3C removal predicted',K.MOLECULAR)
 g.react(key,label+';2hRTdark,5mMEDTAquench',S.CHEMICAL_CONVERSION,((x,R.TEMPLATE),(e,R.CATALYST),(conditions,R.CONDITION)),(y,))
 return cleanup(g,key+'_phenol',y,'Source acid phenol/chloroform RNA purification after demethylation; no assumed loss-free recovery')

def graph_for_group(group):
 if group not in GROUPS:raise ValueError('Explicit source Ctrl/HAC/DM-HAC group required')
 g=ProtocolGraph('HAC-seq — source '+group+' library arm',authority=SOURCE,representation='published_instance')
 x=g.input('ribominus_rna','Source independently prepared MCF7 RiboMinus totalRNA stock for this biological replicate; upstreamTRIzol/RiboMinus kit10ugtotalRNAinput and recovery not molecularlyassembled. Not one stock silently shared by all three group graphs',K.MOLECULAR,primary=True)
 if group in {'Ctrl','DM-HAC'}:
  x=demethylate(g,'pre_hac_alkb',x,'Source pre-HAC AlkBcontrol;Methods77 totalRNA scale5ug, actual recovered stock/transfer relation unassigned')
 used,rest=split(g,'fragment_allocation',x,('Source2ugRiboMinusRNA with/withoutpriorDM enters fragmentation','Untransferred RiboMinus/pre-DM RNA stock; recovered amount/fate unassigned'));g.mark_terminal(Terminal(rest.id,Disposition.UNRESOLVED,'No perfect recovery or automatic whole-stock consumption.'))
 x=op(g,'fragment',used,'Source10ulAM8740RNAfragmentation reagents70C10min;20mMEDTApH8quench; source~200nt size only nominal',S.FRAGMENTATION,'Source chemicalRNAfragmentation conditions; complete initialend inventoryunassembled')
 x=cleanup(g,'fragment_ethanol',x,'Source ethanol precipitation afterfragmentation')
 x=op(g,'dephosphorylate',x,'Source Antarcticphosphatase followingmanufacturer protocol; before HAC',S.DEPHOSPHORYLATION,'SourceAntarcticphosphatase,fullconditionsunassembled')
 x=op(g,'rephosphorylate',x,'Source T4PNK followingmanufacturer protocol; sourceend-repaired RNA5P/3OH,exactsequence/endpopulationsunassembled',S.PHOSPHORYLATION,'SourceT4polynucleotidekinase; before chemicalcleavage,not after')
 x=cleanup(g,'end_repair_clean',x,'SourceZymoRNA-CleanandConcentrator5 afterinitialendprocessing')
 if group in {'HAC','DM-HAC'}:
  x=op(g,'hydrazine',x,'Source25ulicecold10%hydrazine/3MNaCl,onice4h; selectivechemistryreported,notcompleteconversion/absolute specificity imposed',S.CHEMICAL_CONVERSION,'Sourcehydrazine/highsalt conditions')
  x=cleanup(g,'hydrazine_ethanol',x,'Sourceethanolprecipitationstops hydrazinereaction;RNA pellet retained')
  x=op(g,'aniline',x,'Source100ulH2O:glacialacid:aniline7:3:1,RTdark2h;RNA classes includeuncleaved plus3prime5P/3OHfragments and5primefragmentswithdamaged3ends; exactblockedterminalchemistry notassigned',S.CHEMICAL_CONVERSION,'Sourceanilinecleavageconditions; component classes remain together in same tube')
  x=cleanup(g,'aniline_ethanol',x,'Sourceethanolprecipitationafteraniline,redissolveH2O; no perfect removal of blocked5primefragments')
 if group=='HAC':x=demethylate(g,'post_hac_alkb',x,'SourceHAC namedgroup post-cleavageAlkB;Fig2A onlyHAC hasthisstage')
 used,rest=split(g,'library_allocation',x,('Source100ngRNAfragments tolibrarypreparation','Untransferred treatedRNAstock; amount/fateunassigned'));g.mark_terminal(Terminal(rest.id,Disposition.UNRESOLVED,'100ng isactual libraryallocation,notentire2ugfragmentationinput.'))
 stocks=[g.input(key,label,kind) for key,label,kind in (
 ('kit_3_adapter','Source3prime adapter stock; literal body/modifications/ligation ends unassembled',K.ADAPTER),
 ('kit_5_adapter','Source5prime adapter stock; literal body/modifications/ligation ends unassembled',K.ADAPTER),
 ('kit_rt_primer','Source kit RT primer, literal identity/site unassembled',K.OLIGO),
 ('kit_pcr_forward','Source kit PCR forward primer, literal identity/index unassembled',K.OLIGO),
 ('kit_pcr_reverse','Source kit PCR reverse primer, literal identity/index unassembled',K.OLIGO))]
 enzymes=g.input('kit_enzymes','Historical kit ligation/RT/PCR enzyme inventory; stage-specific additions and workups unassembled, not one simultaneous mixture',K.CATALYST)
 x=g.product('kit_preparation','Source recovered cDNA library preparation; kit molecular populations, blocked-RNA/off-product physical fates and exact complete state unassembled',K.MOLECULAR)
 g.react('historical_library','SourceNEBNextSmallRNALibraryPrepSet adapterligation/RT/PCR compound boundary; sourceRT60C1h',S.MAGIC,((used,R.TEMPLATE),*((a,R.STRUCTURAL)for a in stocks),(enzymes,R.CATALYST)),(x,),reason='Historical kitrevision,modifiedoligos,stagegeometry,cycles/workups notassembled; compoundinventory is not a claim of simultaneous addition or perfect blocked-fragment removal.')
 return finish(g,x,'SourceIlluminaNextSeq500;exactrunlength/indexcycles/wholeoligosandprimerlandinggeometryunassembled. Authors report3primecleaved/full-lengthRNAendscompatiblewithlibraryligationand5primecleaved3endblocked. Failedligationdoesnotphysicallydestroyall5primefragments; sourcekit workupandactualoffproducts notmolecularlyassembled. RNAcleavage ratio/sitecalling iscomputational,notliteralRNAremoval. AlkBnotassumedperfect; sourceAlkB-resistanttRNA sites andfalsepositiveslimitquantitation.')
def protocol_graph():return graph_for_group('HAC')
def selftest():
 from checks import Check,Source,have
 from storage import sources_root
 c=Check();graphs={k:graph_for_group(k)for k in GROUPS};h=graphs['HAC'];ctrl=graphs['Ctrl'];dm=graphs['DM-HAC']
 c('SourceHAC post-DM followschemicalcleavage',list(h.steps).index('aniline')<list(h.steps).index('post_hac_alkb'),True)
 c('Sourcecontrolpre-DM precedesfragmentation',all(list(graphs[k].steps).index('pre_hac_alkb')<list(graphs[k].steps).index('fragment')for k in ('Ctrl','DM-HAC')),True)
 c('CtrlhasnoHAC andDM-HAC hasnoinventedsecondAlkB',('hydrazine'in ctrl.steps,'post_hac_alkb'in dm.steps),(False,False))
 c('Endprocessing precedesnewcleavage ends',list(h.steps).index('rephosphorylate')<list(h.steps).index('hydrazine'),True)
 c('BothactualRNAstocksremainphysical',all(k in h.terminals for k in ('fragment_allocation_1','library_allocation_1')),True)
 c('Exactlibrarygeometry remains unassembled',all(bool(g.publication_problems())for g in graphs.values()),True)
 src=Source(sources_root()/'hac-seq__10.1093+nar+gkaa1186'/'nar_gkaa1186_PMC7969016.xml.txt','PublishedMethods85/Fig2A p6/Results103')
 if have(c,src,label='Publishedmain'):
  t=''.join(src.path.read_text().split());c('Source distinguishesblockedandcompatiblecleavageends','5′cleavedfragmentshavea3′damagedend'in t and '3′cleavedfragmentscontainthecorrect5′phosphateand3′OHgroups'in t,True)
  c('SourceRT60C1h islocalkitmodification','RNAswereincubatedat60°Cfor1h'in t,True)
 c.report()
if __name__=='__main__':selftest()
