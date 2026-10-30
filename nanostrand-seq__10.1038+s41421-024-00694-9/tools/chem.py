"""Source-owned private partial NanoStrand preparation; geometry remains unassembled."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from nextera import ME_RC,ADAPTOR_S5
from illumina import TRUSEQ_READ1
TITLE='NanoStrand-seq selected GM12878 route'
AUTHORITY='Published PMC11231365 main76/172–194; publishedSI FigS1caption and wholeTableS1workbook context inspected,3originals. No historical ONT kit guide review.'
INITIAL='Source GM12878 cultured cells for RO3306 synchronization and24hBrdU exposure; mouse concentrations separate'
CELL_INPUT=True
READOUT='Source ONT PromethION48/R9.4.1 SQK-LSK109, native kit finish applied to amplified parental-template derivatives. Distinct 24nt first/second cell barcode identities known but actual codes unassigned due apparent TableS1 conflict; no UMI assigned. PCR products are not original native duplex or methylation-retention claims.'
DEFERRED='Historical VazymeTn5/SQK-LSK109 guides and complete source-specific primer/end geometry unassembled. ApparentTableS1 conflict: all40 Bar1–24/49–64 sequences repeat AATGATACGGCGACCACCGAGATCT while main188/190 require distinct24nt codes. Main176 explicitly5P/3NH2 ME; TableS1lists5P only. Actual code assignment leftopen; no accepted absence claim.'
AMP=TRUSEQ_READ1[:24]
TAG_ARM=ADAPTOR_S5[:21]
OLIGOS={
 'Source ME short5P3NH2':{'segments':lambda:[Segment('ME reverse',ME_RC)],'mod5':'P','mod3':'NH2','chemistry':{'polymer':'DNA','five_prime':'P','three_prime':'NH2'}},
 'Source S5 long':{'segments':lambda:[Segment('S5 adapter',ADAPTOR_S5)],'chemistry':{'polymer':'DNA'}},
 'Source first strand tagging primer':{'segments':lambda:[Segment('amplification arm',AMP),Segment('first cellbarcode','N'*24,placeholder=True,feature=feature('nanostrand_first','cell_barcode','unknown')),Segment('source tagging arm',TAG_ARM)],'chemistry':{'polymer':'DNA','barcode_width':24,'barcode_code':'unassigned;TableS1apparentconflict'}},
 'Source second strand tagging primer':{'segments':lambda:[Segment('amplification arm',AMP),Segment('second cellbarcode','N'*24,placeholder=True,feature=feature('nanostrand_second','cell_barcode','unknown')),Segment('source tagging arm',TAG_ARM)],'chemistry':{'polymer':'DNA','barcode_width':24,'barcode_code':'unassigned;TableS1apparentconflict'}},
 'Source single amplification primer':{'segments':lambda:[Segment('source amplification primer',AMP)],'chemistry':{'polymer':'DNA'}}}

from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,SourceItem,SequencingPlan,Terminal,Disposition
from chemdraw import Scene,Segment,feature
from adapter_ligation import AdapterLigationGeometry
from pcr import PCRGeometry

def scene(label):
    item=OLIGOS.get(label)
    if item is None:return None
    sc=Scene();sc.strand('oligo',item['segments'](),mod5=item.get('mod5',''),mod3=item.get('mod3',''));return sc

def chemistry(label):return OLIGOS.get(label,{}).get('chemistry')
def unknown_scene():
    sc=Scene();sc.strand('unassembled',[Segment('unassembled geometry','[source geometry unassembled]',placeholder=True)]);return sc

def protocol_graph():
    g=ProtocolGraph(TITLE,authority=AUTHORITY,representation='published_instance')
    x=g.input('input',INITIAL,K.BIOLOGICAL if CELL_INPUT else K.MOLECULAR,primary=True)
    def stock(key,label):return g.input(key,label,K.OLIGO,state=scene(label),chemical_state=chemistry(label))
    def op(key,label,kind,reagent,product,oligos=(),symmetric=False):
        nonlocal x
        enzyme=kind in {S.DIGESTION,S.END_REPAIR,S.EXTENSION,S.TAILING,S.LIGATION,S.PCR,S.DEPROTEINIZATION,S.TAGMENTATION,S.TRANSCRIPTION}
        c=g.input(key+'_stock',reagent,K.CATALYST if enzyme else K.CONDITION)
        ins=[(x,R.TEMPLATE),(c,R.CATALYST if enzyme else R.CONDITION)]
        for i,o in enumerate(oligos):ins.append((stock(key+'_oligo_'+str(i),o),R.STRUCTURAL))
        kw={}
        if kind==S.MAGIC:kw['reason']='Source-specific compound transition or molecular geometry remains unassembled in this selected model.'
        if kind==S.GENERIC:kw['changes_molecule']=False
        if kind==S.LIGATION:kw['geometry']=AdapterLigationGeometry('unavailable','unavailable',(unknown_scene(),),'Actual source ligation; full source-specific substrate/junction geometry unassembled, not a source absence claim.')
        if kind==S.PCR:kw['geometry']=PCRGeometry.unavailable(unknown_scene(),reason='Source amplification known; whole source-specific primer landing and product geometry unassembled.',symmetric_species=symmetric)
        y=g.product(key+'_product',product,K.MOLECULAR);g.react(key,label,kind,ins,(y,),**kw);x=y
    def clean(key,label):
        nonlocal x
        c=g.input(key+'_condition',label,K.CONDITION);y=g.product(key+'_kept','Physically recovered source preparation; exact recovery/carryover unassigned',K.MOLECULAR);w=g.product(key+'_waste','Physically separated source workup fractions',K.WASTE)
        g.react(key,label,S.CLEANUP,((x,R.TEMPLATE),(c,R.CONDITION)),(y,w));g.discard(w,'Source physical workup, not perfect molecular exclusion.');x=y
    def split(key,label,rest_label,disposition=Disposition.UNRESOLVED):
        nonlocal x
        y=g.product(key+'_used',label,K.MOLECULAR);r=g.product(key+'_rest',rest_label,K.QC if disposition==Disposition.QC else K.MOLECULAR);g.react(key,label,S.SPLIT,((x,R.TEMPLATE),),(y,r));g.mark_terminal(Terminal(r.id,disposition,rest_label));x=y
    def qc(key,label):split(key,'Production preparation retained after separate analytical aliquot',label,Disposition.QC)
    def external_anneal(key,primer,label,enzyme,programme,product,rt=False):
        nonlocal x
        p=stock(key+'_primer',primer);c=g.input(key+'_anneal_condition',label,K.CONDITION)
        t=g.product(key+'_template','Same-tube template component after annealing; external primer tracked separately, not an isolated fraction',K.MOLECULAR)
        retained=g.product(key+'_retained_primer','Same-tube co-retained external primer component, not a second stock addition or intramolecular geometry',K.OLIGO,state=scene(primer),chemical_state=chemistry(primer))
        g.react(key+'_anneal',label,S.THERMAL_TREATMENT,((x,R.TEMPLATE),(p,R.STRUCTURAL),(c,R.CONDITION)),(t,retained))
        e=g.input(key+'_enzyme',enzyme,K.CATALYST);y=g.product(key+'_product',product,K.MOLECULAR)
        g.react(key,programme,S.REVERSE_TRANSCRIPTION if rt else S.EXTENSION,((t,R.TEMPLATE),(retained,R.STRUCTURAL),(e,R.CATALYST)),(y,));x=y

    # Independent stock preparation; the nucleus receives one loaded-stock aliquot.
    me=stock('stock_me','Source ME short5P3NH2');s5=stock('stock_s5','Source S5 long')
    anneal_condition=g.input('stock_anneal_condition','Source oligos dissolved100uM in PCRgradewater,annealed50uM; thermal programme not assigned',K.CONDITION)
    duplex=g.product('stock_adapter_duplex','Source annealed ME/S5 adapter stock; exact loading geometry unassembled',K.ADAPTER)
    g.react('stock_adapter_anneal','Source independent ME/S5 annealing before transposon loading',S.HYBRIDIZATION,((s5,R.TEMPLATE),(me,R.STRUCTURAL),(anneal_condition,R.CONDITION)),(duplex,))
    tn5=g.input('stock_tn5','Source10ul Vazyme TagmentEnzyme2ugul for50ul loading reaction',K.CATALYST)
    coupling=g.input('stock_coupling','Source7ul50uM adapter/33ulCouplingBuffer,55C10min',K.CONDITION)
    loaded=g.product('loaded_stock','Source50ul adapter-loaded Tn5 preparation; precise stoichiometry/unbound companion species unassembled',K.MOLECULAR)
    g.react('stock_load','Source independent transposon complex assembly55C10min',S.MAGIC,((duplex,R.TEMPLATE),(tn5,R.CATALYST),(coupling,R.CONDITION)),(loaded,),reason='Actual source enzyme/adapter assembly; exact protein-adapter stoichiometry and molecular loading states unassembled, no invented free-oligo addition to nuclei.')
    dilution=g.input('stock_dilution_buffer','Source1:5000 stock dilution in50percentglycerol/50mMHEPESpH7.2/100mMNaCl/0.1mMEDTA/1mMDTT/0.1percentTriton;single-usealiquots/-80C',K.CONDITION)
    diluted=g.product('diluted_loaded_stock','Source diluted loaded-transposon preparation',K.MOLECULAR)
    g.react('stock_dilute','Source1:5000 dilution',S.GENERIC,((loaded,R.TEMPLATE),(dilution,R.CONDITION)),(diluted,),changes_molecule=False)
    used=g.product('loaded_stock_aliquot','Source1ul diluted loaded-transposon stock for chosen nucleus',K.MOLECULAR);rest=g.product('loaded_stock_remainder','Remaining diluted transposon stock; actual remaining volume/fate unassigned',K.MOLECULAR)
    g.react('stock_aliquot','Allocate source1ul diluted complex stock',S.SPLIT,((diluted,R.TEMPLATE),),(used,rest));g.mark_terminal(Terminal(rest.id,Disposition.UNRESOLVED,'Stock remainder is not silently added to the selected nucleus.'))
    bound_enzyme=g.product('loaded_tn5_component','Enzyme component of the same chosen loaded-stock aliquot; not a new Tn5 addition or physical fraction',K.CATALYST)
    bound_adapter=g.product('loaded_adapter_component','ME/S5 adapter component of that same loaded-stock aliquot; not newly added free oligos or a physical fraction',K.ADAPTER)
    g.react('loaded_stock_components','Non-transforming component representation of one loaded aliquot; companions/endgeometry remain unassembled',S.GENERIC,((used,R.TEMPLATE),),(bound_enzyme,bound_adapter),changes_molecule=False)
    op('synchronize','Source10uMRO3306 4h at50–60percentconfluency',S.GENERIC,'Source RO3306 treatment','Source synchronized cells')
    op('brdu','Source release into40uMBrdU24h',S.MAGIC,'Source in-cell nucleotide incorporation; one-cycle labeling classes unassembled','Source GM12878 BrdU-labeled daughter/other cell classes; incomplete labeling retained')
    clean('cellwash','Source twicecoldPBSwash')
    op('nuclei','Source nuclei isolation100mMTrispH7.4/150mMNaCl/1mMCa/0.5mMMg/0.1percentNP40/2percentBSA10minice',S.PERMEABILIZATION,'Source nuclei isolationbuffer','Source nuclei suspension and residual cell material')
    op('stain_sort','SourceHoechst33258/PI10ugml each; choose half-HoechstG1 nuclei',S.MAGIC,'Source staining/gating/physical sorting; nonselectedfractionunassigned','Source chosen single nucleus and unselected nuclei classes')
    split('chosen_nucleus','Source selected single nucleus into2.5ullysisbuffer','Unselected/other collected nuclei outside chosen preparation')
    op('lysis','SourceQiagenprotease1mgml/0.3percentTriton/20mMKCl/2mMTE50C3h70C30min;optional-80Cstorage',S.DEPROTEINIZATION,'Source Qiagenprotease','Source nucleus lysate with DNA/RNA/protein remnants')
    mix=g.input('tagment_condition','Source2.5ul transpositionmix:1ul5xTAPS_PEG8K/1ul chosen loadedstock/0.5ulwater;55C10min',K.CONDITION)
    tagged=g.product('tagment_product','Source singleS5-loaded tagmentation products with residual untagged DNA, protein/adapter companion species and reagents',K.MOLECULAR)
    g.react('tagment','Source tagmentation with prepared1ul diluted loaded-stock aliquot; no fresh free ME/S5 oligos',S.TAGMENTATION,((x,R.TEMPLATE),(bound_enzyme,R.CATALYST),(bound_adapter,R.STRUCTURAL),(mix,R.CONDITION)),(tagged,));x=tagged
    op('tn5_stop','Source1ul0.25percentSDS55C10min',S.QUENCHING,'Source SDS','Source dissociated/denatured transposase preparation')
    op('sds_quench','Source1.25ul20percentTween20 roomtemperature5min',S.QUENCHING,'Source Tween20','Source SDS-quenched mixture')
    op('gap','Source0.5ulBst2.0/1ul10xIsothermal/0.5ul10mMdNTP/0.75ulwater65C30min',S.EXTENSION,'Source Bst2.0 DNA polymerase','Source gap-filled blunt-end fragments; exact product states unassembled')
    clean('gap_clean','Source0.8xAMPure/elute9ulwater')
    op('uv_nick','Source1ulHoechst33258(1mgml),365nm3000Jm2',S.CHEMICAL_CONVERSION,'Source Hoechst-sensitized UVA treatment','Source BrdU-associated nicked daughter-strand preparation; not complete physical daughterDNA removal')
    op('tag1','Source13.5ulMix1,200nM firstprimer;94C2min62C5min68C10min',S.EXTENSION,'Source0.2ulTksGflex/10ul2xbuffer','Source first-strand-tagged derivative classes plus nicked labeled strands/remnants',('Source first strand tagging primer',))
    op('exo1','Source2ulExoI37C30min80C20min',S.DIGESTION,'Source ExoI','Source primer-digested firststrand mixture; remnants retained beforecleanup')
    clean('tag1_clean','Source0.8xAMPure/elute10ulwater')
    op('tag2','Source13.5ulMix2,200nM secondprimer;94C2min62C5min68C10min',S.EXTENSION,'Source separate0.2ulTksGflex/10ul2xbuffer','Source second-strand-tagged derivative classes; first/secondcodesdistinctinmain, valuesunassigned',('Source second strand tagging primer',))
    op('exo2','Source2ulExoI37C30min80C20min',S.DIGESTION,'Source freshExoI addition','Source primer-digested taggedproduct mixture')
    clean('tag2_clean','Source0.8xAMPure/elute11ulwater')
    op('pcr','Source25ul;94C3min21cycles[98C15s62C30s68C8min],68C10min4Chold',S.PCR,'Source1.5ul20uM amplificationprimer/0.3ul1.25UulTksGflex/12.5ul2xbuffer','Source amplified parental-template derivative library classes; PCRbias/backgroundunassigned',('Source single amplification primer',),symmetric=True)
    clean('pcr_clean','Source0.8xAMPure')
    qc('quant','Separate source QubitdsDNAHS analytical aliquot')
    others=g.input('other_cell_amplicons','Other separately prepared source cells with distinct first/secondbarcodepairs',K.MOLECULAR,primary=True);pool=g.product('pool_product','Source equalamount pooled cellamplicons',K.MOLECULAR);g.react('pool','Source equalamount physical pooling',S.POOL,((x,R.TEMPLATE),(others,R.TEMPLATE)),(pool,));x=pool
    clean('long_clean1','Source0.7xAMPure first long-fragment selection')
    clean('long_clean2','Source0.7xAMPure second long-fragment selection')
    qc('size_qc','Separate source representative Bioanalyzer sizeQC; not recoveryinto production')
    split('ont_aliquot','Source~1ug pooled amplicon input for ONT kit','Unselected pooled amplicon mass/fate unassigned')
    op('ont_kit','Source SQK-LSK109 adaptor library construction',S.MAGIC,'Historical ONTkit substrate/end/adapter internals unreviewed','Source ONTadapter-bearing amplified long-fragment preparation')
    library=g.product('library','Source sequencing preparation; full source-specific library geometry unassembled',K.LIBRARY)
    g.react('loading','Physical source loading preparation',S.GENERIC,((x,R.TEMPLATE),),(library,),changes_molecule=False)
    g.sequence(library,SequencingPlan(limitation=READOUT))
    g.add_source_item(SourceItem('pending',DEFERRED,'Selected source scope only; incomplete linked-source/whole molecular review, no accepted absence claim.'))
    return g

def selftest():
    from checks import Check,Source,have
    from storage import sources_root
    c=Check();g=protocol_graph()
    c('Tagment consumes produced loaded-stock components',[(i.material,i.role)for i in g.steps['tagment'].inputs if i.role in {R.CATALYST,R.STRUCTURAL}],[('loaded_tn5_component',R.CATALYST),('loaded_adapter_component',R.STRUCTURAL)])
    c('Raw loading stocks added only once',tuple(sum(i.material==key for st in g.steps.values()for i in st.inputs)for key in ('stock_me','stock_s5','stock_tn5')),(1,1,1))
    c('Prepared stock remainder explicit',g.terminals['loaded_stock_remainder'].disposition,Disposition.UNRESOLVED)
    c('Source S5 whole ordering',ADAPTOR_S5,'TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG')
    c('BrdU nicking follows gap-fill cleanup',g.steps['uv_nick'].inputs[0].material,'gap_clean_kept')
    c('First strand tag precedes ExoI',g.steps['exo1'].inputs[0].material,'tag1_product')
    c('Second strand tag follows first cleanup',g.steps['tag2'].inputs[0].material,'tag1_clean_kept')
    c('Separate second ExoI addition',g.steps['exo2'].inputs[0].material,'tag2_product')
    c('PCR uses cleaned twice-tagged preparation',g.steps['pcr'].inputs[0].material,'tag2_clean_kept')
    c('One amplification-primer species',sum(p.role==R.STRUCTURAL for p in g.steps['pcr'].inputs),1)
    c('First whole primer pattern transcription',''.join(s.top for s in scene('Source first strand tagging primer').strands['oligo'].segs),'ACACTCTTTCCCTACACGACGCTC'+'N'*24+'TCGTCGGCAGCGTCAGATGTG')
    c('Second whole primer pattern transcription',''.join(s.top for s in scene('Source second strand tagging primer').strands['oligo'].segs),'ACACTCTTTCCCTACACGACGCTC'+'N'*24+'TCGTCGGCAGCGTCAGATGTG')
    c('Source ME3prime NH2 retained',chemistry('Source ME short5P3NH2')['three_prime'],'NH2')
    c('Long-fragment selection after pool',g.steps['long_clean1'].inputs[0].material,'pool_product')
    c('Unselected pool not silently consumed',g.terminals['ont_aliquot_rest'].disposition,Disposition.UNRESOLVED)
    src=Source(sources_root()/Path(__file__).resolve().parents[1].name/'s41421-024-00694-9_supp_41421_2024_694_MOESM2_ESM.xlsx','Published supplementaryTableS1 at doi10.1038/s41421-024-00694-9')
    if have(c,src,label='original workbook barcode context'):
        import zipfile,xml.etree.ElementTree as E
        ns={'m':'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
        with zipfile.ZipFile(src.path)as z:
            shared=[''.join(x.itertext())for x in E.fromstring(z.read('xl/sharedStrings.xml'))]
            sh=E.fromstring(z.read('xl/worksheets/sheet2.xml'));rows=sh.findall('m:sheetData/m:row',ns)
            values=[shared[int(row.find("m:c[@r='B"+row.attrib['r']+"']/m:v",ns).text)]for row in rows[7:]]
            c('Original TableS1 sequence column whole barcode rows',values,['AATGATACGGCGACCACCGAGATCT']*40)
            c('Original first whole ordering template',shared[int(rows[4].find("m:c[@r='B5']/m:v",ns).text)],'ACACTCTTTCCCTACACGACGCTC[Bar1-24]TCGTCGGCAGCGTCAGATGTG')
            c('Original second whole ordering template',shared[int(rows[5].find("m:c[@r='B6']/m:v",ns).text)],'ACACTCTTTCCCTACACGACGCTC[Bar49-64]TCGTCGGCAGCGTCAGATGTG')
    c('Private partial graph remains incomplete',bool(g.publication_problems()),True)
    c.report()
if __name__=='__main__':selftest()
