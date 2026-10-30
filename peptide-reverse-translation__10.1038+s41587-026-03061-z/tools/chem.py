from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Scene,Segment,revcomp,feature
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,SourceItem,SequencingPlan,Terminal,Disposition
from adapter_ligation import AdapterLigationGeometry
from pcr import PCRGeometry
class B:
    def __init__(self,title,authority,start,kind=K.BIOLOGICAL):
        self.g=ProtocolGraph(title,authority=authority,representation='published_instance');self.x=self.g.input('source',start,kind,primary=True)
    def stock(self,key,label,kind=K.CATALYST,role=R.CATALYST,**kw):return(self.g.input(key,label,kind,**kw),role)
    def cond(self,key,label):return self.stock(key,label,K.CONDITION,R.CONDITION)
    def op(self,key,label,kind,inputs=(),geometry=None,product='Source mixture including incomplete/unreacted/off-target classes; occupancy/recovery unassigned',reason=None):
        y=self.g.product(key+'_product',product,K.MOLECULAR);kw={}
        if kind==S.MAGIC:kw['reason']=reason or 'Positive source operation; complete occupied molecular inputs/products and historical formulation unassembled, not a missing-information conclusion.'
        if kind==S.GENERIC:kw['changes_molecule']=False
        if geometry is not None:kw['geometry']=geometry
        self.g.react(key,label,kind,((self.x,R.TEMPLATE),)+tuple(inputs),(y,),**kw);self.x=y;return y
    def clean(self,key,label,extra=()):
        y=self.g.product(key+'_kept','Actual source retained physical fraction; yield/carryover unassigned',K.MOLECULAR);w=self.g.product(key+'_waste','Actual separated source workup fraction',K.WASTE)
        self.g.react(key,label,S.CLEANUP,((self.x,R.TEMPLATE),self.cond(key+'_condition',label))+tuple(extra),(y,w));self.g.discard(w,'Actual workup fraction excluded; no perfect specificity or molecular yield claim.');self.x=y
    def split(self,key,label,selected,remainder,disposition=Disposition.UNRESOLVED):
        y=self.g.product(key+'_selected',selected,K.MOLECULAR);r=self.g.product(key+'_remainder',remainder,K.MOLECULAR)
        self.g.react(key,label,S.SPLIT,((self.x,R.TEMPLATE),),(y,r));self.g.mark_terminal(Terminal(r.id,disposition,remainder));self.x=y;return y

def unknown_scene():
    s=Scene();s.strand('unassembled',[Segment('source products unassembled','[unassembled]',placeholder=True)]);return s

ORDERS={'ODN-S24': "5'-/AmC6//Sp18//Sp18//Sp18/-d(Cc7Gc7G Tc7Ac7G c7GCc7A c7GTC TC)",
 'ODN-S28': "5'-d(GTA ATA AGA GAC TGC CTA CCG)",
 'ODN-S29': "5'-/DBCO//Sp18//BiodT//Sp18//Sp18/-d(c7Gc7GC c7Ac7AC c7Gc7GC c7ATT TTC)",
 'ODN-S31': "5'-/5TCO//Sp18//Sp18//Sp18/-d(TAC TGT ACC CTC TGT GCG TAT CGG TAG)",
 'ODN-S33': "5'-/5TCO//Sp18//Sp18//Sp18/-d(TAC TGT ACC CTC TGT GCG ATA CGG TAG)",
 'ODN-S41': "5'-2-O3PO-d(TTc7A TTc7A (N)30 Tc7Gc7G CTC c7Gc7Ac7A c7Ac7AT c7GCC c7GTT "
            'c7GCC)-/Sp18//Sp18//Sp18//DBCO/',
 'ODN-S42': "5'-/AmC6//Sp18//Sp18//Sp18/-d(c7A c7AC c7GTC CCT c7G c7AC c7ACT)",
 'ODN-S43': "5'-2-O3PO-d(CCc7A TTT CCT TTc7A Cc7AC)-/Sp18//Sp18//Sp18//DBCO/",
 'ODN-S44': "5'-d(GGA AAT GGA GTG TCA GGG ACG TT)",
 'ODN-S45': "5'-/DBCO//Sp18//Sp18//Sp18/-d(c7GTc7G Tc7Ac7A c7Ac7Gc7G c7Ac7Ac7A Tc7Gc7G c7Ac7GT c7GTC "
            'c7Ac7Gc7G c7Gc7AC c7GTT)',
 'ODN-S46': "5'-d(TCG TCG GCA GCG TCA GAT GTG TAT AAG AGA CAG NNC ATG GCA ACG GCA TTT TCG AG)",
 'ODN-S47': "5'-d(TCG TCG GCA GCG TCA GAT GTG TAT AAG AGA CAG NNN CAT GGC AAC GGC ATT TTC GAG)",
 'ODN-S48': "5'-d(TCG TCG GCA GCG TCA GAT GTG TAT AAG AGA CAG NNN NAT CGG CAA CGG CAT TTT CGA G)",
 'ODN-S49': "5'-d(TCG TCG GCA GCG TCA GAT GTG TAT AAG AGA CAG NNN NNA TCG GCA ACG GCA TTT TCG AG)",
 'ODN-S50': "5'-d(TCG TCG GCA GCG TCA GAT GTG TAT AAG AGA CAG NNN NNN TCA GGC AAC GGC ATT TTC GAG)",
 'ODN-S51': "5'-d(TCG TCG GCA GCG TCA GAT GTG TAT AAG AGA CAG NNN NNN NTC AGG CAA CGG CAT TTT CGA G)",
 'ODN-S52': "5'-d(TCG TCG GCA GCG TCA GAT GTG TAT AAG AGA CAG NNN NNN NNT ACG GCA ACG GCA TTT TCG AG)",
 'ODN-S53': "5'-d(TCG TCG GCA GCG TCA GAT GTG TAT AAG AGA CAG NNN NNN NNN TAC GGC AAC GGC ATT TTC GAG)",
 'ODN-S54': "5'-d(GTC TCG TGG GCT CGG AGA TGT GTA TAA GAG ACA GNN NNT ACT GTA CCC TCT GTG CG)"}
SOURCE_OLIGO_EXCEPTIONS={x:'Original Table2 c7-deazapurine notation and/or (N)30 expansion is not faithfully normalized by generic source text scanner; model preserves raw whole ordering and actual original ZIP/XML fixture. No altered nucleotide identity or source-absence exemption.' for x in ('molecule_barcode','carrier_43','carrier_45')}
AUTHORITY='Publishedmainabstract/changehistory25March2026; publishedSI1Fig4/6/13–17andTable7; originalTable2MOESM3. Threeoriginals; paywalledmainchemistrybodyunread. SelectedFig5four-cyclepartialscope, nocaption-to-completerecipeinference.'
def oligo(name):
    import re
    raw=ORDERS[name];body=raw[raw.index('-d(')+3:raw.rindex(')')];body=body.replace('(N)30','N'*30).replace('c7','');bases=''.join(body.split())
    if name=='ODN-S41':
        left,tail=bases.split('N'*30);segs=[Segment('Source modified DNA prefix',left),Segment('Peptide molecule UMI','N'*30,feature=feature('peptide_molecule','umi','random',note='Identifies originating peptide molecule; not amino-acid identity')),Segment('Source modified DNA suffix',tail)]
    elif name in ('ODN-S31','ODN-S33'):
        segs=[Segment('Common DNA prefix',bases[:18]),Segment('Source antibody tag DNA',bases[18:21],feature=feature('antibody_'+name.lower().replace('-','_'),'feature_barcode','fixed',note='Physical antibody-associated DNA tag; Table7 Fig5 reported identity is complementary orientation; no full reporter/read geometry claimed')),Segment('Source PEA suffix',bases[21:])]
    else:segs=[Segment(name,bases)]
    sc=Scene();sc.strand(name,segs)
    if name=='ODN-S29':
        sc.strand(name+'_biotin_dt',[Segment('Internal biotinylated deoxythymidine','T')]);sc.note(name+'_biotin_dt','Same covalent source oligo: 5DBCO → Sp18 → BiodT → Sp18 → Sp18 → shown DNA body. Separate DNA tracts depict nonnucleotide linker boundaries, not free oligos or a direct dT-to-body phosphodiester bond.')
    sc.note(name,'Source whole ordering: '+raw);return sc

def oligo_input(b,key,name):
    return b.stock(key,name+' source whole modified ordering',K.OLIGO,R.STRUCTURAL,state=oligo(name),chemical_state={'source_ordering':ORDERS[name],'modifications':'c7 denotes deazapurine where printed; Sp18/DBCO/TCO/BiodT/AmC6 and source 2-O3PO retained verbatim; not normal DNA substitutions','whole_conjugate':'Peptide/protein and linker chemistry represented separately; exact occupied covalent complex unassembled',**({'physical_nucleotide_order':'T'+''.join(seg.top for seg in oligo(name).strands[name].segs),'internal_biotin_nucleotide_index_1based':1,'nonnucleotide_boundaries':'5DBCO → Sp18 → biotin-dT nucleotide1 → Sp18 → Sp18 → DNA body; displayed as separate covalently linker-connected DNA tracts, not free oligos'}if name=='ODN-S29'else{})})

def protocol_graph():
    b=B('Peptide reverse translation: selected Fig5 four-cycle source-qualified partial preparation',AUTHORITY,'Actual Fig5 FWFWGGGX peptide with GGGGGGX carrier preparation; amino-acid polymers, NOT nucleotide sequences');g=b.g
    parts=[oligo_input(b,'anchor','ODN-S24'),oligo_input(b,'molecule_barcode','ODN-S41')]
    parts += [oligo_input(b,'carrier_'+str(i),'ODN-S'+str(i))for i in range(42,46)]
    b.op('prepared_beads','Source peptide/barcoding DNA immobilized preparation; detailed mainbody construction unassembled',S.MAGIC,tuple(parts)+(b.stock('carboxyl_beads','Actualsourcecarboxylicacidmagneticbeads',K.SELECTION,R.SELECTION),),product='Peptide-containingbeadpreparation with associated modifiedDNA barcode/carrierclasses; attachment/occupancyunassembled',reason='SourceFig5/Table7 peptide UMI and Table2S41/carrier inventory retained. Full main preparation unavailable in selected publishedpreview, so no guessed attachment chemistry/stockamounts.')
    releases=[]
    for i in range(4):
        k='cycle_'+str(i+1)
        b.op(k+'_tag','ModifiedEdman: introduceazidePITC then proximityDNAconjugation',S.MAGIC,(b.stock(k+'_pitc','Actualazide-modifiedPITC reagent',K.MOLECULAR,R.STRUCTURAL),oligo_input(b,k+'_primer','ODN-S29'),oligo_input(b,k+'_splint','ODN-S28')),reason='SourcepositiveDNA-encodedEdman operation; completecovalentcomplex, reactionconditions, copiedbarcodeextensionenzyme/dose and occupiedstrandgeometryunassembled. DNAstability/polymerasescreens not imported as selectedrecipe.')
        retained=g.product(k+'_retained','Actual bead-retained shortenedpeptide/preparation fornextcycle; residueclasses/recoveryunassigned',K.MOLECULAR)
        released=g.product(k+'_released','Actual released DNA-barcoded PTC amino-acid preparation, including source incompleteextension/A-tail/offproductclasses',K.MOLECULAR)
        g.react(k+'_cleave','ModifiedEdman cleavage: releasedDNA-amino-acidfraction versus retainedshortenedpeptide',S.MAGIC,((b.x,R.TEMPLATE),b.cond(k+'_cleavage_conditions','Actualsourcecleavageprogramme, completechemistryconditionsunassembled')),(retained,released),reason='Publishedabstract/SIFig6 demonstrate release/retention and truncated/A-tailedproducts; no perfectcleavage/specificity claim.');b.x=retained;releases.append(released)
    g.mark_terminal(Terminal(b.x.id,Disposition.UNRESOLVED,'Actual retainedpeptide/carrier/beadpreparation afterfourcycles; not pooledintoreporterDNA.'))
    outputs=[]
    for i,start in enumerate(releases):
        k='read_'+str(i+1);b.x=start
        b.op(k+'_capture','Capture released biotin-bearing DNA/aminoacid preparation on streptavidinbeads',S.MAGIC,(b.stock(k+'_sa','Actualstreptavidinbeads',K.SELECTION,R.SELECTION),),reason='SIFig16 demonstrates bead-immobilizedDNA-barcodedPTC substrates. Selectedsourcefractionwash/elution/occupancy details unassembled, not assumingperfectcapture.')
        antibody_stock=g.input(k+'_antibody_conjugates','Actual prepared anti-PTC-F/W DNA-conjugated antibody stock; no fresh free-oligo addition to PEA',K.MOLECULAR)
        antibody_products=[];antibody_parts=[]
        for suffix,name in (('tag_f','ODN-S33'),('tag_w','ODN-S31')):
            x=g.product(k+'_'+suffix,'DNA component of actual prepared antibody conjugate',K.OLIGO,state=oligo(name),chemical_state={'source_ordering':ORDERS[name],'physical_form':'Component of same source prepared antibody conjugate, not fresh free oligo','modifications':'Whole TCO/Sp18 ordering retained; complete covalent conjugate occupancy unassembled'})
            antibody_products.append(x);antibody_parts.append((x,R.STRUCTURAL))
        proteins=g.product(k+'_antibodies','Protein components of same actual anti-PTC-F/W conjugate stock; not DNA',K.MOLECULAR);antibody_products.append(proteins);antibody_parts.append((proteins,R.STRUCTURAL))
        g.react(k+'_antibody_components','Expose copresent components of one actual prepared conjugate stock; no physical dissociation',S.GENERIC,((antibody_stock,R.TEMPLATE),),tuple(antibody_products),changes_molecule=False)
        antibody_parts=tuple(antibody_parts)
        b.op(k+'_pea','Antibody-associatedproximityextension writes amino-acid identity',S.MAGIC,antibody_parts+(b.stock(k+'_extension','ActualPEApolymerase/NT/buffer aggregate supplies; formulation, packaging and stock counts unassigned',K.MOLECULAR,R.STRUCTURAL),),reason='SourceSIFig13andTable7 tagconjugation/identity scope; fullantibodyconjugateoccupancy and selectedPEAenzyme/wholejunction unassembled. Table2oligos are components, not freeoligos substituted for antibodies.')
        b.op(k+'_pcr','Sourcecycle-specific PCR-amplifiableDNAreporter preparation',S.MAGIC,tuple(oligo_input(b,k+'_fp_'+str(j),'ODN-S'+str(46+2*i+j))for j in range(2))+(oligo_input(b,k+'_rp','ODN-S54'),b.stock(k+'_pcr_stock','ActualPCRcatalyst/NT/buffer aggregate supplies; formulation, packaging and stock counts unassigned',K.MOLECULAR,R.STRUCTURAL)),reason='SourceTable2 assigns two staggeredforwardorders percycle plusRP. Their exact physicalallocation (mixed/separatereactions), completeamplification and finalIllumina construction unassembled; Nstaggerwindows NOTUMIs.');outputs.append(b.x)
    chosen=[]
    for i,x in enumerate(outputs):
        b.x=x;b.split('pool_'+str(i),'Source reporter preparation selectedfor sequencing; amount unassigned','Selectedreporter aliquot','Untransferredreporter stock');chosen.append(b.x)
    y=g.product('reporter_pool','LogicalcollectionofselectedcycleDNAreporteraliquots; physicalpooling/wholelibrary/readlayoutunassembled',K.MOLECULAR)
    g.react('pool','Opaque finishing interface for selected cycle reporters; physical pooling unassigned',S.MAGIC,tuple((x,R.TEMPLATE)for x in chosen),(y,),reason='Publishedabstract establishes DNA reporter NGS readout; whether/how cycle aliquots physically mix, pool amounts and remaining construction are unassembled. This interface is not a physical pooling claim.');b.x=y
    b.op('library_load','SourceIlluminaDNAreadout',S.MAGIC,(b.cond('sequencing_condition','IlluminaNGS; runlengths/indexorientationunassigned'),),reason='PositivepublishedNGSreporterreadout, notcomplete library readiness.');g.mark_terminal(Terminal(b.x.id,Disposition.UNRESOLVED,'Source-positiveDNA reporter readout; full library/primers/readlengthsunassembled.'))
    g.add_source_item(SourceItem('published_scope',AUTHORITY,'Partial positive component/operation scope, not source absence or library readiness.'))
    return g
def selftest():
    from checks import Check
    import xml.etree.ElementTree as ET,zipfile
    from storage import sources_root
    c=Check();g=protocol_graph()
    c('Four retained peptide fractions are separate from released reporters',all('cycle_'+str(i)+'_retained' in g.materials and 'cycle_'+str(i)+'_released' in g.materials for i in range(1,5)),True)
    c('Peptide source has no DNA scene',g.materials['source'].state,None)
    c('Each next Edman cycle consumes retained peptide fraction',all(any(p.material=='cycle_'+str(i)+'_retained' for p in g.steps['cycle_'+str(i+1)+'_tag'].inputs)for i in range(1,4)),True)
    c('PEA input follows its own physically released cycle fraction',all(any(p.material=='cycle_'+str(i)+'_released'for p in g.steps['read_'+str(i)+'_capture'].inputs)for i in range(1,5)),True)
    c('Reporters retain an untransferred stock',all('pool_'+str(i)+'_remainder' in g.materials for i in range(4)),True)
    c('Antibody tags are physical DNA components, not peptide sequences',all(g.materials['read_'+str(i)+'_tag_f'].kind==K.OLIGO for i in range(1,5)),True)
    c('Each PEA consumes components from one prepared antibody conjugate stock',all(g.steps['read_'+str(i)+'_antibody_components'].inputs[0].material=='read_'+str(i)+'_antibody_conjugates' and set(g.steps['read_'+str(i)+'_antibody_components'].outputs)=={'read_'+str(i)+'_tag_f','read_'+str(i)+'_tag_w','read_'+str(i)+'_antibodies'}for i in range(1,5)),True)
    c('BiodT physical nucleotide survives separate linker-connected tract',oligo('ODN-S29').strands['ODN-S29_biotin_dt'].segs[0].top,'T')
    c('Biotin is attached to physical nucleotide1, not a nonnucleotide replacement',g.materials['cycle_1_primer'].chemical_state['internal_biotin_nucleotide_index_1based'],1)
    c('S29 physical nucleotide ordering includes BiodT before all source DNA-body nucleotides',g.materials['cycle_1_primer'].chemical_state['physical_nucleotide_order'],'TGGCAACGGCATTTTC')
    c('N30 is source peptide molecule UMI',oligo('ODN-S41').strands['ODN-S41'].segs[1].feature.role,'umi')
    c('Staggered PCR windows are not additional UMI features',all(not seg.feature for i in range(46,55) for seg in oligo('ODN-S'+str(i)).strands['ODN-S'+str(i)].segs),True)
    path=sources_root()/Path(__file__).resolve().parents[1].name/'s41587-026-03061-z_41587_2026_3061_MOESM3_ESM.xlsx'
    if path.exists():
        ns={'s':'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
        with zipfile.ZipFile(path) as z:
            strings=[''.join(x.itertext())for x in ET.fromstring(z.read('xl/sharedStrings.xml')).findall('s:si',ns)]
            root=ET.fromstring(z.read('xl/worksheets/sheet1.xml'));rows={}
            for row in root.findall('.//s:row',ns):
                vals=[]
                for cell in row.findall('s:c',ns):
                    v=cell.find('s:v',ns);s=v.text if v is not None else '';vals.append(strings[int(s)]if cell.get('t')=='s'else s)
                if len(vals)>=2:rows[vals[0].strip()]=vals[1]
            for name,value in ORDERS.items():c('Original Table2 whole modified ordering '+name,rows.get(name),value)
    else:print('SKIP Original Table2 workbook absent; source fixtures skipped')
    c.report()
if __name__=='__main__':selftest()
