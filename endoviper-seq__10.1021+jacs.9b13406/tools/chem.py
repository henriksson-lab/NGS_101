"""Private source-qualified partial preparation; no complete molecular/absence claim."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Scene,Segment,feature
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,SourceItem,SequencingPlan,Terminal,Disposition
class B:
    def __init__(self,title,authority,label,kind=K.MOLECULAR):
        self.g=ProtocolGraph(title,authority=authority,representation='published_instance');self.x=self.g.input('source',label,kind,primary=True)
    def stock(self,key,label,kind=K.MOLECULAR,role=R.STRUCTURAL,**kw):return(self.g.input(key,label,kind,**kw),role)
    def cat(self,key,label):return self.stock(key,label,K.CATALYST,R.CATALYST)
    def cond(self,key,label):return self.stock(key,label,K.CONDITION,R.CONDITION)
    def op(self,key,label,kind=S.MAGIC,inputs=(),product='Actual source preparation; whole occupied substrates/products/side products and conversion unassembled',**kw):
        y=self.g.product(key+'_product',product,K.MOLECULAR)
        if kind==S.MAGIC:kw.setdefault('reason','Positive source operation; whole occupied molecules, end distributions, detailed molecular products and conversion unassembled, not a source-absence claim.')
        if kind==S.GENERIC:kw['changes_molecule']=False
        self.g.react(key,label,kind,((self.x,R.TEMPLATE),)+tuple(inputs),(y,),**kw);self.x=y;return y
    def clean(self,key,label):
        y=self.g.product(key+'_kept','Actual retained physical source fraction, carryover/recovery unassigned',K.MOLECULAR);w=self.g.product(key+'_waste','Actual excluded workup fraction',K.WASTE)
        self.g.react(key,label,S.CLEANUP,((self.x,R.TEMPLATE),self.cond(key+'_condition',label)),(y,w));self.g.discard(w,'Actual source workup fraction excluded; no quantitative purification claim.');self.x=y;return y
    def aliquot(self,key,label,kept='Selected source aliquot',rest='Untransferred stock, amount/fate unassigned',qc=False):
        y=self.g.product(key+'_selected',kept,K.MOLECULAR);r=self.g.product(key+'_remainder',rest,K.MOLECULAR)
        self.g.react(key,label,S.SPLIT,((self.x,R.TEMPLATE),),(y,r))
        if qc:self.g.mark_terminal(Terminal(y.id,Disposition.QC,kept));self.x=r;return r
        self.g.mark_terminal(Terminal(r.id,Disposition.UNRESOLVED,rest));self.x=y;return y

    def finish(self):
        y=self.g.product('library','Source selected EndoVIPER SMARTer PicoV2 library; whole molecular state unassembled',K.LIBRARY)
        self.g.react('library_load','Actual selected sequencing loading stock',S.GENERIC,((self.x,R.TEMPLATE),),(y,),changes_molecule=False)
        self.g.sequence(y,SequencingPlan(read_lengths={'Read 1':150,'Read 2':150},limitation='Associated AM Methods96: NextSeq550 paired150; standard8-base i5/i7 values/kit oligos/read landing unassembled. HEADCROP3 of read2 is computational processing, not physical RNA trimming.'))
        return y

AUTHORITY='Associated AM NIHMS1587417/PMC7286354 main76/96/100, not final published equivalence. Synthetic validation76/92 and current kit recipes not imported.'
def protocol_graph():
    b=B('EndoVIPER selected human brain mRNA affinity preparation',AUTHORITY,'Source2ug human brain mRNA; biological extraction upstream outside selected scope');g=b.g
    g.add_source_item(SourceItem('library_input_wording','Associated manuscript Results52 describes approximately4ng for library preparation, whereas detailed Methods96 specifies8ng of each RNAseq/EndoVIPER replicate.','Selected model uses the detailed Methods8ng allocation; actual used value and final published wording unverified. This is an associated-manuscript notation difference, not a claimed final-publication error.'))
    b.op('fragment','Source NEBNext Magnesium RNA Fragmentation Module94C1min',inputs=(b.stock('fragment_reagents','Source fragmentation reagent inventory; packaging and molecular product ends unassembled'),))
    b.clean('fragment_precip','Source ethanol precipitation after fragmentation')
    b.op('glyoxal','Source6percent glyoxal/50percent DMSO in100ul50C1h',inputs=(b.stock('glyoxal_reagent','Actual source glyoxal chemical contributor; adduct occupancy unassembled'),b.cond('glyoxal_conditions','Source50percent DMSO/100ul/50C1h')))
    b.clean('glyoxal_precip','Source ethanol precipitation after glyoxalation')
    b.aliquot('nanodrop','Actual NanoDrop analytical quantification, volume unassigned',qc=True)
    parts=[g.product(k,v,K.MOLECULAR) for k,v in [('control1','Source500ng RNAseq comparator1 in30ul water, frozen at-80C for later use'),('control2','Source500ng RNAseq comparator2 in30ul water, frozen at-80C for later use'),('selected','Source500ng EndoVIPER replicate1'),('replicate2','Source500ng parallel EndoVIPER replicate2'),('remaining','Unallocated glyoxalated RNA stock; recovery and amount unassigned')]]
    g.react('allocations','Source reported500ng to each of four tubes; quantitative recovery not asserted',S.SPLIT,((b.x,R.TEMPLATE),),tuple(parts))
    for z in parts[:2]:g.mark_terminal(Terminal(z.id,Disposition.STORED,z.label+'; actual saved stock, downstream preparation outside selected trace.'))
    for z in parts[3:]:g.mark_terminal(Terminal(z.id,Disposition.UNRESOLVED,z.label+'; downstream preparation outside selected replicate, no assumed pooling.'))
    b.x=parts[2]
    # Source buffer contains BSA; component decomposition is inventory accounting, not a physical separation.
    buf=g.input('ev_buffer','Source1x EV binding buffer with30ug/ml BSA; preparation packaging unassigned',K.MOLECULAR)
    carrier=g.product('buffer_bsa','Same binding-buffer BSA component',K.MOLECULAR);condition=g.product('buffer_conditions','Same buffer19mMTris/100mMNaCl/1mMCaCl2/15uMEDTA/150uMDTT/.025percentTriton/7percentglycerol pH7.4',K.CONDITION)
    g.react('binding_buffer_components','Copresent source buffer components, not physical fractions or extra BSA addition',S.GENERIC,((buf,R.TEMPLATE),),(carrier,condition),changes_molecule=False)
    b.op('bind','Source500ng/250ul100nM MBP-eEndoV+120U RNasin RT30min; binding not source nuclease cleavage',inputs=(b.stock('endov','Actual source MBP-eEndoV recognition protein100nM',K.SELECTION,R.SELECTION),(carrier,R.STRUCTURAL),(condition,R.CONDITION),b.cond('rnasin','Actual120U RNasin Plus protective reagent')))
    beadstock=g.input('bead_slurry','Source300ul anti-MBP magnetic bead slurry; allocation across replicates unassigned',K.MOLECULAR)
    washed=g.product('washed_beads','Source washed anti-MBP bead preparation',K.MOLECULAR);wash=g.product('bead_wash','Actual excluded bead-wash fraction',K.WASTE)
    g.react('wash_beads','Source extensive EV-buffer wash before binding',S.CLEANUP,((beadstock,R.TEMPLATE),b.cond('ev_wash','Source19mMTris/100mMNaCl/1mMCaCl2/7percentglycerol pH7.4')),(washed,wash));g.discard(wash,'Source prewash excluded.')
    bead=g.product('selected_beads','Selected bead allocation for represented replicate; volume unassigned',K.SELECTION);beadrest=g.product('bead_remainder','Other/unallocated washed-bead stock',K.MOLECULAR)
    g.react('bead_allocation','Bead stock selected allocation; per-replicate dose not assigned',S.SPLIT,((washed,R.TEMPLATE),),(bead,beadrest));g.mark_terminal(Terminal(beadrest.id,Disposition.UNRESOLVED,'Stock allocation for parallel sample or remainder unassembled.'))
    bound=g.product('bound','Actual source bead-associated RNA/protein mixture',K.MOLECULAR);unbound=g.product('unbound','Actual discarded capture supernatant',K.WASTE)
    g.react('capture','Source anti-MBP capture RT2h end-over-end rotation',S.AFFINITY_CAPTURE,((b.x,R.TEMPLATE),(bead,R.SELECTION)),(bound,unbound));g.discard(unbound,'Source explicitly discards capture supernatant; no complete specificity/yield claim.');b.x=bound
    b.clean('bound_wash','Source three500ul EV washes')
    eluate=g.product('eluate','Actual95C10min200ul EV eluate SUPERNATANT',K.MOLECULAR);residual=g.product('residual_beads','Residual beads after collection',K.MOLECULAR)
    g.react('elute','Source0.5MTEAA pH8.6/47.5percentformamide/.01percentSDS95C10min',S.CLEANUP,((b.x,R.TEMPLATE),b.cond('elution_buffer','Actual200ul EV elution conditions')),(eluate,residual));g.mark_terminal(Terminal(residual.id,Disposition.UNRESOLVED,'Residual bead fraction fate/remaining RNA unassigned.'));b.x=eluate
    b.clean('spin_filter','Source0.22uM microfuge spin filter removes residual beads; unit retained as printed')
    b.clean('monarch','Source Monarch RNA Cleanup, water elution')
    b.op('deprotect','Source65C2h100ul50percentDMSO/phosphate saline glyoxal removal',inputs=(b.cond('deprotect_buffer','137mMNaCl/2.7mMKCl/8mMNa2HPO4/2mMKH2PO4pH7.4/50percentDMSO'),),product='Source deprotected RNA preparation; full adduct removal and products unassembled')
    b.clean('deprotect_precip','Source ethanol precipitation and water resuspension')
    b.aliquot('bioanalyzer','Actual Bioanalyzer6000RNA Pico analytical allocation',qc=True)
    b.aliquot('kit_input','Source8ng selected enriched RNA for library preparation','Actual8ng EndoVIPER RNA','Untransferred enriched RNA stock')
    kitstock=g.input('kit_inventory','Aggregate historical PicoV2 kit supplies; actual packaging/stock counts and internal order unassigned',K.MOLECULAR)
    i5=g.product('kit_i5','Source8-base i5 sample-index oligo component; whole adapter/index body unassembled',K.OLIGO,chemical_state={'sample_index':feature('endoviper.i5','sample_index','unknown',note='Source8-base i5; value/coordinate/whole oligo unassembled')})
    i7=g.product('kit_i7','Source8-base i7 sample-index oligo component; whole adapter/index body unassembled',K.OLIGO,chemical_state={'sample_index':feature('endoviper.i7','sample_index','unknown',note='Source8-base i7; value/coordinate/whole oligo unassembled')})
    other=g.product('kit_other','Other kit molecular supplies; whole identities and internal reactions unassembled',K.MOLECULAR)
    g.react('kit_components','Inventory decomposition only; not evidence of one prepared mix or fresh duplicate supplies',S.GENERIC,((kitstock,R.TEMPLATE),),(i5,i7,other),changes_molecule=False)
    b.op('kit','Source SMARTer Stranded Total RNASeq Kitv2 Pico Input with8basei5/i7; whole kit chemistry unassembled',inputs=((i5,R.STRUCTURAL),(i7,R.STRUCTURAL),(other,R.STRUCTURAL)))
    b.aliquot('loading','Source selected sequencing loading allocation; mass/volume unassigned')
    b.finish();return g

def selftest():
    from checks import Check,Source,have
    from storage import sources_root
    c=Check();g=protocol_graph()
    c('Native affinity route not EndoV cleavage',g.steps['bind'].kind,S.MAGIC)
    c('Capture consumes selected500ng',g.steps['bind'].inputs[0].material,'selected')
    c('Actual unbound capture fraction discarded',g.terminals['unbound'].disposition,Disposition.DISCARD)
    c('Elution supernatant retained',g.steps['spin_filter'].inputs[0].material,'eluate')
    c('Fresh library input selected after QC',g.steps['kit'].inputs[0].material,'kit_input_selected')
    c('Buffer BSA supplied as molecule',g.materials['buffer_bsa'].kind,K.MOLECULAR)
    c('BSA not fresh duplicate stock','buffer_bsa' in g.steps['binding_buffer_components'].outputs,True)
    c('Only source kit index oligos, no analytical primers',sorted(m.id for m in g.materials.values() if m.kind==K.OLIGO),['kit_i5','kit_i7'])
    c('Both actual kit index components consumed',tuple(i.material for i in g.steps['kit'].inputs[1:3]),('kit_i5','kit_i7'))
    c('Finished loading remainder retained','loading_remainder' in g.terminals,True)
    c('Frozen comparator actual storage fates',tuple(g.terminals[k].disposition for k in ('control1','control2')),(Disposition.STORED,Disposition.STORED))
    p=sources_root()/Path(__file__).resolve().parents[1].name/'jacs.9b13406_PMC7286354.xml.txt'
    if have(c,Source(p,'https://doi.org/10.1021/jacs.9b13406')):
        t=''.join(p.read_text().split());c('Actual source glyoxal deprotection','65°Cfor2hours'in t,True);c('Actual source8ng kit allocation','8ngofeachRNA-seqandEndoVIPER-seq'in t,True)
    c('Private incomplete model',bool(g.publication_problems()),True);c.report()
if __name__=='__main__':selftest()
