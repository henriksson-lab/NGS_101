"""Private workflow-only m6A-label-seq model, published SI FigS4 scope."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,SourceItem,SequencingPlan,Terminal,Disposition
TITLE='m6A-label-seq: published FigS4 workflow boundary'
AUTHORITY='Published SI31 FigS4 originalPDFp5/legend207; published main abstract661 only (subscription preview). Core Methods not available in selected final main. No analytical probe/validation primer is assigned to production library.'
def protocol_graph():
    g=ProtocolGraph(TITLE,authority=AUTHORITY,representation='published_instance');x=g.input('labelled_mrna','Source a6A-labelled cellular mRNA after native metabolic allyl transfer; upstream cellular extraction/formulation unassembled',K.MOLECULAR,primary=True)
    c=g.input('fragment_condition','SourceRNAfragmentation, chemistry/conditions unassembled',K.CONDITION);frag=g.product('fragments','Sourcefragmented a6A-labelled mRNA',K.MOLECULAR);g.react('fragment','SourceFigS4 RNAfragmentation',S.MAGIC,((x,R.TEMPLATE),(c,R.CONDITION)),(frag,),reason='Figure workflow positive; complete fragmentation chemistry unassembled.')
    saved=g.product('saved_input','Source saved input portion before IP; amount and subsequent full path unassembled',K.MOLECULAR);ip=g.product('ip_input','Source remaining fragments for IP; allocation unknown',K.MOLECULAR);g.react('save_input','SourceFigS4 Save some as Input',S.SPLIT,((frag,R.TEMPLATE),),(saved,ip));g.mark_terminal(Terminal(saved.id,Disposition.UNRESOLVED,'Source input exists; no guessed allocation or library recipe.'))
    ab=g.input('a6a_antibody','Source anti-a6A antibody',K.SELECTION);bound=g.product('bound','Sourcea6A-RIP preparation; full binding/wash geometry unassembled',K.MOLECULAR);g.react('ip','Source a6A immunoprecipitation',S.BINDING,((ip,R.TEMPLATE),(ab,R.SELECTION)),(bound,))
    eluent=g.input('a6atp','Sourcea6ATP(4) eluent, concentration unassembled',K.MOLECULAR);elu=g.product('eluate','SourceelutedIP RNA; a6ATP carryover/recovery unassigned',K.MOLECULAR);res=g.product('ip_residual','Unassembled residual immunoprecipitation fraction; physical fate not specified by workflow figure',K.MOLECULAR);g.react('elute','Sourcea6ATPcompetitiveelution',S.MAGIC,((bound,R.TEMPLATE),(eluent,R.STRUCTURAL)),(elu,res),reason='SourceFigS4 elution positive; detailed conditions/recovery unassembled.');g.mark_terminal(Terminal(res.id,Disposition.UNRESOLVED,'Workflow arrow does not establish residual fraction disposal.'))
    untreated=g.product('untreated','Sourcea6A-RIP comparator preparation; actual parallel-stock/split allocation unassembled',K.MOLECULAR);treated=g.product('treated','Source preparation for iodination/cyclization; actual parallel-stock/split allocation unassembled',K.MOLECULAR);g.react('comparison_split','SourceFigS4 untreatedIP versus treatedIP comparisons; actual parallel-stock/split allocation unassembled',S.MAGIC,((elu,R.TEMPLATE),),(untreated,treated),reason='Figure establishes comparison preparations, not a documented physical split recipe; common preparation provenance retained without assigned allocation.')
    iodine=g.input('iodine','Source iodination reagent molecular contribution; exact formulation unassembled',K.MOLECULAR);cond=g.input('iodination_conditions','Source iodination-induced cyclization; complete recipe unassembled',K.CONDITION);cyc=g.product('cyc_rna','Sourcecyc-A treatedRNA and residual a6A/unlabelled classes; conversion yield unassigned',K.MOLECULAR);g.react('cyclize','Source iodination-induced cyclization',S.CHEMICAL_CONVERSION,((treated,R.TEMPLATE),(iodine,R.STRUCTURAL),(cond,R.CONDITION)),(cyc,))
    for arm,start in [('ip',untreated),('treated',cyc)]:
        prep=g.input(arm+'_prep','Source RT/library preparation reagents; complete molecular components and fresh-versus-retained stock identities unassembled',K.MOLECULAR)
        lib=g.product(arm+'_library','Source'+arm+'cDNA library; complete reagents/oligos/products unassembled',K.LIBRARY);g.react(arm+'_library_prep','Source RT and library construction',S.MAGIC,((start,R.TEMPLATE),(prep,R.STRUCTURAL)),(lib,),reason='Published FigS4 explicitly RT/library construction; core Methods inaccessible in selected main. Analytical probe/validation oligos not substituted.');g.sequence(lib,SequencingPlan(limitation='SourceFigS4 high-throughput sequencing. Platform, primers, index identities, full molecular products/read layout unassembled; mutation readout not exact biochemical conversion claim.'))
    g.add_source_item(SourceItem('pending','Actual final main Methods needed; SI31 validation oligos and analytical RT screens are not assigned to core library. Workflow-only source scope.','No confirmed absence or complete source review.'));return g

def selftest():
    from checks import Check
    c=Check();g=protocol_graph();c('Savedinputnotinventedlibrary',g.terminals['saved_input'].disposition,Disposition.UNRESOLVED)
    c('Cyclizationhasmolecularcontribution',any(p.material=='iodine'and p.role==R.STRUCTURAL for p in g.steps['cyclize'].inputs),True)
    c('UntreatedIPnotcyclized',g.steps['ip_library_prep'].inputs[0].material,'untreated')
    c('TreatedIPaftercyclization',g.steps['treated_library_prep'].inputs[0].material,'cyc_rna')
    c('Noanalyticalprimerimport',all(m.kind!=K.OLIGO for m in g.materials.values()),True)
    c('Privatepartialworkflow',bool(g.publication_problems()),True);c.report()
if __name__=='__main__':selftest()
