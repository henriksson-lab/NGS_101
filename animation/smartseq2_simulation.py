"""Protocol-owned reaction recipe: all products come from the connectivity engine."""
from dataclasses import asdict
from .simulator import Simulator,Polymerase
from .model import complement


def seq(state,ids):return ''.join(state.atoms[i].base for i in ids)


def simulate(m,body,*,anchor='GA',cut_offsets=(8,32),i5='TAGATCGC',i7='TCGCCTTA'):
    source='SMART-seq2 protocol definitions; selected synthetic transcript;2014 TSO form'
    s=Simulator();s.supply_nt('DNA','ACGT',source=source)
    rt=Polymerase('SuperScript II',('RNA','DNA'),'DNA',terminal_bases=('C',))
    pol=Polymerase('DNA polymerase',('DNA',),'DNA')
    preamp_pol=Polymerase('KAPA HiFi',('DNA',),'DNA')
    rna=s.add('rna',body+complement(anchor[::-1]).replace('T','U')+'A'*m.DT_LEN,sugar='RNA',features=[(0,len(body),'insert'),(len(body),len(body)+len(anchor),'rna.anchor'),(len(body)+len(anchor),len(body)+len(anchor)+m.DT_LEN,'polyA')],source=source)
    dt=s.add('dt',m.OLIGO_DT[:-2]+anchor,features=[(0,len(m.ISPCR),'smart.handle'),(len(m.ISPCR),len(m.ISPCR)+len(m.DT_LINKER),'smart.dt_linker'),(len(m.ISPCR)+len(m.DT_LINKER),len(m.OLIGO_DT)-len(anchor),'polyT'),(len(m.OLIGO_DT)-len(anchor),len(m.OLIGO_DT),'anchor')],source=source)
    out={'reagents':{},'source':source,'simulator':s,'rna_ids':rna,'dt_ids':dt,'rt_primer':s.sequence('dt')}
    out['anchor_before']=s.snapshot()
    out['anchor_bind']=s.anneal('dt','rna',binding=m.DT_LEN+len(anchor),source=source)
    out['rt_before']=s.snapshot()
    out['rt']=s.extend('dt',rt,source=source);out['first_ids']=s.chain('dt');out['first']=s.sequence('dt')
    out['terminal_before']=s.snapshot()
    out['terminal']=s.terminal_add('dt','C',3,rt,source=source+' explicit illustrative terminal extent')
    out['tail_ids']=s.chain('dt');out['tail']=s.sequence('dt')
    tso_seq=''.join(part.top for part in m.drawn_tso());T=len(tso_seq)
    tso=s.add('tso',tso_seq,sugars=[(0,T-3,'DNA'),(T-3,T-1,'RNA'),(T-1,T,'LNA')],source=source)
    out['tso_ids']=tso
    out['switch_before']=s.snapshot()
    out['switch_bind']=s.anneal('dt','tso',binding=3,template_site=T-1,source=source)
    out['switch']=s.extend('dt',rt,source=source);out['full_ids']=s.chain('dt');out['full']=s.sequence('dt')
    out['reagents']['preamp_reverse']=s.add('preamp_reverse',m.ISPCR,source=source)
    out['firstcopy_before']=s.snapshot()
    out['rna_melt']=s.unpair(rna,source='PCR thermal separation')
    out['tso_melt']=s.unpair(tso,source='PCR thermal separation')
    out['firstcopy_bind']=s.anneal('preamp_reverse','dt',binding=len(m.ISPCR),source=source)
    out['firstcopy']=s.extend('preamp_reverse',preamp_pol,source=source)
    out['first_reverse_ids']=s.chain('preamp_reverse');out['first_reverse']=s.sequence('preamp_reverse')
    for name in ['preamp_r2','preamp_f2']:out['reagents'][name]=s.add(name,m.ISPCR,source=source)
    out['pcr_before']=s.snapshot()
    out['preamp_melt']=s.unpair(out['full_ids'],source='next PCR thermal separation')
    owners={}
    for name,template in [('preamp_r2','dt'),('preamp_f2','preamp_reverse')]:
        out[name+'_bind']=s.anneal(name,template,binding=len(m.ISPCR),source=source)
        owners[s.tail(name)]=name;owners[s.tail(template)]=name+'_parent'
    out['preamp_operations']=s.extend_available(preamp_pol,['dt','preamp_reverse','preamp_r2','preamp_f2'],source=source)
    out['pcr_after']=s.snapshot()
    for op in out['preamp_operations']:
        name=owners.get(op.details['starting_3prime'])
        if name is not None:out[name]=op
    for name in ['preamp_r2','preamp_f2']:
        out[name+'_ids']=s.chain(name)
    top=s.chain('dt');bottom=s.chain('preamp_r2')
    out['target_top_ids']=top;out['target_bottom_ids']=bottom
    a=len(dt)+cut_offsets[0];b=len(dt)+cut_offsets[1]
    out['cut_left']=a;out['cut_right']=b
    nx=m.nx;me=nx.ME
    for name,adapter in [('left_top_adapter',nx.ADAPTOR_S5),('left_bottom_adapter',nx.ADAPTOR_S5),('right_top_adapter',nx.ADAPTOR_S7),('right_bottom_adapter',nx.ADAPTOR_S7)]:
        out['reagents'][name]=s.add(name,adapter,features=[(0,len(adapter)-len(me),'nextera.s5'if name.startswith('left')else'nextera.s7'),(len(adapter)-len(me),len(adapter),'nextera.me')],source=source+' canonical loaded adapter')
        short=name+'_nontransfer'
        out['reagents'][short]=s.add(short,complement(me[::-1]),five='unassigned',three='unassigned',source=source+' ME complement; terminal chemistry unassigned')
        s.anneal(short,name,binding=len(me),template_site=len(adapter)-1,source=source)
    out['tag_before']=s.snapshot()
    out['transfer_left']=s.tn5('left_top_adapter','left_bottom_adapter',top[a],bottom[len(bottom)-a-9],stagger=9,complex_id='tn5_left',source=source)
    out['transfer_right']=s.tn5('right_top_adapter','right_bottom_adapter',top[b],bottom[len(bottom)-b-9],stagger=9,complex_id='tn5_right',source=source)
    # Removing protein and melting only the inter-fragment9bp bridges do not
    # destroy the duplexes within the fragments selected to leave the viewport.
    out['release_left']=s.release_protein('tn5_left',source='source transposase-removal workup')
    out['release_right']=s.release_protein('tn5_right',source='source transposase-removal workup')
    out['left_bridge_melt']=s.unpair(top[a:a+9],source='illustrative separation of overlapping adjacent fragments')
    out['right_bridge_melt']=s.unpair(top[b:b+9],source='illustrative separation of overlapping adjacent fragments')
    selected_top=s.chain(top[a]);selected_bottom=s.chain(bottom[len(bottom)-b-9])
    out['tag_top_ids']=selected_top;out['tag_bottom_ids']=selected_bottom
    out['tag_after']=s.snapshot()
    out['repair_before']=s.snapshot()
    for name,anchor_id,short in [('top',top[a],'right_bottom_adapter_nontransfer'),('bottom',bottom[len(bottom)-b-9],'left_top_adapter_nontransfer')]:
        # Ordinary extension naturally stalls at the still-paired ME complement.
        # No expected nine-base product or desired endpoint is passed to it.
        out[name+'_gap']=s.extend(anchor_id,pol,source=source)
        if out[name+'_gap'].details['stop']!='occupied_template':raise ValueError('selected repair did not encounter its ME complement')
        out[name+'_me_release']=s.unpair(s.chain(short),source='ME complement loss during repair; displacement/digestion pathway unassigned')
        out[name+'_adapter_copy']=s.extend(anchor_id,pol,source=source)
    out['repaired_top_ids']=s.chain(top[a]);out['repaired_bottom_ids']=s.chain(bottom[len(bottom)-b-9])
    out['repair_after']=s.snapshot()
    out['repaired_top']=s.sequence(top[a]);out['repaired_bottom']=s.sequence(bottom[len(bottom)-b-9])
    # Index-PCR calls exactly the same extension operation as RT/gap repair.
    # Physical template ends determine first-cycle product lengths automatically.
    F=nx.n5xx_primer(i5);R=nx.n7xx_primer(i7)
    fsymbols={len(m.P5)+j:f'i5:{j}'for j,c in enumerate(i5)if c not in 'ACGT'}
    rsymbols={len(m.P7)+j:f'i7:{j}'for j,c in enumerate(i7)if c not in 'ACGT'}
    out['index_f']=F;out['index_r']=R;out['i5_oligo']=i5;out['i7_oligo']=i7
    for name,primer,symbols in [('index_f1',F,fsymbols),('index_r1',R,rsymbols)]:out['reagents'][name]=s.add(name,primer,symbols=symbols,features=[(0,len(m.P5)if name.startswith('index_f')else len(m.P7),'illumina.p5'if name.startswith('index_f')else'illumina.p7'),(len(m.P5)if name.startswith('index_f')else len(m.P7),len(primer)-(len(nx.S5)if name.startswith('index_f')else len(nx.S7)),'index.i5'if name.startswith('index_f')else'index.i7'),(len(primer)-(len(nx.S5)if name.startswith('index_f')else len(nx.S7)),len(primer),'nextera.s5'if name.startswith('index_f')else'nextera.s7')],source=source)
    out['index_before']=s.snapshot()
    out['index_melt']=s.unpair(out['repaired_top_ids'],source='index-PCR denaturation')
    def amplify(entries,key):
        owners={}
        for name,template,binding in entries:
            out[name+'_bind']=s.anneal(name,template,binding=binding,expand=True,source=source)
            owners[s.tail(name)]=name
            owners[s.tail(template)]=name+'_parent'
        operations=s.extend_available(pol,[name for name,_,_ in entries],source=source)
        out[key+'_operations']=operations;out[key+'_after']=s.snapshot()
        for op in operations:
            # Names are optional compatibility aliases. All discovered operations
            # remain in the reaction list consumed by the projector.
            name=owners.get(op.details['starting_3prime'])
            if name is not None:out[name+'_op']=op
        for name,template,_ in entries:
            out[name+'_ids']=s.chain(name);out[name]=s.sequence(name)
            out[name+'_parent_ids']=s.chain(template)
    amplify([('index_f1',out['repaired_bottom_ids'][0],len(nx.S5)),('index_r1',out['repaired_top_ids'][0],len(nx.S7))],'index')
    for name,primer,symbols in [('index_f2',F,fsymbols),('index_r2',R,rsymbols)]:out['reagents'][name]=s.add(name,primer,symbols=symbols,features=[(0,len(m.P5)if name.startswith('index_f')else len(m.P7),'illumina.p5'if name.startswith('index_f')else'illumina.p7'),(len(m.P5)if name.startswith('index_f')else len(m.P7),len(primer)-(len(nx.S5)if name.startswith('index_f')else len(nx.S7)),'index.i5'if name.startswith('index_f')else'index.i7'),(len(primer)-(len(nx.S5)if name.startswith('index_f')else len(nx.S7)),len(primer),'nextera.s5'if name.startswith('index_f')else'nextera.s7')],source=source)
    out['index_copy_before']=s.snapshot()
    out['index_f1_melt']=s.unpair(out['index_f1_ids'],source='subsequent PCR denaturation')
    out['index_r1_melt']=s.unpair(out['index_r1_ids'],source='subsequent PCR denaturation')
    amplify([('index_f2','index_r1',len(nx.S5)),('index_r2','index_f1',len(nx.S7))],'index_copy')
    # Continue the same operation over the source's five enrichment cycles.
    # Whether a tail can pair and whether a duplex becomes complete are state
    # queries, not cycle-number rules or expected-product length patches.
    out['pcr_cycles']=[]
    for cycle in range(3,6):
        for sense,primer,symbols in [('f',F,fsymbols),('r',R,rsymbols)]:
            name=f'index_{sense}{cycle}';out['reagents'][name]=s.add(name,primer,symbols=symbols,features=[(0,len(m.P5)if name.startswith('index_f')else len(m.P7),'illumina.p5'if name.startswith('index_f')else'illumina.p7'),(len(m.P5)if name.startswith('index_f')else len(m.P7),len(primer)-(len(nx.S5)if name.startswith('index_f')else len(nx.S7)),'index.i5'if name.startswith('index_f')else'index.i7'),(len(primer)-(len(nx.S5)if name.startswith('index_f')else len(nx.S7)),len(primer),'nextera.s5'if name.startswith('index_f')else'nextera.s7')],source=source)
        if cycle==3:out['index_duplex_before']=s.snapshot()
        for old in [f'index_f{cycle-1}',f'index_r{cycle-1}']:
            out[old+'_melt']=s.unpair(s.chain(old),source='PCR denaturation')
        entries=[(f'index_f{cycle}',f'index_r{cycle-1}',len(nx.S5)),(f'index_r{cycle}',f'index_f{cycle-1}',len(nx.S7))]
        amplify(entries,'index_cycle'+str(cycle))
        for sense,template,binding in [('f',f'index_r{cycle-1}',len(nx.S5)),('r',f'index_f{cycle-1}',len(nx.S7))]:
            name=f'index_{sense}{cycle}'
            anneal=out[name+'_bind']
            out[name+'_anneal']=anneal
        out['pcr_cycles'].append({'cycle':cycle,'f_binding':out[f'index_f{cycle}_anneal'].details['binding'],'r_binding':out[f'index_r{cycle}_anneal'].details['binding']})
    # This pair is actually a duplex in the final simulated state: an existing
    # full-length parent and the new complementary strand copied from it.
    final_plus_ids=s.chain('index_f4');final_minus_ids=s.chain('index_r5')
    if any(s.paired.get(a)!=b for a,b in zip(final_plus_ids,reversed(final_minus_ids))):raise ValueError('selected final state is not a duplex')
    out['final_plus_ids']=final_plus_ids;out['final_minus_ids']=final_minus_ids
    out['final_plus']=s.sequence('index_f4');out['final_minus']=s.sequence('index_r5')
    # This sequence is measured from retained target atoms, not provided as a
    # desired insert or built by substituting into an expected product.
    target_set=set(top)
    chosen=[i for i in out['repaired_top_ids']if i in target_set or s.atoms[i].template in set(bottom)]
    out['insert_ids']=tuple(chosen);out['insert']=seq(s.snapshot(),chosen)
    out['final_state']=s.snapshot()
    out['trace']=s.trace();return out
