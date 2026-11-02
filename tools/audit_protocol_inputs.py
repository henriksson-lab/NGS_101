#!/usr/bin/env python3
"""Inventory prepared molecular inputs and verified cross-protocol handoffs.

The text filter only selects candidates for source review. A library/cDNA label
never supplies an upstream producer or authorizes inheriting its chemistry.
"""
import argparse,ast,csv,json,re,subprocess,sys
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/'lib'),str(ROOT/'tools'),str(ROOT)]


def worker(name):
    from selftest_source_oligos import graph_module_paths,load
    from protocol_graph import MaterialKind
    import contextlib,io,re
    from dataclasses import asdict
    paths=graph_module_paths(ROOT/name)
    if len(paths)!=1:return {'protocol':name,'status':'no_unique_graph','inputs':[]}
    with contextlib.redirect_stdout(io.StringIO()):
        module=load(paths[0]);sys.path.insert(0,str(paths[0].parent));g=module.protocol_graph()
    rows=[]
    for key,m in g.materials.items():
        if g.producer_of(key)is not None or m.kind not in {MaterialKind.MOLECULAR,MaterialKind.LIBRARY,MaterialKind.BIOLOGICAL}:continue
        binding=g.imports.get(key)
        origin=g.input_origins.get(key)
        if not binding and not origin and not re.search(r'(?i)\b(?:cDNA|librar(?:y|ies)|amplicon|tagmented|pre.?amplified)\b',m.label):continue
        pooled=bool(re.search(r'(?i)(?:other|independent|PhiX|calibrat|quantif|standard|sequencing.?pool|other.?run)',key+' '+m.label))
        row={'material':key,'label':m.label,'kind':m.kind.value,'primary':m.primary,'source':m.source,
             'status':'validated_handoff'if binding else'source_named_handoff'if origin else'pool_or_control_review'if pooled else'needs_source_review'}
        if binding:
            from dataclasses import asdict
            row['upstream']=asdict(binding.reference)
        elif origin:row['upstream']=asdict(origin)
        rows.append(row)
    exports=[*g.exports.values(),*getattr(module,'ADDITIONAL_MATERIAL_EXPORTS',())]
    for export in exports:
        if export.reference.protocol!=name:raise ValueError('protocol advertises another owner’s export')
        export.graph.require_publishable()
    return {'protocol':name,'status':'inspected_declaration','inputs':rows,'exports':[asdict(e.reference)for e in exports]}


def note_candidates(directory):
    """Literal input declarations only; never execute incomplete chemistry."""
    result=[]
    origins={}
    for path in sorted((ROOT/directory/'tools').glob('*.py')):
        tree=ast.parse(path.read_text(),filename=str(path))
        for call in ast.walk(tree):
            if isinstance(call,ast.Call)and isinstance(call.func,ast.Name)and call.func.id=='declare_input_origin'and len(call.args)>1:
                if isinstance(call.args[1],ast.Constant)and isinstance(call.args[1].value,str):
                    origins[(str(path.relative_to(ROOT)),call.args[1].value)]={x.arg:x.value.value for x in call.keywords if isinstance(x.value,ast.Constant)and isinstance(x.value.value,str)}
            if not isinstance(call,ast.Call)or not isinstance(call.func,ast.Attribute)or call.func.attr!='input':continue
            args={**dict(zip(('id','label'),call.args[:2])),**{x.arg:x.value for x in call.keywords if x.arg}}
            values={k:v.value for k,v in args.items()if isinstance(v,ast.Constant)and isinstance(v.value,str)}
            label=values.get('label','')
            if not re.search(r'(?i)\b(?:cDNA|librar(?:y|ies)|amplicon|tagmented|pre.?amplified)\b',label):continue
            result.append({'material':values.get('id',''),'label':label,'file':str(path.relative_to(ROOT)),'line':call.lineno,'status':'literal_declaration_needs_source_review'})
    for row in result:
        key=(row['file'],row['material'])
        if key in origins:row['declared_origin']=origins[key]
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('protocols',nargs='*');p.add_argument('--worker');p.add_argument('--output',type=Path,default=None);p.add_argument('--jobs',type=int,default=6);p.add_argument('--notes',action='store_true',help='Inventory literal inputs in notes/draft modules without executing incomplete models; separate snapshot');a=p.parse_args()
    if a.worker:print(json.dumps(worker(a.worker)));return
    with (ROOT/'catalogue/ours.tsv').open()as f:rows=list(csv.DictReader(f,delimiter='\t'))
    if a.notes:
        eligible={r['dir']for r in rows if r['status']in {'notes','draft'}}
        names=a.protocols or sorted(eligible)
        if set(names)-eligible:raise ValueError('--notes only scans notes/draft entries')
        results=[{'protocol':name,'inputs':note_candidates(name)}for name in names]
        output=a.output or ROOT/'catalogue/protocol_input_notes_candidates.json'
        output.write_text(json.dumps({'schema':1,'scope':'Static literal input declarations in private notes/draft models; not executed graphs, source review, material inheritance or public-page readiness. Nonliteral and indirect declarations are not covered.','protocols':results},indent=2)+'\n')
        print(str(len(results))+' notes/draft entries inspected for literal inputs; '+str(sum(len(r['inputs'])for r in results))+' literal input candidates; '+str(output.relative_to(ROOT)if output.is_relative_to(ROOT)else output))
        return
    eligible={r['dir']:r for r in rows if r['status']in {'documented','construct','assay'}}
    names=a.protocols or sorted(eligible)
    if set(names)-set(eligible):raise ValueError('only published schematic protocols are inventoried')
    def run(name):
        try:
            result=subprocess.run([sys.executable,__file__,'--worker',name],capture_output=True,text=True,timeout=60,cwd=ROOT)
            if result.returncode:return {'protocol':name,'status':'load_failed','error':result.stderr[-900:],'inputs':[]}
            return json.loads(result.stdout)
        except Exception as e:return {'protocol':name,'status':'load_failed','error':str(e),'inputs':[]}
    with ThreadPoolExecutor(max_workers=a.jobs)as pool:results=list(pool.map(run,names))
    output={'schema':1,'scope':'Declared molecular handoffs plus text-selected review candidates; not ancestry inference, complete source review or absence evidence. Pool/control rows are not automatically protocol extensions.','protocols':results}
    (a.output or ROOT/'catalogue/protocol_input_inventory.json').write_text(json.dumps(output,indent=2)+'\n')
    found=[(r['protocol'],x)for r in results for x in r['inputs']]
    print(f'{len(results)} public catalogue entries; '+str(sum(r['status']=='inspected_declaration'for r in results))+' graph declarations loaded; '+str(sum(r['status']=='no_unique_graph'for r in results))+' entries without a unique ProtocolGraph; '+f'{len(found)} text-selected input candidates; '+str(sum(x['status']=='validated_handoff'for _,x in found))+' validated handoffs; '+str(sum(x['status']=='source_named_handoff'for _,x in found))+' source-named handoffs; '+str(sum(r['status']=='load_failed'for r in results))+' load failures')
    for name,x in found:
        if x['status']!='pool_or_control_review':print(name,x['material'],x['status'],x['label'][:120],sep='\t')

if __name__=='__main__':main()
