"""Explicit versioned molecular handoffs, separate from method ancestry.

An import copies an actual produced state from a source-owned graph. It does not
establish that upstream chemistry has been validated or recover undisclosed bases.
"""
from copy import deepcopy
from dataclasses import dataclass,fields,is_dataclass
from enum import Enum
from pathlib import PurePosixPath
from collections.abc import Mapping
from urllib.parse import urlsplit
import re
from protocol_graph import Material,MaterialKind


def _signature(value):
    """Strict declaration identity, including uncertainty, features and geometry."""
    if value is None or isinstance(value,(str,int,float,bool)):return value
    if isinstance(value,Enum):return value.value
    if isinstance(value,Mapping):return tuple(sorted((str(k),_signature(v))for k,v in value.items()))
    if isinstance(value,(tuple,list)):return tuple(_signature(x)for x in value)
    if isinstance(value,(set,frozenset)):return tuple(sorted(_signature(x)for x in value))
    if is_dataclass(value):return(type(value).__name__,tuple((f.name,_signature(getattr(value,f.name)))for f in fields(value)))
    if hasattr(value,'__dict__'):return(type(value).__name__,_signature(vars(value)))
    raise TypeError('Unsupported exported material state: '+type(value).__name__)


def material_signature(material):
    return _signature((material.kind,material.state,material.chemical_state,material.inferred))


def _known_strands(material):
    from chemdraw import Construct,Scene,revcomp
    state=material.state
    def sequence(segments):return ''.join('?'*max(1,len(s.top))if s.placeholder or s.is_role_token()else s.top.replace('U','T')for s in segments)
    if isinstance(state,Scene):return tuple(sequence(s.segs)for s in state.strands.values())
    if isinstance(state,Construct):
        top=sequence(state)
        if isinstance(material.chemical_state,Mapping)and material.chemical_state.get('strandedness')=='duplex_pattern':
            # Unknown windows remain barriers, never wildcard-compatible sites.
            bottom=''.join('?'if base=='?'else revcomp(base)for base in top[::-1])
            return top,bottom
        return (top,)
    return ()


def _check_sites(material,required_sites):
    from chemdraw import revcomp
    strands=_known_strands(material)
    for name,primer in required_sites:
        if not primer or set(primer.upper())-set('ACGT'):raise ValueError('required primer site must have explicit DNA bases: '+name)
        if not any(revcomp(primer.upper())in strand for strand in strands):raise ValueError('upstream material lacks a known template site for '+name)


@dataclass(frozen=True)
class ExportReference:
    protocol: str
    version: str
    name: str
    material: str
    producer: str
    page: str
    source: str
    limitation: str

    @property
    def href(self):return '../'+self.protocol+'/'+self.page+'#export-'+self.name


@dataclass(frozen=True)
class MaterialExport:
    reference: ExportReference
    graph: object
    signature: object

    @property
    def material(self):return self.graph.materials[self.reference.material]


@dataclass(frozen=True)
class MaterialImport:
    export: MaterialExport
    expected_protocol: str
    expected_version: str
    evidence: str
    required_sites: tuple = ()

    @property
    def reference(self):return self.export.reference


def export_material(graph,name,material,*,protocol,version,page,source,limitation):
    if not re.fullmatch(r'[a-z][a-z0-9_]*',name):raise ValueError('export needs a stable identifier')
    if name in graph.exports:raise ValueError('duplicate material export '+name)
    if not protocol or PurePosixPath(protocol).name!=protocol or protocol in {'.','..'}:raise ValueError('export protocol must be one directory')
    if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*\.html',page):raise ValueError('export page must be a plain HTML filename')
    if not version.strip() or not source.strip() or not limitation.strip():raise ValueError('export needs exact version, evidence and validation scope')
    identifier=material.id if isinstance(material,Material)else material
    if identifier not in graph.materials:raise ValueError('export material does not exist')
    value=graph.materials[identifier];producer=graph.producer_of(identifier)
    if producer is None:raise ValueError('export requires a produced material, not an unprocessed input')
    if value.kind not in {MaterialKind.MOLECULAR,MaterialKind.LIBRARY}:raise ValueError('export must be molecular material')
    if value.state is None and value.chemical_state is None:raise ValueError('export needs a molecular declaration')
    graph.require_publishable()
    ref=ExportReference(protocol,version,name,identifier,producer.id,page,source,limitation)
    result=MaterialExport(ref,graph,material_signature(value));graph.exports[name]=result;return result


def import_material(graph,identifier,export,*,expected_protocol,expected_version,source,label=None,primary=False,requirements=None,required_sites=None):
    ref=export.reference
    if expected_protocol!=ref.protocol or expected_version!=ref.version:raise ValueError('upstream protocol/version does not match the required material')
    if not source.strip():raise ValueError('material handoff needs downstream source evidence')
    export.graph.require_publishable()
    value=export.material
    if export.graph.exports.get(ref.name)is not export or material_signature(value)!=export.signature:raise ValueError('upstream export changed after declaration')
    for key,want in (requirements or {}).items():
        if not isinstance(value.chemical_state,Mapping) or key not in value.chemical_state or value.chemical_state[key]!=want:raise ValueError('upstream material does not satisfy '+key)
    sites=tuple((required_sites or {}).items());_check_sites(value,sites)
    material=Material(identifier,label or value.label,value.kind,state=deepcopy(value.state),chemical_state=deepcopy(value.chemical_state),primary=primary,source=source,inferred=value.inferred)
    graph.add_input(material);graph.imports[identifier]=MaterialImport(export,expected_protocol,expected_version,source,sites)
    return material


def handoff_problems(graph,stack):
    problems=[]
    for name,export in graph.exports.items():
        if export.graph is not graph or export.reference.name!=name:problems.append('invalid material export '+name);continue
        if material_signature(export.material)!=export.signature:problems.append('exported material changed: '+name)
    for key,binding in graph.imports.items():
        ref=binding.reference;export=binding.export
        if key not in graph.materials or graph.producer_of(key)is not None:problems.append('import is not an external input: '+key);continue
        if ref.protocol!=binding.expected_protocol or ref.version!=binding.expected_version:problems.append('material handoff version mismatch: '+key)
        if export.graph.exports.get(ref.name)is not export:problems.append('material export is no longer declared: '+key)
        if material_signature(graph.materials[key])!=export.signature:problems.append('imported molecular declaration changed: '+key)
        if material_signature(export.material)!=export.signature:problems.append('upstream molecular declaration changed: '+key)
        try:_check_sites(export.material,binding.required_sites)
        except ValueError as error:problems.append(str(error))
        for problem in export.graph.publication_problems(_handoff_stack=stack):problems.append('upstream '+ref.protocol+': '+problem)
    return problems


@dataclass(frozen=True)
class SourceNamedInput:
    """A located upstream-method claim, without importing an unverified state."""
    protocol: str
    stage: str
    page: str
    source: str
    limitation: str
    version: str = ''
    title: str = ''
    url: str = ''

    @property
    def href(self):return self.url or '../'+self.protocol+'/'+self.page


def declare_input_origin(graph,material,*,stage,source,limitation,protocol='',page='',version='',title='',url=''):
    identifier=material.id if isinstance(material,Material)else material
    if identifier not in graph.materials or graph.producer_of(identifier)is not None:raise ValueError('upstream reference must target an external material')
    if graph.materials[identifier].kind not in {MaterialKind.MOLECULAR,MaterialKind.LIBRARY}:raise ValueError('upstream molecular reference cannot target a reagent or biological population')
    if url:
        target=urlsplit(url)
        if protocol or page:raise ValueError('upstream reference must choose a local model or an external original')
        if not title.strip() or target.scheme!='https' or not target.hostname or target.username or target.password or any(c.isspace()for c in url):raise ValueError('external upstream reference needs a title and an HTTPS original URL')
    else:
        if not protocol or PurePosixPath(protocol).name!=protocol or protocol in {'.','..'}:raise ValueError('upstream reference needs one protocol directory')
        if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*\.html',page):raise ValueError('upstream reference needs a plain page filename')
    if not stage.strip() or not source.strip() or not limitation.strip():raise ValueError('upstream reference needs stage, evidence and explicit validation limits')
    if identifier in graph.imports or identifier in graph.input_origins:raise ValueError('material already has upstream provenance')
    result=SourceNamedInput(protocol,stage,page,source,limitation,version,title,url);graph.input_origins[identifier]=result;return result
