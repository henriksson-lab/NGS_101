/* Browser host for the existing Python model. No chemistry in JavaScript. */
let runtimePromise;
let queue = Promise.resolve();
async function runtime() {
  if (!runtimePromise) runtimePromise = (async () => {
    const start = performance.now();
    const response = await fetch(new URL('./anneal_python.json',import.meta.url));
    if (!response.ok) throw new Error('The Python model could not be loaded.');
    const bundle = await response.json();
    if (!/^[0-9]+\.[0-9]+\.[0-9]+$/.test(bundle.runtime_version) || !['cdn','local'].includes(bundle.runtime_source)) throw new Error('Unsupported browser runtime configuration.');
    const indexURL = bundle.runtime_source === 'cdn'
      ? 'https://cdn.jsdelivr.net/pyodide/v'+bundle.runtime_version+'/full/'
      : new URL('./pyodide/',import.meta.url).href;
    let py;
    try {
      const {loadPyodide} = await import(new URL('pyodide.mjs',indexURL).href);
      py = await loadPyodide({indexURL,enableRunUntilComplete:false});
    } catch (error) {
      throw new Error('The pairing engine could not load. Check your connection and try again.');
    }
    const runtimeMs = performance.now()-start;
    if (py.version !== bundle.runtime_version) throw new Error('The browser runtime version does not match its model bundle.');
    for (const name of ['protocol_graph','chemdraw','topology','annealing']) {
      if (typeof bundle.modules[name] !== 'string') throw new Error('Incomplete Python model bundle.');
      py.FS.writeFile('/home/pyodide/'+name+'.py',bundle.modules[name]);
    }
    py.runPython('import annealing, json');
    return {py,timing:{runtimeMs,modelMs:performance.now()-start-runtimeMs}};
  })().catch(error => { runtimePromise = undefined; throw error; });
  return runtimePromise;
}
async function execute(message) {
  const {id,operation,payload={}} = message;
  let started = false;
  try {
    const {py,timing} = await runtime(), start = performance.now();
    started = true;
    if (operation === 'ready') { self.postMessage({id,result:timing}); return; }
    if (!['search','render'].includes(operation)) throw new Error('Unknown annealing operation.');
    py.globals.set('_anneal_input',JSON.stringify({operation,payload}));
    const result = py.runPython(`
_anneal_request = json.loads(_anneal_input)
_anneal_payload = _anneal_request['payload']
try:
    if _anneal_request['operation'] == 'search':
        _anneal_result = annealing.search(_anneal_payload['first'], _anneal_payload['second'],
                                        _anneal_payload.get('mode', 'both'), _anneal_payload.get('minimum', 8))
    else:
        _anneal_result = annealing.render(_anneal_payload['first'], _anneal_payload['second'], _anneal_payload['candidate'])
except (ValueError, KeyError, TypeError) as _anneal_error:
    _anneal_result = {'_anneal_error': str(_anneal_error)}
json.dumps(_anneal_result)
`);
    const parsed = JSON.parse(result);
    if (parsed._anneal_error) throw new Error(parsed._anneal_error);
    self.postMessage({id,result:parsed,timing:{...timing,operationMs:performance.now()-start}});
  } catch (error) { self.postMessage({id,error:String(error.message || error),startupError:!started}); }
}
self.onmessage = ({data}) => { queue = queue.then(()=>execute(data)).catch(error=>self.postMessage({id:data.id,error:String(error)})); };
