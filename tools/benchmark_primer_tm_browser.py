#!/usr/bin/env python3
"""Optional actual-Chrome checks: python3 tools/benchmark_primer_tm_browser.py --site _site."""
import argparse,base64,functools,http.server,json,os,socket,struct,subprocess,tempfile,threading,time,urllib.request
from pathlib import Path
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--site',type=Path,default=Path(__file__).resolve().parents[1]/'_site')
parser.add_argument('--screenshot',type=Path)
parser.add_argument('--report',type=Path,default=Path('/tmp/primer-tm-browser.json'))
args=parser.parse_args()
class WS:
 def __init__(self,url):
  from urllib.parse import urlparse
  u=urlparse(url);self.s=socket.create_connection((u.hostname,u.port),timeout=5);self.buf=b''
  key=base64.b64encode(os.urandom(16)).decode()
  self.s.sendall(('GET '+u.path+' HTTP/1.1\r\nHost: '+u.netloc+'\r\nUpgrade: websocket\r\nConnection: Upgrade\r\nSec-WebSocket-Key: '+key+'\r\nSec-WebSocket-Version: 13\r\n\r\n').encode())
  data=b''
  while b'\r\n\r\n' not in data:data+=self.s.recv(4096)
  header,self.buf=data.split(b'\r\n\r\n',1)
  assert b' 101 ' in header,header
 def read(self,n):
  while len(self.buf)<n:
   chunk=self.s.recv(max(n-len(self.buf),4096))
   if not chunk:raise ConnectionError('Browser debugger disconnected')
   self.buf+=chunk
  result,self.buf=self.buf[:n],self.buf[n:];return result
 def recv(self):
  a,b=self.read(2);length=b&127
  if length==126:length=struct.unpack('!H',self.read(2))[0]
  elif length==127:length=struct.unpack('!Q',self.read(8))[0]
  mask=self.read(4) if b&128 else None
  data=self.read(length)
  if mask:data=bytes(v^mask[i%4] for i,v in enumerate(data))
  return json.loads(data)
 def send(self,value):
  data=json.dumps(value).encode();length=len(data);header=bytes([129,128|length]) if length<126 else bytes([129,254])+struct.pack('!H',length)
  mask=os.urandom(4);self.s.sendall(header+mask+bytes(v^mask[i%4] for i,v in enumerate(data)))
 def call(self,method,params=None):
  self.ident=getattr(self,'ident',0)+1;ident=self.ident;self.send({'id':ident,'method':method,'params':params or {}})
  while True:
   reply=self.recv()
   if reply.get('id')==ident:return reply

class Handler(http.server.SimpleHTTPRequestHandler):
 def log_message(self,*args):pass
 def do_GET(self):
  if not self.path.startswith('/NGS_101/'):
   self.send_error(404);return
  self.path=self.path[len('/NGS_101'):];super().do_GET()
server=http.server.ThreadingHTTPServer(('127.0.0.1',0),functools.partial(Handler,directory=str(args.site.resolve())))
threading.Thread(target=server.serve_forever,daemon=True).start();checks=[]
with tempfile.TemporaryDirectory(prefix='chem-tm-browser-') as profile:
 browser=subprocess.Popen(['google-chrome','--headless','--no-sandbox','--disable-gpu','--disable-dev-shm-usage','--user-data-dir='+profile,'--remote-debugging-port=0','--window-size=1440,1100','about:blank'],stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
 try:
  portfile=Path(profile)/'DevToolsActivePort';deadline=time.monotonic()+25
  while not portfile.exists():
   if time.monotonic()>deadline:raise TimeoutError('Chrome startup')
   time.sleep(.05)
  base='http://127.0.0.1:'+portfile.read_text().splitlines()[0]
  tab=json.load(urllib.request.urlopen(urllib.request.Request(base+'/json/new?about:blank',method='PUT')));ws=WS(tab['webSocketDebuggerUrl']);ws.call('Runtime.enable');ws.call('Page.enable')
  ws.call('Page.addScriptToEvaluateOnNewDocument',{'source':"window.__errors=[];window.addEventListener('error',e=>window.__errors.push(e.message));window.addEventListener('unhandledrejection',e=>window.__errors.push(String(e.reason)));window.__blobs=[];const create=URL.createObjectURL.bind(URL);URL.createObjectURL=b=>{window.__blobs.push(b);return create(b)};const click=HTMLAnchorElement.prototype.click;HTMLAnchorElement.prototype.click=function(){if(this.download)return;return click.call(this)};window.__requests=0;window.fetch=()=>{window.__requests++;throw new Error('Unexpected remote calculation')};"})
  def ev(js,wait=False):
   r=ws.call('Runtime.evaluate',{'expression':js,'returnByValue':True,'awaitPromise':wait})
   if 'exceptionDetails' in r.get('result',{}):raise RuntimeError(r)
   return r.get('result',{}).get('result',{}).get('value')
  def wait(js):
   deadline=time.monotonic()+20
   while not ev(js):
    if time.monotonic()>deadline:raise TimeoutError(js)
    time.sleep(.04)
  def check(label,js):
   assert ev(js),label+' '+str(ev('window.__errors'));checks.append(label)
  def setv(name,value):ev("(()=>{const e=document.getElementById('tm-"+name+"');e.value="+json.dumps(value)+";e.dispatchEvent(new Event('input',{bubbles:true}));return true})()")
  def click(name):ev("document.getElementById('tm-"+name+"').click();true")
  ws.call('Page.navigate',{'url':'http://127.0.0.1:'+str(server.server_port)+'/NGS_101/tm.html'});wait("document.readyState==='complete' && !!window.PrimerTm")
  check('tab navigation works under Pages prefix',"!!document.querySelector('nav a[href=\"tm.html\"]') && [...document.querySelectorAll('nav a')].every(a=>new URL(a.href).pathname.startsWith('/NGS_101/'))")
  check('Roche concentrations and monovalent assumption are explicit',"document.getElementById('tm-mg').value==='2.5' && document.getElementById('tm-dntp').value==='0.3' && document.querySelector('.tm-assumption').textContent.includes('not Roche')")
  click('example');check('example produces two primer results and delta Tm',"!document.getElementById('tm-results').hidden && document.querySelectorAll('.tm-value').length===2 && document.getElementById('tm-difference').textContent.includes('ΔTm')")
  check('per-dNTP concentration is correctly converted to total',"document.getElementById('tm-conditions-summary').textContent.includes('1.200 mM total')")
  original=ev("document.querySelector('.tm-value strong').textContent")
  if args.screenshot:
   ev("document.getElementById('tm-results').scrollIntoView({block:'nearest'});true");args.screenshot.write_bytes(base64.b64decode(ws.call('Page.captureScreenshot',{'format':'png'})['result']['data']))
  click('export');tsv=ev('window.__blobs.at(-1).text()',True);assert 'dntp_each_mM' in tsv and 'dntp_total_mM' in tsv and 'assumed_monovalent_mM' in tsv and '1.2' in tsv;checks.append('TSV includes results, concentration units and assumptions')
  setv('first','GGGGGGACGTCAGTGCATGCTAGTCA');setv('tail-first','6');check('excluded tail preserves binding Tm',"document.querySelector('.tm-value strong').textContent==="+json.dumps(original)+" && document.querySelector('.tm-value code').textContent==='ACGTCAGTGCATGCTAGTCA'")
  setv('first','NNNN');check('invalid primer clears stale results',"document.getElementById('tm-results').hidden && document.getElementById('tm-status').textContent.length>0")
  click('example');setv('second','');check('one-primer mode is supported',"document.querySelectorAll('.tm-value').length===1 && document.getElementById('tm-difference').textContent===''")
  setv('monovalent','150');check('changing assumed salt updates calculated Tm',"document.querySelector('.tm-value strong').textContent!=="+json.dumps(original))
  setv('mg','1.5');check('modified disclosed concentrations are labelled custom',"document.getElementById('tm-conditions-summary').textContent.startsWith('Custom conditions')")
  click('reset');check('reset restores preset concentrations and result',"document.getElementById('tm-mg').value==='2.5' && document.getElementById('tm-dntp').value==='0.3' && document.querySelector('.tm-value strong').textContent==="+json.dumps(original))
  setv('mg','0');setv('monovalent','0');check('no-ion conditions fail without stale output',"document.getElementById('tm-results').hidden && document.getElementById('tm-status').textContent.includes('needs')")
  click('reset');setv('dntp','-1');check('negative dNTP fails explicitly',"document.getElementById('tm-results').hidden && document.getElementById('tm-status').textContent.length>0")
  click('reset');setv('tail-first','999');check('invalid tail fails explicitly',"document.getElementById('tm-results').hidden && document.getElementById('tm-status').textContent.includes('tail')")
  click('example');setv('first','GAATTC');check('self-complementarity and short-sequence scope are reported',"document.getElementById('tm-details').textContent.includes('Self-complementary') && document.getElementById('tm-warnings').textContent.includes('length')")
  check('manufacturer Ta guidance is separate from computed Tm',"document.querySelector('.tm-method').textContent.includes('65°C') && document.querySelector('.tm-method').textContent.includes('no fixed Tm-to-Ta offset')")
  check('calculation makes no fetch requests and raises no uncaught errors',"window.__requests===0 && window.__errors.length===0")
  ws.call('Emulation.setDeviceMetricsOverride',{'width':390,'height':844,'deviceScaleFactor':1,'mobile':True});check('mobile form fits without document-wide horizontal overflow',"document.documentElement.scrollWidth<=window.innerWidth+1")
  click('clear');check('clear removes primers, output and errors',"document.getElementById('tm-first').value==='' && document.getElementById('tm-results').hidden && document.getElementById('tm-status').textContent===''")
  args.report.write_text(json.dumps({'passed':len(checks),'checks':checks},indent=2)+'\n');print(args.report.read_text())
 finally:
  browser.terminate();browser.wait(timeout=10);server.shutdown()
