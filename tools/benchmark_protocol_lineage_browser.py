#!/usr/bin/env python3
"""Optional Chrome checks for the full relationship graph and local view."""
import argparse,base64,functools,http.server,json,os,socket,struct,subprocess,tempfile,threading,time,urllib.request
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--site',type=Path,default=ROOT/'_site')
parser.add_argument('--report',type=Path,default=Path('/tmp/protocol-lineage-browser.json'))
parser.add_argument('--screenshot',type=Path)
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
  self.path=self.path[len('/NGS_101'):]
  if self.path=='/assets/protocol_lineage.json' and self.server.fail_once:
   self.server.fail_once=False;self.send_error(503);return
  if self.path=='/lineage-test/local.html':
   import sys
   sys.path.insert(0,str(ROOT/'lib'))
   from protocol_lineage import local_view
   payload=json.loads((args.site/'assets/protocol_lineage.json').read_text())
   self.send_response(200);self.send_header('Content-Type','text/html; charset=utf-8');self.end_headers()
   self.wfile.write(('<!doctype html><h1>Local view fixture</h1>'+local_view('share-seq__10.1016+j.cell.2020.09.056',payload)).encode());return
  super().do_GET()
server=http.server.ThreadingHTTPServer(('127.0.0.1',0),functools.partial(Handler,directory=str(args.site.resolve())));server.fail_once=True
threading.Thread(target=server.serve_forever,daemon=True).start();checks=[]
with tempfile.TemporaryDirectory(prefix='chem-lineage-browser-') as profile:
 browser=subprocess.Popen(['google-chrome','--headless','--no-sandbox','--disable-gpu','--disable-dev-shm-usage','--user-data-dir='+profile,'--remote-debugging-port=0','--window-size=1600,1100','about:blank'],stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
 try:
  portfile=Path(profile)/'DevToolsActivePort';deadline=time.monotonic()+25
  while not portfile.exists():
   if time.monotonic()>deadline:raise TimeoutError('Chrome startup')
   time.sleep(.05)
  base='http://127.0.0.1:'+portfile.read_text().splitlines()[0]
  tab=json.load(urllib.request.urlopen(urllib.request.Request(base+'/json/new?about:blank',method='PUT')));ws=WS(tab['webSocketDebuggerUrl']);ws.call('Runtime.enable');ws.call('Page.enable')
  ws.call('Page.addScriptToEvaluateOnNewDocument',{'source':"window.__errors=[];window.addEventListener('error',e=>window.__errors.push(e.message));window.addEventListener('unhandledrejection',e=>window.__errors.push(String(e.reason)));window.__blobs=[];const create=URL.createObjectURL.bind(URL);URL.createObjectURL=b=>{window.__blobs.push(b);return create(b)};const click=HTMLAnchorElement.prototype.click;HTMLAnchorElement.prototype.click=function(){if(this.download)return;return click.call(this)};"})
  def ev(js,wait=False):
   r=ws.call('Runtime.evaluate',{'expression':js,'returnByValue':True,'awaitPromise':wait})
   if 'exceptionDetails' in r.get('result',{}):raise RuntimeError(r)
   return r.get('result',{}).get('result',{}).get('value')
  def wait(js):
   deadline=time.monotonic()+20
   while not ev(js):
    if time.monotonic()>deadline:raise TimeoutError(js+'; state='+str(ev("({url:location.href,status:document.getElementById('lineage-status')?.textContent,errors:window.__errors,ready:document.readyState})")))
    time.sleep(.04)
  def check(label,js):
   assert ev(js),label+' '+str(ev('window.__errors'));checks.append(label)
  def click(name):ev("document.getElementById('lineage-"+name+"').click();true")
  def toggle(name,value):ev("(()=>{const e=document.getElementById('lineage-"+name+"');e.checked="+json.dumps(value)+";e.dispatchEvent(new Event('change',{bubbles:true}));return true})()")
  page='http://127.0.0.1:'+str(server.server_port)+'/NGS_101/lineage.html'
  ws.call('Page.navigate',{'url':page});wait("document.getElementById('lineage-status')?.textContent.includes('Reload to retry')")
  checks.append('data loading failure is explicit and retryable by reload')
  ws.call('Page.reload',{});wait("document.querySelectorAll('.lineage-node').length>0")
  check('shared navigation works under Pages prefix',"!!document.querySelector('nav a[href=\"lineage.html\"]') && [...document.querySelectorAll('nav a')].every(a=>new URL(a.href).pathname.startsWith('/NGS_101/'))")
  initial=ev("document.querySelectorAll('[data-relation]').length")
  check('curated graph defaults to documented and inferred ancestry only',"document.querySelectorAll('[data-relation]').length>=8 && document.querySelectorAll('[data-similarity]').length===0")
  check('graph exposes multiple SHARE-seq parents',"document.querySelectorAll('[data-protocol=\"share-seq__10.1016+j.cell.2020.09.056\"]').length===1")
  if args.screenshot:args.screenshot.write_bytes(base64.b64decode(ws.call('Page.captureScreenshot',{'format':'png'})['result']['data']))
  toggle('inferred',False);check('inference filter removes only inferred links',"document.querySelectorAll('[data-relation]').length==="+str(initial-1));toggle('inferred',True)
  ev("document.querySelector('[data-relation]').dispatchEvent(new MouseEvent('click',{bubbles:true}));true")
  check('edge selection shows evidence, scope and limits',"document.getElementById('lineage-detail').textContent.includes('Scope') && document.getElementById('lineage-detail').textContent.includes('Limits') && !!document.querySelector('#lineage-detail a[href^=\"https://doi.org/\"]')")
  ev("document.getElementById('lineage-search').value='SHARE-seq';true");click('find')
  check('search focuses selected protocol and both source parents',"document.querySelectorAll('.lineage-node').length===3 && document.getElementById('lineage-detail').textContent.includes('SPLiT-seq') && document.getElementById('lineage-detail').textContent.includes('Paired-seq')")
  before=ev("document.querySelector('#lineage-canvas svg>g').getAttribute('transform')");click('plus')
  check('zoom controls alter SVG transform',"document.querySelector('#lineage-canvas svg>g').getAttribute('transform')!=="+json.dumps(before))
  ev("document.getElementById('lineage-canvas').focus();true");before=ev("document.querySelector('#lineage-canvas svg>g').getAttribute('transform')")
  ws.call('Input.dispatchKeyEvent',{'type':'keyDown','key':'ArrowRight'});ws.call('Input.dispatchKeyEvent',{'type':'keyUp','key':'ArrowRight'})
  check('graph supports keyboard panning',"document.querySelector('#lineage-canvas svg>g').getAttribute('transform')!=="+json.dumps(before))
  click('fit');click('export');export=ev('window.__blobs.at(-1).text()',True)
  import xml.etree.ElementTree as ET
  assert ET.fromstring(export).tag=='{http://www.w3.org/2000/svg}svg';checks.append('export is a standalone valid SVG')
  click('all');toggle('similar',True)
  check('similarity remains a separate optional undirected layer',"document.querySelectorAll('[data-similarity]').length>0 && document.querySelectorAll('[data-relation]').length==="+str(initial))
  ev("document.querySelector('[data-similarity]').dispatchEvent(new MouseEvent('click',{bubbles:true}));true")
  check('similarity selection shows matched sequences and non-ancestry scope',"document.getElementById('lineage-detail').textContent.includes('not a curated ancestry link') && !!document.querySelector('#lineage-detail code')")
  toggle('similar',False);toggle('isolated',True)
  check('unlinked public protocols can be included without inventing edges',"document.querySelectorAll('.lineage-node').length>100 && document.querySelectorAll('[data-relation]').length==="+str(initial))
  ws.call('Page.navigate',{'url':page.replace('lineage.html','lineage-test/local.html')});wait("!!document.querySelector('.protocol-lineage-local')")
  ev("document.querySelector('details').open=true;true")
  check('local graph renders source-qualified multi-parent view',"document.querySelectorAll('.protocol-lineage-local svg a').length===3 && document.body.textContent.includes('Solid: source-documented')")
  href=ev("document.querySelector('.protocol-lineage-local p a').href");ws.call('Page.navigate',{'url':href});wait("document.querySelectorAll('.lineage-node').length===3")
  check('local link opens focused full graph',"document.getElementById('lineage-detail').textContent.includes('SHARE-seq') && location.search.includes('protocol=share-seq')")
  check('no uncaught browser errors',"window.__errors.length===0")
  ws.call('Emulation.setDeviceMetricsOverride',{'width':390,'height':844,'deviceScaleFactor':1,'mobile':True});check('mobile graph layout avoids document overflow',"document.documentElement.scrollWidth<=window.innerWidth+1")
  args.report.write_text(json.dumps({'passed':len(checks),'checks':checks},indent=2)+'\n');print(args.report.read_text())
 finally:
  browser.terminate();browser.wait(timeout=10);server.shutdown()
