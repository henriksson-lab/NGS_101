#!/usr/bin/env python3
"""Opt-in real-browser compatibility/latency benchmark for the annealing page.

No packages are installed. Uses a temporary browser profile, localhost static server
and a minimal Chrome DevTools transport. The site must already be built; no runtime packages are installed by this tool.
"""
import argparse
from pathlib import Path
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--site','--prototype',dest='prototype',type=Path,required=True)
parser.add_argument('--browser',default='google-chrome')
parser.add_argument('--report',type=Path,default=Path('/tmp/chem-anneal-browser-report.json'))
parser.add_argument('--screenshot',type=Path,default=Path('/tmp/chem-anneal-browser.png'))
parser.add_argument('--url-prefix', default='', help='Serve under a prefix such as /NGS_101 to verify subdirectory hosting.')
args=parser.parse_args()
prefix='/' + args.url_prefix.strip('/') if args.url_prefix.strip('/') else ''
if not (args.prototype/'anneal.html').is_file():parser.error('Build the site or standalone prototype first.')
import base64,functools,http.server,json,os,socket,struct,subprocess,tempfile,threading,time,urllib.request
from pathlib import Path
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
class Quiet(http.server.SimpleHTTPRequestHandler):
 def log_message(self,*args):pass
 def do_GET(self):
  if prefix:
   if not self.path.startswith(prefix+'/'):self.send_error(404);return
   self.path=self.path[len(prefix):]
  if self.path == '/assets/anneal_python.json' and getattr(self.server,'fail_bundle_once',False):
   self.server.fail_bundle_once=False;self.send_error(503,'Intentional transient prototype-test failure');return
  super().do_GET()
import xml.etree.ElementTree as ET
handler=functools.partial(Quiet,directory=str(args.prototype))
server=http.server.ThreadingHTTPServer(('127.0.0.1',0),handler);threading.Thread(target=server.serve_forever,daemon=True).start()
with tempfile.TemporaryDirectory(prefix='chem-anneal-e2e-') as profile:
 browser=subprocess.Popen([args.browser,'--headless','--no-sandbox','--disable-gpu','--disable-dev-shm-usage','--user-data-dir='+profile,'--remote-debugging-port=0','--window-size=1440,1200','about:blank'],stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
 try:
  deadline=time.monotonic()+30;portfile=Path(profile)/'DevToolsActivePort'
  while not portfile.exists():
   if browser.poll() is not None:raise RuntimeError('Chrome exited')
   if time.monotonic()>deadline:raise TimeoutError('No debugger port')
   time.sleep(.05)
  port=int(portfile.read_text().splitlines()[0]);base='http://127.0.0.1:'+str(port)
  tab=json.load(urllib.request.urlopen(urllib.request.Request(base+'/json/new?about:blank',method='PUT')));ws=WS(tab['webSocketDebuggerUrl'])
  ws.call('Runtime.enable');ws.call('Page.enable')
  ws.call('Page.addScriptToEvaluateOnNewDocument',{'source':"""
window.__workerLog=[];window.__beats=0;window.__maxGap=0;let last=performance.now();setInterval(()=>{const now=performance.now();window.__beats++;window.__maxGap=Math.max(window.__maxGap,now-last);last=now;},25);
const NativeWorker=window.Worker;window.Worker=class extends NativeWorker{constructor(...args){super(...args);this.addEventListener('message',({data})=>{if(data.timing)window.__workerLog.push(data.timing);});}};
"""})
  def evaluate(js):
   reply=ws.call('Runtime.evaluate',{'expression':js,'returnByValue':True})
   if 'exceptionDetails' in reply.get('result',{}):raise RuntimeError(reply)
   return reply.get('result',{}).get('result',{}).get('value')
  def wait(js):
   deadline=time.monotonic()+40
   while not evaluate(js):
    if time.monotonic()>deadline:raise TimeoutError(evaluate("document.getElementById('anneal-status')?.textContent"))
    time.sleep(.05)
  ws.call('Page.navigate',{'url':'http://127.0.0.1:'+str(server.server_port)+prefix+'/anneal.html'})
  wait("!!document.getElementById('anneal-example') && document.readyState==='complete'")
  # Cancel a calculation before startup, then recover without a page reload.
  evaluate("document.getElementById('anneal-example').click();document.getElementById('anneal-form').requestSubmit();document.getElementById('anneal-cancel').click();true")
  wait("document.getElementById('anneal-status').textContent==='Cancelled.' && !document.getElementById('anneal-search').disabled")
  # A transient source-bundle failure must not poison the worker's cached startup.
  server.fail_bundle_once=True
  evaluate("document.getElementById('anneal-form').requestSubmit();true")
  wait("document.getElementById('anneal-status').textContent.includes('The Python model could not be loaded')")
  evaluate("document.getElementById('anneal-form').requestSubmit();true")
  wait("!!document.querySelector('#anneal-diagram svg')")
  first=json.loads(evaluate("JSON.stringify({count:document.getElementById('anneal-count').textContent,summary:document.getElementById('anneal-summary').textContent,timings:window.__workerLog,beats:window.__beats,maxGap:window.__maxGap,svg:document.querySelector('#anneal-diagram svg').outerHTML})"))
  assert 'arm bases' in first['summary'],first
  canvas=json.loads(evaluate("JSON.stringify((()=>{const s=document.querySelector('#anneal-diagram svg'),b=s.getBBox(),v=s.viewBox.baseVal;return {width:s.width.baseVal.value,height:s.height.baseVal.value,contentWidth:b.width,contentHeight:b.height,left:b.x-v.x,top:b.y-v.y,right:v.x+v.width-b.x-b.width,bottom:v.y+v.height-b.y-b.height};})())"))
  assert canvas['height'] < 250 and canvas['width'] < 620,canvas
  assert all(0 < canvas[edge] <= 12 for edge in ('left','top','right','bottom')),canvas
  assert evaluate("document.getElementById('anneal-coordinate-text').textContent.includes('Target gap: 11–14 (4 bases)')")
  if evaluate("!!document.querySelector('nav.site-nav')"):
   assert evaluate("[...document.querySelectorAll('nav.site-nav a')].every(a=>new URL(a.href).pathname.startsWith("+json.dumps(prefix+'/')+"))")
  assert ET.fromstring(first['svg']).tag=='{http://www.w3.org/2000/svg}svg'
  evaluate("const select=document.getElementById('anneal-candidate');select.value=[...select.options].find(o=>o.textContent.startsWith('Linear')).value;select.dispatchEvent(new Event('change'));true")
  wait("document.getElementById('anneal-summary').textContent.includes('exact base pairs')")
  warm=json.loads(evaluate("JSON.stringify(window.__workerLog)"))
  # Unknown backbone bases must remain unpaired in the shared Python topology.
  evaluate("document.getElementById('anneal-example').click();const probe=document.getElementById('anneal-second');probe.value=probe.value.slice(0,8)+'NNNN'+probe.value.slice(-8);probe.dispatchEvent(new Event('input'));document.getElementById('anneal-form').requestSubmit();true")
  wait("document.getElementById('anneal-summary').textContent.includes('4 backbone bases')")
  # Editing between queries retains the runtime while invalidating the previous output.
  evaluate("document.getElementById('anneal-first').value='AC-GT';document.getElementById('anneal-first').dispatchEvent(new Event('input'));document.getElementById('anneal-form').requestSubmit();true")
  wait("document.getElementById('anneal-status').textContent.includes('Use A/C/G/T')")
  assert evaluate("!document.getElementById('anneal-status').textContent.includes('Traceback')")
  evaluate("document.getElementById('anneal-mode').value='linear';document.getElementById('anneal-example').click();document.getElementById('anneal-form').requestSubmit();true")
  wait("document.getElementById('anneal-summary').textContent==='12 exact base pairs'")
  # Exercise real SVG download via the UI, intercepting the local link only.
  evaluate("window.__download='';window.__exportText='';const createURL=URL.createObjectURL;URL.createObjectURL=function(blob){blob.text().then(text=>window.__exportText=text);return createURL.call(this,blob);};const click=HTMLAnchorElement.prototype.click;HTMLAnchorElement.prototype.click=function(){window.__download=this.download;};document.getElementById('anneal-download').click();HTMLAnchorElement.prototype.click=click;URL.createObjectURL=createURL;true")
  assert evaluate('window.__download')=='annealing.svg'
  wait("window.__exportText.includes('viewBox')")
  exported=ET.fromstring(evaluate('window.__exportText'))
  visible=json.loads(evaluate("JSON.stringify([...document.querySelector('#anneal-diagram svg').attributes].filter(a=>['width','height','viewBox'].includes(a.name)).map(a=>[a.name,a.value]))"))
  assert all(exported.attrib.get(name)==value for name,value in visible),'Export canvas differs from the compact display'
  report={'browser':json.load(urllib.request.urlopen(base+'/json/version'))['Browser'],'cold':{k:v for k,v in first.items() if k!='svg'},'warm_timings':warm,'final_summary':evaluate("document.getElementById('anneal-summary').textContent"),'svg_download':True,'invalid_input':True,'unknown_backbone':True,'cancellation':True,'runtime_retry':True,'url_prefix':prefix,'compact_canvas':canvas}
  args.report.write_text(json.dumps(report,indent=2));print(json.dumps(report))
  shot=ws.call('Page.captureScreenshot',{'format':'png'})['result']['data'];args.screenshot.write_bytes(base64.b64decode(shot))
 finally:
  browser.terminate();browser.wait(timeout=5);server.shutdown()
