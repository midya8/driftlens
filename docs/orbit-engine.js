'use strict';
(function(root){
const gcd=(a,b)=>{a=a<0n?-a:a;while(b){[a,b]=[b,a%b];}return a;};
class Q{
 constructor(n,d=1n){if(!d)throw Error('Zero denominator');if(d<0n){n=-n;d=-d;}const g=gcd(n,d);this.n=n/g;this.d=d/g;}
 static parse(x,bounded=false){if(typeof x!=='string'||x.length>(bounded?16:2500)||! /^-?(?:0|[1-9]\d*)(?:\/[1-9]\d*)?$/.test(x))throw Error('Canonical rational string required');const p=x.split('/'),q=new Q(BigInt(p[0]),p.length===1?1n:BigInt(p[1]));if(q.toString()!==x||(bounded&&((q.n<0n?-q.n:q.n)>1000000n||q.d>1000000n)))throw Error('Bounded reduced rational required');return q;}
 add(x){return new Q(this.n*x.d+x.n*this.d,this.d*x.d);} sub(x){return new Q(this.n*x.d-x.n*this.d,this.d*x.d);} mul(x){return new Q(this.n*x.n,this.d*x.d);} div(x){return new Q(this.n*x.d,this.d*x.n);} cmp(x){const d=this.n*x.d-x.n*this.d;return d<0n?-1:d>0n?1:0;} toString(){return String(this.n)+(this.d===1n?'':'/'+this.d);}
}
const q=x=>new Q(BigInt(x)),sum=a=>a.reduce((s,x)=>s.add(x),q(0)),keys=(x,w)=>x&&typeof x==='object'&&!Array.isArray(x)&&Object.keys(x).sort().join('|')===w.slice().sort().join('|');
const scope='Complete finite supplied-Gram minimum-p global randomization test. Exact rational arithmetic and all ties retained. Requires a fixed label-independent kernel family, prespecified blocks and joint invariance under allowed reallocations. PSD is finite geometry, not a characteristic population kernel or validated sampling assumption. One global null decision only, not strong component-wise inference, causality, model-quality, equivalence, automatic bandwidth choice or continuous monitoring.';
function validate(r){
 if(!keys(r,['row_names','reference','blocks','kernels','alpha']))throw Error('Request fields row_names/reference/blocks/kernels/alpha required');
 const name=x=>typeof x==='string'&&x.trim()&&[...x].length<=40&&![...x].some(c=>c.codePointAt(0)>=0xD800&&c.codePointAt(0)<=0xDFFF),names=r.row_names;
 if(!Array.isArray(names)||names.length<4||names.length>10||!names.every(name)||new Set(names).size!==names.length)throw Error('Four to ten distinct bounded row names required');
 const n=names.length,left=r.reference;if(!Array.isArray(left)||left.length<2||left.length>n-2||left.some((v,j)=>!Number.isInteger(v)||v<0||v>=n||(j&&v<=left[j-1])))throw Error('Sorted distinct integer reference indices; each side >=2');
 if(r.blocks!==null&&(!Array.isArray(r.blocks)||r.blocks.length!==n||!r.blocks.every(name)))throw Error('Null or aligned bounded block labels required');
 if(!Array.isArray(r.kernels)||r.kernels.length<1||r.kernels.length>6)throw Error('One to six supplied kernels required');const seen=new Set();
 const matrices=r.kernels.map(k=>{if(!keys(k,['name','gram'])||!name(k.name)||seen.has(k.name))throw Error('Distinct kernel names/gram required');seen.add(k.name);if(!Array.isArray(k.gram)||k.gram.length!==n||k.gram.some(row=>!Array.isArray(row)||row.length!==n))throw Error('Aligned square Gram required');const a=k.gram.map(row=>row.map(x=>Q.parse(x,true)));for(let i=0;i<n;i++)for(let j=0;j<n;j++)if(a[i][j].cmp(a[j][i]))throw Error('Exactly symmetric Gram required');return a;});
 const alpha=Q.parse(r.alpha,true);if(alpha.cmp(q(0))<=0||alpha.cmp(q(1))>=0)throw Error('Alpha strictly between zero and one');return matrices;
}
function psd(a){ // Congruence Schur reduction, not the native LDL recurrence.
 const n=a.length,m=a.map(row=>row.slice()),lower=Array.from({length:n},(_,i)=>Array.from({length:n},(_,j)=>q(i===j?1:0))),d=[];
 for(let k=0;k<n;k++){const pivot=m[k][k];if(pivot.cmp(q(0))<0)throw Error('Gram is not PSD');d.push(pivot);
  for(let i=k+1;i<n;i++){if(pivot.n===0n){if(m[i][k].n)throw Error('Zero PSD pivot has nonzero residual');}else lower[i][k]=m[i][k].div(pivot);}
  if(pivot.n)for(let i=k+1;i<n;i++)for(let j=k+1;j<n;j++)m[i][j]=m[i][j].sub(lower[i][k].mul(pivot).mul(lower[j][k]));
 }
 return {lower:lower.map(row=>row.map(String)),diagonal:d.map(String),rank:d.filter(x=>x.n>0n).length};
}
function inspect(request){
 const matrices=validate(request),r=JSON.parse(JSON.stringify(request)),n=r.row_names.length,certificates=matrices.map(psd),labels=r.blocks||Array(n).fill('pooled'),allowed=[];
 for(let mask=0;mask<(1<<n);mask++){const row=Array.from({length:n},(_,i)=>i).filter(j=>mask&(1<<j));if(row.length!==r.reference.length)continue;let good=true;for(const label of new Set(labels))if(row.filter(j=>labels[j]===label).length!==r.reference.filter(j=>labels[j]===label).length)good=false;if(good)allowed.push(row);}
 allowed.sort((a,b)=>{for(let j=0;j<a.length;j++){if(a[j]!==b[j])return a[j]-b[j];}return 0;});const total=allowed.length;
 const scores=allowed.map(left=>{const membership=new Set(left),right=Array.from({length:n},(_,j)=>j).filter(j=>!membership.has(j)),nx=left.length,ny=right.length;return matrices.map(a=>{let value=q(0);for(let i=0;i<n;i++)for(let j=0;j<n;j++){if(i===j)continue;const same=membership.has(i)===membership.has(j);const den=same?(membership.has(i)?nx*(nx-1):ny*(ny-1)):nx*ny;value=same?value.add(a[i][j].div(q(den))):value.sub(a[i][j].div(q(den)));}return value;});});
 const ranks=scores.map(row=>row.map((v,k)=>scores.filter(s=>s[k].cmp(v)>=0).length)),minima=ranks.map(row=>Math.min(...row)),alpha=Q.parse(r.alpha),family=matrices.length;
 const rows=allowed.map((left,i)=>{const raw=ranks[i].map(x=>new Q(BigInt(x),BigInt(total))),order=Array.from({length:family},(_,j)=>j).sort((a,b)=>raw[a].cmp(raw[b])||a-b),holm=Array(family);let prior=q(0);
 order.forEach((k,j)=>{let v=raw[k].mul(q(family-j));if(v.cmp(q(1))>0)v=q(1);if(v.cmp(prior)<0)v=prior;prior=v;holm[k]=v;});const tail=Array.from({length:total},(_,j)=>j).filter(j=>minima[j]<=minima[i]),globalP=new Q(BigInt(tail.length),BigInt(total));
 return {reference:left,current:Array.from({length:n},(_,j)=>j).filter(j=>!left.includes(j)),mmd2_unbiased:scores[i].map(String),tail_counts:ranks[i],raw_p:raw.map(String),minimum_tail_count:minima[i],selected_kernel_indices:ranks[i].map((v,k)=>v===minima[i]?k:-1).filter(k=>k>=0),holm_p:holm.map(String),global_tail_indices:tail,global_p:String(globalP),reject_global:globalP.cmp(alpha)<=0,reject_uncorrected_minimum:new Q(BigInt(minima[i]),BigInt(total)).cmp(alpha)<=0,reject_holm_any:holm.some(v=>v.cmp(alpha)<=0)};});
 const countKeys=['reject_global','reject_uncorrected_minimum','reject_holm_any'],counts=Object.fromEntries(countKeys.map(k=>[k,rows.filter(row=>row[k]).length]));
 return {schema:'driftlens-rational-orbit-v1',request:r,psd_certificates:certificates,allocation_count:total,observed_allocation_index:allowed.findIndex(row=>JSON.stringify(row)===JSON.stringify(r.reference)),allocations:rows,uniform_orbit_rejection_counts:counts,uniform_orbit_rejection_rates:Object.fromEntries(countKeys.map(k=>[k,String(new Q(BigInt(counts[k]),BigInt(total)))])),global_superuniform_at_alpha:new Q(BigInt(counts.reject_global),BigInt(total)).cmp(alpha)<=0,scope};
}
function parse(text){
 if(typeof text!=='string'||new TextEncoder().encode(text).length>65536)throw Error('Input exceeds65536bytes');
 let i=0;const ws=()=>{while(i<text.length&&/[\x20\t\r\n]/.test(text[i]))i++;};
 function value(path){ws();const c=text[i];if(c==='"'){const start=i++;let escaped=false;for(;i<text.length;i++){const ch=text[i];if(!escaped&&ch==='"'){i++;return JSON.parse(text.slice(start,i));}if(!escaped&&ch==='\\')escaped=true;else escaped=false;}throw Error('Unterminated JSON string');}
 if(c==='['){i++;ws();const a=[];if(text[i]===']'){i++;return a;}for(;;){a.push(value([...path,a.length]));ws();if(text[i]===']'){i++;return a;}if(text[i++]!==',')throw Error('Expected array comma');}}
 if(c==='{'){i++;ws();const o={},keys=new Set();if(text[i]==='}'){i++;return o;}for(;;){ws();if(text[i]!=='"')throw Error('Object key must be a string');const k=value(path);if(keys.has(k))throw Error('Duplicate JSON key: '+k);keys.add(k);ws();if(text[i++]!==':')throw Error('Expected colon');Object.defineProperty(o,k,{value:value([...path,k]),writable:true,enumerable:true,configurable:true});ws();if(text[i]==='}'){i++;return o;}if(text[i++]!==',')throw Error('Expected object comma');}}
 for(const [token,result]of [['true',true],['false',false],['null',null]])if(text.startsWith(token,i)){i+=token.length;return result;}
 const token=text.slice(i).match(/^-?(?:0|[1-9]\d*)(?:\.\d+)?(?:[eE][+-]?\d+)?/);if(!token)throw Error('Invalid JSON value');i+=token[0].length;const result=Number(token[0]);if(!Number.isFinite(result))throw Error('Nonfinite JSON number');if(path.at(-2)==='reference'&&!/^-?(?:0|[1-9]\d*)$/.test(token[0]))throw Error('Reference indices require integer JSON tokens');return result;
 }
 const r=value([]);ws();if(i!==text.length)throw Error('Trailing JSON content');return r;
}

const api={inspect,validate,parse,Q};if(typeof module!=='undefined')module.exports=api;root.DRIFT_ORBIT=api;
})(typeof globalThis!=='undefined'?globalThis:this);
