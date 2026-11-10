'use strict';
const test=require('node:test'),assert=require('node:assert/strict'),cp=require('node:child_process');const I=require('../docs/orbit-engine.js');
const script=`import json,random
from fractions import Fraction as F
from driftlens.orbit_controls import definitions
from driftlens.orbit import audit_orbit
r=random.Random(125315);requests=[x['request']for x in definitions()]
for n in range(4,11):
 for seed in range(12):
  kernels=[]
  for k in range(1+seed%6):
   v=[[r.randrange(-3,4)for j in range(n)]for q in range(1+seed%3)]
   kernels.append(dict(name=str(k),gram=[[str(sum(F(row[i]*row[j],3)for row in v)+F((i==j)*(seed%2)))for j in range(n)]for i in range(n)]))
  requests.append(dict(row_names=[str(j)for j in range(n)],reference=list(range(2+seed%(n-3))),blocks=[str(j%3)for j in range(n)]if seed%3==0 else None,kernels=kernels,alpha='1/5'))
print(json.dumps([audit_orbit(x)for x in requests]))`;
test('full independent BigInt orbit reports exactly equal native Fraction reports and native-verify actual serialization',()=>{
 const reports=JSON.parse(cp.execFileSync('python',['-c',script],{maxBuffer:256*1024*1024})),actual=[],counts={reports:0,allocations:0,kernelScores:0,tailMembers:0};
 for(const expected of reports){const a=I.inspect(expected.request);assert.deepEqual(a,expected);actual.push(a);counts.reports++;counts.allocations+=a.allocation_count;counts.kernelScores+=a.allocation_count*a.request.kernels.length;counts.tailMembers+=a.allocations.reduce((s,x)=>s+x.global_tail_indices.length,0);}
 const output=cp.execFileSync('python',['-c',`import json,sys
from driftlens.orbit_verify import verify_orbit
rs=json.load(sys.stdin)
for r in rs:assert verify_orbit(r)
print(len(rs))`],{input:JSON.stringify(actual),maxBuffer:256*1024*1024});assert.equal(Number(output.toString()),reports.length);assert.equal(reports.length,95);console.log(JSON.stringify(counts));
});
test('singular zero pivots, indefinite matrices, schema, rationals, nonmutation and Unicode',()=>{
 require('../docs/orbit-data.js');const r=structuredClone(ORBIT_CONTROLS[0].request),before=structuredClone(r);I.inspect(r);assert.deepEqual(r,before);
 for(const c of ORBIT_CONTROLS)I.inspect(c.request);
 for(const patch of [{alpha:'2/4'},{alpha:'1.0'},{alpha:'-0'},{alpha:'1/0'},{alpha:'1000001/1000002'},{reference:[true,1]},{extra:1},{row_names:['\ud800','b','c','d','e','f']}])assert.throws(()=>I.inspect({...r,...patch}));
 const bad=structuredClone(ORBIT_CONTROLS[2].request);bad.kernels[0].gram[0][0]='-1';assert.throws(()=>I.inspect(bad));bad.kernels[0].gram=Array.from({length:4},(_,i)=>Array.from({length:4},(_,j)=>(i===0&&j===1)||(i===1&&j===0)?'1':'0'));assert.throws(()=>I.inspect(bad));
 const unicode=structuredClone(r);unicode.row_names[0]='α';unicode.row_names[1]='雪😀';assert.equal(I.inspect(unicode).request.row_names[1],'雪😀');
 const Q=I.Q;assert.equal(String(Q.parse('1/2').add(Q.parse('-1/3'))),'1/6');assert.equal(String(Q.parse('-4/3').div(Q.parse('2/9'))),'-6');
});

test('strict duplicate/size/nonfinite/index token parser',()=>{
 require('../docs/orbit-data.js');const text=JSON.stringify(ORBIT_CONTROLS[0].request);assert.deepEqual(I.parse(text),ORBIT_CONTROLS[0].request);for(const bad of ['{"alpha":"1/2","alpha":"1/3"}','NaN','[1e309]',' '.repeat(65537),text.replace('"reference":[0,2,3]','"reference":[0e0,2,3]')])assert.throws(()=>I.validate(I.parse(bad)));
 const x=I.parse('{"__proto__":{"polluted":1}}');assert.equal(Object.hasOwn(x,'__proto__'),true);assert.equal({}.polluted,undefined);assert.throws(()=>I.inspect(x));
});
