import unittest,random,copy,itertools,json
from fractions import Fraction as F
from driftlens.orbit import audit_orbit
from driftlens.orbit_verify import verify_orbit


def request(vectors,reference=None,blocks=None,alpha='1/5'):
    n=len(vectors[0]);return dict(row_names=[str(j)for j in range(n)],reference=reference or list(range(n//2)),blocks=blocks,
        kernels=[dict(name='view'+str(k),gram=[[str(F(x)*F(y))for y in row]for x in row])for k,row in enumerate(vectors)],alpha=alpha)


class OrbitTests(unittest.TestCase):
    def test_exact_selection_penalty_and_holm_comparison(self):
        r=request([[0,2,0,1,-1,-3],[-1,1,-3,1,0,1],[1,3,-3,1,2,3]],[0,2,3]);a=audit_orbit(r);self.assertTrue(verify_orbit(a));row=a['allocations'][a['observed_allocation_index']]
        self.assertEqual(row['mmd2_unbiased'],['-11/9','4/3','64/9']);self.assertEqual(row['tail_counts'],[14,8,2]);self.assertEqual(row['global_p'],'1/5');self.assertEqual(row['holm_p'],['4/5','4/5','3/10']);self.assertEqual(a['uniform_orbit_rejection_counts'],dict(reject_global=4,reject_uncorrected_minimum=8,reject_holm_any=0));self.assertEqual(row['global_tail_indices'],[4,5,14,15])
    def test_one_kernel_duplicates_ties_and_zero_psd(self):
        for vectors in [[[0,1,2,3]],[[1,1,1,1]],[[0,0,0,0]]]:
            r=request(vectors,alpha='1/3');a=audit_orbit(r);self.assertTrue(verify_orbit(a))
            for row in a['allocations']:self.assertEqual(row['raw_p'][0],row['global_p'])
            r['kernels']*=3
            for k,v in enumerate(r['kernels']):r['kernels'][k]=dict(name=str(k),gram=v['gram'])
            b=audit_orbit(r);self.assertTrue(verify_orbit(b));self.assertEqual([x['global_p']for x in a['allocations']],[x['global_p']for x in b['allocations']])
    def test_restrictions_frozen_and_group_swap(self):
        r=request([[0,0,1,1,2,2]],[0,2,4],['a','a','b','b','c','c']);a=audit_orbit(r);self.assertEqual(a['allocation_count'],8);self.assertTrue(verify_orbit(a));self.assertTrue(all(x['global_p']=='1'for x in a['allocations']))
        r['blocks']=['left','right','left','right','left','right'];a=audit_orbit(r);self.assertEqual(a['allocation_count'],1);self.assertEqual(a['allocations'][0]['global_p'],'1')
        r=request([[0,1,2,3,4,5]],[0,1]);a=audit_orbit(r);r['reference']=[2,3,4,5];b=audit_orbit(r);self.assertTrue(verify_orbit(b));lookup={tuple(x['current']):x for x in b['allocations']}
        for row in a['allocations']:self.assertEqual(row['mmd2_unbiased'],lookup[tuple(row['reference'])]['mmd2_unbiased']);self.assertEqual(row['global_p'],lookup[tuple(row['reference'])]['global_p'])
    def test_fraction_psd_full_rank_and_singular(self):
        r=request([['1/2','-2/3','1','0']]);r['kernels'][0]['gram']=[[str(F(i==j)+F(r['kernels'][0]['gram'][i][j]))for j in range(4)]for i in range(4)];a=audit_orbit(r);self.assertTrue(verify_orbit(a));self.assertEqual(a['psd_certificates'][0]['rank'],4)
        r['kernels'][0]['gram']=[['0','1','0','0'],['1','0','0','0'],['0','0','1','0'],['0','0','0','1']]
        with self.assertRaises(ValueError):audit_orbit(r)
        r['kernels'][0]['gram']=[['-1'if i==j else '0'for j in range(4)]for i in range(4)]
        with self.assertRaises(ValueError):audit_orbit(r)
    def test_seeded_full_orbits_all_thresholds(self):
        rng=random.Random(91537);cases=0;allocations=0
        for n in range(4,9):
            for seed in range(10):
                vectors=[[rng.randrange(-3,4)for _ in range(n)]for k in range(1+seed%3)];r=request(vectors,alpha='2/5');
                for k,kernel in enumerate(r['kernels']):
                    features=[[rng.randrange(-2,3)for _ in range(n)]for q in range(seed%4)]
                    kernel['gram']=[[str(F(kernel['gram'][i][j])+sum(F(v[i]*v[j],3)for v in features)+F((i==j)*(seed%2)))for j in range(n)]for i in range(n)]
                a=audit_orbit(r);self.assertTrue(verify_orbit(a));cases+=1;allocations+=a['allocation_count'];rows=a['allocations'];total=len(rows)
                for threshold in range(1,total+1):
                    count=sum(F(row['global_p'])<=F(threshold,total)for row in rows);self.assertLessEqual(count,threshold)
                for row in rows:self.assertGreaterEqual(F(row['global_p']),min(map(F,row['raw_p'])));self.assertLessEqual(F(row['global_p']),min(map(F,row['holm_p'])))
        print('independent rational orbit replay cases',cases,'allocations',allocations)
    def test_every_report_primitive_individually_mutated(self):
        a=audit_orbit(request([[0,2,0,1,-1,-3],[-1,1,-3,1,0,1]],[0,2,3]));self.assertTrue(verify_orbit(a));paths=[]
        def visit(x,path):
            if type(x)is dict:
                for k,v in x.items():
                    if k!='request':visit(v,path+[k])
            elif type(x)is list:
                for j,v in enumerate(x):visit(v,path+[j])
            else:paths.append(path)
        visit(a,[])
        for path in paths:
            b=copy.deepcopy(a);parent=b
            for k in path[:-1]:parent=parent[k]
            x=parent[path[-1]];parent[path[-1]]=not x if type(x)is bool else x+1 if type(x)is int else x+'!'
            with self.assertRaises(ValueError,msg=str(path)):verify_orbit(b)
        print('orbit non-request primitive mutations rejected',len(paths))
    def test_schema_and_metadata_types_are_exact(self):
        r=request([[0,1,2,3]]);before=copy.deepcopy(r);a=audit_orbit(r);self.assertEqual(r,before)
        for patch in [dict(reference=[True,1]),dict(alpha='2/4'),dict(alpha='-0'),dict(alpha='1.0'),dict(alpha='1'),dict(blocks=['a']),dict(row_names=['a','a','b','c']),dict(extra=1)]:
            with self.assertRaises(ValueError):audit_orbit({**r,**patch})
        for entry in [True,1.0,'01','1/0','1000001','1/1000001','1e0']:
            bad=copy.deepcopy(r);bad['kernels'][0]['gram'][0][0]=entry
            with self.assertRaises(ValueError):audit_orbit(bad)
        for key in ['allocation_count','observed_allocation_index']:
            bad=copy.deepcopy(a);bad[key]=float(bad[key])
            with self.assertRaises(ValueError):verify_orbit(bad)
        for key in ['request','psd_certificates','allocations','scope']:
            bad=copy.deepcopy(a);del bad[key]
            with self.assertRaises(ValueError):verify_orbit(bad)
        bad=copy.deepcopy(a);bad['extra']=0
        with self.assertRaises(ValueError):verify_orbit(bad)
    def test_maximum_boundary_six_kernels_252_allocations_and_request_size(self):
        vectors=[[((j+k)%7)-3 for j in range(10)]for k in range(6)];r=request(vectors,alpha='1/7');a=audit_orbit(r);self.assertTrue(verify_orbit(a));self.assertEqual(a['allocation_count'],252);self.assertEqual(len(a['allocations']),252);self.assertLess(len(json.dumps(r).encode()),65536);self.assertLess(len(json.dumps(a).encode()),8*1024*1024)
    def test_label_reordering_and_kernel_permutation_invariance(self):
        r=request([[0,2,0,1,-1,-3],[-1,1,-3,1,0,1],[1,3,-3,1,2,3]],[0,2,3]);a=audit_orbit(r);p=[3,1,5,0,4,2];b=copy.deepcopy(r);b['row_names']=[r['row_names'][j]for j in p];b['reference']=sorted(p.index(j)for j in r['reference'])
        for k,matrix in enumerate(r['kernels']):b['kernels'][k]['gram']=[[matrix['gram'][i][j]for j in p]for i in p]
        out=audit_orbit(b);self.assertTrue(verify_orbit(out));lookup={tuple(sorted(p[j]for j in row['reference'])):row for row in out['allocations']}
        for row in a['allocations']:self.assertEqual(row['raw_p'],lookup[tuple(row['reference'])]['raw_p']);self.assertEqual(row['global_p'],lookup[tuple(row['reference'])]['global_p'])
        b=copy.deepcopy(r);b['kernels']=b['kernels'][::-1];c=audit_orbit(b)
        for x,y in zip(a['allocations'],c['allocations']):self.assertEqual(x['global_p'],y['global_p']);self.assertEqual(x['raw_p'],y['raw_p'][::-1])
    def test_psd_acceptance_matches_all_principal_minors_exhaustively(self):
        def determinant(a):
            n=len(a)
            if n==0:return F(1)
            return sum((-1)**j*a[0][j]*determinant([row[:j]+row[j+1:]for row in a[1:]])for j in range(n))
        r=request([[0,0,0,0]]);cases=0
        for diagonal in itertools.product(range(-1,2),repeat=3):
            for off in itertools.product(range(-1,2),repeat=3):
                m=[[F(0)for _ in range(4)]for _ in range(4)]
                for j,x in enumerate(diagonal):m[j][j]=F(x)
                for (i,j),x in zip([(0,1),(0,2),(1,2)],off):m[i][j]=m[j][i]=F(x)
                positive=all(determinant([[m[i][j]for j in indices]for i in indices])>=0 for k in range(1,5)for indices in itertools.combinations(range(4),k));r['kernels'][0]['gram']=[[str(x)for x in row]for row in m]
                if positive:self.assertTrue(verify_orbit(audit_orbit(r)))
                else:
                    with self.assertRaises(ValueError):audit_orbit(r)
                cases+=1
        self.assertEqual(cases,729);print('PSD accept/reject matched all principal-minor oracle matrices',cases)
