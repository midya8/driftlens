"""Original supplied finite geometries, selected development controls, not power tests."""
import copy


def _request(vectors,reference=None,blocks=None,alpha='1/5'):
    n=len(vectors[0]);return dict(row_names=[chr(97+j)for j in range(n)],reference=reference or list(range(n//2)),blocks=blocks,
        kernels=[dict(name='view '+str(k+1),gram=[[str(x*y)for y in v]for x in v])for k,v in enumerate(vectors)],alpha=alpha)


def definitions():
    selection=_request([[0,2,0,1,-1,-3],[-1,1,-3,1,0,1],[1,3,-3,1,2,3]],[0,2,3])
    duplicate=_request([[0,0,0,1,1,1]]*3,alpha='1/10')
    all_ties=_request([[1,1,1,1]],alpha='1/3')
    zero=_request([[0,0,0,0]],alpha='1/3')
    blocked=_request([[0,0,1,1,2,2]],[0,2,4],['site a','site a','site b','site b','site c','site c'])
    mixture=_request([[0,0,1,1,2,2]],[0,1,2],alpha='1/5')
    restricted=copy.deepcopy(mixture);restricted['blocks']=['low','low','middle','middle','high','high']
    frozen=copy.deepcopy(mixture);frozen['blocks']=['reference'if j in mixture['reference']else 'current'for j in range(6)]
    unequal=_request([[-2,-1,0,1,2,3],[1,-1,1,-1,1,-1]],[0,1],alpha='1/3')
    nonlinear=_request([[0,1,-1,2,-2,1],[0,1,1,4,4,1]],[0,2,4])
    fullrank=_request([[0,1,2,3]],alpha='1/3');fullrank['kernels'][0]['gram']=[[str(i*j+(i==j))for j in range(4)]for i in range(4)]
    return [dict(id=id,title=title,note=note,request=r)for id,title,note,r in [
        ('selection','Choosing a view changes the test','Three supplied rank-one views: raw minimum rejects 8/20 allowed allocations at alpha 1/5, global recalibration 4/20, Holm 0/20. Selected in development, not a power benchmark.',selection),
        ('duplicate','Three copies, one piece of evidence','Duplicated geometry leaves the global p unchanged; Holm counts the prespecified family. Not three independent replications.',duplicate),
        ('all-ties','Ties stay in the tail','Constant geometry has no label contrast. Every tied score and every allocation remains in the tail.',all_ties),
        ('zero','Zero geometry is still PSD','Rank zero is allowed. No nonzero pivot is fabricated and all global p-values equal one.',zero),
        ('blocked','Same features inside each site','Within-site exchanges preserve each side count. They ask a conditional null question, not a causal adjustment.',blocked),
        ('mixture','Unrestricted label exchanges','The entire pooled allocation orbit, with the same score definition as the restricted control.',mixture),
        ('restricted','Restrict the same mixture','The same observed rows and Gram matrix, but only within-stratum labels may move. A different null orbit, not an automatic repair.',restricted),
        ('frozen','No label can move','All blocks contain one side only. The sole allocation gives p=1, explicitly uninformative.',frozen),
        ('unequal','Unequal groups, same rule','Two reference and four current rows. All off-diagonal normalizers change with the group size.',unequal),
        ('nonlinear','Two fixed feature maps','Linear and squared coordinates are supplied as different finite Gram views. They were not chosen by an automatic bandwidth tuner.',nonlinear),
        ('full-rank','A full-rank finite geometry','Identity plus a rank-one view has a nonsingular exact LDL certificate. PSD alone does not validate a population null.',fullrank)]]
