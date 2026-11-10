"""Complete rational finite-Gram randomization orbit; original small implementation."""
from fractions import Fraction as F
from itertools import combinations,product
import copy
from .orbit_input import validate,SCOPE


def _ldl(a):
    n=len(a);l=[[F(i==j)for j in range(n)]for i in range(n)];d=[]
    for j in range(n):
        pivot=a[j][j]-sum(l[j][k]**2*d[k]for k in range(j))
        if pivot<0:raise ValueError('Gram matrix is not positive semidefinite')
        d.append(pivot)
        for i in range(j+1,n):
            residual=a[i][j]-sum(l[i][k]*l[j][k]*d[k]for k in range(j))
            if pivot==0:
                if residual!=0:raise ValueError('zero Gram pivot has nonzero residual column')
            else:l[i][j]=residual/pivot
    return {'lower':[[str(x)for x in row]for row in l],'diagonal':[str(x)for x in d],'rank':sum(x>0 for x in d)}


def _allocations(request):
    n=len(request['row_names']);observed=set(request['reference']);labels=request['blocks']or['pooled']*n;groups=[]
    for label in dict.fromkeys(labels):
        row=[j for j,x in enumerate(labels)if x==label];count=len(observed.intersection(row));groups.append((row,count))
    return sorted(sorted(j for block in choice for j in block)for choice in product(*(combinations(row,k)for row,k in groups)))


def _score(a,left):
    n=len(a);seen=set(left);right=[i for i in range(n)if i not in seen];x=len(left);y=len(right)
    within_x=sum(a[i][j]for i in left for j in left if i!=j)
    within_y=sum(a[i][j]for i in right for j in right if i!=j)
    cross=sum(a[i][j]for i in left for j in right)
    return within_x/(x*(x-1))+within_y/(y*(y-1))-2*cross/(x*y)


def audit_orbit(request):
    matrices=validate(request);certificates=[_ldl(a)for a in matrices];allocations=_allocations(request);total=len(allocations)
    scores=[[_score(a,row)for row in allocations]for a in matrices]
    ranks=[[sum(x>=value for x in col)for value in col]for col in scores]
    minima=[min(col[j]for col in ranks)for j in range(total)]
    tails=[[j for j,value in enumerate(minima)if value<=threshold]for threshold in minima]
    alpha=F(request['alpha']);rows=[];family=len(matrices)
    for i,left in enumerate(allocations):
        raw=[F(col[i],total)for col in ranks];order=sorted(range(family),key=lambda k:(raw[k],k));holm=[None]*family;previous=F(0)
        for step,k in enumerate(order):previous=max(previous,min(F(1),(family-step)*raw[k]));holm[k]=previous
        rows.append({'reference':left,'current':[j for j in range(len(request['row_names']))if j not in left],
            'mmd2_unbiased':[str(col[i])for col in scores],'tail_counts':[col[i]for col in ranks],
            'raw_p':[str(x)for x in raw],'minimum_tail_count':minima[i],
            'selected_kernel_indices':[k for k,col in enumerate(ranks)if col[i]==minima[i]],
            'holm_p':[str(x)for x in holm],'global_tail_indices':tails[i],
            'global_p':str(F(len(tails[i]),total)),'reject_global':F(len(tails[i]),total)<=alpha,
            'reject_uncorrected_minimum':F(minima[i],total)<=alpha,'reject_holm_any':min(holm)<=alpha})
    counts={key:sum(row[key]for row in rows)for key in ['reject_global','reject_uncorrected_minimum','reject_holm_any']}
    return {'schema':'driftlens-rational-orbit-v1','request':copy.deepcopy(request),'psd_certificates':certificates,
        'allocation_count':total,'observed_allocation_index':allocations.index(request['reference']),
        'allocations':rows,'uniform_orbit_rejection_counts':counts,
        'uniform_orbit_rejection_rates':{k:str(F(v,total))for k,v in counts.items()},
        'global_superuniform_at_alpha':F(counts['reject_global'],total)<=alpha,'scope':SCOPE}
