"""Independent binary-mask orbit replay, all ordered pair sums and PSD reconstruction.

Does not import the orbit generator, its combinatorial enumeration or LDL routine.
"""
from fractions import Fraction as F
from .orbit_input import validate,SCOPE


def _q(value):
    if type(value)is not str or len(value)>2500:raise ValueError('report rationals must be strings')
    try:x=F(value)
    except(ZeroDivisionError,ValueError)as e:raise ValueError('invalid rational')from e
    if str(x)!=value or max(x.numerator.bit_length(),x.denominator.bit_length())>4096:raise ValueError('invalid canonical report rational')
    return x


def _same(actual,expected):
    if type(actual)is not type(expected):return False
    if type(expected)is dict:return set(actual)==set(expected)and all(_same(actual[k],v)for k,v in expected.items())
    if type(expected)is list:return len(actual)==len(expected)and all(_same(a,b)for a,b in zip(actual,expected))
    return actual==expected


def verify_orbit(report):
    try:
        if type(report)is not dict:raise ValueError('report object required')
        r=report['request'];matrices=validate(r);n=len(r['row_names']);left=r['reference'];labels=r['blocks']or['pooled']*n
        certs=report['psd_certificates']
        if type(certs)is not list or len(certs)!=len(matrices):raise ValueError('PSD certificate count')
        for a,cert in zip(matrices,certs):
            if type(cert)is not dict or set(cert)!={'lower','diagonal','rank'}:raise ValueError('PSD certificate keys')
            lower=cert['lower'];diagonal=cert['diagonal']
            if type(lower)is not list or len(lower)!=n or any(type(row)is not list or len(row)!=n for row in lower)or type(diagonal)is not list or len(diagonal)!=n:raise ValueError('PSD dimensions')
            l=[[_q(x)for x in row]for row in lower];d=[_q(x)for x in diagonal]
            if any(x<0 for x in d)or type(cert['rank'])is not int or cert['rank']!=sum(x>0 for x in d):raise ValueError('PSD nonnegative diagonal/rank')
            if any(l[i][j]!=(1 if i==j else 0)for i in range(n)for j in range(i,n)):raise ValueError('unit lower triangular required')
            if any(a[i][j]!=sum(l[i][k]*d[k]*l[j][k]for k in range(n))for i in range(n)for j in range(n)):raise ValueError('PSD reconstruction differs')
            # Canonical zero-pivot column, while valid LDL factors can be nonunique.
            if any(l[i][j]!=0 for j in range(n)if d[j]==0 for i in range(j+1,n)):raise ValueError('noncanonical zero-pivot column')
        allowed=[]
        for mask in range(1<<n):
            row=[j for j in range(n)if mask>>j&1]
            if len(row)!=len(left):continue
            if all(sum(labels[j]==label for j in row)==sum(labels[j]==label for j in left)for label in set(labels)):allowed.append(row)
        allowed.sort();total=len(allowed);statistics=[]
        for row in allowed:
            same=set(row);other=set(range(n))-same;values=[]
            for a in matrices:
                value=F(0)
                for i in range(n):
                    for j in range(n):
                        if i==j:continue
                        if i in same and j in same:value+=a[i][j]/(len(same)*(len(same)-1))
                        elif i in other and j in other:value+=a[i][j]/(len(other)*(len(other)-1))
                        else:value-=a[i][j]/(len(same)*len(other))
                values.append(value)
            statistics.append(values)
        rank_rows=[[sum(s[k]>=row[k]for s in statistics)for k in range(len(matrices))]for row in statistics]
        q=[min(row)for row in rank_rows];expected_rows=[];alpha=F(r['alpha'])
        for i,row in enumerate(allowed):
            raw=[F(v,total)for v in rank_rows[i]];holm=[]
            for k,p in enumerate(raw):
                # Equivalent closed max over all ordered ranks up to this tied p.
                ordered=sorted(raw);holm.append(str(min(F(1),max((len(raw)-j)*v for j,v in enumerate(ordered)if v<=p))))
            tail=[j for j in range(total)if q[j]<=q[i]];global_p=F(len(tail),total)
            expected_rows.append({'reference':row,'current':sorted(set(range(n))-set(row)),
                'mmd2_unbiased':[str(x)for x in statistics[i]],'tail_counts':rank_rows[i],'raw_p':[str(x)for x in raw],
                'minimum_tail_count':q[i],'selected_kernel_indices':[k for k,x in enumerate(rank_rows[i])if x==q[i]],
                'holm_p':holm,'global_tail_indices':tail,'global_p':str(global_p),'reject_global':global_p<=alpha,
                'reject_uncorrected_minimum':F(q[i],total)<=alpha,'reject_holm_any':min(map(F,holm))<=alpha})
        counts={key:sum(row[key]for row in expected_rows)for key in ['reject_global','reject_uncorrected_minimum','reject_holm_any']}
        expected={'schema':'driftlens-rational-orbit-v1','request':r,'psd_certificates':certs,'allocation_count':total,
            'observed_allocation_index':allowed.index(left),'allocations':expected_rows,'uniform_orbit_rejection_counts':counts,
            'uniform_orbit_rejection_rates':{k:str(F(v,total))for k,v in counts.items()},'global_superuniform_at_alpha':F(counts['reject_global'],total)<=alpha,'scope':SCOPE}
        if not _same(report,expected):raise ValueError('report differs from independent finite orbit replay')
        return True
    except(KeyError,TypeError,IndexError,OverflowError,ZeroDivisionError)as e:raise ValueError('invalid orbit report')from e
