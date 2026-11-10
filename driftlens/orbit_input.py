"""Strict finite rational-Gram request contract, separate from floating RBF APIs."""
from fractions import Fraction
import re

SCOPE=('Complete finite supplied-Gram minimum-p global randomization test. Exact rational '
 'arithmetic and all ties retained. Requires a fixed label-independent kernel family, '
 'prespecified blocks and joint invariance under allowed reallocations. PSD is finite '
 'geometry, not a characteristic population kernel or validated sampling assumption. '
 'One global null decision only, not strong component-wise inference, causality, '
 'model-quality, equivalence, automatic bandwidth choice or continuous monitoring.')


def rational(value):
    if type(value)is not str or not re.fullmatch(r'-?(?:0|[1-9][0-9]*)(?:/[1-9][0-9]*)?',value)or len(value)>16:
        raise ValueError('rational entries must be bounded canonical fraction strings')
    x=Fraction(value)
    if abs(x.numerator)>1000000 or x.denominator>1000000 or str(x)!=value:
        raise ValueError('rational entry must be reduced and numerator/denominator <=1000000')
    return x


def text(value):
    return type(value)is str and bool(value.strip())and len(value)<=40 and not any(0xD800<=ord(c)<=0xDFFF for c in value)


def validate(request):
    if type(request)is not dict or set(request)!={'row_names','reference','blocks','kernels','alpha'}:
        raise ValueError('request keys must be row_names/reference/blocks/kernels/alpha')
    names=request['row_names']
    if type(names)is not list or not 4<=len(names)<=10 or not all(map(text,names))or len(set(names))!=len(names):
        raise ValueError('four to ten distinct bounded row names required')
    n=len(names);left=request['reference']
    if type(left)is not list or not 2<=len(left)<=n-2 or any(type(j)is not int or not 0<=j<n for j in left)or left!=sorted(set(left)):
        raise ValueError('reference must be sorted distinct integer indices with each side >=2')
    blocks=request['blocks']
    if blocks is not None and(type(blocks)is not list or len(blocks)!=n or not all(map(text,blocks))):
        raise ValueError('blocks must be null or aligned bounded string labels')
    family=request['kernels']
    if type(family)is not list or not 1<=len(family)<=6:raise ValueError('one to six kernels required')
    seen=set();matrices=[]
    for kernel in family:
        if type(kernel)is not dict or set(kernel)!={'name','gram'}or not text(kernel['name'])or kernel['name']in seen:
            raise ValueError('kernels need distinct names and gram only')
        seen.add(kernel['name']);m=kernel['gram']
        if type(m)is not list or len(m)!=n or any(type(row)is not list or len(row)!=n for row in m):raise ValueError('aligned square Gram matrix required')
        a=[[rational(x)for x in row]for row in m]
        if any(a[i][j]!=a[j][i]for i in range(n)for j in range(n)):raise ValueError('Gram matrix must be exactly symmetric')
        matrices.append(a)
    alpha=rational(request['alpha'])
    if not 0<alpha<1:raise ValueError('alpha must be strictly between zero and one')
    return matrices
