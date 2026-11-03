"""Smaller affine-wreath graphs, preserving the original exact graph interface.
For primitive dimension d>=2, replace the old dense affine relation gadgets
by all affine hyperplanes. Dimension-one prime degrees retain old midpoints.
Positive certificates require only the generator-audited inclusion into Aut.
"""
from itertools import combinations,product
from math import prod
from pathlib import Path
import sys
STRUCT=Path(__file__).resolve().parent.parent/'structural_reduction'
sys.path.insert(0,str(STRUCT))
from affine_graphs import Graph,factor,verify_permutation

def add_word(g,w,offset,coordinates,component_index):
 def rec(v,off,isroot=False):
  node=g.vertex(10+component_index if isroot else 2)
  if not v:
   g.edge(node,off);g.tree.append({'vertex':node,'leaves':(off,),'children':(off,),'q':1,'coords':None});return node
  q=v[-1];m=prod(v[:-1]);children=[]
  for i in range(q):
   child=off+i if m==1 else rec(v[:-1],off+i*m)
   children.append(child);g.edge(node,child)
  coords=None
  if q>4:
   p,d=factor(q);coords=[tuple(x)for x in coordinates[str(q)]]
   assert len(coords)==q and set(coords)==set(product(range(p),repeat=d))
   if d>=2:
    normals=[]
    for a in product(range(p),repeat=d):
     if any(a) and next(x for x in a if x)==1:normals.append(a)
    assert len(normals)==(q-1)//(p-1)
    seen=set()
    for a in normals:
     for b in range(p):
      ends=tuple(sorted(children[i]for i in range(q)if sum(a[k]*coords[i][k]for k in range(d))%p==b))
      assert len(ends)==q//p and ends not in seen;seen.add(ends)
      rv=g.vertex(5);g.edge(rv,node)
      for child in ends:g.edge(rv,child)
      g.relations.append({'kind':'plane','geometric_kind':'affine_hyperplane','vertex':rv,'parent':node,'ends':ends})
    assert len(seen)==p*(q-1)//(p-1)
   else:
    # Full one-dimensional AGL needs midpoint relations; singleton
    # hyperplanes would give the larger symmetric group and no useful test.
    assert p==q and p%2==1;index={x:i for i,x in enumerate(coords)};half=pow(2,-1,p)
    for i,j in combinations(range(q),2):
     mid=((coords[i][0]+coords[j][0])*half%p,);k=index[mid]
     rv=g.vertex(3);cv=g.vertex(4)
     g.edge(rv,node);g.edge(rv,children[i]);g.edge(rv,children[j]);g.edge(rv,cv);g.edge(cv,children[k])
     g.relations.append({'kind':'mid','vertex':rv,'port':cv,'parent':node,'ends':tuple(sorted((children[i],children[j]))),'center':children[k]})
  g.tree.append({'vertex':node,'leaves':tuple(range(off,off+q*m)),'children':tuple(children),'q':q,'coords':coords})
  return node
 root=rec(tuple(w),offset,True);g.roots.append(root)

def cover_graph(words,cover,coordinates):
 n=sum(prod(words[i])for i in cover);g=Graph(n);offset=0
 for component,i in enumerate(cover):
  add_word(g,words[i],offset,coordinates,component);offset+=prod(words[i])
 return g
