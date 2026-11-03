"""Colored-graph representation of full affine-wreath cover groups.

Leaf vertices are 0..n-1. Internal vertex colors distinguish rooted-tree
structure from affine relation gadgets. All component roots receive distinct
colors, fixing the direct-product factors. See GRAPH_ENCODING.md.

No automorphism implementation is embedded here: the output is input to
nauty/Traces or an independent graph verifier. Explicit GAP group generators
can be checked to preserve the encoding with --verify-generators.
"""
from pathlib import Path
from itertools import combinations
from math import prod
import argparse,json,re

class Graph:
 def __init__(self,n):
  self.n=n;self.colors=[0]*n;self.edges=set();self.tree=[];self.relations=[];self.roots=[]
 def vertex(self,color):
  x=len(self.colors);self.colors.append(color);return x
 def edge(self,x,y):
  assert x!=y
  self.edges.add((min(x,y),max(x,y)))
 def export(self):
  return {'degree':self.n,'vertex_count':len(self.colors),'colors':self.colors,
          'edges':[list(e) for e in sorted(self.edges)],'leaf_vertices':list(range(self.n)),
          'root_vertices':self.roots}

def factor(q):
 p=next(i for i in range(2,q+1) if q%i==0);d=0;x=q
 while x%p==0:x//=p;d+=1
 assert x==1
 return p,d

def add_word(g,w,offset,coordinates,component_index):
 def rec(v,off,isroot=False):
  color=10+component_index if isroot else 2
  node=g.vertex(color)
  if not v:
   g.edge(node,off)
   g.tree.append({'vertex':node,'leaves':(off,),'children':(off,), 'q':1,'coords':None})
   return node
  q=v[-1];m=prod(v[:-1]);children=[]
  for i in range(q):
   child=off+i if m==1 else rec(v[:-1],off+i*m)
   children.append(child);g.edge(node,child)
  coords=None
  if q>4:
   p,d=factor(q);coords=[tuple(x) for x in coordinates[str(q)]]
   assert len(coords)==q and len(set(coords))==q
   index={x:i for i,x in enumerate(coords)}
   if p!=2:
    half=pow(2,-1,p)
    for i,j in combinations(range(q),2):
     mid=tuple(((coords[i][k]+coords[j][k])*half)%p for k in range(d))
     k=index[mid]
     rv=g.vertex(3);cv=g.vertex(4)
     g.edge(rv,node);g.edge(rv,children[i]);g.edge(rv,children[j]);g.edge(rv,cv);g.edge(cv,children[k])
     g.relations.append({'kind':'mid','vertex':rv,'port':cv,'parent':node,
                         'ends':tuple(sorted((children[i],children[j]))),'center':children[k]})
   else:
    for i,j,k in combinations(range(q),3):
     fourth=tuple(coords[i][h]^coords[j][h]^coords[k][h] for h in range(d))
     ell=index[fourth]
     if ell<=k:continue
     rv=g.vertex(5);g.edge(rv,node)
     for c in(i,j,k,ell):g.edge(rv,children[c])
     g.relations.append({'kind':'plane','vertex':rv,'parent':node,
                         'ends':tuple(sorted(children[c] for c in(i,j,k,ell)))})
  g.tree.append({'vertex':node,'leaves':tuple(range(off,off+q*m)),
                 'children':tuple(children),'q':q,'coords':coords})
  return node
 root=rec(tuple(w),offset,True);g.roots.append(root)

def cover_graph(words,cover,coordinates):
 n=sum(prod(words[i]) for i in cover);g=Graph(n);offset=0
 for component,i in enumerate(cover):
  add_word(g,words[i],offset,coordinates,component);offset+=prod(words[i])
 return g

def verify_permutation(g,perm):
 """Extend a one-based image-list permutation and verify graph automorphism."""
 p=[i-1 for i in perm];n=g.n
 if len(p)!=n or sorted(p)!=list(range(n)):return False,'invalid permutation'
 mapping={i:p[i] for i in range(n)}
 nodes={(g.colors[t['vertex']],frozenset(t['leaves'])):t['vertex'] for t in g.tree}
 for t in g.tree:
  key=(g.colors[t['vertex']],frozenset(p[x] for x in t['leaves']))
  if key not in nodes:return False,'tree block not preserved'
  mapping[t['vertex']]=nodes[key]
 keys={}
 for r in g.relations:
  if r['kind']=='mid':key=('mid',r['parent'],r['ends'],r['center'])
  else:key=('plane',r['parent'],r['ends'])
  assert key not in keys
  keys[key]=r
 for r in g.relations:
  if r['kind']=='mid':key=('mid',mapping[r['parent']],tuple(sorted(mapping[x] for x in r['ends'])),mapping[r['center']])
  else:key=('plane',mapping[r['parent']],tuple(sorted(mapping[x] for x in r['ends'])))
  if key not in keys:return False,'affine relation not preserved'
  t=keys[key];mapping[r['vertex']]=t['vertex']
  if r['kind']=='mid':mapping[r['port']]=t['port']
 if len(mapping)!=len(g.colors) or len(set(mapping.values()))!=len(mapping):return False,'extension not bijective'
 if any(g.colors[i]!=g.colors[mapping[i]] for i in mapping):return False,'color changed'
 if any(tuple(sorted((mapping[u],mapping[v]))) not in g.edges for u,v in g.edges):return False,'edge changed'
 return True,'verified'

def parse_gap_assignment(text,name):
 s=re.search(r'\b'+re.escape(name)+r'\s*:=\s*(.*?);',text,re.S)
 if not s:raise ValueError('missing '+name)
 value=re.sub(r'\[\s*(-?\d+)\s*\.\.\s*(-?\d+)\s*\]',lambda m:json.dumps(list(range(int(m.group(1)),int(m.group(2))+1))),s.group(1))
 return json.loads(value)

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--template',required=True);ap.add_argument('--coordinates',required=True)
 ap.add_argument('--degree',type=int);ap.add_argument('--mode',choices=['covers','words'],default='words')
 ap.add_argument('--output',required=True);ap.add_argument('--verify-generators')
 args=ap.parse_args();data=json.load(open(args.template));coords=json.load(open(args.coordinates));words=data['words']
 covers=data['covers_0based'][str(args.degree)] if args.degree else None
 if args.mode=='covers':
  assert covers is not None
  graphs=[cover_graph(words,c,coords) for c in covers]
  output={'degree':args.degree,'covers_0based':covers,'graphs':[g.export() for g in graphs]}
 else:
  graphs=[cover_graph(words,[i],coords) for i in range(len(words))]
  output={'max_degree':data['max_degree'],'words':words,'word_graphs':[g.export() for g in graphs],
           'covers_0based':data['covers_0based']}
 if args.verify_generators:
  assert covers is not None
  text=Path(args.verify_generators).read_text();degree=parse_gap_assignment(text,'ExplicitDegree')
  covers1=parse_gap_assignment(text,'ExplicitCoverWords');gens=parse_gap_assignment(text,'ExplicitGenerators')
  assert degree==args.degree and covers1==[[i+1 for i in c] for c in covers]
  if args.mode!='covers':graphs=[cover_graph(words,c,coords) for c in covers]
  checks=0
  for index,(g,gg) in enumerate(zip(graphs,gens)):
   for gen in gg:
    ok,why=verify_permutation(g,gen)
    if not ok:raise ValueError(f'cover {index+1}: {why}')
    checks+=1
  output['generator_audit']={'groups_checked':len(graphs),'generators_checked':checks,'status':'every generator extends to a color-preserving graph automorphism'}
 Path(args.output).write_text(json.dumps(output,separators=(',',':')))
 print({'output':args.output,'mode':args.mode,'graphs':len(graphs),'max_vertices':max(len(g.colors) for g in graphs),
        'max_edges':max(len(g.edges) for g in graphs),'generator_audit':output.get('generator_audit')})
if __name__=='__main__':main()
