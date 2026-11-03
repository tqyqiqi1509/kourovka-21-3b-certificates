// Exact coloured-graph certificate checker for affine-envelope intersections.
// The graph encoding is independently supplied and generator-audited.
// Accept only when every common point vertex is fixed by the full colour-
// preserving graph automorphism group. No group-order floating point is used
// in the acceptance condition.
#include <algorithm>
#include <chrono>
#include <fstream>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>
#include <cstdint>
#include <random>
extern "C" int affine_graph_nauty_version(void);
extern "C" int affine_graph_leaf_fix(int,size_t,size_t*,int*,int*,int*,int,int,int*,unsigned long*);
using namespace std;
struct Graph {int nv=0; vector<int> color; vector<pair<int,int>> edges; bool odd=false; vector<int> core;};
struct Family {int n=0; vector<Graph> g;};
Family read_family(const string& path){
 ifstream f(path);string magic;f>>magic; if(magic!="KOUROVKA_GRAPH_COVERS_V1")throw runtime_error("bad cover format");
 Family a;int k;f>>a.n>>k; if(a.n<1||k<1)throw runtime_error("bad family size");
 a.g.resize(k);
 for(int i=0;i<k;i++){
  int id,ne,odd;auto &g=a.g[i];f>>id>>g.nv>>ne>>odd;
  if(id!=i+1||g.nv<a.n||ne<0)throw runtime_error("bad graph header");g.odd=odd;
  g.color.resize(g.nv);for(auto &c:g.color)f>>c;
  for(int v=0;v<g.nv;v++)if((v<a.n)!=(g.color[v]==0))throw runtime_error("point colour mismatch");
  for(int e=0;e<ne;e++){int u,v;f>>u>>v;if(u<0||v<0||u>=g.nv||v>=g.nv||u==v)throw runtime_error("bad edge");if(u>v)swap(u,v);g.edges.emplace_back(u,v);}
  sort(g.edges.begin(),g.edges.end());if(adjacent_find(g.edges.begin(),g.edges.end())!=g.edges.end())throw runtime_error("duplicate edge");
  vector<vector<int>> ln(a.n);for(auto [u,v]:g.edges){if(u<a.n)ln[u].push_back(v);if(v<a.n)ln[v].push_back(u);}
  vector<int>keys;for(int v=0;v<a.n;v++)keys.push_back(ln[v].size()==1?ln[v][0]:g.nv+v);
  auto uniq=keys;sort(uniq.begin(),uniq.end());uniq.erase(unique(uniq.begin(),uniq.end()),uniq.end());
  for(int key:keys)g.core.push_back(lower_bound(uniq.begin(),uniq.end(),key)-uniq.begin());
 }
 if(!f)throw runtime_error("truncated graph family");return a;
}
vector<vector<int>> read_pool(const string&path,int n){
 ifstream f(path);int nn,k;f>>nn>>k;if(nn!=n||k<1)throw runtime_error("bad pool header");
 vector<vector<int>>p(k,vector<int>(n));
 for(int c=0;c<k;c++){int id;f>>id;if(id!=c+1)throw runtime_error("pool index");for(int &v:p[c]){f>>v;--v;}auto a=p[c];sort(a.begin(),a.end());for(int i=0;i<n;i++)if(a[i]!=i)throw runtime_error("invalid pool permutation");}
 if(!f)throw runtime_error("truncated pool");return p;
}
int parity(const vector<int>&p){vector<char>s(p.size(),0);int cyc=0;for(int i=0;i<(int)p.size();i++)if(!s[i]){++cyc;int j=i;do{s[j]=1;j=p[j];}while(!s[j]);}return ((p.size()-cyc)&1)?-1:1;}
struct CheckResult {bool fixed=false;int nv=0,orbits=0;unsigned long nodes=0;};
CheckResult check(const Graph&a,const Graph&b,const vector<int>&x,int n,bool traces){
 const int nv=a.nv+b.nv-n;
 vector<int>colors=a.color;colors.resize(nv);int coff=*max_element(a.color.begin(),a.color.end())+1;
 for(int v=n;v<b.nv;v++)colors[a.nv+v-n]=coff+b.color[v];
 auto mapb=[&](int v){return v<n?x[v]:a.nv+v-n;};
 vector<int>degrees(nv,0);for(auto [u,v]:a.edges){++degrees[u];++degrees[v];}
 for(auto [u,v]:b.edges){++degrees[mapb(u)];++degrees[mapb(v)];}
 vector<size_t>offs(nv);size_t nde=0;for(int i=0;i<nv;i++){offs[i]=nde;nde+=degrees[i];}
 vector<int>neigh(nde);auto at=offs;
 for(auto [u,v]:a.edges){neigh[at[u]++]=v;neigh[at[v]++]=u;}
 for(auto [u0,v0]:b.edges){int u=mapb(u0),v=mapb(v0);neigh[at[u]++]=v;neigh[at[v]++]=u;}
 for(int i=0;i<nv;i++)sort(neigh.begin()+offs[i],neigh.begin()+offs[i]+degrees[i]);
 CheckResult out;out.nv=nv;
 int status=affine_graph_leaf_fix(nv,nde,offs.data(),degrees.data(),neigh.data(),colors.data(),n,traces?1:0,&out.orbits,&out.nodes);
 if(status<0)throw runtime_error("graph automorphism engine error");out.fixed=(status==1);return out;
}

bool core_disjoint(const Graph&a,const Graph&b,const vector<int>&x,vector<unsigned>&mark,unsigned&stamp){
 if(++stamp==0){fill(mark.begin(),mark.end(),0);stamp=1;}
 const int n=x.size();
 for(int i=0;i<n;i++){int key=b.core[i]*n+a.core[x[i]];if(mark[key]==stamp)return false;mark[key]=stamp;}
 return true;
}
int run_search(const string&file,const string&prefix,const string&engine,int maxcand,uint64_t seed){
 bool traces=(engine=="traces");if(!traces&&engine!="nauty")throw runtime_error("bad engine");
 auto f=read_family(file);int n=f.n,k=f.g.size();
 ofstream cert(prefix+"_"+to_string(n)+"_certificates.jsonl"),poolout(prefix+"_"+to_string(n)+"_permutation_pool.jsonl"),fails(prefix+"_"+to_string(n)+"_unresolved.jsonl");
 vector<vector<int>>pool;mt19937_64 rng(seed+n);
 vector<pair<int,int>>pairs;for(int i=0;i<k;i++)pairs.emplace_back(i,i);for(int i=0;i<k;i++)for(int j=i+1;j<k;j++)pairs.emplace_back(i,j);
 vector<unsigned>marks(n*n,0);unsigned stamp=0;long long successes=0,unresolved=0,proposals=0,autcalls=0;
 int all_even=0;for(auto &g:f.g)if(!g.odd)++all_even;
 auto start=chrono::steady_clock::now(),last=start;
 cerr<<"GRAPH_SEARCH_START degree="<<n<<" covers="<<k<<" pairs="<<pairs.size()<<" all_even="<<all_even<<" engine="<<engine<<" seed="<<(seed+n)<<"\n";
 for(size_t z=0;z<pairs.size();z++){
  auto [i,j]=pairs[z];vector<int>parities{1};if(!f.g[i].odd&&!f.g[j].odd)parities.push_back(-1);
  for(int p:parities){bool good=false;int tried=0,called=0,c=0;
   for(c=0;c<maxcand;c++){
    ++tried;++proposals;
    if(c==(int)pool.size()){
     vector<int>x(n);iota(x.begin(),x.end(),0);shuffle(x.begin(),x.end(),rng);if(parity(x)==-1)for(int &v:x){if(v==0)v=1;else if(v==1)v=0;}
     pool.push_back(x);poolout<<"{\"index\":"<<(c+1)<<",\"x\":[";
     for(int a=0;a<n;a++){if(a)poolout<<",";poolout<<(x[a]+1);}poolout<<"]}\n";poolout.flush();
    }
    auto x=pool[c];if(p==-1)for(int &v:x){if(v==0)v=1;else if(v==1)v=0;}
    if(!core_disjoint(f.g[i],f.g[j],x,marks,stamp))continue;
    ++called;++autcalls;
    auto result=check(f.g[i],f.g[j],x,n,traces);
    if(result.fixed){good=true;break;}
   }
   if(good){++successes;cert<<"{\"n\":"<<n<<",\"i\":"<<(i+1)<<",\"j\":"<<(j+1)<<",\"parity\":"<<p<<",\"candidate\":"<<(c+1)<<",\"proposals\":"<<tried<<",\"automorphism_calls\":"<<called<<"}\n";}
   else{++unresolved;fails<<"{\"n\":"<<n<<",\"i\":"<<(i+1)<<",\"j\":"<<(j+1)<<",\"parity\":"<<p<<",\"proposals\":"<<tried<<"}\n";cerr<<"GRAPH_UNRESOLVED "<<n<<" "<<(i+1)<<" "<<(j+1)<<" "<<p<<"\n";}
  }
  auto now=chrono::steady_clock::now();if(chrono::duration<double>(now-last).count()>=15){
   cert.flush();fails.flush();cerr<<"GRAPH_SEARCH_PROGRESS degree="<<n<<" pairs_done="<<(z+1)<<"/"<<pairs.size()<<" certificates="<<successes<<" unresolved="<<unresolved<<" proposals="<<proposals<<" aut_calls="<<autcalls<<" seconds="<<chrono::duration<double>(now-start).count()<<"\n";last=now;
  }
 }
 double sec=chrono::duration<double>(chrono::steady_clock::now()-start).count();
 cerr<<"GRAPH_SEARCH_DONE degree="<<n<<" pairs="<<pairs.size()<<" certificates="<<successes<<" unresolved="<<unresolved<<" proposals="<<proposals<<" aut_calls="<<autcalls<<" pool="<<pool.size()<<" seconds="<<sec<<"\n";
 return unresolved?1:0;
}
int main(int argc,char**argv){
 try{
  if(argc!=7){cerr<<"Usage: graph_engine check COVERS POOL QUERIES OUTPUT ENGINE; or search COVERS PREFIX ENGINE MAX_CANDIDATES SEED\n";return 2;}
  if(string(argv[1])=="search")return run_search(argv[2],argv[3],argv[4],stoi(argv[5]),stoull(argv[6]));
  if(string(argv[1])!="check")throw runtime_error("only check mode implemented");
  bool traces=string(argv[6])=="traces";if(!traces&&string(argv[6])!="nauty")throw runtime_error("bad engine");
  auto f=read_family(argv[2]);auto pool=read_pool(argv[3],f.n);
  ifstream q(argv[4]);long long count;q>>count;if(count<0)throw runtime_error("bad query count");
  ofstream out(argv[5]);auto start=chrono::steady_clock::now();long long pass=0,fail=0;int maxnv=0;
  for(long long t=0;t<count;t++){
   int i,j,p,c;q>>i>>j>>p>>c;if(!q||i<1||j<i||j>(int)f.g.size()||c<1||c>(int)pool.size()||(p!=1&&p!=-1))throw runtime_error("bad query");
   auto x=pool[c-1];if(parity(x)!=p){for(int &v:x){if(v==0)v=1;else if(v==1)v=0;}}
   if(parity(x)!=p)throw runtime_error("parity correction failed");
   auto z=check(f.g[i-1],f.g[j-1],x,f.n,traces);maxnv=max(maxnv,z.nv);
   if(z.fixed)++pass;else{++fail;out<<"FAIL "<<i<<" "<<j<<" "<<p<<" "<<c<<" nv="<<z.nv<<" orbits="<<z.orbits<<"\n";}
   if((t+1)%10000==0){double sec=chrono::duration<double>(chrono::steady_clock::now()-start).count();cerr<<"GRAPH_REPLAY progress="<<(t+1)<<"/"<<count<<" pass="<<pass<<" fail="<<fail<<" seconds="<<sec<<"\n";}
  }
  double sec=chrono::duration<double>(chrono::steady_clock::now()-start).count();
  out<<"GRAPH_REPLAY_DONE engine="<<argv[6]<<" degree="<<f.n<<" queries="<<count<<" pass="<<pass<<" fail="<<fail<<" max_vertices="<<maxnv<<" seconds="<<sec<<" nauty_version="<<affine_graph_nauty_version()<<"\n";
  cerr<<"GRAPH_REPLAY_DONE queries="<<count<<" pass="<<pass<<" fail="<<fail<<" seconds="<<sec<<"\n";return fail?1:0;
 }catch(const exception&e){cerr<<"ERROR "<<e.what()<<"\n";return 2;}
}
