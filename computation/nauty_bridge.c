/* C ABI bridge to the unmodified nauty 2.8.8 library. */
#include <stdlib.h>
#include <stddef.h>
#include "nausparse.h"
#include "traces.h"
static int *sort_colours;
static int cmp_vertex(const void *aa,const void *bb){int a=*(const int*)aa,b=*(const int*)bb;if(sort_colours[a]!=sort_colours[b])return sort_colours[a]<sort_colours[b]?-1:1;return (a>b)-(a<b);}
int affine_graph_nauty_version(void){return NAUTYVERSIONID;}
int affine_graph_leaf_fix(int nv,size_t nde,size_t *v,int *d,int *e,int *col,int n,int use_traces,int *norb,unsigned long *nodes){
 int *lab=malloc((size_t)nv*sizeof(int)),*ptn=malloc((size_t)nv*sizeof(int)),*orb=malloc((size_t)nv*sizeof(int));
 unsigned char *seen=calloc((size_t)nv,1);if(!lab||!ptn||!orb||!seen)abort();
 sparsegraph sg={nde,v,nv,d,e,NULL,(size_t)nv,(size_t)nv,nde,0};
 for(int i=0;i<nv;i++)lab[i]=i;
 sort_colours=col;qsort(lab,(size_t)nv,sizeof(int),cmp_vertex);
 for(int i=0;i<nv;i++)ptn[i]=(i<nv-1&&col[lab[i]]==col[lab[i+1]])?1:0;
 SG_DECL(canon);int err=0;
 if(use_traces){DEFAULTOPTIONS_TRACES(opt);opt.defaultptn=FALSE;opt.getcanon=FALSE;TracesStats st={0};Traces(&sg,lab,ptn,orb,&opt,&st,&canon);/* interrupted counts normal refinement pruning, not process interruption. */err=st.errstatus;*norb=st.numorbits;*nodes=st.numnodes;}
 else{DEFAULTOPTIONS_SPARSEGRAPH(opt);opt.defaultptn=FALSE;opt.getcanon=FALSE;statsblk st={0};sparsenauty(&sg,lab,ptn,orb,&opt,&st,&canon);err=st.errstatus;*norb=st.numorbits;*nodes=st.numnodes;}
 SG_FREE(canon);int result=1;
 if(err)result=-1;
 else for(int i=0;i<n;i++){if(orb[i]<0||orb[i]>=nv){result=-2;break;}if(seen[orb[i]]){result=0;break;}seen[orb[i]]=1;}
 free(lab);free(ptn);free(orb);free(seen);return result;
}
