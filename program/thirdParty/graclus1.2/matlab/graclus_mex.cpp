/*
 * Copyright 2005
 *
 * graclus_mex.c
 *
 * This file contains the driving routine for multi-level kernel k-means
 *
 * Started 12/2004
 * Yuqiang Guan
 *
 * Modified 4/2008
 * Brian Kulis
 *
 */

#include <metis.h>
#include "mex.h"

int boundary_points;
int spectral_initialization;
int cutType; /*cut type, default is normalized cut */
int memory_saving; /* forbid using local search or empty cluster removing */
/*char mlwkkm_fname[256]; */ /*used to store coarsest file*/

/*************************************************************************
* multi-level weighted kernel k-means main function
**************************************************************************/
void mexFunction(int nlhs, mxArray *plhs[],
		 int nrhs, const mxArray* prhs[])
{
  mwSize *jc, *ir;
  double *pr, *out, *out2;
  int edgenum, curredge, numedges, currval;

  int readew, readvw, fmt, ncon, tmp;
  idxtype *xadj, *adjncy, *vwgt, *adjwgt;
  int i, nparts=1, options[11];
  idxtype *part;
  float rubvec[MAXNCON], lbvec[MAXNCON];
  float result;
  GraphType graph;
  char clusterIDFile[256], *program_name;
  int numflag = 0, wgtflag = 0, edgecut, chain_length;
  timer TOTALTmr, METISTmr, IOTmr;
  int no_args = 1, clusteringEva =0, levels;
  idxtype *t1, *t2;

  ir = mxGetIr(prhs[0]);
  jc = mxGetJc(prhs[0]);
  pr = mxGetPr(prhs[0]);
  edgenum = (int) (mxGetScalar(prhs[1]));
  nparts = (int) (mxGetScalar(prhs[2]));
  cutType = (int) (mxGetScalar(prhs[3]));
  chain_length = (int) (mxGetScalar(prhs[4]));
  spectral_initialization = (int) (mxGetScalar(prhs[5]));
  
  memory_saving = 0;
  boundary_points = 0;

  cleartimer(TOTALTmr);
  cleartimer(METISTmr);
  cleartimer(IOTmr);

  starttimer(TOTALTmr);
  starttimer(IOTmr);

  /*printf("edgenum = %d, edgenum2 = %d\n", edgenum, mxGetNumberOfElements(prhs[0]));*/
   readew=1;
   InitGraph(&graph);

   graph.nvtxs = mxGetM(prhs[0]);
   graph.nedges = edgenum;
   xadj = graph.xadj = idxsmalloc(graph.nvtxs+1, 0, "graclus_mex: xadj");
   adjncy = graph.adjncy = idxmalloc(graph.nedges, "graclus_mex: adjncy");
   adjwgt = graph.adjwgt = (readew ? idxmalloc(graph.nedges, "graclus_mex: adjwgt") : NULL);

   ncon = 0;
   fmt = 1;
   ncon = graph.ncon = (ncon == 0 ? 1 : ncon);
   readew = (fmt%10 > 0);
   readvw = ((fmt/10)%10 > 0);
   wgtflag = 0;
   if (readew)
     wgtflag += 1;
   if (readvw)
     wgtflag += 2;

   for(i = 0; i < edgenum; i++)
   {
      adjncy[i] = (idxtype)ir[i];
      adjwgt[i] = (idxtype) pr[i];
   }
   for(i = 0; i <= graph.nvtxs; i++)
      xadj[i] = (idxtype) jc[i];

   /*levels = 5*nparts;*/
   levels = amax((graph.nvtxs)/(40*log2_metis(nparts)), 5*(nparts));
   /*printf("Will coarsen until %d nodes...\n", levels);*/

  if (graph.nvtxs <= 0) {
    printf("Empty graph. Nothing to do.\n");
    exit(0);
  }
  stoptimer(IOTmr);
  printf("\n----------------------------------------------------------------------\n");
  printf("%s", MLKKMTITLE);
  printf("Graph Information:\n");
  printf("  #Vertices: %d, #Edges: %d, ", graph.nvtxs, graph.nedges/2);
  if (graph.ncon > 1)
    printf("  Balancing Constraints: %d\n", graph.ncon);

  part = idxmalloc(graph.nvtxs, "main: part");
  options[0] = 0;
  starttimer(METISTmr);
    printf("#Clusters: %d\n", nparts);
    printf("#local search steps: %d\n", chain_length);
    if (graph.ncon == 1) 
    {
      /*METIS_PartGraphKway(&graph.nvtxs, graph.xadj, graph.adjncy, graph.vwgt, graph.adjwgt, 
	&wgtflag, &numflag, &nparts, options, &edgecut, part); 
      */
      MLKKM_PartGraphKway(&graph.nvtxs, graph.xadj, graph.adjncy, graph.vwgt, graph.adjwgt, 
			  &wgtflag, &numflag, &nparts, &chain_length, options, &edgecut, part, levels);
    
      /* ends */
    }
    else {
      for (i=0; i<graph.ncon; i++)
	rubvec[i] = HORIZONTAL_IMBALANCE;
      /*
	METIS_mCPartGraphKway(&graph.nvtxs, &graph.ncon, graph.xadj, graph.adjncy, graph.vwgt, 
	graph.adjwgt, &wgtflag, &numflag, &nparts, rubvec, options, &edgecut, part);
      */
    }
  stoptimer(METISTmr);
  ComputePartitionBalance(&graph, nparts, part, lbvec);
  if (cutType == 0){
    result = ComputeNCut(&graph, part, nparts);
    printf("\nNormalized-Cut... \n   Cut value: %7f, Balance: ", result);
  }
  else{
    result = ComputeRAsso(&graph, part, nparts);
    printf("\nRatio Association...  \n  Association value: %7f, Balance: ", result);
  }
    
  for (i=0; i<graph.ncon; i++)
    printf("%5.2f ", lbvec[i]);
  printf("\n");

/*  if(clusteringEva  ==0){
    starttimer(IOTmr);
    WritePartition(filename, part, graph.nvtxs, nparts); 
    stoptimer(IOTmr);
    } */
  stoptimer(TOTALTmr);
  
  printf("\nTiming Information:\n");
  printf("  I/O:          \t\t %7.3f\n", gettimer(IOTmr));
  printf("  Clustering:   \t\t %7.3f   (Graclus time)\n", gettimer(METISTmr));
  printf("  Total:        \t\t %7.3f\n", gettimer(TOTALTmr));
  printf("----------------------------------------------------------------------\n");
  plhs[0] = mxCreateDoubleMatrix(graph.nvtxs,1,mxREAL);
  out = mxGetPr(plhs[0]);
  for (i = 0; i < graph.nvtxs; i++)
    out[i] = (double) part[i]+1;
  plhs[1] = mxCreateDoubleMatrix(1,1,mxREAL);
  out2 = mxGetPr(plhs[1]);
  out2[0] = result;
  GKfree((void **) &graph.xadj, (void **) &graph.adjncy, (void **) &graph.vwgt, (void **) &graph.adjwgt, (void **) &part, LTERM);
}  


