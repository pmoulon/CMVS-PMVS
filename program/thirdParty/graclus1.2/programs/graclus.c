/*
 * Copyright 2005
 *
 * mlkkm.c
 *
 * This file contains the driving routine for multi-level kernel k-means
 *
 * Started 12/2004
 * Yuqiang Guan
 *
 * Modified 9/2006
 * Brian Kulis
 *
 */

#include <metis.h>
//#include <io.h>

int boundary_points;
int spectral_initialization;
int cutType; //cut type, default is normalized cut
int memory_saving; // forbid using local search or empty cluster removing
//char mlwkkm_fname[256]; //used to store coarsest file

/*************************************************************************
* multi-level weighted kernel k-means main function
**************************************************************************/
main(int argc, char *argv[])
{
  int i, nparts=1, options[11];
  idxtype *part;
  float rubvec[MAXNCON], lbvec[MAXNCON];
  float result;
  GraphType graph;
  char filename[256], clusterIDFile[256], *program_name;
  int numflag = 0, wgtflag = 0, edgecut, chain_length;
  timer TOTALTmr, METISTmr, IOTmr;
  int no_args = 1, clusteringEva =0, levels = 0;

  cutType = 0;
  memory_saving =0;
  chain_length = 0;
  boundary_points = 0;
  spectral_initialization = 0;
  program_name = argv[0];

  for (argv++; *argv != NULL; argv++){
    if ((*argv)[0] == '-')
      switch ((*argv)[1])
	{
	case 'e':
	case 'E':
	  strcpy(clusterIDFile, *(++argv));
	  clusteringEva = 1;
	  break;
	case 'o':
	case 'O':
	  if(strcmp(*(++argv), "ncut") == 0)
	    cutType = NCUT;
	  else if(strcmp(*argv, "rassoc") == 0)
	    cutType = RASSO;
	  else{
	    printf("Invalid option %s\n", *argv);
	    print_help(program_name);
	    exit(0);
	  }
	  break;
	case 'l':
	case 'L':
	  chain_length = atoi(*(++argv));
	  break;
	case 'm':
	case 'M':
	  memory_saving = 1;
	  break;
	//case 'c':
	//case 'C':
	//  levels = atoi(*(++argv));
	//  break;
        case 'b':
        case 'B':
	  boundary_points = 1;
	  break;
        //case 's':
	//case 'S':
        //  spectral_initialization = 1;
        //  break;
	default:
	  printf("Invalid switch %s\n", *argv);
	  print_help(program_name);
	  exit(0);
	}
    else{
      strcpy(filename, *argv);
      if(clusteringEva ==0)
	nparts = atoi(*(++argv));
      no_args = 0;
    }
  }
  if(no_args >0){
    print_help(program_name);
    exit(1);
  }
  /*
    if (argc <3 || argc>5) {
    //printf("Usage: %s <GraphFile> <Nparts> [<Cut Type> <chain length>]\n",argv[0]);
    print_help(argv[0]);
    exit(0);
    }
    else if(argc ==3){
    strcpy(filename, argv[1]);
    nparts = atoi(argv[2]);
    }
    else if(argc==4){
    strcpy(filename, argv[1]);
    nparts = atoi(argv[2]);
    cutType = atoi(argv[3]);
    }
    else{
    strcpy(filename, argv[1]);
    nparts = atoi(argv[2]);
    cutType = atoi(argv[3]);
    chain_length = atoi(argv[4]);
    }
  */

  
  if ((clusteringEva == 0) && (nparts < 2)) {
    printf("The number of partitions should be greater than 1!\n");
    exit(0);
  }
  
  //extractfilename(filename, mlwkkm_fname);
  //strcat(mlwkkm_fname, argv[2]);
  //strcat(mlwkkm_fname, *argv);



  cleartimer(TOTALTmr);
  cleartimer(METISTmr);
  cleartimer(IOTmr);

  starttimer(TOTALTmr);
  starttimer(IOTmr);
  ReadGraph(&graph, filename, &wgtflag);
  if (graph.nvtxs <= 0) {
    printf("Empty graph. Nothing to do.\n");
    exit(0);
  }
  stoptimer(IOTmr);
  if(levels == 0)
	  levels = amax((graph.nvtxs)/(40*log2_metis(nparts)), 20*(nparts));
  //levels = graph.nvtxs/128;
  printf("\n----------------------------------------------------------------------\n");
  printf("%s", MLKKMTITLE);
  printf("Graph Information:\n");
  printf("  Name: %s, \n  #Vertices: %d, #Edges: %d, ", filename, graph.nvtxs, graph.nedges/2);
  if (graph.ncon > 1)
    printf("  Balancing Constraints: %d\n", graph.ncon);

  part = idxmalloc(graph.nvtxs, "main: part");
  options[0] = 0;
  starttimer(METISTmr);
  if(clusteringEva >0){
	  //nparts = readClustering(clusterIDFile, part, graph.nvtxs);
    printf("#Clusters: %d\n\nClustering file:\n  %s\n", nparts, clusterIDFile);
  }
  else{
    printf("#Clusters: %d\n", nparts);
    if (graph.ncon == 1) {
      // modification 
      /*METIS_PartGraphKway(&graph.nvtxs, graph.xadj, graph.adjncy, graph.vwgt, graph.adjwgt, 
	&wgtflag, &numflag, &nparts, options, &edgecut, part); 
      */
      MLKKM_PartGraphKway(&graph.nvtxs, graph.xadj, graph.adjncy, graph.vwgt, graph.adjwgt, 
			  &wgtflag, &numflag, &nparts, &chain_length, options, &edgecut, part,levels);
    
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
  }
  stoptimer(METISTmr);
  ComputePartitionBalance(&graph, nparts, part, lbvec);
 
  if (cutType == NCUT){
    result = ComputeNCut(&graph, part, nparts);
    printf("\nNormalized-Cut... \n   Cut value: %7f, Balance: ", result);
  }
  else{
    result = ComputeRAsso(&graph, part, nparts);
    printf("\nRatio Association...  \n  Association value: %7f, Balance: ", result);
  }
    
  //printf("  %d-way Edge-Cut: %7d\n", nparts, ComputeCut(&graph, part));
  for (i=0; i<graph.ncon; i++)
    printf("%5.2f ", lbvec[i]);
  printf("\n");

  if(clusteringEva  ==0){
    starttimer(IOTmr);
    WritePartition(filename, part, graph.nvtxs, nparts); 
    stoptimer(IOTmr);
  }
  stoptimer(TOTALTmr);
  

  printf("\nTiming Information:\n");
  printf("  I/O:          \t\t %7.3f\n", gettimer(IOTmr));
  printf("  Clustering:   \t\t %7.3f   (Graclus time)\n", gettimer(METISTmr));
  printf("  Total:        \t\t %7.3f\n", gettimer(TOTALTmr));
  printf("----------------------------------------------------------------------\n");
  FILE *fp;
  /*fp = fopen("results.txt", "a");
  fprintf(fp,"%7f %7.3f\n", result, gettimer(METISTmr));
  fclose(fp);*/

  GKfree((void **) &graph.xadj, (void **) &graph.adjncy, (void **) &graph.vwgt, (void **) &graph.adjwgt, (void **) &part, LTERM);
}  


