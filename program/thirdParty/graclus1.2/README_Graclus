/*      Graclus -- Efficient graph clustering software for
	normalized cut and ratio association on undirected graphs.
        Copyright(c) 2008 Brian Kulis, Yuqiang Guan (version 1.2)
	kulis@cs.utexas.edu
	http://www.cs/utexas.edu/users/kulis
*/

*To compile graclus: 
 make

NOTE: If you are running on a 64-bit machine, change the COPTIONS line in
Makefile.in to COPTIONS = -DNUMBITS=64

In Windows, we had no trouble compiling with CYGWIN.  Getting the Matlab
interface to work was a little bit more difficult.  Our precompiled mex files
are included in this distribution.  For information about how to compile the
mex in Windows, see the end of this document.

*To remove .o and executable files:
 make clean

*After compiling the program, you can compile the Matlab interface.  In the
 matlab directory, run graclus_compile_script.m.  NOTE that the compiler used
 by mex in Matlab must be the same as the compiler used to compile graclus.
 The Matlab interface includes the option of using spectral methods at the
 coarsest level.  This is typically slower than the region growing method but
 often produces better results, especially for image data.  The Matlab
 interface also makes it easier to integrate with existing Matlab code.
 NOTE: for 64-bit machines, change the -DNUMBITS=32 to -DNUMBITS=64 in each
 line. 

*To cluster a graph into a given number of clusters:
 graclus [options] <GraphFile> <Nparts> 

 options: -o ncut | rassoc
		 ncut --- normalized cut (default)
                 rassoc --- ratio association

	 -l number_of_local_search_steps (default is 0)
	 -b use only boundary points (default is to use all points) 
 <GraphFile> contains matrix in adjacent list format, see an example below
 <Nparts> is number of clusters.

*To compute objective function value for a given graph and clustering:
 graclus [options] -e <clusteringFile> <GraphFile> 

 options: -o ncut | rassoc
		 ncut --- normalized cut (default)
                 rassoc --- ratio association

*Here is an example graph and its representation. When all the edges have the 
 same weight,

 1 ----- 2
 |	 |
 |	 |
 |	 |
 3 ----- 5
  \     /
   \   /
    \ /
     4

 the matrix representation is

 5 6          	<--- # of nodes and edges
 2 3	     	<--- nodes adjacent to 1
 1 5		.
 1 4 5		.
 3 5		.
 2 3 4	     	<--- nodes adjacent to 5

 Otherwise, if edge weights (must be integer values) are different,
      10	
  1 ----- 2
  |	  |
 9|	  |6
  |   7   |
  3 ----- 5
   \     /
  11\   /28
     \ /    
      4

 then the matrix representation becomes

 5 6 1		<--- # of nodes and edges and format
 2 10 3 9	<--- nodes adjacent to 1 and corresponding edge weight
 1 10 5 6	.
 1 9 4 11 5 7	.
 3 11 5 28	.
 2 6 3 7 4 28	<--- nodes adjacent to 5 and corresponding edge weight

*Output:
 The output of the program is a file contains cluster membership in one column,
 which is generated in the current directory.
 
 For example:
 0		<--- node 1 belongs to cluster 0
 0		.
 1		.
 1		.
 1		<--- node 5 belongs to cluster 1

*Some command examples:
 graclus mygraph 100		    <--- 100 clusters, normalized cut, no local search
 graclus -o ncut mygraph 100	    <--- 100 clusters, normalized cut, no local search
 graclus -o rassoc -b mygraph 100      <--- 100 clusters, ratio association, no local search, only cluster using boundary points
 graclus -o ncut -l 20 mygraph 100  <--- 100 clusters, normalized cut, number of local search steps= 20  
 graclus -e cluster_file mygraph    <--- compute normalized cut value for the given clustering of mygraph
 graclus -o rassoc -e cluster_file mygraph <--- compute ratio association value for the given clustering of mygraph
 
*Matlab interface:

function [partition, obj] = graclus(G, k, cutType, l, spectral)
%% Graclus Matlab interface
% G: Graph adjacency matrix
% k: number of clusters
% cutType: 0 for NCut, 1 for RAssoc (default is NCut)
% l: number of local search steps (default is 0)
% spectral: 1 for spectral clustering at coarsest level, 0
% otherwise (default is not to use spectral clustering)

Getting the Matlab Interface to Work in Windows
===============================================

Below is what I did to get the Matlab interface working in Windows using
CYGWIN.  I'm sure there is a better way to do this, but what I did seemed to
work.  I cannot guarantee that this will work in all cases.

The simplest thing to do is to try the version of graclus_mex.mexw32 that is
included in the release.  If this works and you do not need to modify
Graclus, then nothing more needs to be done.  

Otherwise, take the following steps:

1. Change the CC line in Makefile.in from the main directory to the one that
is commented.  ( CC = g++ -mno-cygwin -mwindows )

2. Change graclus_compile_script.m in the matlab directory to the commented
version.  

3. Download gnumex and run it.  

4. Change the mexopts.bat file so that the C compiler is g++, not gcc.
Change the Makefile in the metisLib directory so that the INCLUDES line is
INCLUDES = -I. -I/usr/include/mingw
and the LD line to
LD = $(CC) -L. -L/usr/lib
Similarly, change the Makefile in the multilevelLib directory so that
INCLUDES = -I../metisLib -I/usr/include/mingw
LD = $(CC) -L. -L/usr/lib
And the Makefile in the programs directory so that
INCLUDES = -I/usr/include/mingw -I../multilevelLib -I../metisLib
LD = $(CC) $(LDOPTIONS) -L/usr/lib L. -L..

5. Now run make in the main graclus directory (you probably want to run make
clean and make realclean first to remove any old files).  It should compile
successfully. 

6. In the main graclus directory, run cp libmetis.a metis.lib and cp
libmultilevel.a multilevel.lib.  Then in Matlab, run the
graclus_compile_script.  This should compile graclus_mex.mexw32 correctly.

This seems to get things working in Matlab for me.
