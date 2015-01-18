Authors : 
[main author] 		Yasutaka Furukawa furukawa[AT]cs.washington.edu
[main author] 		Jean Ponce Jean.Ponce[AT]ens.fr
[Windows Porting] 	Pierre moulon pmoulon[AT]gmail.com

Special thanks to ASTRE Henri for the PTHREAD 64 bits lib and dll : http://www.visual-experiments.com/

Date : 13 July 2011

--------------------
- Web ressources : - 
--------------------
[CMVS] http://grail.cs.washington.edu/software/cmvs/
[windows porting] http://francemapping.free.fr/Portfolio/Prog3D/CMVS.html


--------------------
-- Compilation  --
--------------------

Windows => Use precompiled binary, or compile it with VS2008/2010 (Express or pro, Pro will allow you to enable Opemp in CMVS)
        => Use CMake GUI in order to generate the Visual Studio project file (in ./program you will find the main CMakeLists.txt).

Linux => use makefile in program/main.

 Or use CMake build system :
=> Install the following libraries : jpeg boost boost-graph
 $ mkdir OutputLinux
 $ cd OutputLinux
 $ cmake . ..
 $ make
 => That's all. Openmp is not activated yet. Add openmp in the cmvs link option and define the _OPENMP cxx flags


--------------------
----  Notes :   ----
--------------------
To test it : ./WinX-VS2010/Readme.txt "Usage of the binaries"

What have been done on native Yasutaka Furukawa source code :

- Create the CMake build system for CMVS/PMVS2 and required library.

- Optimize a little bit JPEG image loading.

- Update CMVS source code in order to compile.

- Add changes from Nghia Ho http://nghiaho.com/
  - memoize pow(2,X)
  - Change GSL simplex to lmfit.

- Replaced GSL simplex/lmfit with nlopt optimizer

- Replaced image loading routines with CImg. Now PPMs are supported properly, with optional support for PNG and TIFF

- Replaced BLAS/LAPACK with Eigen

- Updated internal jpeg library and miniBoost

- CMake-system now supports system boost, jpeg and other libraries if available. 

- Replaced pthread with tinycthread to get rid of pthread.dll on Windows
