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

- See https://github.com/pmoulon/CMVS-PMVS/BUILD.md

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
