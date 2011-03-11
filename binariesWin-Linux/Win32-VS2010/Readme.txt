Author :
[Windows porting]	Pierre moulon 	pmoulon[AT]gmail.com

Date : 10 April 2010

Those binaries are binaries compiled under VS2010.
To launch them you require the VS2010 redistribuable.
http://www.microsoft.com/downloads/en/details.aspx?FamilyID=a7b7a05e-6de6-4d3a-a423-37bf0912db84

Windows porting of original GPL sources of CMVS
http://grail.cs.washington.edu/software/cmvs/

-------------------------
-- Usage of the binaries --
-------------------------

Command specification :
http://grail.cs.washington.edu/software/cmvs/documentation.html
A .bat sample is also given in ./PreCompiledBinary

CMVS prefix maximage[=100] CPU[=4]

prefix : (directory path of bundler.out). Ended with /
maximage : the max image allowed in a given cluster (to change according your ram size that will require PMVS for take in charge the given cluster)
CPU : CPU core of the machine

genOption path

path : Path where CMVS output data have been exported

Go to the path you have specified and lauch the pmvs.bat files.
If you have message 'invalid input parameter' adjust the path in the .bat file.
"pmvs path option-X"


PMVS2.exe "theoptionfilepath/" option.txt
If you have DLL error message you have to install VS2008 redistribuables files.
http://www.microsoft.com/downloads/details.aspx?FamilyID=9B2DA534-3E03-4391-8A4D-074B9F2BC1BF&displaylang=en

// dataset could be download from the PMVS original website:
http://grail.cs.washington.edu/software/cmvs/
Search : (CMVS linux binaries and source codes (two samples dataset are included))


--------------------------------------------------
----- Notes about modifications of the code : ----
--------------------------------------------------

IN PMVS :
- Static array with dynamic value are not accepted by microsoft VS 2008 compiler.

- Use jpeg/gsl/pthread/blas/lapack/f2c precompiled library.

- iterator.begin()-1 could crash on windows stl.

- Optimize a little bit the JPEG loading.

Comments "// Pierre Moulon" or // PM in the code show some of the modifications.

IN CMVS :
- Static array with dynamic value are not accepted by microsoft VS 2008 compiler.

- Fix for M_PI on windows.

- Fix out of size array access.

