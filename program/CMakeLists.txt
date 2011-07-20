# CMake powered build system for CMVS
# Create by Pierre Moulon
# First draft the 4 june 2010
# Second draft 24 sept 2010
# Third draft (linux) 12 oct 2010 

CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

PROJECT(CMVS-PMVS2 C CXX)

IF(WIN32)
SET(ADDITIONAL_WIN_INCLUDES 
      ./thirdParty/miniBoost
      ./thirdParty/pthread/include
      ./thirdParty/jpeg)
ENDIF(WIN32)

INCLUDE_DIRECTORIES(
./program/base/image
./program/base/numeric
./program/base/pmvs

./thirdParty/
./thirdParty/gsl-1.13/include
#http://eris.liralab.it/iCub/downloads/packages/windows/msvc9/
./thirdParty/graclus1.2/metisLib
./thirdParty/lmfit-3.2/lib
${ADDITIONAL_WIN_INCLUDES}
)

#if you build on win64 you have to set to true the following variable
SET(WIN64_PTHREAD false)

IF(WIN32)
IF(WIN64_PTHREAD)
  SET(PTHREAD_LIB ${CMAKE_SOURCE_DIR}/thirdParty/pthreadVC2.lib)
ELSE(WIN64_PTHREAD)
  SET(PTHREAD_LIB ${CMAKE_SOURCE_DIR}/thirdParty/pthread.lib)
ENDIF(WIN64_PTHREAD)

SET(PRECOMPILED_LIBS ${PTHREAD_LIB})

ELSE(WIN32)
  FIND_PACKAGE(Boost REQUIRED)
  SET(PRECOMPILED_LIBS jpeg lapack gsl blas pthread)
ENDIF (WIN32)

ADD_SUBDIRECTORY(base)
ADD_SUBDIRECTORY(main)
ADD_SUBDIRECTORY(thirdParty)

