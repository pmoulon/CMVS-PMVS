
ADD_EXECUTABLE(pmvs2 pmvs2.cc)

IF(WIN32)
TARGET_LINK_LIBRARIES(pmvs2 pmvs_lib image_lib numeric_lib ${PRECOMPILED_LIBS} jpeg gsl cblas lapack lmfit)
ELSE(WIN32)
TARGET_LINK_LIBRARIES(pmvs2 pmvs_lib image_lib numeric_lib ${PRECOMPILED_LIBS} lmfit)
ENDIF(WIN32)


ADD_EXECUTABLE(genOption genOption.cc)

ADD_EXECUTABLE(cmvs cmvs.cc)
IF(WIN32)
TARGET_LINK_LIBRARIES(cmvs graclus_lib cmvs_lib image_lib ${PRECOMPILED_LIBS} jpeg)
ELSE(WIN32)
TARGET_LINK_LIBRARIES(cmvs cmvs_lib graclus_lib image_lib ${PRECOMPILED_LIBS})
ENDIF(WIN32)
