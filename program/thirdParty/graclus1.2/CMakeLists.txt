
PROJECT(graclus_lib CXX)
FILE(GLOB LIBRARY_FILES_C "metisLib/*.c")
FILE(GLOB LIBRARY_FILES_H "metisLib/*.h")  

ADD_LIBRARY(graclus_lib
${LIBRARY_FILES_C}
${LIBRARY_FILES_H}
./multilevelLib/wkkm.c
./multilevelLib/mlkkm.c
)
