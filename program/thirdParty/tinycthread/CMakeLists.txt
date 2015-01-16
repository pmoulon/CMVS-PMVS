PROJECT(tinycthread C CXX)

ADD_LIBRARY(tinycthread STATIC tinycthread.c tinycthread.h rwmutex.cpp rwmutex.h)
