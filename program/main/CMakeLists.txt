
ADD_EXECUTABLE(pmvs2 pmvs2.cc)

TARGET_LINK_LIBRARIES(pmvs2 pmvs_lib image_lib numeric_lib ${PMVS_LIBRARIES} nlopt tinycthread)

ADD_EXECUTABLE(genOption genOption.cc)

ADD_EXECUTABLE(cmvs cmvs.cc)
TARGET_LINK_LIBRARIES(cmvs cmvs_lib image_lib graclus_lib ${PMVS_LIBRARIES} tinycthread)

INSTALL(TARGETS genOption DESTINATION bin/)
INSTALL(TARGETS cmvs DESTINATION bin/)
INSTALL(TARGETS pmvs2 DESTINATION bin/)
