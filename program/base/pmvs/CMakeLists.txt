ADD_LIBRARY(pmvs_lib
	detectFeatures.h
	detectFeatures.cc
	detector.h
	detector.cc
dog.h
dog.cc
expand.h
expand.cc
filter.h
filter.cc
findMatch.h
findMatch.cc
harris.h
harris.cc
optim.h
optim.cc
option.h
option.cc
patch.h
patch.cc
patchOrganizerS.h
patchOrganizerS.cc
point.h
point.cc
seed.h
seed.cc
)

set_target_properties(pmvs_lib PROPERTIES
    CXX_STANDARD 17
    CXX_EXTENSIONS OFF)
