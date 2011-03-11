
INCLUDE_DIRECTORIES("../metisLib")
# list sources files
SET(SRC wkkm.c mlkkm.c)


ADD_LIBRARY(multilevelLib
${SRC}
)