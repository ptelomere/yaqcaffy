# ========================================================
# Create a new environment containing the probenames
# required by the yaqc functions. The data is stored in
# the data directory of the package as tab delimited file.
# --------------------------------------------------------

.yaqcEnv <- new.env(parent=emptyenv(), hash=TRUE)

data(morespikes, envir = .yaqcEnv, package="yaqcaffy")



