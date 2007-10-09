# ==========================================================================
# YAQCStats class: holds the data necessary for the quality control
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setClass("YAQCStats", 
    contains="QCStats",
    representation(
        average.noise="numeric",
        morespikes="matrix",
        gcos.probes="matrix",
        bio.calls="matrix",
        log="logical"      
    )
)
