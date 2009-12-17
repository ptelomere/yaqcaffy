##========================================
## Generic methods for YAQCStats objects
## - - - - - - - - - - - - - - - - - - - -
setGeneric("avns",                function(object)     standardGeneric("avns"))
setGeneric("moreSpikeInProbes",   function(object)     standardGeneric("moreSpikeInProbes"))
setGeneric("gcosProbes",          function(object)     standardGeneric("gcosProbes"))
setGeneric("bioCalls",            function(object)     standardGeneric("bioCalls"))
setGeneric("arrays",              function(object)     standardGeneric("arrays"))
setGeneric("isLog",               function(object)     standardGeneric("isLog"))
setGeneric("yaqc",                function(object,...) standardGeneric("yaqc"))
setGeneric("getYaqcControlProbes",function(object)     standardGeneric("getYaqcControlProbes"))
setGeneric("objectVersion",       function(object)     standardGeneric("objectVersion"))


##================================================
## Generic methods for YaqcControlProbes objects
## - - - - - - - - - - - - - - - - - - - - - - - -
setGeneric("bio",   function(object)       standardGeneric("bio"))
setGeneric("spk",   function(object)       standardGeneric("spk"))
setGeneric("deg",   function(object)       standardGeneric("deg"))
setGeneric("info",  function(object)       standardGeneric("info"))
setGeneric("info<-",function(object,value) standardGeneric("info<-"))
