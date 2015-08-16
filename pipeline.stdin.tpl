@< path("ElasticConstantsFirstStrainGradientCubic_TD") >@
@< TESTNAME >@
@< MODELNAME >@
fcc
@< query({
"flat"      : "on", 
"database"  : "data", 
"fields"    : {"_id": 0, "meta.runner._id": 1, "a.source-value": 1}, 
"limit"     : 1, 
"query"     : {"meta.runner._id": {"$regex":"LatticeConstantCubicEnergy"},"meta.subject._id": MODELNAME,"short-name.source-value": "fcc"},
"project"   : ["a.source-value"]
}) >@
@< query({
"flat"      : "on", 
"database"  : "data", 
"fields"    : {"_id": 0, "meta.runner._id": 1, "species.source-value": 1}, 
"limit"     : 1, 
"query"     : {"meta.runner._id": {"$regex":"LatticeConstantCubicEnergy"},"meta.subject._id": MODELNAME,"short-name.source-value": "fcc"},
"project"   : ["species.source-value"]
}) >@
