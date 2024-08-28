"""
Functions for SISTR simulations
"""

def RunSimulation(sval=0):
	afreqs = [] # TODO
	res = {}
	res["afreqs"] = afreqs
	res["afreqs_string"] = ",".join([str(item) for item in afreqs])
	return res