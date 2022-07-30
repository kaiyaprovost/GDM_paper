import sys
import glob
import os

species = "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/BBS_data/SpeciesList_delimit.txt"
# Seq,AOU,ORDER,Family,Genus,Species
state = "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/BBS_data/NotSouthwest.csv"
#RouteDataID,CountryNum,StateNum,Route,RPID,Year,AOU,Count10,Count20,Count30,Count40,Count50,StopTotal,SpeciesTotal
routes = "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/BBS_data/routes_delimit.csv"
# CountryNum,StateNum,Route,RouteName,Active,Latitude,Longitude,Stratum,BCR,RouteTypeID,RouteTypeDetailID
strata = "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/BBS_data/BBSStrata_delimit.txt"
# Stratum,StratumName

## AOU links species with state
## CountryNum, Statenum, Roure link state with routes 
## stratum links routes with strata

with open(species,"r") as speciesfile:
	specieslines = speciesfile.readlines()

with open(state,"r") as statefile:
	statelines = statefile.readlines()

with open(routes,"r") as routesfile:
	routeslines = routesfile.readlines()

with open(strata,"r") as stratafile:
	stratalines = stratafile.readlines()

def lines2dict(lines,header,keynums,valuenums,delim):
	## defaults to taking 0th column as key, col 1:2 as values
	dicti = {}
	for i in range(len(lines)):
		line = lines[i].upper()
		split = line.strip().split(delim)
		keys = [split[k].lstrip("0") for k in keynums]
		values = [split[j].lstrip("0") for j in valuenums]
		if header == True and i == 0:
			print("Key: "+" ".join(keys))
			print("Values: "+" ".join(values))
		else:
			dicti[tuple(keys)] = (values)
	return(dicti)

## delimit into dictionaries
speciesdict = lines2dict(specieslines,header=False,keynums=[1],valuenums=[0,2,3,4,5],delim=",")
statedict = lines2dict(statelines,header=False,keynums=[1,2,3,5,6],valuenums=[0,4,7,8,9,10,11,12,13],delim=",")
routesdict = lines2dict(routeslines,header=False,keynums=[0,1,2],valuenums=[3,4,5,6,7,8,9,10],delim=",")
stratadict = lines2dict(stratalines,header=False,keynums=[0],valuenums=[1],delim=",")

## combine these dictionaries accordingly 
## AOU links species with state
## CountryNum, Statenum, Roure link state with routes 
## stratum links routes with strata

## first merge species with state by AOU
for statekey,statedata in statedict.items():
	#print(statedata)
	AOUnum = statekey[4].lstrip("0")
	AOUtup = tuple([AOUnum])
	speciesdata = speciesdict[AOUtup]
	newstatedata = speciesdata+statedata
	statedict[statekey] = newstatedata

## then merge routes with strata by Stratum
for routeskey,routesdata in routesdict.items():
	#print(routesdata)
	Stratanum = routesdata[4] 
	Stratatup = tuple([Stratanum])
	stratadata = stratadict[Stratatup]
	newroutesdata = routesdata[0:4]+stratadata+routesdata[5:]
	routesdict[routeskey] = newroutesdata

## then finally merge states with routes
for statekey,statedata in statedict.items():
	#print(statedata)
	routeskey = statekey[0:3]
	routesdata = routesdict[routeskey]
	newstatedata = statedata+routesdata
	statedict[statekey] = newstatedata

## write out the dictionary
outfilestr = state+"_fulldata.csv"
with open(outfilestr,"a") as outfile:
	for statekey,statedata in statedict.items():
		towrite = ",".join(statekey)+","+",".join(statedata)+"\n"
		outfile.write(towrite)



