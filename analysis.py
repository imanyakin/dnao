import numpy as np, re,sys
import matplotlib.pyplot as plt
from dls.model import ExperimentModel
from dls.pipelines import cumulant_pipeline, basic_contin_pipeline
from nplab import datafile as df
FILEPATH = "/home/ilya/workspace/dnao/data/extracted_file.h5"

file = df.DataFile(FILEPATH,"r")

def make_model(temperature):
	refractive_index = 1.33 #refractive index of water
	viscosity = 1.002*1e-3 #water viscosity at 20C
	return ExperimentModel(
		T=temperature+273.15,
		n=refractive_index,
		wavelength=633e-9,
		viscosity=viscosity,
		theta=(170*np.pi)/180.0
	)

def analyse(run):

	#get temperature of current run:
	attr_keys = run.attrs.keys()
	temperature_pattern = re.compile("Temperature.*")
	temperature_key = [k for k in attr_keys if bool(temperature_pattern.match(k))][0]
	temperature = float(run.attrs[temperature_key])

	measured_countrate_pattern = re.compile("Mean Count Rate.*")
	mean_countrate_key = [k for k in attr_keys if bool(measured_countrate_pattern.match(k))][0]
	mean_count_rate = float(run.attrs[mean_countrate_key])

	#This is the model we will use for fitting data
	model = make_model(temperature)


	#extract autocorrelation and times:
	data = np.asarray(run)
	times = np.asarray(data[0],dtype=float)*1e-6 #we know times are given in microseconds
	acs = np.asarray(data[1],dtype=float) + 1 #add one - we assume the "Correlation Data" is g2

	intercept = acs[0]	

	ipsd = None
	radii = np.linspace(1e-9,300e-9,300)

	try:
		contin_pipeline = basic_contin_pipeline(g2_autocorrelation=np.asarray(acs),
		    times=np.asarray(times),
		    scattering_angle=model.theta,
		    temperature=model.T,
		    wavelength=model.wavelength,
		    solvent_viscosity=model.viscosity,
		    solvent_refractive_index=model.n,
		    min_radius = 1e-9,
		    max_radius = 300e-9,
		    radii_steps = 300,
		    alpha = 0.1,
		    time_truncation_limit=6e-4,
		    debug=False
		)
		ipsd = contin_pipeline["ipsd"]
	except Exception,e:
		print e

	print times
	outp = {
		"temperature":temperature,
		"acs":acs,
		"times":times,
		"intercept":intercept,
		"ipsd":ipsd,
		"radii":radii,
		"model": model,
		"mean_count_rate":mean_count_rate
	}

	return outp

plot_columns = 2
analysis_struct = [
	{"acss":[],"times":None,"runs":[],"intercepts":[],"count_rates":[],"ipsds":[],"radii":None,"temperature":None,"model":None},
	{"acss":[],"times":None,"runs":[],"intercepts":[],"count_rates":[],"ipsds":[],"radii":None,"temperature":None,"model":None}
]

for k in file["DLS_Zetasizer"].keys()[0:100]:

	run = file["DLS_Zetasizer"][k] 
	result = analyse(run)
	print "-"*100,
	print result
	print "-"*100
	temperature = result["temperature"]


	if int(temperature) == 25:
		plot_column = 0
	elif int(temperature) == 40:
		plot_column = 1

	analysis_struct[plot_column].update({"runs": analysis_struct[plot_column]["runs"] + [run]})
	analysis_struct[plot_column]["intercepts"].append(result["intercept"])
	analysis_struct[plot_column]["count_rates"].append(result["mean_count_rate"])

	analysis_struct[plot_column].update({"acss" : analysis_struct[plot_column]["acss"] + [ result["acs"]]})
	if result["ipsd"] is not None:

		analysis_struct[plot_column].update({"ipsds" :  analysis_struct[plot_column]["ipsds"] + [result["ipsd"]] })

	
	analysis_struct[plot_column].update({"times": result["times"]})
	analysis_struct[plot_column]["radii"] = result["radii"]
	analysis_struct[plot_column]["model"] = result["model"]

print "-"*50,len(analysis_struct[0]["ipsds"])
print "-"*50,len(analysis_struct[1]["ipsds"])
for i in [0,1]:
	analysis_struct[i].update({"averaged_acs":  np.mean(analysis_struct[i]["acss"],axis=0) })
	analysis_struct[i].update({"averaged_ipsd" : np.mean(analysis_struct[i]["ipsds"],axis=0)})

	model = analysis_struct[i]["model"] 

	times = analysis_struct[i]["times"]
	acs = analysis_struct[i]["averaged_acs"]


	print "times.shape",times.shape
	print "acs.shape",acs.shape
	print acs
	
	contin_pipeline = basic_contin_pipeline(g2_autocorrelation=np.asarray(acs),
	    times=np.asarray(times),
	    scattering_angle=model.theta,
	    temperature=model.T,
	    wavelength=model.wavelength,
	    solvent_viscosity=model.viscosity,
	    solvent_refractive_index=model.n,
	    min_radius = 1e-9,
	    max_radius = 300e-9,
	    radii_steps = 600,
	    alpha = 0.1,
	    time_truncation_limit=6e-4,
	    debug=False
	)

	analysis_struct[i]["averaged_acs_ipsd"] = contin_pipeline["ipsd"]
	analysis_struct[i]["averaged_acs_ipsd_fit"] = contin_pipeline["acs_fit"]


	
fig,axarr = plt.subplots(7,2,figsize=(16*2,16*7))
for i in [0,1]:

	#count rate distribution
	axarr[0][i].hist(analysis_struct[i]["count_rates"],bins=100)
	axarr[0][i].set_title("Mean Count Rate (True, not derived)")
	axarr[0][i].set_xlabel("Count rate (kcps)")
	axarr[0][i].set_ylabel("Occurenace frequency")

	#histogram of intercepts
	axarr[1][i].hist(analysis_struct[i]["intercepts"],bins=np.linspace(1.1,2.6,50))
	axarr[1][i].set_title("Intecepts (T={})".format(analysis_struct[i]["temperature"]))
	axarr[1][i].set_xlabel("$g^{(2)}$(0)-1")
	axarr[1][i].set_ylabel("Occurenace frequency")


	#individual autocorrelations
	for acs in analysis_struct[i]["acss"]:
		axarr[2][i].semilogx(analysis_struct[i]["times"],acs)
	axarr[2][i].set_title("Individual Autocorrelations (T={})".format(analysis_struct[i]["temperature"]))
	axarr[2][i].set_xlabel("Time [s]")
	axarr[2][i].set_ylabel("$g^{(2)}(t)$")
	
	#individual IPSDs
	for ipsd in analysis_struct[i]["ipsds"]:
		axarr[3][i].plot(analysis_struct[i]["radii"],ipsd)
	axarr[3][i].set_title("Individual IPSDs (T={}) [alpha:0.1,truncation_time:6e-4]".format(analysis_struct[i]["temperature"]))
	axarr[3][i].set_xlabel("Hydrodynamic radius [m]")
	axarr[3][i].set_ylabel("IPSD")
	

	#averaged individual IPSD:
	axarr[4][i].plot(analysis_struct[i]["radii"],analysis_struct[i]["averaged_ipsd"])
	axarr[4][i].set_title("Averaged IPSD (T={})[alpha:0.1,truncation_time:6e-4]".format(analysis_struct[i]["temperature"]))
	axarr[4][i].set_xlabel("Hydrodynamic radius [m]")
	axarr[4][i].set_ylabel("IPSD")
	


	#averaged autocorrelation	
	axarr[5][i].semilogx(analysis_struct[i]["times"],analysis_struct[i]["averaged_acs"])
	axarr[5][i].set_title("Averaged ACS (T={})".format(analysis_struct[i]["temperature"]))
	axarr[5][i].set_xlabel("Time [s]")
	axarr[5][i].set_ylabel("$g^{(2)}(t)$")

	contin_pipeline = basic_contin_pipeline(g2_autocorrelation=analysis_struct[i]["averaged_acs"],
	    times=np.asarray(times),
	    scattering_angle=model.theta,
	    temperature=model.T,
	    wavelength=model.wavelength,
	    solvent_viscosity=model.viscosity,
	    solvent_refractive_index=model.n,	    min_radius = 1e-9,
	    max_radius = 300e-9,
	    radii_steps = 600,
	    alpha = 0.1,
	    time_truncation_limit=6e-4,
	    debug=False
	)
	ipsd =  contin_pipeline["ipsd"]
	radii =  np.linspace(1e-9,300e-9,600)
	axarr[6][i].plot(radii,ipsd)
	axarr[6][i].set_title("Averaged ACS IPSD (T={})[alpha:0.1,truncation_time:6e-4]".format(analysis_struct[i]["temperature"]))
	axarr[6][i].set_xlabel("Hydrodynamic radius [m]")
	axarr[6][i].set_ylabel("IPSD")
	

plt.tight_layout()
plt.legend()
plt.savefig("analysis.png")



