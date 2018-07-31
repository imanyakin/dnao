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

	#This is the model we will use for fitting data
	model = make_model(temperature)


	#extract autocorrelation and times:
	data = np.asarray(run)
	times = np.asarray(data[0],dtype=float)*1e-6 #we know times are given in microseconds
	acs = np.asarray(data[1],dtype=float) + 1 #add one - we assume the "Correlation Data" is g2

	intercept = acs[0]	

	
	
	
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
	    alpha = 0.0,
	    time_truncation_limit=1e-3,
	    debug=False
	)
	ipsd = contin_pipeline["ipsd"]
	radii = np.linspace(1e-9,300e-9,300)

	outp = {
		"temperature":temperature,
		"acs":acs,
		"times":times,
		"intercept":intercept,
		"ipsd":ipsd,
		"radii":radii
	}

	return outp

for k in file["DLS_Zetasizer"].keys():

	run = file["DLS_Zetasizer"][k] 
	results = analyse(run)
	data = np.asarray(run)

	

	#get temperature of current run:
	attr_keys = run.attrs.keys()
	temperature_pattern = re.compile("Temperature.*")
	temperature_key = [k for k in attr_keys if bool(temperature_pattern.match(k))][0]
	temperature = int(float(run.attrs[temperature_key]))
	
	type_index = int(temperature == 25) #separate into 40oC and 25oC on temperature

	times = np.asarray(data[0],dtype=float)*1e-6 #we know times are given in microseconds
	acs = np.asarray(data[1],dtype=float) + 1

	analysis_struct[type_index]["intercepts"].append(acs[0])
	analysis_struct[type_index]["acs"] = analysis_struct[type_index]["acs"] + [acs]
	analysis_struct[type_index]["times"] = times

	model = make_model(temperature)

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
		    alpha = 0.0,
		    time_truncation_limit=1e-3,
		    debug=False
		)
		ipsd = contin_pipeline["ipsd"]
		analysis_struct[type_index]["ipsds"].append([ipsd])
		analysis_struct[type_index]["contin_radii"] = np.linspace(1e-9,300e-9,300)
	except:
		print "Run:", k, "FAILED"

	
fig,axarr = plt.subplots(4,2,figsize=(3*4,4*2))
for i in [0,1]:
	analysis_struct[i]["mean_acs"] = np.mean(np.asarray([a for a in analysis_struct[i]["acs"] if a[0] < 1.4]),axis=0)


	axarr[0][i].hist(analysis_struct[i]["intercepts"],bins=np.linspace(1.1,2.6,50))
	axarr[0][i].set_title("Autocorrelation Intecepts (T={})".format(analysis_struct[i]["temperature"]))
	
	for acs in analysis_struct[i]["acs"]:
		axarr[1][i].semilogx(analysis_struct[i]["times"],acs)
		axarr[1][i].set_title("Autocorrelations (T={})".format(analysis_struct[i]["temperature"]))

	axarr[2][i].semilogx(analysis_struct[i]["times"],analysis_struct[i]["mean_acs"])
	axarr[2][i].set_title("Mean Autocorrelation (T={})".format(analysis_struct[i]["temperature"]))

	contin_pipeline = basic_contin_pipeline(g2_autocorrelation=analysis_struct[i]["mean_acs"],
	    times=np.asarray(times),
	    scattering_angle=model.theta,
	    temperature=model.T,
	    wavelength=model.wavelength,
	    solvent_viscosity=model.viscosity,
	    solvent_refractive_index=model.n,
	    min_radius = 1e-9,
	    max_radius = 300e-9,
	    radii_steps = 300,
	    alpha = 0,
	    time_truncation_limit=1e-2,
	    debug=False
	)
	ipsd =  contin_pipeline["ipsd"]
	radii =  np.linspace(1e-9,300e-9,300)
	axarr[3][i].plot(radii,ipsd)
	axarr[3][i].set_title("Mean ACS IPSD (T={})".format(analysis_struct[i]["temperature"]))

	

plt.tight_layout()
plt.legend()
plt.show()



