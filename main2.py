import numpy as np 
import csv,sys,re 
import matplotlib.pyplot as plt 


from dls.utils import truncate_dataset
from dls.algorithms.convex.contin import contin_g1_2,contin_g1, contin
from dls.algorithms.convex.generate_data import make_g1_matrix
from dls.algorithms.convex.test_contin import L2_matrix

from dls.algorithms.cumulants.cumulants import cumulants_fitting, g2_2nd_order,g2_4th_order
from dls.model import ExperimentModel, DEFAULT_TEMP, DEFAULT_REFRACTIVE_INDEX,DEFAULT_WAVELENGTH, DEFAULT_VISCOSITY, DEFAULT_THETA


FILES = ["/home/ilya/workspace/dnao/single_measurements_reduced.csv","/home/ilya/workspace/dnao/single_measurements2_reduced.csv"]

print "------------------------------STARTING-------"

#Keys to ignore printing
nKOI = [
				"Cumulants Residuals.*",
				"Residuals.*",
				"Oversize By Number.*",
				"Undersize By Number.*",
				"Oversize By Intensity.*",
				"Undersize By Intensity.*",
				"Oversize By Volume.*",
				"Undersize By Volume.*",
				"Diffusions.*",
				"Mark-Houwink Mol..*",
				"Relaxation Times.*",
				"Cumulants.*",
				"Distribution Fit.*",
				"Numbers.*",
				"Volumes.*",
				"Intensities.*",
				"Sizes.*"
				]


def get_acs_from_row(keys,row):

	time_pattern = re.compile("Correlation Delay Times.*")
	acs_pattern = re.compile("Correlation Data.*")
		
	times = [-1]*193
	acs = [-1]*193
	for i,key in enumerate(keys):
		# print key
		# assert(len(keys) == len(row))
		for r in row:
			value = row[i]
			if bool(time_pattern.match(key))==True:
				# print "KEY",key
				replaced = (key.split("[")[1]).split("]")[0]
				
				index = int(replaced)
				times[index] = float(value)

			elif bool(acs_pattern.match(key))==True:
				# print "KEY",key
				index = int((key.split("[")[1]).split("]")[0])
				acs[index] = float(value)
				
	return np.asarray(times[1:])*1e-6,np.asarray(acs[1:])	

def f(arg):
	i,rows,keys = arg[0],arg[1],arg[2]
	data = rows[i]
	times, acs = get_acs_from_row(keys[0:len(data)],data)
	print "i:",i
	return (times,acs)

def get_correlations_from_file(filepath,min_index,max_index,temperature):

	with open(filepath, 'r') as csvfile:
		dialect = csv.Sniffer().sniff(csvfile.read(1024), delimiters=";,\t")
		csvfile.seek(0)
		reader = csv.reader(csvfile, dialect)
		rows = [r for r in reader]
		keys = rows[0]

		rows = rows[1:len(rows)]
		
		tss = []
		acss = []
		
		
		for i,r in enumerate(rows[min_index:max_index]):
			print "i:",i
			ts, acs = get_acs_from_row(keys,r)
			tss = tss + [np.asarray(ts)]
			acss = acss + [np.asarray(acs)]
			fig, ax = plt.subplots(1,figsize=(4*3,3*3))
			ax.set_title("{0}run{1}".format(str(temperature),i))
			ax.semilogx(tss[-1],acss[-1],'x-',label=filepath.split("/")[-1]+"_images_{0}.png".format(str(i)))
			ax.set_xlabel("Correlation Time x 1e-6")
			ax.set_xlabel("Correlation Data")
			subfolders = temperature.split("_")
			image_name = "/home/ilya/workspace/dnao/images/{0}/{1}/{2}/run_{3}.png".format(subfolders[0],subfolders[1],subfolders[2],i%51)
			fig.savefig(image_name)
			plt.close(fig)
		intercepts = [a[0] for a in acss]


		mean_acs = np.mean(acss,axis=0)
		return tss, acss, intercepts, mean_acs

# datasets = 
def f(arg):
	file = arg[0]
	min_index = arg[1]
	max_index = arg[2]
	temp = arg[3]
	return get_correlations_from_file(file,min_index,max_index,temp)


def contin_analysis(g12s, times, radii, thetas,wavelength,alpha):

    y = np.asarray(g12s)

    N = len(radii)
    M = len(thetas)
    T = len(times)

    #make g1 matrix for the expected radial ranges
    g1 = make_g1_matrix(times=np.asarray(times), thetas=np.asarray(thetas), radii = radii,wavelength=wavelength, verbose = False)

    #make square of g1 matrix
    A = g1**2

    R = L2_matrix(radii.shape[0])
    r = np.zeros(R.shape[0])
    Me = np.identity(y.shape[1])

    #solve problem
    solution = contin(A=A,y=y,alpha=alpha,R=R,r=r,Me=Me,debug=False)
    
    #return outputs
    output_weights = np.asarray([solution[i,0] for i in range(len(radii))])

    return output_weights    

def simple_analysis(times, acs,theta=np.pi*(173.0/180.0),alpha=0.3,wavelength=633e-9,threshold_active=False, temperature=temp_deg):

    radii  = np.linspace(1e-9,3e-7,300)
    
    #Pull data from file
    AU_MODEL = ExperimentModel(T=temp_deg+273.15,n=2.1,wavelength=wavelength,viscocity=DEFAULT_VISCOSITY,theta=theta,verbose=False)
  
    #truncate dataset
    truncated = truncate_dataset(times=times,values = acs, limit = 1e-2)
    times, acs = truncated["times"], np.asarray(truncated["values"])
    cumulants_parameters,is4thOrder = cumulants_fitting(times=times,g2=acs,debug=True)

    #cumulant fitting
    if is4thOrder == True:
        [B,beta,Gamma,mu2,_,_]  = cumulants_parameters
        acsm = np.asarray([g2_4th_order(t, B,beta,Gamma,mu2,mu3,k4) for t in times])
    elif is4thOrder == False:
        [B,beta,Gamma,mu2]  = cumulants_parameters
        acsm = np.asarray([g2_2nd_order(t,B,beta,Gamma,mu2) for t in times])

    radius = AU_MODEL.calc_HydroRadius(Gamma)
    diameter = 2.0*radius 

    # return diameter
    print "PARAMETERS [g2(0): {0:.3g}, beta: {1}], Baseline:{2}".format(acs[0], beta,B)
    g12 = (np.asarray(acs)-B)/beta

    weights = contin_analysis(g12s=np.asarray([g12]),times=times,radii=radii, thetas=np.asarray([theta]),wavelength=wavelength,alpha=alpha)

    return radii,weights
    




from  multiprocessing import Pool
import sys
pool = Pool(6)

labels = ["25deg_file1_series1_","40deg_file1_series1_","25deg_file2_series1_","40deg_file2_series1_","25deg_file2_series2_","40deg_file2_series2_"]

jobs = [(FILES[1],0,10,labels[1])]#,(FILES[0],0,50,labels[1]),(FILES[1],0,50,labels[2]),(FILES[1],50,100,labels[3]),(FILES[1],100,150,labels[4]),(FILES[1],150,200,labels[5])]
results = pool.map(f,jobs)

for (i,(r,label)) in enumerate(zip(results,labels)):

	fig, [ax1,ax2] = plt.subplots(2,figsize=(4*4*2,3*4))
		
	accs = 0
	counts = 0
	for j,(ts,acs) in enumerate(zip(r[0],r[1])):
		acs = acs + 1
		
		if type(accs) == int:
			accs = acs
			 
		else:
			accs = accs + acs
		counts = counts + 1
		radii,weights = simple_analysis(times=ts,acs=acs)
		


		ax1.semilogx(ts,acs,label="Run:{}".format(j))
		ax2.plot(radii*1e9, weights,label="Run:{}".format(j))

		ax1.set_xlabel("Time [s]")
		ax1.set_ylabel("Data + 1 = g2?")

		ax2.set_xlabel("Radius [nm]")
		ax2.set_ylabel("Intensity PSD")

	accs = accs/float(counts)
	radii,weights = simple_analysis(times=ts,acs=accs)

	fig2, ax = plt.subplots(1)
	ax.plot(radii*1e9,weights, label="Averaged ACS IPSD")

	ax.set_xlabel("Radius [nm]")
	ax.set_ylabel("Intensity PSD")

ax1.legend()
ax2.legend()
ax.legend()
plt.show()
