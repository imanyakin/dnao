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





from  multiprocessing import Pool
import sys
pool = Pool(6)

labels = [("25deg_file1_series1_",25),("40deg_file1_series1_",40),("25deg_file2_series1_",25),("40deg_file2_series1_",40),("25deg_file2_series2_",25),("40deg_file2_series2_",40)]

jobs = [(FILES[0],0,50,labels[0][0]),(FILES[0],50,100,labels[1][0]),(FILES[1],0,50,labels[2][0]),(FILES[1],50,100,labels[3][0]),(FILES[1],100,150,labels[4][0]),(FILES[1],150,200,labels[5][0])]
# results = pool.map(f,jobs)

for job in jobs:
	file = job[0]
	min_index = job[1]
	max_index = job[2]
	label = job[3]
	title = "{0}, {1} to {2}, label: {3}".format(file.replace(".csv","").split("/")[-1],min_index,max_index,label)
	results = pool.map(f,[job])

	for i,(r,(label,temp)) in enumerate(zip(results,labels)):

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
			try:
				
				radii,weights = simple_analysis(times=ts,acs=acs,temperature=temp)
				


				ax1.semilogx(ts,acs,label="Run:{}".format(j))
				ax2.plot(radii*1e9, weights,label="Run:{}".format(j))

				ax1.set_xlabel("Time [s]")
				ax1.set_ylabel("Data + 1 = g2?")

				ax2.set_xlabel("Radius [nm]")
				ax2.set_ylabel("Intensity PSD")
				ax1.set_title(title)
				ax2.set_title(title)
			except:
				pass
		accs = accs/float(counts)
		radii,weights = simple_analysis(times=ts,acs=accs,temperature=temp)

		fig2, ax = plt.subplots(1)
		ax.set_title(title)
		ax.plot(radii*1e9,weights, label="Averaged ACS IPSD")

		ax.set_xlabel("Radius [nm]")
		ax.set_ylabel("Intensity PSD")

	ax1.legend()
	ax2.legend()
	ax.legend()
	plt.show()
