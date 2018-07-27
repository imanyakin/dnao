import numpy as np 
import csv,sys,re 
import matplotlib.pyplot as plt 

FILES = ["/home/ilya/workspace/dnao/raw_data_test_open.csv","/home/ilya/workspace/dnao/raw_data_test_closed.csv"]

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
def get_acs_from_file(filepath):
	with open(filepath,'r') as f:

		reader = csv.reader(f)
		rows = []
		for row in reader:
			r = str(row).replace("\\t","||").split("||")
			rows = rows + [r]


		time_pattern = re.compile("Correlation Delay Times.*")
		acs_pattern = re.compile("Correlation Data.*")
		times = [-1]*193
		acs = [-1]*193
		metadata_outp = {}
		for i in range(1,min(len(rows[1]),len(rows[0]))):
			key = rows[0][i]
			value = rows[1][i]

			if bool(time_pattern.match(key))==True:
				# print "KEY",key
				replaced = (key.split("[")[1]).split("]")[0]
				
				index = int(replaced)
				times[index] = float(value)

			elif bool(acs_pattern.match(key))==True:
				# print "KEY",key
				index = int((key.split("[")[1]).split("]")[0])
				acs[index] = float(value)
				
			else:
				
				ign = [re.compile(k) for k in nKOI]
				if any([bool(i.match(key))==True for i in ign]) == False:
					# print key
					pass
				
				# KOI = ["Intecept","Measured Intercept","Signal To Noise Ratio","Derived Count Rate (kcps)","Mean Count Rate (kcps)","Number of CONTIN Runs"]
				metadata = [
				"Intercept",
				"Mean Count Rate (kcps)",
				"Derived Count Rate (kcps)",
				"Attenuation Factor",
				"Measurement Status",
				"Duration Used (s)",
				"Attenuator",
				"Measurement Position (mm)",
				"Data Filtering Factor",
				"Probability to Reject",
				"Size Merit",
				"Measured Size Baseline",
				"Concentration (%)",
				"Cell Compensation",
				"Temperature (\\xb0C)"

				]
		
				
				if key in metadata:
					metadata_outp.update({key:value})
		
		return np.asarray(times[1:])/1e-6,np.asarray(acs[1:]), metadata_outp
		

print "------{}-----".format(FILES[0])
t1,acs1,meta= get_acs_from_file(FILES[0])
for it in meta.items():
	print it
t2,acs2,meta2= get_acs_from_file(FILES[1])

print "------{}-----".format(FILES[1])
for it in meta2.items():
	print it
plt.semilogx(t1,acs1,label=FILES[0].split("/")[-1])
# plt.semilogx(t2,acs2,label=FILES[1].split("/")[-1])

plt.xlabel("Time $\\tau $ [$\\mu s]$")
plt.ylabel("$g^{(1)}(\\tau)$")

plt.legend()
# plt.show()