import numpy as np 
import csv,sys,re 
import matplotlib.pyplot as plt 

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

def get_correlations_from_file_by_header(filepath,min_index,max_index):

	def extract_by_indices(iterable, indices):
		return [iterable[i] for i in indices]


	def extract_by_koi(header,dataset,koi):

		patterns = [re.compile(k) for k in koi]
		matched_headers = []
		for h in header:
			if np.sum([bool(p.match(h)) for p in patterns]) >= 1:
				matched_headers = matched_headers + [h]
			else:
				if h == 'Derived Count Rate (kcps)':
					print np.sum([bool(p.match(h)) for p in [re.compile("Derived Count Rate.*")]])
		
		dataset_indices = [header.index(m) for m in matched_headers]

		outp = []
		for data in dataset:
			outp_data = [None]*len(dataset_indices)
			for i in range(len(dataset_indices)):
				if dataset_indices[i] < len(data):
					outp_data[i] = data[dataset_indices[i]]
			outp = outp + [outp_data]
		return matched_headers,outp

	with open(filepath, 'r') as csvfile:

		KOI = [
			"Correlation Data.*",
			"Correlation Delay Times.*",
			"Intercept.*",
			"Mean Count Rate.*",
			"Derived Count Rate.*",
			"Temperature.*"
		]

		reader = csv.reader(csvfile,delimiter="\t")

		header = next(reader)
		rows = list(reader) #get all rows of file
		print min_index,max_index
		rows = [rows[i] for i in range(min_index,max_index)] #extract rows we want
		
		all_headers,all_data = extract_by_koi(header,rows,KOI)

		acs_koi = [KOI[0]]
		times_koi = [KOI[1]]
		metadata_koi = KOI[2:]

		acs_header,acs_data = extract_by_koi(all_headers,all_data,acs_koi)
		times_header,times_data = extract_by_koi(all_headers,all_data,times_koi)
		metadata_header,metadata = extract_by_koi(all_headers,all_data,metadata_koi)
		
		assert(len(acs_data)==len(times_data)==len(metadata))

		dataset_size = len(acs_data)

		outp = []
		for i in range(dataset_size):
			outp = outp + [{
				"acs":acs_data[i],
				"acs_header":acs_header,
				"times_header":times_header,
				"times":times_data[i],
				"metadata":metadata[i],
				"metadata_header":metadata_header
			}]
		
		return outp 


def main():
	FILES = ["/home/ilya/workspace/dnao/data/single_measurements.csv","/home/ilya/workspace/dnao/data/single_measurements2.csv"]

	#Put files you want to load here:
	FILES = ["/home/ilya/workspace/dnao/data/single_measurements.csv","/home/ilya/workspace/dnao/data/single_measurements2.csv"]

	#make separate jobs, specify the temperature:
	jobs = [
	{"file":FILES[0],"start":0,"stop":100},
	{"file":FILES[1],"start":0,"stop":200}]

	from nplab import datafile as df 
	datafile = df.DataFile("/home/ilya/workspace/dnao/data/extracted_file.h5","w")
	group = datafile.create_group("DLS_Zetasizer")
	for i,job in enumerate(jobs):
		print "JOB:", i
		dataset = get_correlations_from_file_by_header(job["file"],job["start"],job["stop"])
		
		for d in dataset:
			
			attrs  =  {}
			for (header,value) in zip(d["metadata_header"],d["metadata"]):
				attrs.update({header:value})

			attrs.update({"acs_header":d["acs_header"]})
			attrs.update({"times_header":d["times_header"]})
			group.create_dataset(name="run_%d",data=[d["times"],d["acs"]],attrs=attrs)

main()