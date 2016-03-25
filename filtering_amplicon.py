import pysam

class FILTER_AMP:

	def bed_pro(self,file_name):
		f = open(file_name,"r")
		result_arr = []
		for x in f.readlines():
			temp = x.split("\t")
			if len(temp)<=4:
				result_arr.append([temp[0],temp[1],temp[2],temp[3].strip()])
			elif len(temp) > 4:
				result_arr.append([temp[0],temp[1],temp[2],temp[3].strip(),temp[4].strip()])
		f.close()

		return result_arr

	def string_count(self,dat):
		count = 0
		dele = ''
		del_count = 0
		for a in dat:
			if del_count ==0:
				if a.isalpha()==True:
					count += 1
				elif a=='^':
					del_count = 1
					count += 1
			elif del_count ==1:
				if a.isalpha()==True:
					count = count
				elif a.isdigit()==True:
					del_count = 0

		return count

	def distance_measure(self,read, bed_arr):
		result_arr = []
		#if read.is_unmapped==False:
		strand = read.is_reverse
		if strand==False:
			strand = 'F'
		else:
			strand = 'R'
		ref_st = read.reference_start
		ref_end = read.reference_end
		tag = read.get_tag('MD')

		for y in bed_arr:
			if strand == 'F':
				result_arr.append(ref_st-int(y[1]))

			else:
				result_arr.append(int(y[2])-ref_end)

		result_edit_arr = [None]*len(result_arr)
		for y in range(len(result_arr)):
			result_edit_arr[y] = abs(int(result_arr[y]))

		min_value = min(result_edit_arr)
		min_value_index = 0

		if result_edit_arr.count(min_value) == 2:####### count2
			#min_value_index = result_edit_arr.index(-int(min_value))
			min_value_index = result_edit_arr.index(int(min_value))
		else:
			min_value_index = result_edit_arr.index(int(min_value))



		if result_arr[min_value_index] > 0:
			if self.string_count(tag) > 10:
				result_edit_arr[min_value_index] = 100000
				min_value = min(result_edit_arr)
				min_value_index = result_edit_arr.index(min_value)

		result = bed_arr[min_value_index]

		if result[4]=='T':
			return 'T'
		else:
			return 'B'

	def write_proc(self):

		#bed = [['chr7',55242328,55242488],['chr7',140453119,140453273],['chr7',55259386,55259560]]
		read_name_arr = []

		bed_tarr=[]
		bed_earr = []

		for x in self.bed:
			if x[4]=='T':
				bed_tarr.append([x[0],x[1],x[2],x[3]])


		for x in bed_tarr:
			temp = []
			for bc in self.bed:
				if x[3]==bc[3]:
					temp.append(bc)

			for read in self.bamfile.fetch(x[0],int(x[1]),int(x[2])):

				if read.is_unmapped==False:
					result = self.distance_measure(read, temp)
					if result=='T':
						read_name_arr.append(read.query_name)


		for read in self.bamfile.fetch():
			temp = read.cigar
			count = 0

			if read.is_unmapped==False:
				for x in temp:
					if int(x[0])==0 or int(x[0])==1 or int(x[0])==2:
						count += int(x[1])
						#print str(count)
				if count <= 40:
					read_name_arr.append(read.query_name)

		for read in self.bamfile.fetch():
			if read.query_name not in read_name_arr:
				self.rw_bamfile.write(read)



	def __init__(self,fn, bed_fn):
		file_bam = fn.replace(".bam",".aedit.bam")
		self.bamfile = pysam.AlignmentFile(fn,"rb")
		self.rw_bamfile = pysam.AlignmentFile(file_bam,"wb", template=self.bamfile)
		self.bed = self.bed_pro(bed_fn)




	def __del__(self):
		self.bamfile.close()
		self.rw_bamfile.close()

