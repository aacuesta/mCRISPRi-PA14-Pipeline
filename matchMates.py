import sys, simplesam, csv

#Installed simplesam parser from here: https://github.com/mdshw5/simplesam
#Able to install using "pip install simplesam"

read1path = sys.argv[1]
read2path = sys.argv[2]
read1fasta = sys.argv[3]
outputpath = sys.argv[4]

counts = {}
readoneUnaligned = 0
readtwoUnaligned = 0
bothUnaligned = 0

with open(read1path, 'r') as samfile1:
	with open(read2path, 'r') as samfile2:
		sam1 = simplesam.Reader(samfile1)
		sam2 = simplesam.Reader(samfile2)
		r1headers = sam1.seqs
		r2headers = sam2.seqs
		for read in sam1:
			read1_id = read.qname
			aligned_name_1 = read.rname

			read2 = next(sam2)
			read2_id = read2.qname
			aligned_name_2 = read2.rname

			if read1_id == read2_id:
				if (aligned_name_1 != "*") and (aligned_name_2 != "*"):

					pairing = aligned_name_1 + "_" + aligned_name_2
					try:
						counts[pairing] += 1
					except:
						counts[pairing] = 1
				elif (aligned_name_1 == "*") and (aligned_name_2 == "*"):
					bothUnaligned +=1
				elif (aligned_name_1 == "*") and (aligned_name_2 != "*"):
					readoneUnaligned +=1
				elif (aligned_name_1 != "*") and (aligned_name_2 == "*"):
					readtwoUnaligned += 1
			else:
				print("Out of register!")
				print(read1_id)
				print(read2_id)

blocklength = 2
with open(read1fasta, 'r') as fasta1:
	line1 = []
	for each in fasta1:
		line1.append(each.strip())
		if len(line1) == 2:
			promgeneP1 = line1[0][1:] + "_P1"
			promgeneP2 = line1[0][1:] + "_P2"
			promgeneP3 = line1[0][1:] + "_P3"
			if promgeneP1 not in counts.keys():
				counts[promgeneP1] = 0
			if promgeneP2 not in counts.keys():
				counts[promgeneP2] = 0
			if promgeneP3 not in counts.keys():
				counts[promgeneP3] = 0
			line1 = []


with open(outputpath, 'w') as countsfile:
	csvwriter = csv.writer(countsfile)
	for target in counts.keys():
		csvwriter.writerow([target, counts[target]])

	#for header1 in r1headers:
	#	for header2 in r2headers:
	#		fullname = header1 + "_" + header2
	#		if fullname not in counts.keys():
	#			csvwriter.writerow([fullname, 0])
	csvwriter.writerow(["Only Read 1 unaligned", readoneUnaligned])
	csvwriter.writerow(["Only Read 2 unaligned", readtwoUnaligned])
	csvwriter.writerow(["Both Read 1 and 2 unaligned", bothUnaligned])
