import os
import fileinput
import csv
with open("/pine/scr/r/u/rujin/10XGenomics/breast_tissue_A_2k/output/breast_tissue_A_2k_per_cell_summary_metrics.csv") as f:
    reader = csv.reader(f)
    next(reader) # skip header
    data = []
    for r in reader:
        # get the barcode list from the 1st column
        bclist=r[0]
        data.append(bclist)
# total count of reads
tcount=0
bcount=0
for line in fileinput.input():
    tcount+=1
    line=line.strip()
    tags = line.split()[11:]
    tags_sort = sorted(tags)
    for tag in tags_sort:
	tag_split=tag.split(':')
        if 'CB' in tag_split:
	    barcode=tag_split[-1]
	    if barcode in data:
            	bcount+=1
            	print(str(tcount)+" "+str(bcount))
	    	directory = "../align/"+barcode[:-2]
            	outfile = directory+"/"+barcode+".sam"
            	try:
                    f=open(outfile,"a+")
		    print >> f,line
                    f.close()
            	except:
                    os.mkdir(directory)
                    f=open(outfile,"a+")
		    print >> f,line
                    f.close()
	    	break
pct=float(bcount)/float(tcount)
print(pct)

