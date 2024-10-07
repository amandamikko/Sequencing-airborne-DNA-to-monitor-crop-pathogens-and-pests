import os 
import sys
filename1 = sys.argv[1]
filename2 = sys.argv[2]
file_name_start = filename1.split("/")[-1][:3]
out1 = sys.argv[3][:-3]
out2 = sys.argv[4][:-3]
#Unzip files using sys
gzip_d = "gzip -d {f}".format(f=filename1)
os.system(gzip_d)
gzip_d = "gzip -d {f}".format(f=filename2)
os.system(gzip_d)
without_gzip1 = filename1[:-3]
without_gzip2 = filename2[:-3]

#Update read id's so that there are unique ids for forward and revere reads
#Adds /1 to forward reads, and /2 to reverse reads 
nr_contig = -1
nr_reads = -1
current_contig = ""
#Go through all readids in the original file and write a modified readid to the file specified by out1
#the moddified read is unique for each read and makes it possible to extract reads after for example classification 
with open(out1, "w") as out_f:
    with open(without_gzip1, "r") as wgzip:
        for line in wgzip:
            if line[0:4] == "@" + file_name_start:
                contig = line.split()[0].split("-")[-1]
                if contig != current_contig:
                    nr_contig += 1
                    nr_reads = 0
                    current_contig = contig
                else:
                    nr_reads += 1
                id = "_{n}_{c}/1".format(n=nr_reads, c=nr_contig)

                out_line = contig + id + " \n"
                out_f.write("@" + out_line)
            else:
                out_f.write(line)

nr_contig = -1
nr_reads = -1
current_contig = ""
#Go through all readids in the original file and write a modified readid to the file specified by out2
with open(out2, "w") as out_f:
    with open(without_gzip2, "r") as wgzip:
        for line in wgzip:
            if line[0:4] == "@" + file_name_start:
                contig = line.split()[0].split("-")[-1]
                if contig != current_contig:
                    nr_contig += 1
                    nr_reads = 0
                    current_contig = contig
                else:
                    nr_reads += 1
                id = "_{n}_{c}/2".format(n=nr_reads, c=nr_contig)

                out_line = contig + id + " \n"
                out_f.write("@" + out_line)
            else:
                out_f.write(line)
                
#Zip original files and created files
gzip = "gzip {f}".format(f=without_gzip1)
os.system(gzip)
gzip = "gzip {f}".format(f=without_gzip2)
os.system(gzip)
gzip = "gzip {f}".format(f=out1)
os.system(gzip)
gzip = "gzip {f}".format(f=out2)
os.system(gzip)