# Author: Jacqueline Robinson
# Adapted by: Joseph Curti (jcurti3@g.ucla.edu)

# Script to count number of called genotypes and number of heterozygotes per sample in 
# sliding windows.
# Input file is a single- or multi-sample VCF file that has been filtered (passing sites 
# have "PASS" in the FILTER column) and compressed with gzip/bgzip.
#
# Usage: 
# python ./SlidingWindowHet.py [vcf] [window size] [step size] [chromosome/scaffold name]
#
# Windows will be non-overlapping if step size == window size.
#
# Example: 
# python ./SlidingWindowHet.py input.vcf.gz 100000 10000 chr01

import sys
import pysam
import os
import gzip


# Open input file and make sure the VCF file is indexed (if not, create index)
filename = sys.argv[1]
VCF = gzip.open(filename, 'r')

if not os.path.exists("%s.tbi" % filename):
    pysam.tabix_index(filename, preset="vcf")
parsevcf = pysam.Tabixfile(filename)


# Set variables
window_size = int(sys.argv[2])
step_size = int(sys.argv[3])
scaff = sys.argv[4]
scaff_size = {'JALIRH010000001.1':121926015,'JALIRH010000002.1':80829437,'JALIRH010000003.1':68643687,'JALIRH010000004.1':53416606,'JALIRH010000005.1':43834372,'JALIRH010000006.1':42153796,'JALIRH010000007.1':30073597,'JALIRH010000008.1':22756518,'JALIRH010000010.1':20566091,'JALIRH010000011.1':19876581,'JALIRH010000012.1':19408658,'JALIRH010000013.1':18906995,'JALIRH010000014.1':18310531,'JALIRH010000015.1':15888862,'JALIRH010000016.1':15455407,'JALIRH010000017.1':14490932,'JALIRH010000018.1':14030056,'JALIRH010000019.1':13691472,'JALIRH010000020.1':13641242,'JALIRH010000021.1':12744333,'JALIRH010000022.1':12664268,'JALIRH010000024.1':11909977,'JALIRH010000026.1':11083105,'JALIRH010000027.1':10594039,'JALIRH010000028.1':10313052,'JALIRH010000029.1':10126570,'JALIRH010000030.1':9342254,'JALIRH010000032.1':6589112,'JALIRH010000033.1':6552069,'JALIRH010000034.1':6413858,'JALIRH010000035.1':6246747,'JALIRH010000036.1':6170350,'JALIRH010000037.1':5953211,'JALIRH010000038.1':5893186,'JALIRH010000039.1':5486612,'JALIRH010000040.1':5478314,'JALIRH010000042.1':4869296,'JALIRH010000043.1':4864553,'JALIRH010000044.1':4848042,'JALIRH010000045.1':4406244,'JALIRH010000046.1':4372448,'JALIRH010000047.1':4200305,'JALIRH010000048.1':4089568,'JALIRH010000050.1':3968916,'JALIRH010000051.1':3874674,'JALIRH010000053.1':3618797,'JALIRH010000054.1':3601799,'JALIRH010000055.1':3477356,'JALIRH010000056.1':3430521,'JALIRH010000057.1':3317826,'JALIRH010000058.1':3254601,'JALIRH010000059.1':3225810,'JALIRH010000060.1':3225353,'JALIRH010000061.1':3025772,'JALIRH010000063.1':3006878,'JALIRH010000064.1':3006518,'JALIRH010000065.1':2968981,'JALIRH010000067.1':2566882,'JALIRH010000068.1':2524086,'JALIRH010000069.1':2517896,'JALIRH010000070.1':2491622,'JALIRH010000071.1':2472232,'JALIRH010000072.1':2470316,'JALIRH010000073.1':2445880,'JALIRH010000074.1':2266200,'JALIRH010000075.1':2237334,'JALIRH010000076.1':2126800,'JALIRH010000077.1':2121156,'JALIRH010000078.1':2046405,'JALIRH010000079.1':2033782,'JALIRH010000080.1':1941001,'JALIRH010000081.1':1806888,'JALIRH010000082.1':1806726,'JALIRH010000083.1':1777743,'JALIRH010000084.1':1776689,'JALIRH010000086.1':1716458,'JALIRH010000088.1':1648492,'JALIRH010000089.1':1645390,'JALIRH010000090.1':1580646,'JALIRH010000091.1':1538330,'JALIRH010000092.1':1530834,'JALIRH010000093.1':1519857,'JALIRH010000097.1':1465439,'JALIRH010000101.1':1397041,'JALIRH010000102.1':1393513,'JALIRH010000104.1':1309764,'JALIRH010000106.1':1186712,'JALIRH010000108.1':1122020,'JALIRH010000109.1':1120247,'JALIRH010000110.1':1119225,'JALIRH010000111.1':1118613,'JALIRH010000112.1':1107949,'JALIRH010000113.1':1099624,'JALIRH010000114.1':1099513,'JALIRH010000115.1':1045722,'JALIRH010000116.1':1038886,'JALIRH010000117.1':1023579,'JALIRH010000118.1':1022764,'JALIRH010000119.1':1004856} 


# Get list of samples from VCF file header
samples=[]
for line in VCF:
    if line.startswith('##'):
        pass
    else:
        for i in line.split()[9:]: samples.append(i)
        break


# Get start and end positions of scaffold
for line in VCF:
    if line[0] != '#':
        start_pos = int(line.strip().split()[1])
        end_pos = int(scaff_size[scaff])
        break


# Create output file
output = open(filename + '_het_%swin_%sstep.txt' % (window_size, step_size), 'w')
output.write('scaff\twindow_start\tsites_total\tcalls_%s\thets_%s\n' % ('\tcalls_'.join(samples), '\thets_'.join(samples)) )


# Fetch a region, ignore sites that fail filters, tally genotype calls and heterozygotes        
def snp_cal(scaff,window_start,window_end):
    print("%s:%s" % (scaff,window_start))
    rows = tuple(parsevcf.fetch(region="%s:%s-%s" % (scaff, window_start, window_end), parser=pysam.asTuple()))    
    sites_total=0
    calls=[0]*len(samples)
    hets=[0]*len(samples)
    for line in rows:
        if line[6]!="PASS": continue
        sites_total+=1
        for i in range(0,len(samples)):
            if line[i+9][:1]=='.': continue
            calls[i]+=1
            GT=line[i+9].split(':')[0]
            if '/' in GT: sp='/'
            if '|' in GT: sp='|'
            if GT.split(sp)[0]!=GT.split(sp)[1]: hets[i]+=1
    output.write('%s\t%s\t%s\t%s\t%s\n' % (scaff,window_start,sites_total,'\t'.join(map(str,calls)),'\t'.join(map(str,hets))) )


# Initialize window start and end coordinates
window_start = start_pos
window_end = start_pos+window_size-1


# Calculate stats for window, update window start and end positions, 
# repeat to end of scaffold
while window_end <= end_pos:    
    if window_end < end_pos:
        snp_cal(scaff,window_start,window_end)
        window_start = window_start + step_size
        window_end = window_start + window_size - 1
    else:
        snp_cal(scaff,window_start,window_end)
        break    
else:
    window_end = end_pos
    snp_cal(scaff,window_start,window_end)


# Close files and exit
VCF.close()
output.close()

exit()
