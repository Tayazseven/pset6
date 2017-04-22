import pysam
import gffutils as gfu

bam_in = pysam.AlignmentFile('SRR1175.sorted.bam') #opening the alignment file


datbase = gfu.create_db('saccharomyces_cerevisiae.gff', dbfn='sacc.db')# since I am going to run this on server
#I need to create database again
gen_name=[]
gen_fpkm = [] #dictionary to keep gene id's and fpkms
reads = 0
for gene in datbase.features_of_type('mRNA'):
    g_name = gene.attributes['Name'][0]
    gen_name.append(gene.attributes['Name'][0])
    try:
        fetch_gene = bam_in.count(reference=gene.chrom, start=gene.start, end=gene.stop)
    except ValueError:
        continue
    gen_fpkm.append(len([x for x in fetch_gene])) #putting the values into the dictionary
    reads += 1
fpkm =[]

for i in range(len(gen_fpkm)):
    fpkm.append(gen_fpkm[i]/(len(gen_fpkm) * reads) *1e9)

for x in fpkm:
    print ('%s %d' % (str(gen_name), (fpkm[x])))