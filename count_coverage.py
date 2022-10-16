import subprocess as sp

stepsize = 10000

fh = open("seqs.txt")
for x in fh:
    y = x.strip().split('\t')
    myname = y[0]
    mylength = int(y[1])
    for start in range(0,mylength,stepsize):
        proc = sp.Popen(['samtools','coverage','-r','{0}:{1}-{2}'.format(myname,start,start+stepsize),"Paspalum_genome_sorted.bam"],stdout=sp.PIPE)
        proc.wait()
        abc = proc.stdout.read().split("\n")[1].split("\t")
        print(",".join([myname,str(start),str(start+stepsize),abc[5],abc[6]]))
