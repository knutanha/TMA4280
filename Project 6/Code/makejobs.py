from string import Template
import os

textfile = open('tempjob_nojoin', 'r')
template = textfile.read()

# data = [MPI, OMP, SIZE, MIN, SEK, LNODES, PPN, name]

data = [
[6,6,16384,15,0,3,12,'foo'],
[6,6,16384,15,0,3,12,'bar'],
]

workdirectory = "/home/emilys/test"

for i in xrange(0,len(data)):
    name = "%i_%i_%i_%s" % (data[i][2], data[i][0], data[i][1],data[i][7])
    mpipn = "%i" % (data[i][0]/data[i][5])
    omp = "%i" % (data[i][1])
    size = "%i" % (data[i][2])
    minutes = "%i" % (data[i][3])
    secs = "%i" % (data[i][4])
    lnodes = "%i" % (data[i][5])
    ppn = "%i" % (data[i][6])
    test = Template(template)
    a = test.safe_substitute(NAME=name,MIN=minutes, SEK=secs, LNODES=lnodes, PPN=ppn, OMP=omp,MPIPN=mpipn,SIZE=size,WORKDIR=workdirectory)
    name = name + '.sh'
    newfile = open(name,'w')
    newfile.write(a)
    newfile.close()
    os.system("qsub " + name)
    
os.system("qstat")

