#!/usr/bin/env python
#coding:utf-8
from collections import OrderedDict,defaultdict
import argparse,multiprocessing,gzip,bz2file,os,gc,datetime
from os.path import abspath,join
from itertools import count

def read_fq_list(filename):
    if filename.endswith(".gz"):
        h = gzip.open(filename)
    elif filename.endswith(".bz2"):
        h = bz2file.open(filename)
    else:
        h = open(filename)
    seq = []
    for n,line in enumerate(h):
        line = line.strip()
        if n % 4 ==0:
            if seq:
                yield seq
            seq = [line,]
        else:
            seq.append(line)
    h.close()
    yield seq

def sum_by_pos(seqs,outfile1,outfile2,phread=33):
    seqtype = "ATGCN"
    pos_info = {}
    qual_pos = defaultdict(int)
    total,reads,q20,q30,q40 = 0,0,0,0,0
    for num,oneseq in enumerate(seqs):
        if num % 2000000 == 0:
            print "---(%s) Process id %d read 2000000 sequences -----"%(datetime.datetime.now().strftime("%x %X"),os.getpid())
        reads += 1
        oneseq[-1] = [ord(x)-phread for x in oneseq[-1]]
        total += len(oneseq[-1])
        for index in xrange(len(oneseq[1])):
            if oneseq[-1][index] >= 20: q20+=1
            if oneseq[-1][index] >= 30: q30+=1
            if oneseq[-1][index] >= 40: q40+=1
            #qual_pos.setdefault(index,[]).append(oneseq[-1][index])
            qual_pos[index] += oneseq[-1][index]
            pos_info.setdefault(oneseq[1][index],defaultdict(int))[index] += 1
    for i in seqtype:
        pos_info.setdefault(i,defaultdict(int))
    print "---(%s) Process id %d read total %d sequences -----"%(datetime.datetime.now().strftime("%x %X"),os.getpid(),num+1)
    print "Start write %s and %s file ... (%s)"%(outfile1,outfile2,datetime.datetime.now().strftime("%x %X"))
    with open(outfile1,"w") as fo1,open(outfile2,"w") as fo2:
        #fo2.write("# seq_num (%d), total_base(%d), Q20(%.2f), Q30(%.2f), Q40(%.2f)\n"%(reads,total,q20/total*100,q30/total*100,q40/total*100))
        fo2.write("# seq_num (%d), total_base(%d), Q20(%d,%.2f), Q30(%d,%.2f), Q40(%d,%.2f)\n"%(reads,total,q20,float(q20)/total*100,q30,float(q30)/total*100,q40,float(q40)/total*100))
        fo2.write("#position\tavg_qual\terror_rate\n")
        for i in count(0):
            A_count = pos_info["A"].get(i,0)
            T_count = pos_info["T"].get(i,0)
            G_count = pos_info["G"].get(i,0)
            C_count = pos_info["C"].get(i,0)
            N_count = pos_info["N"].get(i,0)
            all_count = sum([A_count,T_count,G_count,C_count,N_count])
            if all_count == 0:break
            avgq = float(qual_pos[i])/all_count
            errate = pow(10,-(avgq)/10.0)*100
            fo1.write("\t".join([str(i),"A",str(A_count),"%.2f"%(float(A_count)/all_count*100),"T",str(T_count),"%.2f"%(float(T_count)/all_count*100),"G",str(G_count),"%.2f"%(float(G_count)/all_count*100),"C",str(C_count),"%.2f"%(float(C_count)/all_count*100),"N",str(N_count),"%.2f"%(float(N_count)/all_count*100)]) +"\n")            
            fo2.write("\t".join([str(i),"%.6f"%avgq,"%.6f"%errate]) + "\n")
    print "Write %s and %s files Done! (%s)"%(outfile1,outfile2,datetime.datetime.now().strftime("%x %X"))
    
def average(seq): 
	return float(sum(seq)) / len(seq)
    
def parseArg():
    parser = argparse.ArgumentParser(description="This Script is used to summarize the fastq seq file, not less than 4 cpus avaliable on you work machine")
    group = parser.add_mutually_exclusive_group()
    group.add_argument("-pe", "--PE",type=str ,help="you input paired-end sequence fastq files, first will be considered as R1 and second will be considered as R2, *.gz or .bz2 can be allowed",nargs=2)
    group.add_argument("-se", "--SE", type=str,help="you input sigle-end sequence fastq file, *.gz or .bz2 can be allowed")
    parser.add_argument("-phread","--phread",type=int,help="PHRED quality type, 33 or 64, default:33",default=33)
    parser.add_argument("-o","--out_dir",type=str,help="the output directory",required=True)
    return parser.parse_args()
    
def main():
    args = parseArg()
    if args.SE:
        seqfile = abspath(args.SE)
        if not os.path.isfile(seqfile):
            print "Error: %s file exists, exit!"
            sys.exit(1)
        print "Start time: %s"%datetime.datetime.now().strftime("%x %X")
        seqlist = read_fq_list(seqfile)
        sample = os.path.basename(os.path.dirname(seqfile))
        if not os.path.isdir(abspath(args.out_dir)): os.makedirs(abspath(args.out_dir))
        sum_by_pos(seqlist,abspath(join(args.out_dir,sample + ".clean.fq.sum_pos.txt")),abspath(join(args.out_dir,sample + ".clean.fq.sum_avgQ.txt")),args.phread)
    elif args.PE:
        seqfile1,seqfile2 = abspath(args.PE[0]),abspath(args.PE[1])
        sample = os.path.basename(os.path.dirname(seqfile1))
        if not os.path.isfile(seqfile1):
            print "Error: %s file not exists, exit!"%seqfile1
            sys.exit(1)
        if not os.path.isfile(seqfile2):
            print "Error: %s file not exists, exit!"%seqfile2
            sys.exit(1)         
        print "Start time: %s"%datetime.datetime.now().strftime("%x %X")
        if not os.path.isdir(abspath(args.out_dir)): os.makedirs(abspath(args.out_dir))
        r1seqlist = read_fq_list(seqfile1)
        r2seqlist = read_fq_list(seqfile2)
        p1 = multiprocessing.Process(target=sum_by_pos,args=(r1seqlist,abspath(join(args.out_dir,sample +".sum_pos.R1.txt")),abspath(join(args.out_dir,sample +".sum_avgQ.R1.txt")),args.phread));p1.start()
        p2 = multiprocessing.Process(target=sum_by_pos,args=(r2seqlist,abspath(join(args.out_dir,sample +".sum_pos.R2.txt")),abspath(join(args.out_dir,sample +".sum_avgQ.R2.txt")),args.phread));p2.start()
        p1.join();p2.join()
        if os.path.exists(abspath(join(args.out_dir,sample +".sum_pos.R1.txt"))) and os.path.exists(abspath(join(args.out_dir,sample +".sum_pos.R2.txt"))):
            posr1 = [i.split()[1:] for i in open(abspath(join(args.out_dir,sample +".sum_pos.R1.txt"))).readlines()]
            posr2 = [i.split()[1:] for i in open(abspath(join(args.out_dir,sample +".sum_pos.R2.txt"))).readlines()]
            with open(abspath(join(args.out_dir,sample + ".sum_pos.txt")),"w") as fo_pos:
                for line,lines in enumerate(posr1+posr2):
                    fo_pos.write(str(line) + "\t" + "\t".join(lines) + "\n")
        if os.path.exists(abspath(join(args.out_dir,sample +".sum_avgQ.R1.txt"))) and os.path.exists(abspath(join(args.out_dir,sample +".sum_avgQ.R2.txt"))):
            avgQr1 = [i.split()[1:] for i in open(abspath(join(args.out_dir,sample +".sum_avgQ.R1.txt"))).readlines()]
            avgQr2 = [i.split()[1:] for i in open(abspath(join(args.out_dir,sample +".sum_avgQ.R2.txt"))).readlines()]
            with open(abspath(join(args.out_dir,sample + ".sum_avgQ.txt")),"w") as fo_avgQ:
                for line,lines in enumerate(avgQr1[2:] + avgQr2[2:]):
                    fo_avgQ.write(str(line) + "\t" + "\t".join(lines) + "\n")

    else:
        print "Error: please give you fastq file by --SE or --PE"
        sys.exit(1)
      
if __name__ == "__main__":
    main()

