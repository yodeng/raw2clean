#!/usr/bin/env python
#coding:utf-8
from collections import OrderedDict,defaultdict,deque
import argparse,multiprocessing,gzip,bz2file,os,gc,datetime,sys
from os.path import abspath,join
from itertools import count

def get_gen(x,n):
    return (x.next() for _ in xrange(n))

def sliding_window(array,q,k):
    if q <= 0 or k <= 0:
        return array
    if k >= len(array):
        if average(array) < q:
            return []
        else:
            return array
    queue = deque()
    for n,i in enumerate(array):
        if len(queue) < k:
            queue.append(i)
        else:        
            if average(queue) < q:
                return array[:n-1]
            queue.popleft()
            queue.append(i)
    if average(queue) < q:
        return array[:n-1]
    return array

def read_fq_list(filename):
    if filename.endswith(".gz"):
        h = gzip.open(filename)
    elif filename.endswith(".bz2"):
        h = bz2file.open(filename)
    else:
        h = open(filename)
    seq = []
    allseq = []
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
    
def get_low_index(arr,n):
    for i,arr_n in enumerate(arr):
        if arr_n >= n:
            return i
    return len(arr)
    
def get_low_index_r(arr,n):
    for i,arr_n in enumerate(arr[::-1]):
        if arr_n >= n:
            return len(arr) - i - 1
    return -1
                              
def qc_seq(oneseq):
    if not oneseq: return
    oneseq[-1] = [ord(x)-phread for x in oneseq[-1]]
    if trimN:
        length = len(oneseq[-1])
        oneseq[1] = oneseq[1].lstrip("N")
        tl = len(oneseq[1])
        oneseq[1] = oneseq[1].rstrip("N")
        tr = len(oneseq[1])
        oneseq[-1] = oneseq[-1][length-tl:length-tl+tr]
    if trimL is not None:
        il = get_low_index(oneseq[-1],trimL)
        oneseq[1] = oneseq[1][il:]
        oneseq[-1] = oneseq[-1][il:]
    if trimR is not None:
        ir = get_low_index_r(oneseq[-1],trimR)
        oneseq[1] = oneseq[1][:ir+1]
        oneseq[-1] = oneseq[-1][:ir+1]
    if len(oneseq[1]) <= minlen:return
    oneseq[-1] = sliding_window(oneseq[-1],windowavg,windowsize);oneseq[1] = oneseq[1][:len(oneseq[-1])]
    if len(oneseq[1]) <= minlen:return
    if oneseq[1].count("N")/float(len(oneseq[1])) > N_percent:return
    if len(filter(lambda x:x<lowq,oneseq[-1]))/float(len(oneseq[-1])) > low_percent:return
    oneseq[-1] = "".join([chr(i+phread) for i in oneseq[-1]])
    return oneseq
    
def qc_seq_pe(inseq_para):
    if not inseq_para:
        return
    r1,r2 = inseq_para
    result_pe = [1,2]
    if r1[0].split()[0] != r2[0].split()[0]:
        sys.stdout.write("reads not paired ( %s, %s ), please check in you sequence file\n"%(r1[0].lstrip(">"),r2[0].lstrip(">")))
        return [None,None]
    r1[-1] = [ord(x)-phread for x in r1[-1]]
    r2[-1] = [ord(x)-phread for x in r2[-1]]
    if trimN:
        length = len(r1[-1])
        r1[1] = r1[1].lstrip("N")
        tl = len(r1[1])
        r1[1] = r1[1].rstrip("N")
        tr = len(r1[1])
        r1[-1] = r1[-1][length-tl:length-tl+tr]
        length = len(r2[-1])
        r2[1] = r2[1].lstrip("N")
        tl = len(r1[1])
        r2[1] = r2[1].rstrip("N")
        tr = len(r2[1])
        r2[-1] = r2[-1][length-tl:length-tl+tr]        
    if trimL > 0:
        il = get_low_index(r1[-1],trimL)
        r1[1] = r1[1][il:]
        r1[-1] = r1[-1][il:]
        il = get_low_index(r2[-1],trimL)
        r2[1] = r2[1][il:]
        r2[-1] = r2[-1][il:]        
    if trimR > 0:
        ir = get_low_index_r(r1[-1],trimR)
        r1[1] = r1[1][:ir+1]
        r1[-1] = r1[-1][:ir+1]
        ir = get_low_index_r(r1[-1],trimR)
        r2[1] = r2[1][:ir+1]
        r2[-1] = r2[-1][:ir+1]              
    if len(r1[1]) <= minlen:result_pe[0] = None
    if len(r2[1]) <= minlen:result_pe[1] = None      
    if result_pe[0] is not None: r1[-1] = sliding_window(r1[-1],windowavg,windowsize);r1[1] = r1[1][:len(r1[-1])]
    if result_pe[1] is not None: r2[-1] = sliding_window(r2[-1],windowavg,windowsize);r2[1] = r2[1][:len(r2[-1])]
    if len(r1[1]) <= minlen:result_pe[0] = None
    if len(r2[1]) <= minlen:result_pe[1] = None 
    if result_pe[0] is not None:
        if r1[1].count("N")/float(len(r1[1])) > N_percent:result_pe[0] = None
    if result_pe[1] is not None:
        if r2[1].count("N")/float(len(r2[1])) > N_percent:result_pe[1] = None
    if result_pe[0] is not None:
        if len(filter(lambda x:x<lowq,r1[-1]))/float(len(r1[-1])) > low_percent:result_pe[0] = None
    if result_pe[1] is not None:
        if len(filter(lambda x:x<lowq,r2[-1]))/float(len(r2[-1])) > low_percent:result_pe[1] = None
    r1[-1] = "".join([chr(i+phread) for i in r1[-1]])
    r2[-1] = "".join([chr(i+phread) for i in r2[-1]])
    if result_pe[0] is not None:result_pe[0] = r1
    if result_pe[1] is not None:result_pe[-1] = r2
    return result_pe

def average(seq): 
	return float(sum(seq)) / len(seq)
    
def parseArg():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,description= '''This Script is used to do quality control for rawdata fastq sequence.

Version
    Author: Yong Deng,yongdeng@capitalbiotech.com
    Company: Capitalbiotech
    Version: 1.0
    Created: 04/30/2018 10:23:53 AM
''')
    group = parser.add_mutually_exclusive_group()
    group.add_argument("-pe", "--PE",type=str ,help="you input paired-end sequence fastq files, first will be considered as R1 and second will be considered as R2, *.gz or .bz2 can be allowed",nargs=2,metavar = 'reads')
    group.add_argument("-se", "--SE", type=str,help="you input sigle-end sequence fastq file, *.gz or .bz2 can be allowed",metavar="se-reads")
    parser.add_argument("-l","--length",type=int,help="The minmum length of one read to keep, default:50",default=50,metavar = "minlength")
    parser.add_argument("-lq","--left_low_q",type=int,help="If passed, trim the base has qual score less than the given value in left ends of the read. 0 by default means do not trim",default=0,metavar = "q-score")
    parser.add_argument("-rq","--right_low_q",type=int,help="If passed, trim the base has qual score less than the given value in left ends of the read. 0 by default means do not trim",default=0,metavar = "q-score")
    parser.add_argument("-stripN",action="store_true",help="If passed, trim the N base in both ends of the read. Off by default",default=False)
    parser.add_argument("-N","--cutNper",type=float,help="cutoff the reads that N content, default:0.05",default=0.05,metavar="N_percent")
    parser.add_argument("-q","--QualScore",type=int,help="The cut-off value for PHRED quality score for high-quality filtering, default:20",default=20,metavar="low_q_score")
    parser.add_argument("-qc","--Qualcontent",type=float,help="The cut-off value for percentage of read length that should be of given quality, default:0.3",default=0.3,metavar = "bad_q_percent")
    parser.add_argument("-w","--window",type=int,help="if passed, two values must be given, first represent the slide window size, second is the cut-off average quality score. default: 0 0, means do not trim by windows",nargs=2)
    parser.add_argument("-c","--cpus",type=int,help="Number of CPUs to be used, default:8",default=8,metavar="ncpus")
    parser.add_argument("-phread","--phread",type=int,help="PHRED quality type, 33 or 64, default:33",default=33)
    parser.add_argument("-o","--out_dir",type=str,help="the output directory",required=True)
    return parser.parse_args()
    
def main():
  global args
  args=parseArg()
  global minlen,trimN,trimL,trimR,lowq,phread,N_percent,low_percent,windowsize,windowavg
  minlen,trimN,trimL,trimR,lowq,phread,N_percent,low_percent = args.length,args.stripN,args.left_low_q,args.right_low_q,args.QualScore,args.phread,args.cutNper,args.Qualcontent
  windowsize,windowavg = args.window if args.window else (0,0)
  if windowsize <0 or windowavg < 0:
    print "Error: only positive values be allowed of --window"
    sys.exit(1)
  if args.SE:
    seqfile = abspath(args.SE)
    if not os.path.isfile(seqfile):
        sys.stderr.write("Error: %s file exists, exit!\n")
        sys.exit(1)
    with open(join(args.out_dir,os.path.splitext(os.path.basename(seqfile))[0]+".clean.fastq"),"w") as fo:pass
    seqlist = read_fq_list(seqfile)
    sys.stdout.write("----------(%s) Start: QC for raw fastq file, %d cpus ----------\n"%(datetime.datetime.now().strftime("%x %X"),min(args.cpus,multiprocessing.cpu_count()-2)))
    pool = multiprocessing.Pool(min(args.cpus,multiprocessing.cpu_count()-2))
    for seqcycle in count(0):
        per_inseq = get_gen(seqlist,2000000)
        per_inseq_list = list(per_inseq)
        if len(per_inseq_list) == 0:break
        sys.stdout.write("%d\tProcessing %d sequences ... (%s)\n"%(seqcycle+1,len(per_inseq_list),datetime.datetime.now().strftime("%x %X")))
        cleanseqs = pool.imap(qc_seq,per_inseq_list)
        with open(join(args.out_dir,os.path.splitext(os.path.basename(seqfile))[0]+".clean.fastq"),"a+") as fo:
            for r in cleanseqs:
                if r:
                    fo.writelines([i + "\n" for i in r]) 
    sys.stdout.write("----------(%s) End. ----------\n"%datetime.datetime.now().strftime("%x %X"))
    pool.close();pool.join()
  elif args.PE:
    seqfile1,seqfile2 = abspath(args.PE[0]),abspath(args.PE[1])
    if not os.path.isfile(seqfile1):
        sys.stderr.write("Error: %s file not exists, exit!\n"%seqfile1)
        sys.exit(1)
    if not os.path.isfile(seqfile2):
        sys.stderr.write("Error: %s file not exists, exit!\n"%seqfile2)
        sys.exit(1)        
    if not os.path.isdir(abspath(args.out_dir)): os.makedirs(abspath(args.out_dir))
    with open(join(args.out_dir,"R1.clean.fq"),"w") as fp1,open(join(args.out_dir,"R2.clean.fq"),"w") as fp2,open(join(args.out_dir,"R1.clean.unpaired.fq"),"w") as f1,open(join(args.out_dir,"R2.clean.unpaired.fq"),"w") as f2:
        pass
    r1seqlist = read_fq_list(seqfile1)
    r2seqlist = read_fq_list(seqfile2)
    sys.stdout.write("----------(%s) Start: QC for raw fastq file, %d cpus ----------\n"%(datetime.datetime.now().strftime("%x %X"),min(args.cpus,multiprocessing.cpu_count()-2)))
    pool = multiprocessing.Pool(min(args.cpus,multiprocessing.cpu_count()-2))
    inseq = ([r1,r2seqlist.next()] for r1 in r1seqlist)
    for seqcycle in count(0):
        per_inseq = get_gen(inseq,2000000)
        per_inseq_list = list(per_inseq)
        if len(per_inseq_list) == 0:
            break 
        sys.stdout.write("%d\tProcessing %d sequences ... (%s)\n"%(seqcycle+1,len(per_inseq_list),datetime.datetime.now().strftime("%x %X")))
        cleanseqs = pool.imap(qc_seq_pe,per_inseq_list)
        with open(join(args.out_dir,"R1.clean.fq"),"a+") as fp1,open(join(args.out_dir,"R2.clean.fq"),"a+") as fp2,open(join(args.out_dir,"R1.clean.unpaired.fq"),"a+") as f1,open(join(args.out_dir,"R2.clean.unpaired.fq"),"a+") as f2:
            for r12 in cleanseqs:
                if r12[0] and r12[1]:
                    fp1.writelines([i+"\n" for i in r12[0]])
                    fp2.writelines([i+"\n" for i in r12[1]])
                else:
                    if r12[0]:
                        f1.writelines([i+"\n" for i in r12[0]])
                    if r12[1]:
                        f2.writelines([i+"\n" for i in r12[1]])
    sys.stdout.write("----------(%s) End. ----------\n"%datetime.datetime.now().strftime("%x %X"))
    pool.close();pool.join()
  else:
    print "Error: please give you fastq file by --SE or --PE"
    sys.exit(1)

if __name__ == "__main__":
    #Making sure you are running a version of python that works with this script.
    if sys.version_info[0] != 2 or sys.version_info[1] < 7 or sys.version_info[2] < 5:
        print("This script requires Python version 2.7.5 or higher within major version 2")
        sys.exit(1)
    main()
                    
