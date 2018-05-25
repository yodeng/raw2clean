# raw2clean.py

quality countrol for raw fastq files


功能：对已除去接头的原始fastq序列数据进行质控，支持单端测序和双端测序数据。

具体质控标准如下：

	1. reads基于overlap进行分析，检测adapter并去除
    
	2. 截去reads两端的N碱基 （可选）
    
	3. 截去reads两端的低质量碱基，低质量参数可设置 （可选）
    
	4. 按指定滑动窗口对reads两端进行质控，窗口平均质量小于给定参数时，截去两端低质量序列（可选）
    
	5. 去除reads中N碱基占5%的reads （默认5%）
    
	6. 去除reads中低质量碱基大于30%的reads（默认30%）
    
	7. 去除长度小于给定长度的reads （默认50）
    
	8. 去除非配对的reads
    
$ python raw2clean.py  -h
usage: raw2clean.py [-h] [-pe reads reads | -se se-reads] [-l minlength]
                    [-lq q-score] [-rq q-score] [-stripN] [-N N_percent]
                    [-q low_q_score] [-qc bad_q_percent] [-w WINDOW WINDOW]
                    [-c ncpus] [-phread PHREAD] -o OUT_DIR

This Script is used to do quality control for rawdata fastq sequence.

Version
    Author: Yong Deng,yodeng@tju.edu.cn
    Version: 1.0
    Created: 04/30/2018 10:23:53 AM

optional arguments:
  -h, --help            show this help message and exit
  -pe reads reads, --PE reads reads
                        you input paired-end sequence fastq files, first will
                        be considered as R1 and second will be considered as
                        R2, *.gz or .bz2 can be allowed
  -se se-reads, --SE se-reads
                        you input sigle-end sequence fastq file, *.gz or .bz2
                        can be allowed
  -l minlength, --length minlength
                        The minmum length of one read to keep, default:50			
  -lq q-score, --left_low_q q-score
                        If passed, trim the base has qual score less than the
                        given value in left ends of the read. 0 by default
                        means do not trim
  -rq q-score, --right_low_q q-score
                        If passed, trim the base has qual score less than the
                        given value in left ends of the read. 0 by default
                        means do not trim
  -stripN               If passed, trim the N base in both ends of the read.
                        Off by default
  -N N_percent, --cutNper N_percent
                        cutoff the reads that N content, default:0.05
  -q low_q_score, --QualScore low_q_score
                        The cut-off value for PHRED quality score for high-
                        quality filtering, default:20
  -qc bad_q_percent, --Qualcontent bad_q_percent
                        The cut-off value for percentage of read length that
                        should be of given quality, default:0.3
  -w WINDOW WINDOW, --window WINDOW WINDOW
                        if passed, two values must be given, first represent
                        the slide window size, second is the cut-off average
                        quality score. default: 0 0, means do not trim by
                        windows
  -c ncpus, --cpus ncpus
                        Number of CPUs to be used, default:8
  -phread PHREAD, --phread PHREAD
                        PHRED quality type, 33 or 64, default:33
  -o OUT_DIR, --out_dir OUT_DIR
                        the output directory
    



# summary_fastq.py

get a summary along the position of you fastq sequence file


功能：对fastq序列文件进行统计，统计每个位置的碱基和质量分布情况。

$ python summary_fastq.py -h

usage: summary_fastq.py [-h] [-pe PE PE | -se SE] [-phread PHREAD] -o OUT_DIR

This Script is used to summarize the fastq seq file, not less than 4 cpus
avaliable on you work machine

optional arguments:
  -h, --help            show this help message and exit
  -pe PE PE, --PE PE PE
                        you input paired-end sequence fastq files, first will
                        be considered as R1 and second will be considered as
                        R2, *.gz or .bz2 can be allowed
  -se SE, --SE SE       you input sigle-end sequence fastq file, *.gz or .bz2
                        can be allowed
  -phread PHREAD, --phread PHREAD
                        PHRED quality type, 33 or 64, default:33
  -o OUT_DIR, --out_dir OUT_DIR
                        the output directory
