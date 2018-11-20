raw2clean
======================
    quality countrol for raw fastq files without adapter
    
## 安装
        pip install raw2clean
    Or
        git clone https://github.com/yodeng/raw2clean.git
        cd raw2clean
        python setup.py install
    

## 功能：
    对已除去接头的原始fastq序列数据进行质控，支持单端测序和双端测序数据。

## 具体质控标准如下：

* 1. reads基于overlap进行分析，检测只有街头序列的空载reads并去除
* 2. 截去reads两端的N碱基 （可选）
* 3. 截去reads两端的低质量碱基，低质量参数可设置 （可选）
* 4. 按指定滑动窗口对reads两端进行质控，窗口平均质量小于给定参数时，截去两端低质量序列（可选）
* 5. 去除reads中N碱基占5%的reads （默认5%）
* 6. 去除reads中低质量碱基大于30%的reads（默认30%）
* 7. 去除长度小于给定长度的reads （默认50）
* 8. 去除非配对的reads


sumfq
=====
    summary fastq

## 功能：
    统计fastq文件中每个位置的不同碱基的数目以及碱基错误率，支持SE和PE
    
## output:
* *.clean.fq.sum_avgQ.txt
* *.clean.fq.sum_pos.txt
