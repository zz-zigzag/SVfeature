# SVfeature
Program: SVfeature (collect SV features base on bam file)
Version: 1.3.5
Contact: zz_zigzag <zz_zigzag@outlook.com>

Usage:   SVfeature [-i/D/I] [options]

Function mode
         -i        call features of indel(< 50bp)
         -D        call features of deletion(> 50bp)
         -I        call features of insertion(> 50bp)
Required parameters:
         -b FILE   the bam config file
                   Per line: bam_file, read_length, mean 
                   and std_dev of insert size, variant_file
                   Per line of variant_file:
                   if deletion: chr, start, end, isExist(0/1)
                   if insertion: chr, start, length, isExist(0/1)
                   if indel: chr, type(D/I), start, end/length, isExist(0/1)
                   For example: /data/test.bam 500 50 /data/test.v
         -o FILE   Output prefix
Optional parameters:
         -d        Output the detail feature

1. need bai to each bam
2. output format
	 label,feature1,feature2.......
3. need samtools in $PATH
4. compiling used SeqAn. SeqAn is an open source C++ library...	(www.seqan.de) see Makefile.
