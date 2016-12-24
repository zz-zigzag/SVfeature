## Run and Usage
* make
* usage: `./SVfeature -D -b test.list` or `./SVfeature -D -e 0.5 -m 1000 -b test.list`

## Note
* Need bai for each bam
* Output finename is `bam_filename_normalized` and `bam_filename_absolute`. The former is used for training model.
* Output format: feature1, feature2.......
* Need samtools in `$PATH`
* Need [SeqAn](www.seqan.de) included