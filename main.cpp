/*  The MIT License

    Copyright (C) 2015 zz_zigzag <zz_zigzag@outlook.com>

*/

#include <iostream>
#include <vector>
#include <zlib.h>
#include <unistd.h>
#include <ctime>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <seqan/bam_io.h>
#include "indel_list.h"
#include "bam.h"
#include "feature.h"
#include "file.h"

using namespace seqan;
using namespace std;

#define PACKAGE_VERSION "1.3.6"



static int usage(void )
{
    fprintf(stderr, "\n");
    fprintf(stderr, "Program: SVfeature (collect SV features base on bam file)\n");
    fprintf(stderr, "Version: %s\n", PACKAGE_VERSION);
    fprintf(stderr, "Contact: zz_zigzag <zz_zigzag@outlook.com>\n\n");
    fprintf(stderr, "Usage:   SVfeature [-i/D/I] [options]\n\n");
    fprintf(stderr, "Function mode\n");
    fprintf(stderr, "         -i        call features of indel(< 50bp)\n");
    fprintf(stderr, "         -D        call features of deletion(> 50bp)\n");
    fprintf(stderr, "         -I        call features of insertion(> 50bp)\n");
    fprintf(stderr, "Required parameters:\n");
    fprintf(stderr, "         -b FILE   the bam config file\n");
    fprintf(stderr, "                   Per line: bam_file, read_length, mean \n");
    fprintf(stderr, "                   and std_dev of insert size, variant_file\n");
    fprintf(stderr, "                   For example: /data/test.bam 100 500 50 /data/test.v\n");
    fprintf(stderr, "                   Per line of variant_file:\n");
    fprintf(stderr, "                   if deletion: chr, start, end+1\n");
    fprintf(stderr, "                   if insertion: chr, start, length\n");
    fprintf(stderr, "                   if indel: chr, type(D/I), start, end+1/length\n");
    fprintf(stderr, "Optional parameters:\n");
    fprintf(stderr, "         -d        Output more detail feature\n");
    fprintf(stderr, "         -m INT    maximum base deviation of breakpoint [100]\n");
    fprintf(stderr, "         -e FLOAT  maximum percent deviation of breakpoint (-e*sv_length) [0.1]\n");
    fprintf(stderr, "\n");
    return 1;
}



int main(int argc, char* argv[])
{
    time_t tstart = time(NULL);
    int mode = -1;      //deletion or insertion
    bool detail=false;
    double insertSize = -1;
    double insertSizeStd  = -1;
    FILE *fpbam=NULL;
    char c;
    double percent_diff = 0.1;
    int max_diff = 100;

    while((c = getopt(argc, argv,"diDIb:m:e:")) !=EOF )
    {
        switch (c)
        {
            case 'd': detail=true; break;
            case 'D': mode = 0; break;
            case 'I': mode = 1; break;
            case 'i': mode = 2; break;
            case 'b':
                openFile(&fpbam, optarg, "r"); break;
            case 'e': percent_diff = atof(optarg); break;
            case 'm': max_diff = atof(optarg); break;
            default: return usage();break;
        }
    }
    if(mode == -1 \
       || fpbam == NULL)
       return usage();


    char bamFilename[FILENAME_MAX];
    char svFilename[FILENAME_MAX];
    int readLen;
    int indiID = 1;
    while(fscanf(fpbam, "%s%d%lf%lf%s", bamFilename, &readLen, &insertSize, &insertSizeStd, svFilename) != EOF)
    {
        time_t istart = time(NULL);
        fprintf(stderr, "[SVfeature_main] Processing the %d: %s\n", indiID, bamFilename);
        BAM *bam = new BAM(readLen, bamFilename);
        IndelList *indelList = new IndelList(svFilename, mode);
        FILE *fpout, *fpstats;
        char outputFilename[FILENAME_MAX], statsFilename[FILENAME_MAX];
        sprintf(outputFilename, "%s_normalized", bamFilename);
        sprintf(statsFilename, "%s_absolute", bamFilename);
        openFile(&fpout,outputFilename,"w");
        openFile(&fpstats,statsFilename,"w");
        Feature feature(bam, indelList, fpout, fpstats, insertSize, insertSizeStd, detail);
        if(mode == 0)
            feature.deletionFeature(percent_diff, max_diff);
        if(mode == 1)
            feature.insertionFeature();
        if(mode == 2)
            feature.indelFeature();
        delete bam;
        delete indelList;
        fclose(fpout);
        fclose(fpstats);
        time_t iend = time(NULL);
        fprintf(stderr, "[SVfeature_main] The %d: %s finished and time used %ld sec.\n", indiID++, bamFilename, iend-istart);
    }

    time_t tend = time(NULL);
    fprintf(stderr, "[SVfeature_main] Total %d finished and time used %ld sec.\n", indiID - 1, tend-tstart);
    return 0;
}




