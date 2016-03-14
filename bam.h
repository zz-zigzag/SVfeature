#ifndef INCLUDE_BAM_H
#define INCLUDE_BAM_H 1

#include <stdio.h>
#include <seqan/bam_io.h>
#include "file.h"

using namespace seqan;


class BAM
{
public:
    char *fileName;
    int recordCnt;
    int readLen;
    BamFileIn bamFileIn;
    BamIndex<Bai> index;
    BamHeader header;
    std::vector <BamAlignmentRecord *> bamRecord;
public:
    BAM(int len, char *fileName);
    int getIdByPos(int pos, char *chrName);
    int getMateId(const BamAlignmentRecord *);
    double getDepth(char *chrName, int left, int right);
    bool isUniquelyMapped(BamAlignmentRecord *record);
    bool readRecord(char *chrName, int pos, int End);
    void clearRecord();
    ~BAM();
};

BAM::BAM(int len, char *fileName)
{
    readLen = len;
    recordCnt = 0;
    this ->fileName = fileName;
    if(!open(bamFileIn, fileName))
    {
        fprintf(stderr, "[SVfeature_bam] ERROR: Could not open %s.\n", fileName);
        exit(1);
    }
    char bai1[FILENAME_MAX];
    char bai2[FILENAME_MAX];
    sprintf(bai1, "%s.bai", fileName);
    strcpy(bai2, fileName);
    bai2[strlen(fileName) - 1] = 'i';
    if (!open(index, bai1) && !open(index, bai2)) {
        fprintf(stderr, "[SVfeature_bam] ERROR: Could not open BAM index file.\n");
        exit(1);
    }
    readHeader(header, bamFileIn);
}
bool BAM::readRecord(char *chrName, int pos, int posEnd)
{
    clearRecord();
    int chrID = nameToId(contigNamesCache(context(bamFileIn)), chrName);
    //fprintf(stderr, "[SVfeature_bam] %s: %d-%d is reading...", chrName, pos, posEnd);
    bool hasAlignments;
    jumpToRegion(bamFileIn, hasAlignments, chrID, pos, posEnd, index);
    BamAlignmentRecord *record;
    while (!atEnd(bamFileIn))
    {
        record = new BamAlignmentRecord();
        seqan::readRecord(*record, bamFileIn);
        bamRecord.push_back(record);
        if (record ->beginPos >posEnd)
            break;
    }
    recordCnt = bamRecord.size();
    //fprintf(stderr, "OK\n[SVfeature_bam] Total %d reads.\n", recordCnt);
    return hasAlignments;
}
void BAM::clearRecord()
{
    for(int i=0 ; i < recordCnt; i++)
        delete bamRecord[i];
    bamRecord.clear();
    recordCnt = 0;
}


int BAM::getIdByPos(int pos, char *chrName)
{
    pos -= 1;
    int chrID = nameToId(contigNamesCache(context(bamFileIn)), chrName);
    int l = 0;
    int r = recordCnt - 1;
    int m;
    if(bamRecord[l] -> rID > chrID ||\
        (bamRecord[l] -> rID == chrID &&\
         bamRecord[l]->beginPos >= pos)) return l;
    if(bamRecord[r] -> rID < chrID ||\
        (bamRecord[r] -> rID == chrID &&\
         bamRecord[r]->beginPos <= pos)) return r;
    while(l < r - 1)
    {
        m = (l+r)/2;
        if(bamRecord[m] -> rID == chrID &&\
            bamRecord[m]->beginPos == pos)  return m;
        else if(bamRecord[m] -> rID > chrID ||\
            (bamRecord[m] -> rID == chrID &&\
             bamRecord[m]->beginPos > pos)) r = m;
        else l = m;
    }
    return l;
}
int BAM::getMateId(const BamAlignmentRecord *record)
{
    int pos = record->pNext;
    int chrID = record->rID;
    int l = 0;
    int r = recordCnt - 1;
    int m;
    if(bamRecord[l] -> rID > chrID ||\
        (bamRecord[l] -> rID == chrID &&\
         bamRecord[l]->beginPos >= pos)) return -1;
    if(bamRecord[r] -> rID < chrID ||\
        (bamRecord[r] -> rID == chrID &&\
         bamRecord[r]->beginPos <= pos)) return -1;
    while(l < r - 1)
    {
        m = (l+r)/2;
        if(bamRecord[m] -> rID == chrID &&\
            bamRecord[m]->beginPos == pos) {
                int i = 0;
                while (bamRecord[m + i] ->qName !=  record->qName && bamRecord[m + i]->beginPos == pos)
                    i++;
                if (bamRecord[m + i] ->qName ==  record->qName)
                    return m + i;
                i = 1;
                while (bamRecord[m - i] ->qName !=  record->qName && bamRecord[m - i]->beginPos == pos)
                    i++;
                if (bamRecord[m - i] ->qName ==  record->qName)
                    return m - i;
                return -1;
            }
        else if(bamRecord[m] -> rID > chrID ||\
            (bamRecord[m] -> rID == chrID &&\
             bamRecord[m]->beginPos > pos)) r = m;
        else l = m;
    }
    return -1;
}


double BAM::getDepth(char *chrName, int from, int to)
{
    //int regionLen = to - from + 1;
    if(from > to)
    {
        fprintf(stderr, "[SVfeature_bam]    getDepth %s:%d-%d error\n", chrName, from, to);
        exit(1);
    }
    int regionLen = 0;
    int depth;
    long long sumDepth = 0;
    double avgDepth;
    char *sys = (char *)malloc(strlen(fileName) + FILENAME_MAX);
    sprintf(sys, "samtools depth -r %s:%d-%d %s", chrName, from, to, this->fileName);
    FILE *fp = popen(sys, "r");
    while(fscanf(fp, "%*s%*d%d", &depth) !=EOF)
    {
        sumDepth += depth;
        regionLen++;
    }
    if (regionLen == 0)
        return 0.0;
    avgDepth = sumDepth*1.0 / regionLen;
    pclose(fp);
    free(sys);
    return avgDepth;
}

bool BAM::isUniquelyMapped(BamAlignmentRecord *record)
{
    int id;
    char c;
    BamTagsDict tags(record->tags);
    if(findTagKey(id, tags, "XT"))
    {
        extractTagValue(c, tags, id);
        return (c == 'U') ? true : false;
    }
    else
        return true;    //if don't have XT tag, default unique.
}

BAM::~BAM()
{
    clearRecord();
}


#endif // INCLUDE_BAM_H
