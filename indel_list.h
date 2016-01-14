#ifndef INCLUDE_INDEL_LIST_H
#define INCLUDE_INDEL_LIST_H 1

#include <vector>
#include <string.h>
#include <stdlib.h>
#include "indel_record.h"
#include "file.h"


class IndelList     // Indel list
{
public:
    char *fileName;
    IndelList(char* fileName, int mode);

public:
    std::vector <IndelRecord *> indelRecord;
    int indelCnt;
    ~IndelList();
};

IndelList::IndelList(char* fileName, int mode)
{
    this ->fileName = fileName;
    long fileSize = getFileSize(fileName);
    char *buffer = new char[fileSize];
    FILE *fp;
    openFile(&fp, fileName, "r");
    fprintf(stderr, "[SVfeature_indel] %s is reading...", fileName);
    while(fgets(buffer, fileSize, fp)!=NULL)
    {
        if(buffer[0] == '\n') continue;
        indelRecord.push_back(new IndelRecord(buffer, mode));
    }
    indelCnt = indelRecord.size();
    fprintf(stderr, "OK\n[SVfeature_indel] Total %d variants.\n", indelCnt);
    delete buffer;
    fclose(fp);
}

IndelList::~IndelList()
{
    for(int i=0 ; i < indelCnt; i++)
        delete indelRecord[i];
}

#endif // INCLUDE_INDEL_H_
