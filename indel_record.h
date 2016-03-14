#ifndef INCLUDE_INDEL_RECORD_H
#define INCLUDE_INDEL_RECORD_H 1

#include <vector>
#include <stdio.h>
#include "file.h"

class IndelRecord
{
public:
    IndelRecord(char *buffer, int mode);
public:
    std::pair <int, int> breakPoint;
    int isExist, len;
    char type[2];   //'D' or 'I'
    char chrName[20];
};

IndelRecord::IndelRecord(char* buffer, int mode)
{
    int lefPos, rigPos;
    if(mode == 0)
    {
        sscanf(buffer, "%s%d%d", chrName, &lefPos, &rigPos);
        len = rigPos - lefPos - 1;
    }
    else if(mode == 1)
    {
        sscanf(buffer, "%s%d%d", chrName, &lefPos, &len);
        rigPos = lefPos + 1;
    }
    else
    {
        sscanf(buffer, "%s%s%d%d", chrName, type, &lefPos, &rigPos);
        if(type[0] == 'D')
        {
            len = rigPos - lefPos - 1;
        }
        else
        {
            len = rigPos;
            rigPos = lefPos + 1;
        }
    }
    breakPoint = std::make_pair(lefPos, rigPos);
}

#endif
