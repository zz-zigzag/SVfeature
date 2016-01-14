#ifndef INCLUDE_FILE_H
#define INCLUDE_FILE_H

#include <sys/stat.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>

void openFile(FILE **file, const char *filename, const char *type)
{
    if((*file = fopen(filename, type)) == NULL)
    {
       perror(filename);
       exit(1);
    }
    return ;
}

int getFileRowsNumber(const char *filename)
{
    char *sys = (char *)malloc(strlen(filename) + 10);
    sprintf(sys,"wc -l %s", filename);
    FILE *fp = popen(sys, "r");
    int n;
    if(fscanf(fp, "%d", &n) != 1)
       exit(1);
    pclose(fp);
    free(sys);
    return n;
}
int getFileSize(const char *filename)
{
    struct stat fileStat;
    int result = stat(filename, &fileStat);
    if(result != 0)
    {
       perror(filename);
       exit(1);
    }
    return fileStat.st_size;
}

#endif // INCLUDE_FILE_H
