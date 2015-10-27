#ifndef INCLUDE_FEATURE_H
#define INCLUDE_FEATURE_H 1

#include "bam.h"
#include "indel_list.h"
#include "file.h"

using namespace seqan;


const double EXP=1e-5;
#define ISEXP(a) if(a==0)a=EXP
class Feature
{
private:
    BAM *bam;
    IndelList *indelList;
    FILE *fpout, *fpstats;
    int readLen;
    int insertSize;
    int insertSizeStd;
    bool DETAIL;
    int _diff_ ;      // Allow a 1/5 length deviation

public:
    Feature(BAM *b, IndelList *i, FILE *o,FILE *o2, double m, double stdv,bool detail);
    void deletionFeature();
    void insertionFeature();
    void indelFeature();
    void getRegionStats(char *chrName, int start, int end, int &um_ef, int &um_e, int &mm_ef, int &mm_e);

};

Feature::Feature(BAM *b, IndelList *i, FILE *o, FILE *o2, double m, double stdv,bool detail)
{
    bam = b; indelList = i; fpout = o; fpstats=o2; readLen = b ->readLen;
    insertSize = (int)m; insertSizeStd = (int)stdv;
    DETAIL=detail;
}

void Feature::deletionFeature()
{
    fprintf(stderr, "[SVfeature_feature] Collecting features of variants.\n");

    std::vector <IndelRecord *> &indel = indelList->indelRecord;
    std::vector <BamAlignmentRecord *> &bamR = bam->bamRecord;
    int range = insertSize + 3*insertSizeStd;
    for(int i = 0; i < indelList->indelCnt; i++)
    {
        fprintf(stderr, "[SVfeature_feature] the %d variant is collecting...", i + 1);
        fprintf(fpout,"%d\t",indel[i]->isExist);
        fprintf(fpstats,"%d\t",indel[i]->isExist);

        /**
         * um/mm: uniquely mapped/multiple mapped
         * ef/e: error free/error
         * discord/concord: discordant/concordant
         */
        int um_ef_discord, um_ef_concord, um_e_discord, um_e_concord ;
        int mm_ef_discord, mm_ef_concord, mm_e_discord, mm_e_concord ;
        um_ef_discord = um_ef_concord = um_e_discord = um_e_concord = 0;
        mm_ef_discord = mm_ef_concord = mm_e_discord = mm_e_concord = 0;
        double mapQ_um_ef_discord, mapQ_um_ef_concord, mapQ_um_e_discord, mapQ_um_e_concord;
        double mapQ_mm_ef_discord, mapQ_mm_ef_concord, mapQ_mm_e_discord, mapQ_mm_e_concord;
        mapQ_um_ef_discord = mapQ_um_ef_concord = mapQ_um_e_discord = mapQ_um_e_concord = 0;
        mapQ_mm_ef_discord = mapQ_mm_ef_concord = mapQ_mm_e_discord = mapQ_mm_e_concord = 0;
        /**
         * l/r: left/right
         * upAn/downAn: anchor in up/down
         * um/mm: uniquely mapped/multiple mapped
         * clip/full/other:
         */
        int l_upAn_um_clip, l_upAn_um_full, l_upAn_um_other, l_upAn_mm_clip, l_upAn_mm_full, l_upAn_mm_other;
        int l_downAn_um_clip, l_downAn_um_full, l_downAn_um_other, l_downAn_mm_clip, l_downAn_mm_full, l_downAn_mm_other;
        int r_upAn_um_clip, r_upAn_um_full, r_upAn_um_other, r_upAn_mm_clip, r_upAn_mm_full, r_upAn_mm_other;
        int r_downAn_um_clip, r_downAn_um_full, r_downAn_um_other, r_downAn_mm_clip, r_downAn_mm_full, r_downAn_mm_other;
        l_upAn_um_clip = l_upAn_um_full = l_upAn_um_other = l_upAn_mm_clip = l_upAn_mm_full = l_upAn_mm_other = 0;
        l_downAn_um_clip = l_downAn_um_full = l_downAn_um_other = l_downAn_mm_clip = l_downAn_mm_full = l_downAn_mm_other = 0;
        r_upAn_um_clip = r_upAn_um_full = r_upAn_um_other = r_upAn_mm_clip = r_upAn_mm_full = r_upAn_mm_other = 0;
        r_downAn_um_clip = r_downAn_um_full = r_downAn_um_other = r_downAn_mm_clip = r_downAn_mm_full = r_downAn_mm_other = 0;

        int indelLen = indel[i] ->len;
        int brk1 = indel[i] ->breakPoint.first;
        int brk2 = indel[i] ->breakPoint.second;
        int brk0 = brk1 - indelLen + 1;
        if(brk0 < 1) brk0 = 1;
        int brk3 = brk2 + indelLen - 1;
        _diff_ = indelLen/10;
        if(_diff_ > 100)    _diff_ = 100;
        int rightMost = brk2 + range - readLen + _diff_ ;                  //_diff_
        int leftMost = brk1 - range - _diff_ ;                             //_diff_
        if(leftMost < 1)    leftMost = 1;
        bam->readRecord(indel[i] ->chrName, std::min(brk0, leftMost) - range, std::max(brk3, rightMost) + range);
        for(int j = bam->getIdByPos(leftMost, indel[i] ->chrName);\
             j < bam->recordCnt && bamR[j]->beginPos < rightMost; j++)
        {
            if((bamR[j]->flag & 0x4) || (bamR[j]->flag & 0x8))          //the read unmapped or mate unmapped
                continue;
            if(bamR[j]->rNextId != bamR[j]->rID)                        //the read and mate mapped in different chromosome
                continue;
            int pos  = bamR[j]->beginPos + 1;                           //0-based, so + 1
            int mpos = bamR[j]->pNext + 1;
            int mID = bam->getMateId(bamR[j]);
            bool unique_flag, unique_mate_flag;
            if (mID != -1)
                unique_mate_flag = bam->isUniquelyMapped(bamR[mID]);
            unique_flag  = bam->isUniquelyMapped(bamR[j]);

            if(mID != -1 && pos < mpos && pos < brk2 - _diff_ && mpos > brk1 + _diff_)       // _diff_
            {
                if(unique_flag && unique_mate_flag)
                {
                    if(bamR[j]->cigar[0].count == (uint)readLen && bamR[j]->cigar[0].operation == 'M' &&\
                       bamR[mID]->cigar[0].count == (uint)readLen && bamR[mID]->cigar[0].operation == 'M' )
                    {

                        if(bamR[j]->tLen > range)
                        {
                            um_ef_discord++;
                            mapQ_um_ef_discord += (bamR[j]->mapQ+bamR[mID]->mapQ)/2.0;
                        }
                        else
                        {
                            um_ef_concord++;
                            mapQ_um_ef_concord += (bamR[j]->mapQ+bamR[mID]->mapQ)/2.0;
                        }
                    }
                    else
                    {
                        if(bamR[j]->tLen > range)
                        {
                            um_e_discord++;
                            mapQ_um_e_discord += (bamR[j]->mapQ+bamR[mID]->mapQ)/2.0;
                        }
                        else
                        {
                            um_e_concord++;
                            mapQ_um_e_concord += (bamR[j]->mapQ+bamR[mID]->mapQ)/2.0;
                        }
                    }
                }
                else
                {
                    if(bamR[j]->cigar[0].count == (uint)readLen && bamR[j]->cigar[0].operation == 'M' &&\
                       bamR[mID]->cigar[0].count == (uint)readLen && bamR[mID]->cigar[0].operation == 'M' )
                    {

                        if(bamR[j]->tLen > range)
                        {
                            mm_ef_discord++;
                            mapQ_mm_ef_discord += (bamR[j]->mapQ+bamR[mID]->mapQ)/2.0;
                        }
                        else
                        {
                            mm_ef_concord++;
                            mapQ_mm_ef_concord += (bamR[j]->mapQ+bamR[mID]->mapQ)/2.0;
                        }
                    }
                    else
                    {
                        if(bamR[j]->tLen > range)
                        {
                            mm_e_discord++;
                            mapQ_mm_e_discord += (bamR[j]->mapQ+bamR[mID]->mapQ)/2.0;
                        }
                        else
                        {
                            mm_e_concord++;
                            mapQ_mm_e_concord += (bamR[j]->mapQ+bamR[mID]->mapQ)/2.0;
                        }
                    }
                }
            }

            if(pos > brk1-readLen && pos <= brk1 + _diff_)
            {
                if(bamR[j]->tLen > 0)       //anchor in down
                {
                    if(unique_flag)
                    {
                        if(bamR[j]->cigar[0].operation == 'M' && bamR[j]->cigar[0].count == (uint)readLen)
                        {
                            l_upAn_um_full++;
                        }
                        else if(bamR[j]->cigar[bamR[j]->_n_cigar-1].operation == 'S')
                        {
                            l_upAn_um_clip++;
                        }
                        else
                        {
                            l_upAn_um_other++;
                        }
                    }
                    else
                    {
                        if(bamR[j]->cigar[0].operation == 'M' && bamR[j]->cigar[0].count == (uint)readLen)
                        {
                            l_upAn_mm_full++;
                        }
                        else if(bamR[j]->cigar[0].operation == 'S')
                        {
                            l_upAn_mm_clip++;
                        }
                        else
                        {
                            l_upAn_mm_other++;
                        }
                    }
                }
                else
                {
                    if(unique_flag)
                    {
                        if(bamR[j]->cigar[0].operation == 'M' && bamR[j]->cigar[0].count == (uint)readLen)
                        {
                            l_downAn_um_full++;
                        }
                        else if(bamR[j]->cigar[bamR[j]->_n_cigar-1].operation == 'S')
                        {
                            l_downAn_um_clip++;
                        }
                        else
                        {
                            l_downAn_um_other++;
                        }
                    }
                    else
                    {
                        if(bamR[j]->cigar[0].operation == 'M' && bamR[j]->cigar[0].count == (uint)readLen)
                        {
                            l_downAn_mm_full++;
                        }
                        else if(bamR[j]->cigar[0].operation == 'S')
                        {
                            l_downAn_mm_clip++;
                        }
                        else
                        {
                            l_downAn_mm_other++;
                        }
                    }
                }

            }
            else if(pos >= brk2-_diff_ && pos <= brk2 + _diff_)
            {
                if(bamR[j]->tLen > 0)       //anchor in down
                {
                    if(unique_flag)
                    {
                        if(bamR[j]->cigar[0].operation == 'M' && bamR[j]->cigar[0].count == (uint)readLen)
                        {
                            r_upAn_um_full++;
                        }
                        else if(bamR[j]->cigar[0].operation == 'S')
                        {
                            r_upAn_um_clip++;
                        }
                        else
                        {
                            r_upAn_um_other++;
                        }
                    }
                    else
                    {
                        if(bamR[j]->cigar[0].operation == 'M' && bamR[j]->cigar[0].count == (uint)readLen)
                        {
                            r_upAn_mm_full++;
                        }
                        else if(bamR[j]->cigar[0].operation == 'S')
                        {
                            r_upAn_mm_clip++;
                        }
                        else
                        {
                            r_upAn_mm_other++;
                        }
                    }
                }
                else
                {
                    if(unique_flag)
                    {
                        if(bamR[j]->cigar[0].operation == 'M' && bamR[j]->cigar[0].count == (uint)readLen)
                        {
                            r_downAn_um_full++;
                        }
                        else if(bamR[j]->cigar[0].operation == 'S')
                        {
                            r_downAn_um_clip++;
                        }
                        else
                        {
                            r_downAn_um_other++;
                        }
                    }
                    else
                    {
                        if(bamR[j]->cigar[0].operation == 'M' && bamR[j]->cigar[0].count == (uint)readLen)
                        {
                            r_downAn_mm_full++;
                        }
                        else if(bamR[j]->cigar[0].operation == 'S')
                        {
                            r_downAn_mm_clip++;
                        }
                        else
                        {
                            r_downAn_mm_other++;
                        }
                    }
                }
            }

        }
        if(um_ef_discord)   mapQ_um_ef_discord /= um_ef_discord;
        if(um_ef_concord)   mapQ_um_ef_concord/= um_ef_concord;
        if(um_e_discord)   mapQ_um_e_discord /= um_e_discord;
        if(um_e_concord)   mapQ_um_e_concord /= um_e_concord;
        if(mm_ef_discord)   mapQ_mm_ef_discord /= mm_ef_discord;
        if(mm_ef_concord)   mapQ_mm_ef_concord /= mm_ef_concord;
        if(mm_e_discord)   mapQ_mm_e_discord /= mm_e_discord;
        if(mm_e_concord)   mapQ_mm_e_concord /= mm_e_concord;

        //print read pair
        fprintf(fpstats, "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t",um_ef_discord,um_ef_concord,um_e_discord,um_e_concord,\
                mm_ef_discord,mm_ef_concord,mm_e_discord,mm_e_concord);
/*
        fprintf(fpstats, "%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t",mapQ_um_ef_discord,mapQ_um_ef_concord,mapQ_um_e_discord,\
                mapQ_um_e_concord,mapQ_mm_ef_discord,mapQ_mm_ef_concord,mapQ_mm_e_discord,mapQ_mm_e_concord);
*/
        double all;
        all=um_ef_discord+um_ef_concord+um_e_discord+um_e_concord+mm_ef_discord+mm_ef_concord+mm_e_discord+mm_e_concord;
        if(all==0)all=EXP;
        fprintf(fpout, "%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t",um_ef_discord/all,um_ef_concord/all,\
                um_e_discord/all,um_e_concord/all,mm_ef_discord/all,mm_ef_concord/all,\
                mm_e_discord/all,mm_e_concord/all);
        if(DETAIL)
        {
            double um_ef, um_e, mm_ef, mm_e, mm, um;
            double um_dis, um_con, mm_dis, mm_con, dis, con;
            double ef_dis, ef_con, e_dis, e_con, ef, e;
            um_ef=um_ef_concord+um_ef_discord;  ISEXP(um_ef);
            um_e=um_e_concord+um_e_discord;     ISEXP(um_e);
            um_dis=um_ef_discord+um_e_discord;  ISEXP(um_dis);
            um_con=um_ef_concord+um_e_concord;  ISEXP(um_con);

            mm_ef=mm_ef_concord+mm_ef_discord;  ISEXP(mm_ef);
            mm_e=mm_e_concord+mm_e_discord;     ISEXP(mm_e);
            mm_dis=mm_e_discord+mm_ef_discord;  ISEXP(mm_dis);
            mm_con=mm_e_concord+mm_ef_concord;  ISEXP(mm_con);

            ef_dis=um_ef_discord+mm_ef_discord; ISEXP(ef_dis);
            ef_con=um_ef_concord+mm_ef_concord; ISEXP(ef_con);
            e_dis=um_e_discord+mm_e_discord;    ISEXP(e_dis);
            e_con=um_e_concord+mm_e_concord;    ISEXP(e_con);

            //two Equation for each var
            ef=ef_dis+ef_con;   ISEXP(ef);
            e=e_dis+e_con;      ISEXP(e);
            um=um_ef+um_e;      ISEXP(um);
            mm=mm_ef+mm_e;      ISEXP(mm);
            dis=um_dis+mm_dis;  ISEXP(dis);
            con=um_con+mm_con;  ISEXP(con);

            fprintf(fpout, "%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t",\
                    um_ef_concord/um_ef,um_ef_discord/um_ef,\
                    um_e_concord/um_e,um_e_discord/um_e,\
                    um_ef_discord/um_dis,um_e_discord/um_dis,\
                    um_ef_concord/um_con,um_e_concord/um_con);
            fprintf(fpout, "%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t",\
                    mm_ef_concord/mm_ef,mm_ef_discord/mm_ef,\
                    mm_e_concord/mm_e,mm_e_discord/mm_e,\
                    mm_ef_discord/mm_dis,mm_e_discord/mm_dis,\
                    mm_ef_concord/mm_con,mm_e_concord/mm_con);
            fprintf(fpout, "%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t",\
                    um_ef_discord/ef_dis,mm_ef_discord/ef_dis,\
                    um_ef_concord/ef_con,mm_ef_concord/ef_con,\
                    um_e_discord/e_dis,mm_e_discord/e_dis,
                    um_e_concord/e_con,mm_e_concord/e_con);
            //two Equation for each var
            fprintf(fpout, "%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t",\
                    um_dis/um,um_con/um,mm_dis/mm,mm_con/mm,\
                    um_ef/um,um_e/um,mm_ef/mm,mm_e/mm,\
                    ef_dis/ef,ef_con/ef,e_dis/e,e_con/e);
            fprintf(fpout, "%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t",\
                    um_dis/dis,mm_dis/dis,um_con/con,mm_con/con,\
                    um_ef/ef,mm_ef/ef,um_e/e,mm_e/e,\
                    ef_dis/dis,ef_con/con,e_dis/dis,e_con/con);
            fprintf(fpout, "%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t",\
                    ef/all,e/all,um/all,mm/all,dis/all,con/all);
        }

        //print split read
        fprintf(fpstats, "%d\t%d\t%d\t%d\t%d\t%d\t",l_upAn_um_clip,l_upAn_um_full,l_upAn_um_other,\
                l_upAn_mm_clip,l_upAn_mm_full,l_upAn_mm_other);
        fprintf(fpstats, "%d\t%d\t%d\t%d\t%d\t%d\t",l_downAn_um_clip,l_downAn_um_full,l_downAn_um_other,\
                l_downAn_mm_clip,l_downAn_mm_full,l_downAn_mm_other);
        fprintf(fpstats, "%d\t%d\t%d\t%d\t%d\t%d\t",r_upAn_um_clip,r_upAn_um_full,r_upAn_um_other,\
                r_upAn_mm_clip,r_upAn_mm_full,r_upAn_mm_other);
        fprintf(fpstats, "%d\t%d\t%d\t%d\t%d\t%d\t",r_downAn_um_clip,r_downAn_um_full,r_downAn_um_other,\
                r_downAn_mm_clip,r_downAn_mm_full,r_downAn_mm_other);


        double l_all,r_all;
        l_all=l_upAn_um_clip+l_upAn_um_full+l_upAn_um_other+l_upAn_mm_clip+l_upAn_mm_full+l_upAn_mm_other+\
                l_downAn_um_clip+l_downAn_um_full+l_downAn_um_other+l_downAn_mm_clip+l_downAn_mm_full+l_downAn_mm_other;
        r_all=r_upAn_um_clip+r_upAn_um_full+r_upAn_um_other+r_upAn_mm_clip+r_upAn_mm_full+r_upAn_mm_other+\
                r_downAn_um_clip+r_downAn_um_full+r_downAn_um_other+r_downAn_mm_clip+r_downAn_mm_full+r_downAn_mm_other;
        ISEXP(l_all);
        ISEXP(r_all);
        fprintf(fpout, "%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t",\
                l_upAn_um_clip/l_all,l_upAn_um_full/l_all,l_upAn_um_other/l_all,\
                l_upAn_mm_clip/l_all,l_upAn_mm_full/l_all,l_upAn_mm_other/l_all);
        fprintf(fpout, "%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t",\
                l_downAn_um_clip/l_all,l_downAn_um_full/l_all,l_downAn_um_other/l_all,\
                l_downAn_mm_clip/l_all,l_downAn_mm_full/l_all,l_downAn_mm_other/l_all);
        fprintf(fpout, "%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t",\
                r_upAn_um_clip/r_all,r_upAn_um_full/r_all,r_upAn_um_other/r_all,\
                r_upAn_mm_clip/r_all,r_upAn_mm_full/r_all,r_upAn_mm_other/r_all);
        fprintf(fpout, "%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t",\
                r_downAn_um_clip/r_all,r_downAn_um_full/r_all,r_downAn_um_other/r_all,\
                r_downAn_mm_clip/r_all,r_downAn_mm_full/r_all,r_downAn_mm_other/r_all);
        if(DETAIL)
        {
            //left breakpoint
            double l_upAn_um,l_upAn_mm,l_downAn_um,l_downAn_mm;
            double l_upAn_clip,l_upAn_full,l_upAn_other,l_downAn_clip,l_downAn_full,l_downAn_other;
            double l_um_clip,l_um_full,l_um_other,l_mm_clip,l_mm_full,l_mm_other;
            double l_upAn,l_downAn,l_um,l_mm,l_clip,l_full,l_other;
            double l;
            l_upAn_um=l_upAn_um_clip+l_upAn_um_full+l_upAn_um_other;    ISEXP(l_upAn_um);
            l_upAn_mm=l_upAn_mm_clip+l_upAn_mm_full+l_upAn_mm_other;    ISEXP(l_upAn_mm);
            l_downAn_um=l_downAn_um_clip+l_downAn_um_full+l_downAn_um_other;    ISEXP(l_downAn_um);
            l_downAn_mm=l_downAn_mm_clip+l_downAn_mm_full+l_downAn_mm_other;    ISEXP(l_downAn_mm);

            l_upAn_clip=l_upAn_um_clip+l_upAn_mm_clip;  ISEXP(l_upAn_clip);
            l_upAn_full=l_upAn_um_full+l_upAn_mm_full;  ISEXP(l_upAn_full);
            l_upAn_other=l_upAn_um_other+l_upAn_mm_other;   ISEXP(l_upAn_other);
            l_downAn_clip=l_downAn_um_clip+l_downAn_mm_clip;    ISEXP(l_downAn_clip);
            l_downAn_full=l_downAn_um_full+l_downAn_mm_full;    ISEXP(l_downAn_full);
            l_downAn_other=l_downAn_um_other+l_downAn_mm_other; ISEXP(l_downAn_other);

            l_um_clip=l_upAn_um_clip+l_downAn_um_clip;  ISEXP(l_um_clip);
            l_um_full=l_upAn_um_full+l_downAn_um_full;  ISEXP(l_um_full);
            l_um_other=l_upAn_um_other+l_downAn_um_other;   ISEXP(l_um_other);
            l_mm_clip=l_upAn_mm_clip+l_downAn_mm_clip;  ISEXP(l_mm_clip);
            l_mm_full=l_upAn_mm_full+l_downAn_mm_full;  ISEXP(l_mm_full);
            l_mm_other=l_upAn_mm_other+l_downAn_mm_other;   ISEXP(l_mm_other);

            //two Equation for each var
            l_upAn=l_upAn_clip+l_upAn_full+l_upAn_other;    ISEXP(l_upAn);         //l_upAn=l_upAn_um+l_upAn_mm;
            l_downAn=l_downAn_clip+l_downAn_full+l_downAn_other;    ISEXP(l_downAn);//l_downAn=l_downAn_um+l_downAn_mm;
            l_um=l_um_clip+l_um_full+l_um_other;    ISEXP(l_um);
            l_mm=l_mm_clip+l_mm_full+l_mm_other;    ISEXP(l_mm);
            l_clip=l_upAn_clip+l_downAn_clip;   ISEXP(l_clip);
            l_full=l_upAn_full+l_downAn_full;   ISEXP(l_full);
            l_other=l_upAn_other+l_downAn_other;ISEXP(l_other);

            l=l_upAn+l_downAn;  ISEXP(l);

            //  *_*_*_*
            fprintf(fpout,"%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t",\
                    l_upAn_um_clip/l_upAn_um,l_upAn_um_full/l_upAn_um,l_upAn_um_other/l_upAn_um,\
                    l_upAn_mm_clip/l_upAn_mm,l_upAn_mm_full/l_upAn_mm,l_upAn_mm_other/l_upAn_mm,\
                    l_downAn_um_clip/l_downAn_um,l_downAn_um_full/l_downAn_um,l_downAn_um_other/l_downAn_um,\
                    l_downAn_mm_clip/l_downAn_mm,l_downAn_mm_full/l_downAn_mm,l_downAn_mm_other/l_downAn_mm);
            fprintf(fpout,"%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t",\
                    l_upAn_um_clip/l_upAn_clip,l_upAn_mm_clip/l_upAn_clip,\
                    l_upAn_um_full/l_upAn_full,l_upAn_mm_full/l_upAn_full,\
                    l_upAn_um_other/l_upAn_other,l_upAn_mm_other/l_upAn_other,\
                    l_downAn_um_clip/l_downAn_clip,l_downAn_mm_clip/l_downAn_clip,\
                    l_downAn_um_full/l_downAn_full,l_downAn_mm_full/l_downAn_full,\
                    l_downAn_um_other/l_downAn_other,l_downAn_mm_other/l_downAn_other);
            fprintf(fpout,"%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t",\
                    l_upAn_um_clip/l_um_clip,l_downAn_um_clip/l_um_clip,\
                    l_upAn_um_full/l_um_full,l_downAn_um_full/l_um_full,\
                    l_upAn_um_other/l_um_other,l_downAn_um_other/l_um_other,\
                    l_upAn_mm_clip/l_mm_clip,l_downAn_mm_clip/l_mm_clip,\
                    l_upAn_mm_full/l_mm_full,l_downAn_mm_full/l_mm_full,\
                    l_upAn_mm_other/l_mm_other,l_downAn_mm_other/l_mm_other);
            //  *_*_*
            fprintf(fpout,"%.2f\t%.2f\t%.2f\t%.2f\t",\
                    l_upAn_um/l_upAn,l_upAn_mm/l_upAn,\
                    l_downAn_um/l_downAn,l_downAn_mm/l_downAn);
            fprintf(fpout,"%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t",\
                    l_upAn_clip/l_upAn,l_upAn_full/l_upAn,l_upAn_other/l_upAn,\
                    l_downAn_clip/l_downAn,l_downAn_full/l_downAn,l_downAn_other/l_downAn);
            fprintf(fpout,"%.2f\t%.2f\t%.2f\t%.2f\t",\
                    l_upAn_um/l_um,l_downAn_um/l_um,\
                    l_upAn_mm/l_mm,l_downAn_mm/l_mm);
            fprintf(fpout,"%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t",\
                    l_um_clip/l_um,l_um_full/l_um,l_um_other/l_um,\
                    l_mm_clip/l_mm,l_mm_full/l_mm,l_mm_other/l_mm);
            fprintf(fpout,"%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t",\
                    l_upAn_clip/l_clip,l_downAn_clip/l_clip,\
                    l_upAn_full/l_full,l_upAn_full/l_full,\
                    l_upAn_other/l_other,l_downAn_other/l_other);
            fprintf(fpout,"%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t",\
                    l_um_clip/l_clip,l_mm_clip/l_clip,\
                    l_um_full/l_full,l_mm_full/l_full,\
                    l_um_other/l_other,l_mm_other/l_other);
            //  *_*
            fprintf(fpout,"%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t",\
                    l_upAn/l,l_downAn/l,l_um/l,l_mm/l,l_clip/l,l_full/l,l_other/l);

            //right breakpoint
            double r_upAn_um,r_upAn_mm,r_downAn_um,r_downAn_mm;
            double r_upAn_clip,r_upAn_full,r_upAn_other,r_downAn_clip,r_downAn_full,r_downAn_other;
            double r_um_clip,r_um_full,r_um_other,r_mm_clip,r_mm_full,r_mm_other;
            double r_upAn,r_downAn,r_um,r_mm,r_clip,r_full,r_other;
            double r;
            r_upAn_um=r_upAn_um_clip+r_upAn_um_full+r_upAn_um_other;    ISEXP(r_upAn_um);
            r_upAn_mm=r_upAn_mm_clip+r_upAn_mm_full+r_upAn_mm_other;    ISEXP(r_upAn_mm);
            r_downAn_um=r_downAn_um_clip+r_downAn_um_full+r_downAn_um_other;    ISEXP(r_downAn_um);
            r_downAn_mm=r_downAn_mm_clip+r_downAn_mm_full+r_downAn_mm_other;    ISEXP(r_downAn_mm);

            r_upAn_clip=r_upAn_um_clip+r_upAn_mm_clip;  ISEXP(r_upAn_clip);
            r_upAn_full=r_upAn_um_full+r_upAn_mm_full;  ISEXP(r_upAn_full);
            r_upAn_other=r_upAn_um_other+r_upAn_mm_other;   ISEXP(r_upAn_other);
            r_downAn_clip=r_downAn_um_clip+r_downAn_mm_clip;    ISEXP(r_downAn_clip);
            r_downAn_full=r_downAn_um_full+r_downAn_mm_full;    ISEXP(r_downAn_full);
            r_downAn_other=r_downAn_um_other+r_downAn_mm_other; ISEXP(r_downAn_other);

            r_um_clip=r_upAn_um_clip+r_downAn_um_clip;  ISEXP(r_um_clip);
            r_um_full=r_upAn_um_full+r_downAn_um_full;  ISEXP(r_um_full);
            r_um_other=r_upAn_um_other+r_downAn_um_other;   ISEXP(r_um_other);
            r_mm_clip=r_upAn_mm_clip+r_downAn_mm_clip;  ISEXP(r_mm_clip);
            r_mm_full=r_upAn_mm_full+r_downAn_mm_full;  ISEXP(r_mm_full);
            r_mm_other=r_upAn_mm_other+r_downAn_mm_other;   ISEXP(r_mm_other);

            //two Equation for each var
            r_upAn=r_upAn_clip+r_upAn_full+r_upAn_other;    ISEXP(r_upAn);         //r_upAn=r_upAn_um+r_upAn_mm;
            r_downAn=r_downAn_clip+r_downAn_full+r_downAn_other;    ISEXP(r_downAn);//r_downAn=r_downAn_um+r_downAn_mm;
            r_um=r_um_clip+r_um_full+r_um_other;    ISEXP(r_um);
            r_mm=r_mm_clip+r_mm_full+r_mm_other;    ISEXP(r_mm);
            r_clip=r_upAn_clip+r_downAn_clip;   ISEXP(r_clip);
            r_full=r_upAn_full+r_downAn_full;   ISEXP(r_full);
            r_other=r_upAn_other+r_downAn_other;ISEXP(r_other);

            r=r_upAn+r_downAn;  ISEXP(l);

            //  *_*_*_*
            fprintf(fpout,"%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t",\
                    r_upAn_um_clip/r_upAn_um,r_upAn_um_full/r_upAn_um,r_upAn_um_other/r_upAn_um,\
                    r_upAn_mm_clip/r_upAn_mm,r_upAn_mm_full/r_upAn_mm,r_upAn_mm_other/r_upAn_mm,\
                    r_downAn_um_clip/r_downAn_um,r_downAn_um_full/r_downAn_um,r_downAn_um_other/r_downAn_um,\
                    r_downAn_mm_clip/r_downAn_mm,r_downAn_mm_full/r_downAn_mm,r_downAn_mm_other/r_downAn_mm);
            fprintf(fpout,"%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t",\
                    r_upAn_um_clip/r_upAn_clip,r_upAn_mm_clip/r_upAn_clip,\
                    r_upAn_um_full/r_upAn_full,r_upAn_mm_full/r_upAn_full,\
                    r_upAn_um_other/r_upAn_other,r_upAn_mm_other/r_upAn_other,\
                    r_downAn_um_clip/r_downAn_clip,r_downAn_mm_clip/r_downAn_clip,\
                    r_downAn_um_full/r_downAn_full,r_downAn_mm_full/r_downAn_full,\
                    r_downAn_um_other/r_downAn_other,r_downAn_mm_other/r_downAn_other);
            fprintf(fpout,"%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t",\
                    r_upAn_um_clip/r_um_clip,r_downAn_um_clip/r_um_clip,\
                    r_upAn_um_full/r_um_full,r_downAn_um_full/r_um_full,\
                    r_upAn_um_other/r_um_other,r_downAn_um_other/r_um_other,\
                    r_upAn_mm_clip/r_mm_clip,r_downAn_mm_clip/r_mm_clip,\
                    r_upAn_mm_full/r_mm_full,r_downAn_mm_full/r_mm_full,\
                    r_upAn_mm_other/r_mm_other,r_downAn_mm_other/r_mm_other);
            //  *_*_*
            fprintf(fpout,"%.2f\t%.2f\t%.2f\t%.2f\t",\
                    r_upAn_um/r_upAn,r_upAn_mm/r_upAn,\
                    r_downAn_um/r_downAn,r_downAn_mm/r_downAn);
            fprintf(fpout,"%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t",\
                    r_upAn_clip/r_upAn,r_upAn_full/r_upAn,r_upAn_other/r_upAn,\
                    r_downAn_clip/r_downAn,r_downAn_full/r_downAn,r_downAn_other/r_downAn);
            fprintf(fpout,"%.2f\t%.2f\t%.2f\t%.2f\t",\
                    r_upAn_um/r_um,r_downAn_um/r_um,\
                    r_upAn_mm/r_mm,r_downAn_mm/r_mm);
            fprintf(fpout,"%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t",\
                    r_um_clip/r_um,r_um_full/r_um,r_um_other/r_um,\
                    r_mm_clip/r_mm,r_mm_full/r_mm,r_mm_other/r_mm);
            fprintf(fpout,"%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t",\
                    r_upAn_clip/r_clip,r_downAn_clip/r_clip,\
                    r_upAn_full/r_full,r_upAn_full/r_full,\
                    r_upAn_other/r_other,r_downAn_other/r_other);
            fprintf(fpout,"%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t",\
                    r_um_clip/r_clip,r_mm_clip/r_clip,\
                    r_um_full/r_full,r_mm_full/r_full,\
                    r_um_other/r_other,r_mm_other/r_other);
            //  *_*
            fprintf(fpout,"%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t",\
                    r_upAn/r,r_downAn/r,r_um/r,r_mm/r,r_clip/r,r_full/r,r_other/r);
        }



        /**
         *  in_depth, up_depth, down_depth
         */
        double in_depth, up_depth, down_depth;
        in_depth = up_depth = down_depth = 0;
        in_depth = bam->getDepth(indel[i] ->chrName, brk1 + 1, brk2 - 1);
        up_depth = bam->getDepth(indel[i] ->chrName, brk0, brk1);
        down_depth = bam->getDepth(indel[i] ->chrName, brk2, brk3);
        //print read depth
        fprintf(fpstats, "%.2f\t%.2f\t%.2f\t", in_depth,up_depth,down_depth);
        fprintf(fpout, "%.2f\t%.2f\t%.2f\t%.2f\t",up_depth/(in_depth+up_depth),in_depth/(in_depth+up_depth),\
                in_depth/(down_depth+in_depth),down_depth/(down_depth+in_depth));

        /**
         * in_um_ef, in_um_e, in_mm_ef, in_mm_e;
         * up_um_ef, up_um_e, up_mm_ef, up_mm_e;
         * down_um_ef, down_um_e, down_mm_ef, down_mm_e;
         */
        int in_um_ef, in_um_e, in_mm_ef, in_mm_e;
        int up_um_ef, up_um_e, up_mm_ef, up_mm_e;
        int down_um_ef, down_um_e, down_mm_ef, down_mm_e;
        in_um_ef = in_um_e = in_mm_ef = in_mm_e = 0;
        up_um_ef = up_um_e = up_mm_ef = up_mm_e = 0;
        down_um_ef = down_um_e = down_mm_ef = down_mm_e = 0;
        getRegionStats(indel[i] ->chrName, brk1 + 1, brk2 - 1, in_um_ef, in_um_e, in_mm_ef, in_mm_e);
        getRegionStats(indel[i] ->chrName, brk0, brk1, up_um_ef, up_um_e, up_mm_ef, up_mm_e);
        getRegionStats(indel[i] ->chrName, brk2, brk3, down_um_ef, down_um_e, down_mm_ef, down_mm_e);

        fprintf(fpstats, "%d\t%d\t%d\t%d\t",in_um_ef,in_um_e,in_mm_ef,in_mm_e);
        fprintf(fpstats, "%d\t%d\t%d\t%d\t",up_um_ef, up_um_e, up_mm_ef, up_mm_e);
        fprintf(fpstats, "%d\t%d\t%d\t%d\t",down_um_ef, down_um_e, down_mm_ef, down_mm_e);

        all=in_um_ef+in_um_e+in_mm_ef+in_mm_e+\
            up_um_ef+up_um_e+up_mm_ef+up_mm_e+\
            down_um_ef+down_um_e+down_mm_ef+down_mm_e;
        ISEXP(all);
        fprintf(fpout, "%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t",\
                in_um_ef/all,in_um_e/all,in_mm_ef/all,in_mm_e/all,\
                up_um_ef/all,up_um_e/all, up_mm_ef/all, up_mm_e/all,\
                down_um_ef/all,down_um_e/all,down_mm_ef/all,down_mm_e/all);
        if(DETAIL)
        {
            double in_um,in_mm,up_um,up_mm,down_um,down_mm;
            double in_ef,in_e,up_ef,up_e,down_ef,down_e;
            double um_ef,um_e,mm_ef,mm_e;
            double in,up,down,um,mm,ef,e;

            in_um=in_um_e+in_um_ef; ISEXP(in_um);
            in_mm=in_mm_e+in_mm_ef; ISEXP(in_mm);
            up_um=up_um_e+up_um_ef; ISEXP(up_um);
            up_mm=up_mm_e+up_mm_ef; ISEXP(up_mm);
            down_um=down_um_e+down_um_ef;   ISEXP(down_um);
            down_mm=down_mm_e+down_mm_ef;   ISEXP(down_mm);

            in_ef=in_um_ef+in_mm_ef;    ISEXP(in_ef);
            in_e=in_um_e+in_mm_e;       ISEXP(in_e);
            up_ef=up_um_ef+up_mm_ef;    ISEXP(up_ef);
            up_e=up_um_e+up_mm_e;       ISEXP(up_e);
            down_ef=down_um_ef+down_mm_ef;  ISEXP(down_ef);
            down_e=down_um_e+down_mm_e;     ISEXP(down_e);

            um_ef=in_um_ef+up_um_ef+down_um_ef; ISEXP(um_ef);
            um_e=in_um_e+up_um_e+down_um_e;     ISEXP(um_e);
            mm_ef=in_mm_ef+up_mm_ef+down_mm_ef; ISEXP(mm_ef);
            mm_e=in_mm_e+up_mm_e+down_mm_e;     ISEXP(mm_e);

            //two Equation for each var
            in=in_um+in_mm; ISEXP(in);
            up=up_um+up_mm; ISEXP(up);
            down=down_um+down_mm;   ISEXP(down);
            um=in_um+up_um+down_um; ISEXP(um);
            mm=in_mm+up_mm+down_mm; ISEXP(mm);
            ef=in_ef+up_ef+down_ef; ISEXP(ef);
            e=in_e+up_e+down_e;     ISEXP(e);
            //  *_*_*
            fprintf(fpout,"%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t",\
                    in_um_e/in_um,in_um_ef/in_um,\
                    in_mm_e/in_mm,in_mm_ef/in_mm,\
                    up_um_e/up_um,up_um_ef/up_um,\
                    up_mm_e/up_mm,up_mm_ef/up_mm,\
                    down_um_e/down_um,down_um_ef/down_um,\
                    down_mm_e/down_mm,down_mm_ef/down_mm);
            fprintf(fpout,"%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t",\
                    in_um_ef/in_ef,in_mm_ef/in_ef,\
                    in_um_e/in_e,in_mm_e/in_e,\
                    up_um_ef/up_ef,up_mm_ef/up_ef,\
                    up_um_e/up_e,up_mm_e/up_e,\
                    down_um_ef/down_ef,down_mm_ef/down_ef,\
                    down_um_e/down_e,down_mm_e/down_e);
            fprintf(fpout,"%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t",\
                    in_um_ef/um_ef,up_um_ef/um_ef,down_um_ef/um_ef,\
                    in_um_e/um_e,up_um_e/um_e,down_um_e/um_e,\
                    in_mm_ef/mm_ef,up_mm_ef/mm_ef,down_mm_ef/mm_ef,\
                    in_mm_e/mm_e,up_mm_e/mm_e,down_mm_e/mm_e);
            // *_*
            fprintf(fpout,"%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t",\
                    in_um/in,in_mm/in,\
                    in_ef/in,in_e/in,\
                    up_um/up,up_mm/up,\
                    up_ef/up,up_e/up,\
                    down_um/down,down_mm/down,\
                    down_ef/down,down_e/down);
            fprintf(fpout,"%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t",\
                    in_um/um,up_um/um,down_um/um,\
                    um_ef/um,um_e/um,\
                    in_mm/mm,up_mm/mm,down_mm/mm,\
                    mm_ef/mm,mm_e/mm,\
                    in_ef/ef,up_ef/ef,down_ef/ef,\
                    um_ef/ef,mm_ef/ef,\
                    in_e/e,up_e/e,down_e/e,\
                    um_e/e,mm_e/e);
            //  *
            fprintf(fpout,"%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t",in/all,up/all,down/all,um/all,mm/all,ef/all,e/all);
        }

        fprintf(fpstats, "\n");
        fprintf(fpout, "\n");
        fflush(fpout);
        fflush(fpstats);
        fprintf(stderr, "OK\n");
    }
    fprintf(stderr, "[SVfeature_feature] Total %d variants collected.\n", indelList->indelCnt);
}

void Feature::getRegionStats(char *chrName, int start, int end, int &um_ef, int &um_e, int &mm_ef, int &mm_e)
{
    std::vector <BamAlignmentRecord *> &bamR = bam->bamRecord;

    int s = start - _diff_ , e = end + _diff_;
    for(int j = bam->getIdByPos(s, chrName);\
        bamR[j]->beginPos < e && j < bam->recordCnt; j++)
    {
        if(bamR[j]->flag & 0x4)
            continue;
        bool unique_flag = bam->isUniquelyMapped(bamR[j]);
        if(unique_flag)
        {
            if(bamR[j]->cigar[0].count == (uint)readLen && bamR[j]->cigar[0].operation == 'M')
            {
                um_ef++;
            }
            else
            {
                um_e++;
            }
        }
        else
        {
            if(bamR[j]->cigar[0].count == (uint)readLen && bamR[j]->cigar[0].operation == 'M')
            {
                mm_ef++;
            }
            else
            {
                mm_e++;
            }
        }
    }
}



void Feature::insertionFeature()
{

}
void Feature::indelFeature()
{

}










#endif // INCLUDE_FEATURE_H
