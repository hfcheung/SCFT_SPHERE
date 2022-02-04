#ifndef CALC_H
#define CALC_H

#include "struct.h"

void removeAvgFld(SCFT_INFO *scftInfo, double *wA, double *wB);
void calcFreeE(SCFT_INFO *scftInfo);
void updateFldsHist(SCFT_INFO *scftInfo, double *del, double *outs, double *wADiff, double *wBDiff,
                      double *wANew, double *wBNew);
void AndersonMixing(SCFT_INFO *scftInfo, int nRec, double *del, double *outs);
void calcNewFlds(SCFT_INFO *scftInfo, int itr);
void calcJointDistr(SCFT_INFO *scftInfo);

#endif // CALC_H