#ifndef TOPO_H
#define TOPO_H

#include "struct.h"

void checkLenSum(SCFT_INFO* scftInfo);
void calcTotalNs(SCFT_INFO* scftInfo);
void findSolveScheme(SCFT_INFO *scftInfo);
void runSolveScheme(SCFT_INFO *scftInfo);

#endif