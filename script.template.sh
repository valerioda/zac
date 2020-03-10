#!/bin/sh

#source /nfs/gerda2/LNGSMiB/MJsetupMiB.sh
source /afs/lngs.infn.it/experiment/gerda/MJsetupDev.sh
#source  /nfs/gerda5/gerda-sw/bin/swmod-env-cpath.sh gerda@dp-v2.2.1 > /dev/null
#export GA_BASE_DIR=/nfs/gerda5/gerda-data/blind/v02.00/gen
#export MU_CAL=/nfs/gerda5/gerda-data/blind/v02.00/meta/config/_aux/geruncfg


!filtercommand!
!gelatiocommand!
!cmdanalysis!

chmod 666 *.*

