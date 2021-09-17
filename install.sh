#!/bin/bash

. /u/local/Modules/default/init/modules.sh
module load xz

hts_lib=/u/project/pajukant/malvarez/programs/htslib-1.10.2/

export CPATH=${hts_lib}:$CPATH
export LD_LIBRARY_PATH=${hts_lib}:$LD_LIBRARY_PATH
export LIBRARY_PATH=${hts_lib}:$LIBRARY_PATH

make

