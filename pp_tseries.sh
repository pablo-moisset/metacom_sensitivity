#!/bin/sh

#Calls MATLAB function pp_time_series(filename_in, filename_out, p)

echo $0
command="/nfs/MATLAB/R2018a/bin/matlab -nojvm -singleCompThread -r \"pp_time_series($1,$2,$3); quit\""
echo "$command"
/nfs/MATLAB/R2018a/bin/matlab -nojvm -singleCompThread -r "pp_time_series($1,$2,$3); quit"
