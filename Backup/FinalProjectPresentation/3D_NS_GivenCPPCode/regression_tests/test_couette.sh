#! /bin/bash
# Couette flow test case
# Usage: ./test_couette.sh 40
# Returns L1 norm of the difference with theoretical solution

set -e
set -u

# resolution is an input parameter
if [ $# -gt 0 ]; then
    Nx=$1
else
    Nx=20
fi

# Size of the domain and final times
L=2.0
time=2.0
res=$(awk -v L=${L} -v Nx=${Nx} 'BEGIN {print L/Nx}')

m4 -D m_res=${res} -D m_Nx=${Nx} -D m_time=${time} couette.m4.xml > couette.xml
../main couette.xml > couette.log

# Use paraview for post-processing
pvbatch couette.py
awk -v FS=',' 'NR>1{print $6}' couette.csv  | awk -v t=${time} -v L=${L} -v V0=1.0 -v nu=0.02 -f ../scripts/couette.awk > couette.th
awk -v FS=',' 'NR>1{print $6, $2}'  couette.csv > couette.sim

gnuplot couette.gp 2>&1 > couette.gnuplot.log
cp couette.png couette.res${res}.png


thold=1.9
L1=$(paste couette.sim couette.th | awk -v th=${thold} '$1<th{s+=($2-$4)^2; n++} END {print s/n}')

# output resolution and L1 norm
echo ${res} ${L1}