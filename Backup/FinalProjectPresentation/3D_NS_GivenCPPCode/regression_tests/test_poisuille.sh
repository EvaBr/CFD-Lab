#! /bin/bash
# Poisuille flow test case
# Usage: ./test_poisuille.sh 40
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

m4 -D m_res=${res} -D m_Nx=${Nx} -D m_time=${time} poisuille.m4.xml > poisuille.xml
../main poisuille.xml > poisuille.log

# Use paraview for post-processing
pvbatch poisuille.py
awk -v FS=',' 'NR>1{print $6}' poisuille.csv  | awk -v t=${time} -v L=${L} -v F=1.0 -v nu=0.02 -f ../scripts/poisuille.awk > poisuille.th
awk -v FS=',' 'NR>1{print $6, $2}'  poisuille.csv > poisuille.sim

gnuplot poisuille.gp 2>&1 > poisuille.gnuplot.log
cp poisuille.png poisuille.res${res}.png


thold=0.1
L1=$(paste poisuille.sim poisuille.th | awk -v L=${L} -v th=${thold} '$1>th&&$1<L-th{s+=($2-$4)^2; n++} END {print s/n}')

# output resolution and L1 norm
echo ${res} ${L1}