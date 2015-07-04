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

# Size of the domain and the final times
L=1.0
time=8.0
res=$(awk -v L=${L} -v Nx=${Nx} 'BEGIN {print L/Nx}')

m4 -D m_res=${res} -D m_Nx=${Nx} -D m_time=${time} lid.m4.xml > lid.xml

../main lid.xml
pvbatch lid.py

awk -v FS="," 'NR>1{print $5, $6, $2, $3} 1' lid.csv > lid.dat

