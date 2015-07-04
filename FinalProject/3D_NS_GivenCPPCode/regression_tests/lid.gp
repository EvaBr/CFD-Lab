set term x11 1
plot "<awk '$2>0.45&&$2<0.55{print $1, $4}' lid.dat | sort -g" w lp, \
     "ghia.ref.vy" w lp

set term x11 2
plot "<awk '$1>0.45&&$1<0.55{print $2, $3}' lid.dat | sort -g" w lp, \
     "ghia.ref.vx" w lp