set term dumb

plot "couette.th" u 1:2 t "theory", \
     "couette.sim" u 1:2 w p ps 2 pt 6 t "simulation"

call "saver.gp" "couette"