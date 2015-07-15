set term dumb

plot "poisuille.th" u 1:2 t "theory", \
     "poisuille.sim" u 1:2 w p ps 2 pt 6 t "simulation"

call "saver.gp" "poisuille"