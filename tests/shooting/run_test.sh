time ~/bin/ahkab --t-max-nr 100 -o diffamp_tran diffamp_tran.spc
time ~/bin/ahkab --t-max-nr 100 -o diffamp_shooting diffamp_shooting.spc
time ~/bin/ahkab --t-max-nr 100 -o diffamp_bruteforce diffamp_bruteforce.spc
gnuplot diffamp.plot
