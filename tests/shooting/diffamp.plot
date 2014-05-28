set terminal png
set output 'diffamp_vout.png'
set title "Amplificatore differenziale: vout"
plot "diffamp_tran.tran" using 1:($6-$4) title 'TRAN: V2-V1' with linespoints,\
"diffamp_bruteforce.shooting" using 1:($6-$4) title 'BF: V2-V1' with linespoints,\
"diffamp_shooting.shooting" using 1:($6-$4) title 'shooting: V2-V1' with linespoints

set output 'diffamp_vin.png'
set title "Amplificatore differenziale: vin"
plot "diffamp_tran.tran" using 1:($9-$10) title 'TRAN: V30-V31' with linespoints,\
"diffamp_bruteforce.shooting" using 1:($9-$10) title 'BF: V30-V31' with linespoints,\
"diffamp_shooting.shooting" using 1:($9-$10) title 'shooting: V30-V31' with linespoints

set output 'diffamp_v3.png'
set title "Amplificatore differenziale: v3"
plot "diffamp_tran.tran" using 1:3 title 'TRAN: V3' with linespoints,\
"diffamp_bruteforce.shooting" using 1:3 title 'BF: V3' with linespoints,\
"diffamp_shooting.shooting" using 1:3 title 'shooting: V3' with linespoints

