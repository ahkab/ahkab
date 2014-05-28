set terminal png
set output 'nmos_test-vgs_sweep-vgd_-2.png'
set title "VGS Sweep, VGD -2V"
plot './nmos_test-vgs_sweep-vgd_-2.txt' using 1:4 title 'gm@vgd=-2' with linespoints, \
'./nmos_test-vgs_sweep-vgd_-2.txt' using 1:3 title 'I@vgd=-2' w linespoints

set output 'nmos_test-vgs_sweep-vgd_-1.5.png'
set title "VGS Sweep, VGD -1.5V"
plot './nmos_test-vgs_sweep-vgd_-1.5.txt' using 1:4 title 'gm@vgd=-1.5' with linespoints, \
'./nmos_test-vgs_sweep-vgd_-1.5.txt' using 1:3 title 'I@vgd=-1.5' w linespoints

set output 'nmos_test-vgs_sweep-vgd_-1.png'
set title "VGS Sweep, VGD -1V"
plot './nmos_test-vgs_sweep-vgd_-1.txt' using 1:4 title 'gm@vgd=-1' with linespoints, \
'./nmos_test-vgs_sweep-vgd_-1.txt' using 1:3 title 'I@vgd=-1' w linespoints

set output 'nmos_test-vgs_sweep-vgd_-0.5.png'
set title "VGS Sweep, VGD -0.5V"
plot './nmos_test-vgs_sweep-vgd_-0.5.txt' using 1:4 title 'gm@vgd=-0.5' with linespoints, \
'./nmos_test-vgs_sweep-vgd_-0.5.txt' using 1:3 title 'I@vgd=-0.5' w linespoints

set output 'nmos_test-vgs_sweep-vgd_-0.2.png'
set title "VGS Sweep, VGD -0.2V"
plot './nmos_test-vgs_sweep-vgd_-0.2.txt' using 1:4 title 'gm@vgd=-0.2' with linespoints, \
'./nmos_test-vgs_sweep-vgd_-0.2.txt' using 1:3 title 'I@vgd=-0.2' w linespoints

set output 'nmos_test-vgs_sweep-vgd_0.png'
set title "VGS Sweep, VGD 0V"
plot './nmos_test-vgs_sweep-vgd_0.txt' using 1:4 title 'gm@vgd=0' with linespoints, \
'./nmos_test-vgs_sweep-vgd_0.txt' using 1:3 title 'I@vgd=0' w linespoints

set output 'nmos_test-vgs_sweep-vgd_2.png'
set title "VGS Sweep, VGD 2V"
plot './nmos_test-vgs_sweep-vgd_2.txt' using 1:4 title 'gm@vgd=2' with linespoints, \
'./nmos_test-vgs_sweep-vgd_2.txt' using 1:3 title 'I@vgd=2' w linespoints

set output 'nmos_test-vgs_sweep-vgd_1.5.png'
set title "VGS Sweep, VGD 1.5V"
plot './nmos_test-vgs_sweep-vgd_1.5.txt' using 1:4 title 'gm@vgd=1.5' with linespoints, \
'./nmos_test-vgs_sweep-vgd_1.5.txt' using 1:3 title 'I@vgd=1.5' w linespoints

set output 'nmos_test-vgs_sweep-vgd_1.png'
set title "VGS Sweep, VGD 1V"
plot './nmos_test-vgs_sweep-vgd_1.txt' using 1:4 title 'gm@vgd=1' with linespoints, \
'./nmos_test-vgs_sweep-vgd_1.txt' using 1:3 title 'I@vgd=1' w linespoints

set output 'nmos_test-vgs_sweep-vgd_0.5.png'
set title "VGS Sweep, VGD 0.5V"
plot './nmos_test-vgs_sweep-vgd_0.5.txt' using 1:4 title 'gm@vgd=0.5' with linespoints, \
'./nmos_test-vgs_sweep-vgd_0.5.txt' using 1:3 title 'I@vgd=0.5' w linespoints

set output 'nmos_test-vgs_sweep-vgd_0.2.png'
set title "VGS Sweep, VGD 0.2V"
plot './nmos_test-vgs_sweep-vgd_0.2.txt' using 1:4 title 'gm@vgd=0.2' with linespoints, \
'./nmos_test-vgs_sweep-vgd_0.2.txt' using 1:3 title 'I@vgd=0.2' w linespoints

set output 'nmos_test-vgd_sweep-vgs_-2.png'
set title "vgd Sweep, vgs -2V"
plot './nmos_test-vgd_sweep-vgs_-2.txt' using 1:4 title 'go@vgs=-2' with linespoints, \
'./nmos_test-vgd_sweep-vgs_-2.txt' using 1:3 title 'I@vgs=-2' w linespoints

set output 'nmos_test-vgd_sweep-vgs_-1.5.png'
set title "vgd Sweep, vgs -1.5V"
plot './nmos_test-vgd_sweep-vgs_-1.5.txt' using 1:4 title 'go@vgs=-1.5' with linespoints, \
'./nmos_test-vgd_sweep-vgs_-1.5.txt' using 1:3 title 'I@vgs=-1.5' w linespoints

set output 'nmos_test-vgd_sweep-vgs_-1.png'
set title "vgd Sweep, vgs -1V"
plot './nmos_test-vgd_sweep-vgs_-1.txt' using 1:4 title 'go@vgs=-1' with linespoints, \
'./nmos_test-vgd_sweep-vgs_-1.txt' using 1:3 title 'I@vgs=-1' w linespoints

set output 'nmos_test-vgd_sweep-vgs_-0.5.png'
set title "vgd Sweep, vgs -0.5V"
plot './nmos_test-vgd_sweep-vgs_-0.5.txt' using 1:4 title 'go@vgs=-0.5' with linespoints, \
'./nmos_test-vgd_sweep-vgs_-0.5.txt' using 1:3 title 'I@vgs=-0.5' w linespoints

set output 'nmos_test-vgd_sweep-vgs_-0.2.png'
set title "vgd Sweep, vgs -0.2V"
plot './nmos_test-vgd_sweep-vgs_-0.2.txt' using 1:4 title 'go@vgs=-0.2' with linespoints, \
'./nmos_test-vgd_sweep-vgs_-0.2.txt' using 1:3 title 'I@vgs=-0.2' w linespoints

set output 'nmos_test-vgd_sweep-vgs_0.png'
set title "vgd Sweep, vgs 0V"
plot './nmos_test-vgd_sweep-vgs_0.txt' using 1:4 title 'go@vgs=0' with linespoints, \
'./nmos_test-vgd_sweep-vgs_0.txt' using 1:3 title 'I@vgs=0' w linespoints

set output 'nmos_test-vgd_sweep-vgs_2.png'
set title "vgd Sweep, vgs 2V"
plot './nmos_test-vgd_sweep-vgs_2.txt' using 1:4 title 'go@vgs=2' with linespoints, \
'./nmos_test-vgd_sweep-vgs_2.txt' using 1:3 title 'I@vgs=2' w linespoints

set output 'nmos_test-vgd_sweep-vgs_1.5.png'
set title "vgd Sweep, vgs 1.5V"
plot './nmos_test-vgd_sweep-vgs_1.5.txt' using 1:4 title 'go@vgs=1.5' with linespoints, \
'./nmos_test-vgd_sweep-vgs_1.5.txt' using 1:3 title 'I@vgs=1.5' w linespoints

set output 'nmos_test-vgd_sweep-vgs_1.png'
set title "vgd Sweep, vgs 1V"
plot './nmos_test-vgd_sweep-vgs_1.txt' using 1:4 title 'go@vgs=1' with linespoints, \
'./nmos_test-vgd_sweep-vgs_1.txt' using 1:3 title 'I@vgs=1' w linespoints

set output 'nmos_test-vgd_sweep-vgs_0.5.png'
set title "vgd Sweep, vgs 0.5V"
plot './nmos_test-vgd_sweep-vgs_0.5.txt' using 1:4 title 'go@vgs=0.5' with linespoints, \
'./nmos_test-vgd_sweep-vgs_0.5.txt' using 1:3 title 'I@vgs=0.5' w linespoints

set output 'nmos_test-vgd_sweep-vgs_0.2.png'
set title "vgd Sweep, vgs 0.2V"
plot './nmos_test-vgd_sweep-vgs_0.2.txt' using 1:4 title 'go@vgs=0.2' with linespoints, \
'./nmos_test-vgd_sweep-vgs_0.2.txt' using 1:3 title 'I@vgs=0.2' w linespoints

set output 'nmos_test-char.png'
set title 'NMOS CHAR'
plot './nmos_test-char.txt' using 1:2 t 'I@vgs=0.5' with linespoints, \
'./nmos_test-char.txt' using 1:3 t 'I@vgs=1' with linespoints, \
'./nmos_test-char.txt' using 1:4 t 'I@vgs=1.2' with linespoints, \
'./nmos_test-char.txt' using 1:5 t 'I@vgs=1.4' with linespoints, \
'./nmos_test-char.txt' using 1:6 t 'I@vgs=1.6' with linespoints, \
'./nmos_test-char.txt' using 1:7 t 'I@vgs=1.8' with linespoints, \
'./nmos_test-char.txt' using 1:8 t 'I@vgs=2' with linespoints

set output 'pmos_test-vgs_sweep-vgd_-2.png'
set title "VGS Sweep, VGD -2V"
plot './pmos_test-vgs_sweep-vgd_-2.txt' using 1:4 title 'gm@vgd=-2' with linespoints, \
'./pmos_test-vgs_sweep-vgd_-2.txt' using 1:3 title 'I@vgd=-2' w linespoints

set output 'pmos_test-vgs_sweep-vgd_-1.5.png'
set title "VGS Sweep, VGD -1.5V"
plot './pmos_test-vgs_sweep-vgd_-1.5.txt' using 1:4 title 'gm@vgd=-1.5' with linespoints, \
'./pmos_test-vgs_sweep-vgd_-1.5.txt' using 1:3 title 'I@vgd=-1.5' w linespoints

set output 'pmos_test-vgs_sweep-vgd_-1.png'
set title "VGS Sweep, VGD -1V"
plot './pmos_test-vgs_sweep-vgd_-1.txt' using 1:4 title 'gm@vgd=-1' with linespoints, \
'./pmos_test-vgs_sweep-vgd_-1.txt' using 1:3 title 'I@vgd=-1' w linespoints

set output 'pmos_test-vgs_sweep-vgd_-0.5.png'
set title "VGS Sweep, VGD -0.5V"
plot './pmos_test-vgs_sweep-vgd_-0.5.txt' using 1:4 title 'gm@vgd=-0.5' with linespoints, \
'./pmos_test-vgs_sweep-vgd_-0.5.txt' using 1:3 title 'I@vgd=-0.5' w linespoints

set output 'pmos_test-vgs_sweep-vgd_-0.2.png'
set title "VGS Sweep, VGD -0.2V"
plot './pmos_test-vgs_sweep-vgd_-0.2.txt' using 1:4 title 'gm@vgd=-0.2' with linespoints, \
'./pmos_test-vgs_sweep-vgd_-0.2.txt' using 1:3 title 'I@vgd=-0.2' w linespoints

set output 'pmos_test-vgs_sweep-vgd_0.png'
set title "VGS Sweep, VGD 0V"
plot './pmos_test-vgs_sweep-vgd_0.txt' using 1:4 title 'gm@vgd=0' with linespoints, \
'./pmos_test-vgs_sweep-vgd_0.txt' using 1:3 title 'I@vgd=0' w linespoints

set output 'pmos_test-vgs_sweep-vgd_2.png'
set title "VGS Sweep, VGD 2V"
plot './pmos_test-vgs_sweep-vgd_2.txt' using 1:4 title 'gm@vgd=2' with linespoints, \
'./pmos_test-vgs_sweep-vgd_2.txt' using 1:3 title 'I@vgd=2' w linespoints

set output 'pmos_test-vgs_sweep-vgd_1.5.png'
set title "VGS Sweep, VGD 1.5V"
plot './pmos_test-vgs_sweep-vgd_1.5.txt' using 1:4 title 'gm@vgd=1.5' with linespoints, \
'./pmos_test-vgs_sweep-vgd_1.5.txt' using 1:3 title 'I@vgd=1.5' w linespoints

set output 'pmos_test-vgs_sweep-vgd_1.png'
set title "VGS Sweep, VGD 1V"
plot './pmos_test-vgs_sweep-vgd_1.txt' using 1:4 title 'gm@vgd=1' with linespoints, \
'./pmos_test-vgs_sweep-vgd_1.txt' using 1:3 title 'I@vgd=1' w linespoints

set output 'pmos_test-vgs_sweep-vgd_0.5.png'
set title "VGS Sweep, VGD 0.5V"
plot './pmos_test-vgs_sweep-vgd_0.5.txt' using 1:4 title 'gm@vgd=0.5' with linespoints, \
'./pmos_test-vgs_sweep-vgd_0.5.txt' using 1:3 title 'I@vgd=0.5' w linespoints

set output 'pmos_test-vgs_sweep-vgd_0.2.png'
set title "VGS Sweep, VGD 0.2V"
plot './pmos_test-vgs_sweep-vgd_0.2.txt' using 1:4 title 'gm@vgd=0.2' with linespoints, \
'./pmos_test-vgs_sweep-vgd_0.2.txt' using 1:3 title 'I@vgd=0.2' w linespoints

set output 'pmos_test-vgd_sweep-vgs_-2.png'
set title "vgd Sweep, vgs -2V"
plot './pmos_test-vgd_sweep-vgs_-2.txt' using 1:4 title 'go@vgs=-2' with linespoints, \
'./pmos_test-vgd_sweep-vgs_-2.txt' using 1:3 title 'I@vgs=-2' w linespoints

set output 'pmos_test-vgd_sweep-vgs_-1.5.png'
set title "vgd Sweep, vgs -1.5V"
plot './pmos_test-vgd_sweep-vgs_-1.5.txt' using 1:4 title 'go@vgs=-1.5' with linespoints, \
'./pmos_test-vgd_sweep-vgs_-1.5.txt' using 1:3 title 'I@vgs=-1.5' w linespoints

set output 'pmos_test-vgd_sweep-vgs_-1.png'
set title "vgd Sweep, vgs -1V"
plot './pmos_test-vgd_sweep-vgs_-1.txt' using 1:4 title 'go@vgs=-1' with linespoints, \
'./pmos_test-vgd_sweep-vgs_-1.txt' using 1:3 title 'I@vgs=-1' w linespoints

set output 'pmos_test-vgd_sweep-vgs_-0.5.png'
set title "vgd Sweep, vgs -0.5V"
plot './pmos_test-vgd_sweep-vgs_-0.5.txt' using 1:4 title 'go@vgs=-0.5' with linespoints, \
'./pmos_test-vgd_sweep-vgs_-0.5.txt' using 1:3 title 'I@vgs=-0.5' w linespoints

set output 'pmos_test-vgd_sweep-vgs_-0.2.png'
set title "vgd Sweep, vgs -0.2V"
plot './pmos_test-vgd_sweep-vgs_-0.2.txt' using 1:4 title 'go@vgs=-0.2' with linespoints, \
'./pmos_test-vgd_sweep-vgs_-0.2.txt' using 1:3 title 'I@vgs=-0.2' w linespoints

set output 'pmos_test-vgd_sweep-vgs_0.png'
set title "vgd Sweep, vgs 0V"
plot './pmos_test-vgd_sweep-vgs_0.txt' using 1:4 title 'go@vgs=0' with linespoints, \
'./pmos_test-vgd_sweep-vgs_0.txt' using 1:3 title 'I@vgs=0' w linespoints

set output 'pmos_test-vgd_sweep-vgs_2.png'
set title "vgd Sweep, vgs 2V"
plot './pmos_test-vgd_sweep-vgs_2.txt' using 1:4 title 'go@vgs=2' with linespoints, \
'./pmos_test-vgd_sweep-vgs_2.txt' using 1:3 title 'I@vgs=2' w linespoints

set output 'pmos_test-vgd_sweep-vgs_1.5.png'
set title "vgd Sweep, vgs 1.5V"
plot './pmos_test-vgd_sweep-vgs_1.5.txt' using 1:4 title 'go@vgs=1.5' with linespoints, \
'./pmos_test-vgd_sweep-vgs_1.5.txt' using 1:3 title 'I@vgs=1.5' w linespoints

set output 'pmos_test-vgd_sweep-vgs_1.png'
set title "vgd Sweep, vgs 1V"
plot './pmos_test-vgd_sweep-vgs_1.txt' using 1:4 title 'go@vgs=1' with linespoints, \
'./pmos_test-vgd_sweep-vgs_1.txt' using 1:3 title 'I@vgs=1' w linespoints

set output 'pmos_test-vgd_sweep-vgs_0.5.png'
set title "vgd Sweep, vgs 0.5V"
plot './pmos_test-vgd_sweep-vgs_0.5.txt' using 1:4 title 'go@vgs=0.5' with linespoints, \
'./pmos_test-vgd_sweep-vgs_0.5.txt' using 1:3 title 'I@vgs=0.5' w linespoints

set output 'pmos_test-vgd_sweep-vgs_0.2.png'
set title "vgd Sweep, vgs 0.2V"
plot './pmos_test-vgd_sweep-vgs_0.2.txt' using 1:4 title 'go@vgs=0.2' with linespoints, \
'./pmos_test-vgd_sweep-vgs_0.2.txt' using 1:3 title 'I@vgs=0.2' w linespoints

set output 'pmos_test-char.png'
set title 'PMOS CHAR'
plot './pmos_test-char.txt' using 1:2 t 'I@vsg=0.5' with linespoints, \
'./pmos_test-char.txt' using 1:3 t 'I@vsg=1' with linespoints, \
'./pmos_test-char.txt' using 1:4 t 'I@vsg=1.2' with linespoints, \
'./pmos_test-char.txt' using 1:5 t 'I@vsg=1.4' with linespoints, \
'./pmos_test-char.txt' using 1:6 t 'I@vsg=1.6' with linespoints, \
'./pmos_test-char.txt' using 1:7 t 'I@vsg=1.8' with linespoints, \
'./pmos_test-char.txt' using 1:8 t 'I@vsg=2' with linespoints

