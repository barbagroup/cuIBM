reset;
set terminal pdf enhanced color font "Palatino, 16" size 20cm, 15cm;

set output "`echo ${CUIBM_DIR}`/DirectForcingMethod/convergence.pdf";

set xlabel 'Mesh size'
set ylabel 'L2-norm of solution differences'
set log xy

f(x) = a*x+b
fit f(x) "`echo ${CUIBM_DIR}`/DirectForcingMethod/errNorms.txt" u (log10($1)):(log10($2)) via a,b

plot [50:1000] [0.001:1] \
"`echo ${CUIBM_DIR}`/DirectForcingMethod/errNorms.txt" u 1:2 w lp lw 5 pt 13 ps 1.5 lc rgb '#4B5ED7' title 'L2-norm', \
10/x lt 2 lw 3 lc rgb 'grey' title 'First-order convergence', \
1000/(x*x) lt 3 lw 3 lc rgb 'black' title 'Second-order convergence'