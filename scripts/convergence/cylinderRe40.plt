reset;
set terminal pdf enhanced color font "Palatino, 16" size 20cm, 15cm;

set output "`echo ${CUIBM_DIR}`/orderOfConvergence/cylinderRe40/cylRe40Convergence.pdf";

set xlabel 'Mesh size'
set ylabel 'L2-norm of solution differences'
set log xy

f(x) = a*x+b
fit f(x) "`echo ${CUIBM_DIR}`/orderOfConvergence/cylinderRe40/errNorms.txt" u (log10($1)):(log10($2)) via a,b

plot [100:1000] [0.1:5] \
"`echo ${CUIBM_DIR}`/orderOfConvergence/cylinderRe40/errNorms.txt" u (0.5*$1):2 w lp lw 5 pt 13 ps 1.5 lc rgb '#4B5ED7' title 'L2-norm', \
100/x lt 2 lw 3 lc rgb 'grey' title 'First-order convergence'
