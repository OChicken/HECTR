set terminal epslatex standalone color 8
set output 'mpc-tracking.tex'
set multiplot layout 3,1 margins 0.12,0.8,0.1,.96 spacing 0.03,0.03
set grid
unset xlabel
set mxtics 5

set format x " "
unset key
set ylabel "$u_{op}$"
set yrange [-0.5:1]
set ytics -0.5,0.25,1
plot "mpc-tracking-5.txt" using 1:2 with linespoints lw 2 t "test5", \
     "mpc-tracking-6.txt" using 1:2 with linespoints lw 2 t "test6", \
     "mpc-tracking-7.txt" using 1:2 with linespoints lw 2 t "test7", \
     "mpc-tracking-8.txt" using 1:2 with linespoints lw 2 t "test8", \
     "mpc-tracking-9.txt" using 1:2 with linespoints lw 2 t "test9", \
     "mpc-tracking-11.txt" using 1:2 with linespoints lw 2 t "test11", \
     "mpc-tracking-12.txt" using 1:2 with linespoints lw 2 t "test12", \

set format x " "
unset key
set ylabel "$x_0$"
set yrange [-2:1.5]
set ytics -2,0.5,1.5
plot "mpc-tracking-5.txt" using 1:3 with linespoints lw 2 t "test5", \
     "mpc-tracking-6.txt" using 1:3 with linespoints lw 2 t "test6", \
     "mpc-tracking-7.txt" using 1:3 with linespoints lw 2 t "test7", \
     "mpc-tracking-8.txt" using 1:3 with linespoints lw 2 t "test8", \
     "mpc-tracking-9.txt" using 1:3 with linespoints lw 2 t "test9", \
     "mpc-tracking-11.txt" using 1:3 with linespoints lw 2 t "test11", \
     "mpc-tracking-12.txt" using 1:3 with linespoints lw 2 t "test12", \

set format x
set key at screen 1,0.96
set ylabel "$x_1$"
set xlabel "Horizon $n$"
set yrange [-1.25:1.25]
set ytics -1.25,0.5,1.25
plot "mpc-tracking-5.txt" using 1:4 with linespoints lw 2 t "test5", \
     "mpc-tracking-6.txt" using 1:4 with linespoints lw 2 t "test6", \
     "mpc-tracking-7.txt" using 1:4 with linespoints lw 2 t "test7", \
     "mpc-tracking-8.txt" using 1:4 with linespoints lw 2 t "test8", \
     "mpc-tracking-9.txt" using 1:4 with linespoints lw 2 t "test9", \
     "mpc-tracking-11.txt" using 1:4 with linespoints lw 2 t "test11", \
     "mpc-tracking-12.txt" using 1:4 with linespoints lw 2 t "test12", \

unset multiplot


set output

!pdflatex --interaction=batchmode mpc-tracking.tex
!rm *.aux *-inc.eps *converted-to.pdf *.log
