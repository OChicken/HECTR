set terminal epslatex standalone color 8
set output 'inverted-pendulum-mpc-control.tex'
set multiplot layout 3,1 margins 0.12,0.96,0.1,.96 spacing 0.03,0.03
set grid
unset xlabel
set mxtics 5

set format x " "
set ylabel "$u_{op}$"
plot "inverted-pendulum-mpc-control.txt" using 1:2 with linespoints lw 2 t "u"

set format x " "
set key bottom right
set ylabel "pos"
plot "inverted-pendulum-mpc-control.txt" using 1:3 with linespoints lw 2 t "$x$",\
     "inverted-pendulum-mpc-control.txt" using 1:5 with linespoints lw 2 t "$\\theta$"

set format x
set key bottom right
set xlabel "Horizon n"
set ylabel "vel"
plot "inverted-pendulum-mpc-control.txt" using 1:4 with linespoints lw 2 t "$\\dot{x}$",\
     "inverted-pendulum-mpc-control.txt" using 1:6 with linespoints lw 2 t "$\\dot{\\theta}$"

unset multiplot

set output

!pdflatex --interaction=batchmode inverted-pendulum-mpc-control.tex
!rm *.aux *-inc.eps *converted-to.pdf *.log
