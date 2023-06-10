set terminal epslatex standalone color 8
set output 'cstr-ode.tex'
set multiplot layout 2,1 margins 0.12,0.96,0.1,.96 spacing 0.03,0.03
set grid
unset xlabel

set format x " "
set key bottom right
set ylabel "$c_a$ (kmol/m\\textsuperscript{3})"
plot "cstr-ode.txt" using 1:2 with lines lw 2 t "ode45",\
     "cstr-ode.txt" using 1:4 with lines lw 2 t "ode15s"

set format x
set key bottom right
set xlabel "Time t"
set ylabel "$T$ (K)"
plot "cstr-ode.txt" using 1:3 with lines lw 2 t "ode45",\
     "cstr-ode.txt" using 1:5 with lines lw 2 t "ode15s"

unset multiplot

set output

!pdflatex --interaction=batchmode cstr-ode.tex
!rm *.aux *-inc.eps *converted-to.pdf *.log
