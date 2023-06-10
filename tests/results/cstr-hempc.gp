set terminal epslatex standalone color 8

set output 'cstr-hempc-x.tex'
set multiplot layout 3,1 margins 0.12,0.96,0.1,.96 spacing 0.05,0.05
set grid
unset xlabel

set format x " "
set key bottom right
set ylabel "$c_a$ (kmol/m\\textsuperscript{3})"
set yrange [0.87:0.9]
plot "cstr-mpc.bin"   binary format="%int %double %double %double %double %double" using 1:2 with linespoints lw 2 t "plain",\
     "cstr-hempc.bin" binary format="%int %double %double %double %double %double" using 1:2 with linespoints lw 2 t "encrypted"

set format x " "
set key bottom right
set ylabel "$T$ (K)"
set yrange [320:326]
plot "cstr-mpc.bin"   binary format="%int %double %double %double %double %double" using 1:3 with linespoints lw 2 t "plain",\
     "cstr-hempc.bin" binary format="%int %double %double %double %double %double" using 1:3 with linespoints lw 2 t "encrypted"

set format x
set key bottom right
set xlabel "Time t"
set ylabel "$h$ (m)"
set yrange [0.64:0.80]
plot  "cstr-mpc.bin"   binary format="%int %double %double %double %double %double" using 1:4 with linespoints lw 2 t "plain",\
      "cstr-hempc.bin" binary format="%int %double %double %double %double %double" using 1:4 with linespoints lw 2 t "encrypted"

unset multiplot

set output

!pdflatex --interaction=batchmode cstr-hempc-x.tex
!rm *.aux *-inc.eps *converted-to.pdf *.log

set output 'cstr-hempc-u.tex'
set multiplot layout 2,1 margins 0.12,0.96,0.1,.96 spacing 0.05,0.05
set grid
unset xlabel

set format x " "
set key bottom right
set ylabel "$T_c$ (K)"
set yrange [298:301]
plot "cstr-mpc.bin"   binary format="%int %double %double %double %double %double" using 1:5 with steps lw 2 t "plain",\
     "cstr-hempc.bin" binary format="%int %double %double %double %double %double" using 1:5 with steps lw 2 t "encrypted"

set format x
set key bottom right
set xlabel "Time t"
set ylabel "$F$ (m\\textsuperscript{3}/min)"
set yrange [0.095:0.115]
plot "cstr-mpc.bin"   binary format="%int %double %double %double %double %double" using 1:6 with steps lw 2 t "plain",\
     "cstr-hempc.bin" binary format="%int %double %double %double %double %double" using 1:6 with steps lw 2 t "encrypted"

unset multiplot

set output

!pdflatex --interaction=batchmode cstr-hempc-u.tex
!rm *.aux *-inc.eps *converted-to.pdf *.log
