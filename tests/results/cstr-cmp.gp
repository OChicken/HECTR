set terminal epslatex standalone color 8

set output 'cstr-cmp-x.tex'
set multiplot layout 3,1 margins 0.15,0.96,0.1,.96 spacing 0.05,0.05
set grid
unset xlabel

set format x " "
set key bottom right
set ylabel "$||c-c_{\\textrm{he}}||$"
set logscale y
plot "cstr-cmp.bin" binary format="%int %double %double %double %double %double" using 1:2 with linespoints lw 2 t ""

set format x " "
set key bottom right
set ylabel "$||T-T_{\\textrm{he}}||$"
set logscale y
plot "cstr-cmp.bin" binary format="%int %double %double %double %double %double" using 1:3 with linespoints lw 2 t ""

set format x
set key bottom right
set xlabel "time (min)"
set ylabel "$||h-h_{\\textrm{he}}||$"
set logscale y
plot "cstr-cmp.bin" binary format="%int %double %double %double %double %double" using 1:4 with linespoints lw 2 t ""

unset multiplot

set output

!pdflatex --interaction=batchmode cstr-cmp-x.tex
!rm *.aux *-inc.eps *converted-to.pdf *.log

set output 'cstr-cmp-u.tex'
set multiplot layout 2,1 margins 0.15,0.96,0.1,.96 spacing 0.05,0.05
set grid
unset xlabel

set format x " "
set key bottom right
set ylabel "$||T_c-T_{c,\\textrm{he}}||$"
set logscale y
plot "cstr-cmp.bin" binary format="%int %double %double %double %double %double" using 1:5 with steps lw 2 t ""

set format x
set key bottom right
set xlabel "time (min)"
set ylabel "$||F-F_{\\textrm{he}}||$"
set logscale y
plot "cstr-cmp.bin" binary format="%int %double %double %double %double %double" using 1:6 with steps lw 2 t ""

unset multiplot

set output

!pdflatex --interaction=batchmode cstr-cmp-u.tex
!rm *.aux *-inc.eps *converted-to.pdf *.log
