set term pdfcairo enhanced color size 5,5 font 'Times-Roman, 20'
# set term postscript landscape enhanced color 'Times-Roman' 20
 set output 'empty.pdf'
 set noxzeroaxis
# set tics out
 set noxtics
# set key bottom right
# show key
 set nokey
 set grid
 set xrange[    0:  1.]
# set xtics 0,0.5,1
# set mxtics 5
 set yrange[    0.: 5]
 set ytics 0.,1.,10.
 set mytics 2
 set ylabel ' Energy (eV)'
# set xtics ('L'0.00000,'{/Symbol G}'0.26401,'X'0.56887,'U'0.67665,'{/Symbol G}'1.00000)
# set xtics ('L'0.00000,'{/Symbol G}'0.23297,'X'0.50199,'K'0.71467,'{/Symbol G}'1.00000)
# set xtics ('U'0.,'L'0.11618,'{/Symbol G}'0.28047,'X'0.47019,'L'0.63448,'{/Symbol G}'0.79878,'U'1.)
# set xtics ('L'0.,'{/Symbol G}'0.1995,'X'0.42987,'U'0.51132,'{/Symbol G}/K'0.75566,'{/Symbol G}'1.0)
# set xtics ('X'0.0000,'U'0.095931,'L'0.262088,'{/Symbol G}'0.497069,'X'0.768403, ' W ' 0.904069,'K'1.0000)
 set xtics ('{/Symbol G}'0.,'X'0.2228,'W'0.3343,'L' 0.4919,'{/Symbol G}'0.6848,'K,U'0.9212,'X'1.0000)


# set style data  lines
# set style line 1 lt 1 lw 1
# set style line 2 lt 2 lw 1
 set style data points
 set style line 1 pt 7 ps 0.6
 set style line 2 pt 6 ps 0.6

 plot 'empty.out'us 1:5 w p pt 7 ps 0.6, 'empty_vl.out' w l lt 3





# w l --> with lines,  w p --> with points, w lp --> with linespoints
# lc --> linecolor, ls --> linestyle, dt --> dashtype
# pt --> point type, ps --> point size
# to get help: enter interactive mode with > gnuplot, then gnuplot> help 
# or search in google!


set style line 1  lt 1 lw 2 lc rgb "red"     pt 6 ps 0.6 
set style line 2  lt 2 lw 2 lc rgb "blue"    pt 7 ps 0.6 
set style line 3  lt 3 lw 2 lc rgb "black"   pt 5 ps 0.3
set style line 4  lt 4 lw 2 lc rgb "green"   pt 4 ps 0.3
set style line 5  lt 5 lw 2 lc rgb "cyan"    pt 3 ps 0.3
set style line 6  lt 8 lw 2 lc rgb "magenta" pt 2 ps 0.3
set style line 10 lt 7 lw 2 lc rgb "green"   pt 4 ps 0.3


#set style line 1  lt 1 lw 2 pt 6 ps 0.6 lc rgb "black" 
#set style line 2  lt 1 lw 2 pt 7 ps 0.6 lc rgb "black" 

 plot  'empty.out'us 1:5  w p ls 1, 'empty_vl.out' w l lt 3, \




























