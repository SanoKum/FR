set grid

set xlabel "x"
set ylabel "u"

set yrange [-0.1:1.1]

plot 'initVal.dat' w lp lc rgb "red" lw 2 title "Initial"
replot 'endVal.dat' w l lc rgb "green" lw 2 title "Last" 



#plot for[i=10:50:10] sprintf('Val_%06d.dat',i) using 1:2 with lp title sprintf("data %06d",i)
#replot "initVal.dat" u 1:2 w lp

#do for[i=10:110:100]{
#    filename = sprintf('Val_%06d.dat',i)
#    #set key title filename
#    replot filename u 1:2 w l title filename
#}
 