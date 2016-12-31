reset

set size 1.0, 1.0
set term png enhanced size 1024, 1024
unset key

set origin 0.05, 0.0
set size 0.95, 0.95

set xtics font "Monospace, 20" out
set ytics font "Monospace, 20"



set xlabel "Y" font "Monospace, 20"
set ylabel "t_{yy}" font "Monospace, 20"
set xtics ("0.1" 0.1, "0.2" 0.2, "0.3" 0.3, "0.4" 0.4)
#set xtics add ("0.25" 0.25)
set grid xtics lw 5

set output "t_yy.png"
plot "01/stress_line_x0.gpd" u 2:4 w l lw 2, "02/stress_line_x0.gpd" u 2:4 w l lw 2, "03/stress_line_x0.gpd" u 2:4 w l lw 2, "04/stress_line_x0.gpd" u 2:4 w l lw 2

set ylabel "t_{xy}" font "Monospace, 20"
set output "t_xy.png"
plot "01/stress_line_x0.gpd" u 2:3 w l lw 2, "02/stress_line_x0.gpd" u 2:3 w l lw 2, "03/stress_line_x0.gpd" u 2:3 w l lw 2, "04/stress_line_x0.gpd" u 2:3 w l lw 2

#set key
#set xrange [0:3]
#set xtics 0,1.0,5.0
#unset grid

#set xlabel "X" font "Monospace, 20"
#set ylabel "t_{yy}" font "Monospace, 20"
#set output "t_yy_y.png"
#plot "stress_line_y05.gpd" u 1:4 w l lw 2 t "y=0.5", "stress_line_y025.gpd" u 1:4 w l lw 2 t "y=0.25", "stress_line_y0125.gpd" u 1:4 w l lw 2 t "y=0.125"

#set xlabel "X" font "Monospace, 20"
#set ylabel "t_{xy}" font "Monospace, 20"
#set output "t_xy_y.png"
#plot "stress_line_y05.gpd" u 1:3 w l lw 2 t "y=0.5", "stress_line_y025.gpd" u 1:3 w l lw 2 t "y=0.25", "stress_line_y0125.gpd" u 1:3 w l lw 2 t "y=0.125"

