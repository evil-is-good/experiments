reset

set size 1.0, 1.0
set term png enhanced size 1024, 1024
unset key

set origin 0.05, 0.0
set size 0.95, 0.95
#set rotate by 45
#set xtics axis scale 0.5,0 nomirror rotate by -270  offset 0, 0.7

set xtics font "Monospace, 20" out
set ytics font "Monospace, 20"



set xlabel "t_{xx}" font "Monospace, 20"
set ylabel "X" font "Monospace, 20"
#set xtics 0,.2,1.0
#set xtics add ("0.25" 0.25)
#set grid xtics lw 5

set output "test.png"
plot "stress_line_x0.gpd" u 4:2 w p lw 2 #rotate=45 #

