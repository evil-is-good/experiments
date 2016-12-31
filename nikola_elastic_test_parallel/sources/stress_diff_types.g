reset

set size 1.0, 1.0
set term png enhanced size 1024, 1024
#unset key

set origin 0.05, 0.0
set size 0.95, 0.95

set xtics font "Monospace, 20" out
set ytics font "Monospace, 20"



set xlabel "Y" font "Monospace, 20"
set ylabel "t_{yy}" font "Monospace, 20"
set xtics ("0.3" 0.3)
#set xtics add ("0.25" 0.25)
set grid xtics lw 5

set output "t_yy_types.png"
plot "03/stress_line_x0.gpd" u 2:4 w l lw 2 t "E_1=1.0 E_2=1.0 {/Symbol n_1}=0.2 {/Symbol n_2}=0.1", "03_E1=10/stress_line_x0.gpd" u 2:4 w l lw 2 t "E_1=1.0 E_2=10.0 {/Symbol n_1}=0.2 {/Symbol n_2}=0.1", "03_E0=10/stress_line_x0.gpd" u 2:4 w l lw 2 t "E_1=10.0 E_2=1.0 {/Symbol n_1}=0.2 {/Symbol n_2}=0.1" 

set output "t_yy_types_inv.png"
plot "inv_03/stress_line_x0.gpd" u 2:4 w l lw 2 t "E_1=1.0 E_2=1.0 {/Symbol n_1}=0.1 {/Symbol n_2}=0.2", "inv_03_E1=10/stress_line_x0.gpd" u 2:4 w l lw 2 t "E_1=1.0 E_2=10.0 {/Symbol n_1}=0.1 {/Symbol n_2}=0.2", "inv_03_E0=10/stress_line_x0.gpd" u 2:4 w l lw 2 t "E_1=10.0 E_2=1.0 {/Symbol n_1}=0.1 {/Symbol n_2}=0.2"

set ylabel "t_{xy}" font "Monospace, 20"
set output "t_xy_types.png"
plot "03/stress_line_x0.gpd" u 2:3 w l lw 2 t "E_1=1.0 E_2=1.0 {/Symbol n_1}=0.2 {/Symbol n_2}=0.1", "03_E1=10/stress_line_x0.gpd" u 2:3 w l lw 2 t "E_1=1.0 E_2=10.0 {/Symbol n_1}=0.2 {/Symbol n_2}=0.1", "03_E0=10/stress_line_x0.gpd" u 2:3 w l lw 2 t "E_1=10.0 E_2=1.0 {/Symbol n_1}=0.2 {/Symbol n_2}=0.1" 

set output "t_xy_types.png"
plot "inv_03/stress_line_x0.gpd" u 2:3 w l lw 2 t "E_1=1.0 E_2=1.0 {/Symbol n_1}=0.1 {/Symbol n_2}=0.2", "inv_03_E1=10/stress_line_x0.gpd" u 2:3 w l lw 2 t "E_1=1.0 E_2=10.0 {/Symbol n_1}=0.1 {/Symbol n_2}=0.2", "inv_03_E0=10/stress_line_x0.gpd" u 2:3 w l lw 2 t "E_1=10.0 E_2=1.0 {/Symbol n_1}=0.1 {/Symbol n_2}=0.2"

