awk '{if ($1 == 0.0) print $0}' stress_y_1.gpd > stress_line_x0.gpd
awk '{if ($2 == 0.5) print $0}' stress_y_1.gpd > stress_line_y05.gpd
awk '{if ($2 == 0.1) print $0}' stress_y_1.gpd > stress_line_y01.gpd
awk '{if ($2 == 0.2) print $0}' stress_y_1.gpd > stress_line_y02.gpd
awk '{if ($2 == 0.3) print $0}' stress_y_1.gpd > stress_line_y03.gpd
awk '{if ($2 == 0.4) print $0}' stress_y_1.gpd > stress_line_y04.gpd
#awk '{if ($2 == 0.25) print $0}' stress_y_1.gpd > stress_line_y025.gpd
#awk '{if ($2 == 0.125) print $0}' stress_y_1.gpd > stress_line_y0125.gpd
