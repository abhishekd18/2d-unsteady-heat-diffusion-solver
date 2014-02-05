echo "set term post color solid enh
set output \"./verification_par.ps\"
set title 'Temperature distribution along diagonal line'
set xlabel 'diagonal line'
set ylabel 'Temperatture (K)'
set grid x y
set label 't = 0.1 s' at 1,600
set label 't = 0.01 s' at 0.4,500
set label 't = 0.001 s' at 0.15,400
plot \"./serial_0001.csv\" u 3:1 w lines lt 1 lw 2 lc rgb \"green\" title \"Serial\",\
\"./serial_001.csv\" u 3:1 w lines lt 1 lw 2 lc rgb \"green\" title \"\",\
\"./serial_01.csv\" u 3:1 w lines lt 1 lw 2 lc rgb \"green\" title \"\",\
\"./parallel_0001.csv\" u 3:1 w points pt 1 lw 2 lc rgb \"red\" title \"Parallel np=4\",\
\"./parallel_001.csv\" u 3:1 w points pt 1 lw 2 lc rgb \"red\" title \"\",\
\"./parallel_01.csv\" u 3:1 w points pt 1 lw 2 lc rgb \"red\" title \"\"">gnuplot.gp
gnuplot gnuplot.gp
