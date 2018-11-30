set parametric
set pm3d depthorder hidden3d 1
set hidden3d
set style fill transparent solid 0.15
set palette rgbformulae 8, 9, 7
set style fill  transparent solid 0.30 border
set pm3d depthorder border linecolor rgb "#a0a0f0"  linewidth 0.5

R = 1.   # radius of sphere
set urange [-pi/2:pi/2]
set vrange [0:pi]
splot 0.5*cos(u)*cos(v),0.5*cos(u)*sin(v),0.5*sin(u) w surface lc rgb "yellow", 0.6*cos(u)*cos(v),0.6*cos(u)*sin(v),0.6*sin(u)  w surface lc rgb "red", 0.7*cos(u)*cos(v),0.7*cos(u)*sin(v),0.7*sin(u)  w surface lc rgb "blue", 0.8*cos(u)*cos(v),0.8*cos(u)*sin(v),0.8*sin(u)  w surface lc rgb "orange", "short.lst" u 1:2:3 w l , "-" w p 
0 0 0 
e