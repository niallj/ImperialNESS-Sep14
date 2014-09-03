#######################
# set record list
# 1 Bin
# 2 Coord
# 3 Ncount
# 4 density
# 5 vx
# 6 vy
# 7 vz
# 8 atmKE
# 9 atmPE
# 10 atmPxx
# 11 atmPyy
# 12 atmPzz
# 13 atmPxy
# 14 atmPxz
# 15 atmPyz
# 16 J_Kx
# 17 J_Ky
# 18 J_Kz
# 19 J_Ux
# 20 J_Uy
# 21 J_Uz
# 22 J_Cx
# 23 J_Cy
# 24 J_Cz
# 25 J_Qx
# 26 J_Qy
# 27 J_Qz


#######################
# graphic output
# set terminal x11 size 800,640 position 0,0 enhanced font 'Times-Roman,20' persist
# set terminal epslatex size 8.0cm,6.4cm color solid colortext rounded font 'Times-Roman,8' standalone

# set output 'fig_pres.tex'
filename = 'run.2/2.out.prof_fluid'
beg = 1
end = 100


#######################
reset

set format x '%.2f'; set xlabel 'z'
set format y '%.2g'; set ylabel 'Pxy(z)'
set key top right
# set nokey

# xmin = -4.0
# xmax =  4.0
# set xrange [xmin:xmax]; set xtics 2.0; set mxtics 10

# ymin =  0.0
# ymax =  4.0
# set yrange [ymin:ymax]; set ytics 1.0; set mytics 10

# plot sprintf("../pore_W5_1_F0_20/%s", filename) every ::beg::end u 2:13 w lp ps 0.1 lw 2.0 title 'Fe=0.2', \
#      sprintf("../pore_W5_1_F0_40/%s", filename) every ::beg::end u 2:13 w lp ps 0.1 lw 2.0 title 'Fe=0.4'

plot sprintf("../pore_W5_1_F0_20/%s", filename) every ::beg::end u 2:($13/0.2) w lp ps 0.1 lw 2.0 title 'Fe=0.2', \
     sprintf("../pore_W5_1_F0_40/%s", filename) every ::beg::end u 2:($13/0.4) w lp ps 0.1 lw 2.0 title 'Fe=0.4'

