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
set terminal epslatex size 8.0cm,6.4cm color solid colortext rounded font 'Times-Roman,8' standalone

set output 'fig_dens_fit.tex'
filename = '../pore_W10_2_F0_20/run.2/2.out.prof_fluid'
beg = 10
end = 100


#######################
reset

set format x '%.2f'; set xlabel 'z'
set format y '%.2g'; set ylabel 'rho(z)'
set key top center
# set nokey

# xmin = -4.0
# xmax =  4.0
# set xrange [xmin:xmax]; set xtics 2.0; set mxtics 10

# ymin =  0.0
# ymax =  4.0
# set yrange [ymin:ymax]; set ytics 1.0; set mytics 10


#######################
# fit density profile
FIT_LIMIT = 1.e-8
FIT_LAMBDA_FACTOR = 2.0
lx = 10.2
rho(x) = a0 +\
         a1  * cos(1.0*2.0*pi*(x-x0)/lx) + \
         a2  * cos(2.0*2.0*pi*(x-x0)/lx) + \
         a3  * cos(3.0*2.0*pi*(x-x0)/lx) + \
         a4  * cos(4.0*2.0*pi*(x-x0)/lx) + \
         a5  * cos(5.0*2.0*pi*(x-x0)/lx) + \
         a6  * cos(6.0*2.0*pi*(x-x0)/lx) + \
         a7  * cos(7.0*2.0*pi*(x-x0)/lx) + \
         a8  * cos(8.0*2.0*pi*(x-x0)/lx) + \
         a9  * cos(9.0*2.0*pi*(x-x0)/lx) + \
         a10 * cos(10.0*2.0*pi*(x-x0)/lx)

x0  = -1.0
a0  = 2.0
a1  = 1.0
a2  = 1.0
a3  = 1.0
a4  = 1.0
a5  = 1.0
a6  = 1.0
a7  = 1.0
a8  = 1.0
a9  = 1.0
a10 = 1.0

fit rho(x) filename every ::beg::end u 2:4 \
    via x0, a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10

print "density parameters"
print "x0  = ", x0
print "a0  = ", a0
print "a1  = ", a1
print "a2  = ", a2
print "a3  = ", a3
print "a4  = ", a4
print "a5  = ", a5
print "a6  = ", a6
print "a7  = ", a7
print "a8  = ", a8
print "a9  = ", a9
print "a10 = ", a10


#######################
plot filename every ::beg::end u 2:4 w lp ps 1.0 lw 3.0 lc rgb 'dark-gray' title 'rho', \
     rho(x) w l lw 1.0 lc rgb 'red' title 'fit'
