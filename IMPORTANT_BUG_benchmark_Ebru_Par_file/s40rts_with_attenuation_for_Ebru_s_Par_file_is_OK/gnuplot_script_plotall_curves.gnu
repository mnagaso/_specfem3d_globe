#set term x11
set term wxt

set xrange [500:3615]

plot "specfem3d_globe_ifort_mcmodel_medium_s40RTS/AAE.IU.MXE.sem.ascii" w l lc 1, "specfem3d_globe_ifort_s40RTS/AAE.IU.MXE.sem.ascii" w l lc 3
pause -1 "hit key..."
