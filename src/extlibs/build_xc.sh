LIBXC=/home/efefer/mysoftwares/libxc-3.0.0/lib/libxc.so
INC=/home/efefer/mysoftwares/libxc-3.0.0/include

gcc -c -I$INC -fPIC -O3 calc_Vxc_VWN.c
gcc -c -I$INC -fPIC -O3 calc_epsxc_VWN.c

gcc -c -I$INC -fPIC -O3 calc_Vxc_PBE.c
gcc -c -I$INC -fPIC -O3 calc_epsxc_PBE.c

gcc -c -I$INC -fPIC -O3 calc_Vxc_PBE_spinpol.c
gcc -c -I$INC -fPIC -O3 calc_epsxc_PBE_spinpol.c

gcc -c -I$INC -fPIC -O3 calc_epsxc_VWN_spinpol.c
gcc -c -I$INC -fPIC -O3 calc_Vxc_VWN_spinpol.c

gcc -shared -o libxc_interface.so -Wl,--whole-archive \
calc_Vxc_VWN.o calc_epsxc_VWN.o calc_epsxc_VWN_spinpol.o calc_Vxc_VWN_spinpol.o \
calc_Vxc_PBE.o calc_epsxc_PBE.o calc_Vxc_PBE_spinpol.o calc_epsxc_PBE_spinpol.o \
$LIBXC -Wl,--no-whole-archive

