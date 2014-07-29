mkdir -p ./lib
mkdir -p ./objs
( cd ./src; /usr/bin/gmake clean; /usr/bin/gmake )
rm -f ../lib/libfishpack.a ../objs/blktri.o ../objs/cblktri.o ../objs/cmgnbn.o ../objs/comf.o ../objs/fftpack.o ../objs/genbun.o ../objs/gnbnaux.o ../objs/hstcrt.o ../objs/hstcsp.o ../objs/hstcyl.o ../objs/hstplr.o ../objs/hstssp.o ../objs/hw3crt.o ../objs/hwscrt.o ../objs/hwscsp.o ../objs/hwscyl.o ../objs/hwsplr.o ../objs/hwsssp.o ../objs/pois3d.o ../objs/poistg.o ../objs/sepaux.o ../objs/sepeli.o ../objs/sepx4.o 
gfortran  -c blktri.f -o ../objs/blktri.o
gfortran  -c cblktri.f -o ../objs/cblktri.o
gfortran  -c cmgnbn.f -o ../objs/cmgnbn.o
gfortran  -c comf.f -o ../objs/comf.o
gfortran  -c fftpack.f -o ../objs/fftpack.o
gfortran  -c genbun.f -o ../objs/genbun.o
gfortran  -c gnbnaux.f -o ../objs/gnbnaux.o
gfortran  -c hstcrt.f -o ../objs/hstcrt.o
gfortran  -c hstcsp.f -o ../objs/hstcsp.o
gfortran  -c hstcyl.f -o ../objs/hstcyl.o
gfortran  -c hstplr.f -o ../objs/hstplr.o
gfortran  -c hstssp.f -o ../objs/hstssp.o
gfortran  -c hw3crt.f -o ../objs/hw3crt.o
gfortran  -c hwscrt.f -o ../objs/hwscrt.o
gfortran  -c hwscsp.f -o ../objs/hwscsp.o
gfortran  -c hwscyl.f -o ../objs/hwscyl.o
gfortran  -c hwsplr.f -o ../objs/hwsplr.o
gfortran  -c hwsssp.f -o ../objs/hwsssp.o
gfortran  -c pois3d.f -o ../objs/pois3d.o
gfortran  -c poistg.f -o ../objs/poistg.o
gfortran  -c sepaux.f -o ../objs/sepaux.o
gfortran  -c sepeli.f -o ../objs/sepeli.o
gfortran  -c sepx4.f -o ../objs/sepx4.o
/usr/bin/ar -rv ../lib/libfishpack.a ../objs/blktri.o ../objs/cblktri.o ../objs/cmgnbn.o ../objs/comf.o ../objs/fftpack.o ../objs/genbun.o ../objs/gnbnaux.o ../objs/hstcrt.o ../objs/hstcsp.o ../objs/hstcyl.o ../objs/hstplr.o ../objs/hstssp.o ../objs/hw3crt.o ../objs/hwscrt.o ../objs/hwscsp.o ../objs/hwscyl.o ../objs/hwsplr.o ../objs/hwsssp.o ../objs/pois3d.o ../objs/poistg.o ../objs/sepaux.o ../objs/sepeli.o ../objs/sepx4.o 
ar: creating archive ../lib/libfishpack.a
a - ../objs/blktri.o
a - ../objs/cblktri.o
a - ../objs/cmgnbn.o
a - ../objs/comf.o
a - ../objs/fftpack.o
a - ../objs/genbun.o
a - ../objs/gnbnaux.o
a - ../objs/hstcrt.o
a - ../objs/hstcsp.o
a - ../objs/hstcyl.o
a - ../objs/hstplr.o
a - ../objs/hstssp.o
a - ../objs/hw3crt.o
a - ../objs/hwscrt.o
a - ../objs/hwscsp.o
a - ../objs/hwscyl.o
a - ../objs/hwsplr.o
a - ../objs/hwsssp.o
a - ../objs/pois3d.o
a - ../objs/poistg.o
a - ../objs/sepaux.o
a - ../objs/sepeli.o
a - ../objs/sepx4.o
( cd ./test; /usr/bin/gmake clean; /usr/bin/gmake )
rm -f  tblktri.exe tcblktri.exe tcmgnbn.exe tgenbun.exe thstcrt.exe thstcsp.exe thstcyl.exe thstplr.exe thstssp.exe thw3crt.exe thwscrt.exe thwscsp.exe thwscyl.exe thwsplr.exe thwsssp.exe tpois3d.exe tpoistg.exe tsepeli.exe tsepx4.exe
rm -f tblktri.exe
gfortran  tblktri.f -o tblktri.exe -L../lib -l fishpack
./tblktri.exe
1                    SUBROUTINE BLKTRI EXAMPLE

          ** important for TBLKTRI example: to avoid    **
          ** large discretization error, compile it     **
          ** and FISHPACK routines in double precision  **


          THE OUTPUT FROM THE NCAR CONTROL DATA 7600 WAS

                                IERROR = 0
                  DISCRETIZATION ERROR = 1.64478E-05
            REQUIRED LENGTH OF W ARRAY = 823

          THE OUTPUT FROM YOUR COMPUTER IS

                                IERROR = 0
                  DISCRETIZATION ERROR = 0.12737E-01
            REQUIRED LENGTH OF W ARRAY = 823
rm -f tcblktri.exe
gfortran  tcblktri.f -o tcblktri.exe -L../lib -l fishpack
./tcblktri.exe
1                    SUBROUTINE CBLKTR EXAMPLE

          ** important for TCBLKTRI example: to avoid  **
          ** large discretization error, compile it    **
          ** and FISHPACK routines in double precision **


          THE OUTPUT FROM THE NCAR CONTROL DATA 7600 WAS

                                IERROR = 0
                  DISCRETIZATION ERROR = 1.64572E-05
            REQUIRED LENGTH OF W ARRAY = 1123

          THE OUTPUT FROM YOUR COMPUTER IS

                                IERROR = 0
                  DISCRETIZATION ERROR = 0.12737E-01
            REQUIRED LENGTH OF W ARRAY = 1123
rm -f tcmgnbn.exe
gfortran  tcmgnbn.f -o tcmgnbn.exe -L../lib -l fishpack
./tcmgnbn.exe
1                    SUBROUTINE CMGNBN EXAMPLE


          THE OUTPUT FROM THE NCAR CONTROL DATA 7600 WAS

                                IERROR = 0
                  DISCRETIZATION ERROR = 9.16200E-03
            REQUIRED LENGTH OF W ARRAY = 380

          THE OUTPUT FROM YOUR COMPUTER IS

                                IERROR = 0
                  DISCRETIZATION ERROR = 0.91855E-02
            REQUIRED LENGTH OF W ARRAY =380.
rm -f tgenbun.exe
gfortran  tgenbun.f -o tgenbun.exe -L../lib -l fishpack
./tgenbun.exe
1                    SUBROUTINE GENBUN EXAMPLE


          THE OUTPUT FROM THE NCAR CONTROL DATA 7600 WAS

                                IERROR = 0
                  DISCRETIZATION ERROR = 9.64063E-03
            REQUIRED LENGTH OF W ARRAY = 380

          THE OUTPUT FROM YOUR COMPUTER IS

                                IERROR = 0
                  DISCRETIZATION ERROR = 0.96569E-02
            REQUIRED LENGTH OF W ARRAY =380.
rm -f thstcrt.exe
gfortran  thstcrt.f -o thstcrt.exe -L../lib -l fishpack
./thstcrt.exe
1                    SUBROUTINE HSTCRT EXAMPLE


          THE OUTPUT FROM THE NCAR CONTROL DATA 7600 WAS

                                IERROR = 0
                  DISCRETIZATION ERROR = 1.26001E-03
            REQUIRED LENGTH OF W ARRAY = 884

          THE OUTPUT FROM YOUR COMPUTER IS

                                IERROR = 0
                  DISCRETIZATION ERROR = 0.12592E-02
            REQUIRED LENGTH OF W ARRAY =884.
rm -f thstcsp.exe
gfortran  thstcsp.f -o thstcsp.exe -L../lib -l fishpack
./thstcsp.exe
1                    SUBROUTINE HSTCSP EXAMPLE


          THE OUTPUT FROM THE NCAR CONTROL DATA 7600 WAS

                                IERROR = 0
                  DISCRETIZATION ERROR = 5.58432E-03
            REQUIRED LENGTH OF W ARRAY = 583

          THE OUTPUT FROM YOUR COMPUTER IS

                                IERROR = 0
                  DISCRETIZATION ERROR = 0.55845E-02
            REQUIRED LENGTH OF W ARRAY =583.
rm -f thstcyl.exe
gfortran  thstcyl.f -o thstcyl.exe -L../lib -l fishpack
./thstcyl.exe
1                    SUBROUTINE HSTCYL EXAMPLE


          THE OUTPUT FROM THE NCAR CONTROL DATA 7600 WAS

                                IERROR = 0
                                PERTRB =-4.43114E-04
                  DISCRETIZATION ERROR = 7.52796E-05
            REQUIRED LENGTH OF W ARRAY = 958

          THE OUTPUT FROM YOUR COMPUTER IS

                                IERROR = 0
                                PERTRB =-0.44321E-03
                  DISCRETIZATION ERROR = 0.73282E-04
            REQUIRED LENGTH OF W ARRAY =958.
rm -f thstplr.exe
gfortran  thstplr.f -o thstplr.exe -L../lib -l fishpack
./thstplr.exe
1                    SUBROUTINE HSTPLR EXAMPLE


          THE OUTPUT FROM THE NCAR CONTROL DATA 7600 WAS

                                IERROR = 0
                  DISCRETIZATION ERROR = 1.13038E-03
            REQUIRED LENGTH OF W ARRAY = 1042

          THE OUTPUT FROM YOUR COMPUTER IS

                                IERROR = 0
                  DISCRETIZATION ERROR = 0.11299E-02
            REQUIRED LENGTH OF W ARRAY =1042.
rm -f thstssp.exe
gfortran  thstssp.f -o thstssp.exe -L../lib -l fishpack
./thstssp.exe
1                    SUBROUTINE HSTSSP EXAMPLE


          THE OUTPUT FROM THE NCAR CONTROL DATA 7600 WAS

                                IERROR = 0
                                PERTRB = 6.35830E-04
                  DISCRETIZATION ERROR = 3.37523E-03
            REQUIRED LENGTH OF W ARRAY = 540

          THE OUTPUT FROM YOUR COMPUTER IS

                                IERROR = 0
                                PERTRB = 0.63592E-03
                  DISCRETIZATION ERROR = 0.33788E-02
            REQUIRED LENGTH OF W ARRAY =540.
rm -f thw3crt.exe
gfortran  thw3crt.f -o thw3crt.exe -L../lib -l fishpack
./thw3crt.exe
1                    SUBROUTINE HW3CRT EXAMPLE


          THE OUTPUT FROM THE NCAR CONTROL DATA 7600 WAS

                                IERROR = 0
                  DISCRETIZATION ERROR = 9.64802E-03

          THE OUTPUT FROM YOUR COMPUTER IS

                                IERROR = 0
                  DISCRETIZATION ERROR = 0.96481E-02
rm -f thwscrt.exe
gfortran  thwscrt.f -o thwscrt.exe -L../lib -l fishpack
./thwscrt.exe
1                    SUBROUTINE HWSCRT EXAMPLE


          THE OUTPUT FROM THE NCAR CONTROL DATA 7600 WAS

                                IERROR = 0
                  DISCRETIZATION ERROR = 5.36508E-04
            REQUIRED LENGTH OF W ARRAY = 880

          THE OUTPUT FROM YOUR COMPUTER IS

                                IERROR = 0
                  DISCRETIZATION ERROR = 0.51117E-03
            REQUIRED LENGTH OF W ARRAY =880.
rm -f thwscsp.exe
gfortran  thwscsp.f -o thwscsp.exe -L../lib -l fishpack
./thwscsp.exe
1                    SUBROUTINE HWSCSP EXAMPLE 1


          THE OUTPUT FROM THE NCAR CONTROL DATA 7600 WAS

                                IERROR = 0
                  DISCRETIZATION ERROR = 7.99842E-04
            REQUIRED LENGTH OF W ARRAY = 775

          THE OUTPUT FROM YOUR COMPUTER IS

                                IERROR = 0
                  DISCRETIZATION ERROR = 0.79924E-03
            REQUIRED LENGTH OF W ARRAY = 775
1                    SUBROUTINE HWSCSP EXAMPLE 2


          THE OUTPUT FROM THE NCAR CONTROL DATA 7600 WAS

                                IERROR = 0
                  DISCRETIZATION ERROR = 5.86824E-05
            REQUIRED LENGTH OF W ARRAY = 775

          THE OUTPUT FROM YOUR COMPUTER IS

                                IERROR = 0
                  DISCRETIZATION ERROR = 0.60886E-04
            REQUIRED LENGTH OF W ARRAY = 775
rm -f thwscyl.exe
gfortran  thwscyl.f -o thwscyl.exe -L../lib -l fishpack
./thwscyl.exe
1                    SUBROUTINE HWSCYL EXAMPLE


          THE OUTPUT FROM THE NCAR CONTROL DATA 7600 WAS

                                IERROR = 0
                                PERTRB = 2.26734E-04
                  DISCRETIZATION ERROR = 3.73672E-04
            REQUIRED LENGTH OF W ARRAY = 1118

          THE OUTPUT FROM YOUR COMPUTER IS

                                IERROR = 0
                                PERTRB = 0.22698E-03
                  DISCRETIZATION ERROR = 0.36252E-03
            REQUIRED LENGTH OF W ARRAY =1118.
rm -f thwsplr.exe
gfortran  thwsplr.f -o thwsplr.exe -L../lib -l fishpack
./thwsplr.exe
1                    SUBROUTINE HWSPLR EXAMPLE


          THE OUTPUT FROM THE NCAR CONTROL DATA 7600 WAS

                                IERROR = 0
                  DISCRETIZATION ERROR = 6.19134E-04
            REQUIRED LENGTH OF W ARRAY = 882

          THE OUTPUT FROM YOUR COMPUTER IS

                                IERROR = 0
                  DISCRETIZATION ERROR = 0.62260E-03
            REQUIRED LENGTH OF W ARRAY =882.
rm -f thwsssp.exe
gfortran  thwsssp.f -o thwsssp.exe -L../lib -l fishpack
./thwsssp.exe
1                    SUBROUTINE HWSSSP EXAMPLE


          THE OUTPUT FROM THE NCAR CONTROL DATA 7600 WAS

                                IERROR = 0
                  DISCRETIZATION ERROR = 3.38107E-03
            REQUIRED LENGTH OF W ARRAY = 600

          THE OUTPUT FROM YOUR COMPUTER IS

                                IERROR = 0
                  DISCRETIZATION ERROR = 0.33740E-02
            REQUIRED LENGTH OF W ARRAY = 600
rm -f tpois3d.exe
gfortran  tpois3d.f -o tpois3d.exe -L../lib -l fishpack
./tpois3d.exe
1                    SUBROUTINE POIS3D EXAMPLE


          THE OUTPUT FROM THE NCAR CONTROL DATA 7600 WAS

                                IERROR = 0
                  DISCRETIZATION ERROR = 2.93277E-02

          THE OUTPUT FROM YOUR COMPUTER IS

                                IERROR = 0
                  DISCRETIZATION ERROR = 0.29339E-01
rm -f tpoistg.exe
gfortran  tpoistg.f -o tpoistg.exe -L../lib -l fishpack
./tpoistg.exe
1                    SUBROUTINE POISTG EXAMPLE


          THE OUTPUT FROM THE NCAR CONTROL DATA 7600 WAS

                                IERROR = 0
                  DISCRETIZATION ERROR = 5.64171E-04
            REQUIRED LENGTH OF W ARRAY = 560

          THE OUTPUT FROM YOUR COMPUTER IS

                                IERROR = 0
                  DISCRETIZATION ERROR = 0.56316E-03
            REQUIRED LENGTH OF W ARRAY =560.
rm -f tsepeli.exe
gfortran  tsepeli.f -o tsepeli.exe -L../lib -l fishpack
./tsepeli.exe
1                    SUBROUTINE SEPELI EXAMPLE


                    THE OUTPUT FROM THE NCAR CONTROL DATA 7600 WAS

                    IERROR = 0
                    SECOND ORDER DISCRETIZATION ERROR = 9.78910E-05
                    FOURTH ORDER DISCRETIZATION ERROR = 1.47351E-06
                    REQUIRED LENGTH OF W ARRAY = 1118

                    THE OUTPUT FROM YOUR COMPUTER IS

                    IERROR = 0
                    SECOND ORDER DISCRETIZATION ERROR =  0.12445E-03
                    FOURTH ORDER DISCRETIZATION ERROR =  0.29206E-04
                    REQUIRED LENGTH OF W ARRAY = 1118
rm -f tsepx4.exe
gfortran  tsepx4.f -o tsepx4.exe -L../lib -l fishpack
./tsepx4.exe
1                    SUBROUTINE SEPX4  EXAMPLE


                    THE OUTPUT FROM THE NCAR CONTROL DATA 7600 WAS

                    IERROR = 0
                    SECOND ORDER DISCRETIZATION ERROR =  1.5985E-04 
                    FOURTH ORDER DISCRETIZATION ERROR =  1.85749E-06
                    REQUIRED LENGTH OF W ARRAY = 1024

                    THE OUTPUT FROM YOUR COMPUTER IS

                    IERROR = 0
                    SECOND ORDER DISCRETIZATION ERROR =  0.15593E-03
                    FOURTH ORDER DISCRETIZATION ERROR =  0.16212E-04
                    REQUIRED LENGTH OF W ARRAY =  1024
