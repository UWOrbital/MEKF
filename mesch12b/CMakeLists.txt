cmake_minimum_required(VERSION 3.22)
project(mesch12b C)

set(CMAKE_C_STANDARD 99)

# Derived from sources in Makefile produced by ./configure
set(MESCH12B_SOURCES
    # List 1
    src/copy.c
    src/err.c
    src/matrixio.c
    src/memory.c
    src/vecop.c
    src/matop.c
    src/pxop.c
    src/submat.c
    src/init.c
    src/otherio.c
    src/machine.c
    src/matlab.c
    src/ivecop.c
    src/version.c
    src/meminfo.c
    src/memstat.c

    # List 2
    src/lufactor.c
    src/bkpfacto.c
    src/chfactor.c
    src/qrfactor.c
    src/solve.c
    src/hsehldr.c
    src/givens.c
    src/update.c
    src/norm.c
    src/hessen.c
    src/symmeig.c
    src/schur.c
    src/svd.c
    src/fft.c
    src/mfunc.c
    src/bdfactor.c

    # List 3
    src/sparse.c
    src/sprow.c
    src/sparseio.c
    src/spchfctr.c
    src/splufctr.c
    src/spbkp.c
    src/spswap.c
    src/iter0.c
    src/itersym.c
    src/iternsym.c

    # List 4
    src/zmachine.c
    src/zcopy.c
    src/zmatio.c
    src/zmemory.c
    src/zvecop.c
    src/zmatop.c
    src/znorm.c
    src/zfunc.c

    # List 5
    src/zlufctr.c
    src/zsolve.c
    src/zmatlab.c
    src/zhsehldr.c
    src/zqrfctr.c
    src/zgivens.c
    src/zhessen.c
    src/zschur.c
)
add_library(mesch12b STATIC ${MESCH12B_SOURCES})
target_include_directories(mesch12b PUBLIC include/)