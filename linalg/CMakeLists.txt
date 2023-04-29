
# Support Linux from Ubuntu 20.04LTS onward, CentOS 7.x (with EPEL),
# macOS, MSVC 2019 (=Version 16)
cmake_minimum_required(VERSION 3.16)

# set up project
project(linalg VERSION 3.1 DESCRIPTION "Linear Algebra library" LANGUAGES CXX)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# hacks for MSVC to prevent lots of pointless warnings about "unsafe" functions
if(MSVC)
  add_compile_options(/wd4244)
  add_compile_options(/wd4267)
  add_compile_definitions(_CRT_SECURE_NO_WARNINGS)
endif()

add_library(linalg STATIC
  d_lmp_cnjg.cpp
  d_lmp_imag.cpp
  d_lmp_sign.cpp
  dcabs1.cpp
  dcopy.cpp
  dgemm.cpp
  dgemv.cpp
  disnan.cpp
  dlacpy.cpp
  dladiv.cpp
  dlae2.cpp
  dlaed0.cpp
  dlaed1.cpp
  dlaed2.cpp
  dlaed3.cpp
  dlaed4.cpp
  dlaed5.cpp
  dlaed6.cpp
  dlaed7.cpp
  dlaed8.cpp
  dlaed9.cpp
  dlaeda.cpp
  dlaev2.cpp
  dlaev2.cpp
  dlamc3.cpp
  dlamch.cpp
  dlamrg.cpp
  dlanst.cpp
  dlapy2.cpp
  dlapy3.cpp
  dlartg.cpp
  dlascl.cpp
  dlaset.cpp
  dlasr.cpp
  dlasrt.cpp
  dlassq.cpp
  dnrm2.cpp
  drot.cpp
  dscal.cpp
  dstedc.cpp
  dsteqr.cpp
  dsterf.cpp
  dswap.cpp
  dznrm2.cpp
  i_lmp_len.cpp
  i_lmp_nint.cpp
  idamax.cpp
  ieeeck.cpp
  ilaenv.cpp
  ilazlc.cpp
  ilazlr.cpp
  iparmq.cpp
  lsame.cpp
  pow_lmp_di.cpp
  pow_lmp_ii.cpp
  s_lmp_cat.cpp
  s_lmp_cmp.cpp
  s_lmp_copy.cpp
  xerbla.cpp
  z_lmp_abs.cpp
  zaxpy.cpp
  zcopy.cpp
  zdotc.cpp
  zdrot.cpp
  zdscal.cpp
  zgemm.cpp
  zgemv.cpp
  zgerc.cpp
  zheevd.cpp
  zhemv.cpp
  zher2.cpp
  zher2k.cpp
  zhetd2.cpp
  zhetrd.cpp
  zlacgv.cpp
  zlacpy.cpp
  zlacrm.cpp
  zladiv.cpp
  zlaed0.cpp
  zlaed7.cpp
  zlaed8.cpp
  zlanhe.cpp
  zlarf.cpp
  zlarfb.cpp
  zlarfg.cpp
  zlarft.cpp
  zlascl.cpp
  zlaset.cpp
  zlasr.cpp
  zlassq.cpp
  zlatrd.cpp
  zscal.cpp
  zstedc.cpp
  zsteqr.cpp
  zswap.cpp
  ztrmm.cpp
  ztrmv.cpp
  zunm2l.cpp
  zunm2r.cpp
  zunmql.cpp
  zunmqr.cpp
  zunmtr.cpp
)
target_include_directories(linalg PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
