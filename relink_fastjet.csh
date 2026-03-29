#!/bin/tcsh

set OBJDIR=.sl73_x8664_gcc485/obj/StRoot/StPicoD0JetAnaMaker
set LIBDIR=.sl73_x8664_gcc485/lib
set FJ=/gpfs01/star/pwg/lomicond/Ondrej/Jets/Alma9FastJet/fastjet-install

g++ -g -m64 -shared -Wl,-Bdynamic \
  -o ${OBJDIR}/StPicoD0JetAnaMaker.so \
  ${OBJDIR}/StCuts.o \
  ${OBJDIR}/StFJWrapper.o \
  ${OBJDIR}/StHIOverlayAngularities.o \
  ${OBJDIR}/StJet.o \
  ${OBJDIR}/StKaonPion.o \
  ${OBJDIR}/StPicoD0JetAnaMaker.o \
  ${OBJDIR}/StPicoD0JetAnaMaker_Cint.o \
  ${FJ}/lib/libNsubjettiness.a \
  -L${FJ}/lib \
  -lfastjet -lfastjetplugins -lfastjettools -lsiscone -lsiscone_spherical

cp ${OBJDIR}/StPicoD0JetAnaMaker.so ${LIBDIR}/StPicoD0JetAnaMaker.so
cp ${OBJDIR}/StPicoD0JetAnaMaker.so ${LIBDIR}/libStPicoD0JetAnaMaker.so
