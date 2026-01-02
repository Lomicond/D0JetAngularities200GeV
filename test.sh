# stejný include set jako cons (zkopíruj z logu)
INC="-IStRoot/StPicoD0JetAnaMaker -I. -IStRoot -I.sl73_gcc485/include \
-I/afs/rhic.bnl.gov/star/packages/SL22c \
-I/afs/rhic.bnl.gov/star/packages/SL22c/StRoot \
-I/afs/rhic.bnl.gov/star/packages/SL22c/.sl73_gcc485/include \
-I/afs/rhic.bnl.gov/star/ROOT/5.34.38/.sl73_gcc485/rootdeb/include"


# testuj po jednom
rootcint -f /tmp/testCint.cxx -c -DROOT_CINT -D__ROOT__ "-IStRoot/StPicoD0JetAnaMaker -I. -IStRoot -I.sl73_gcc485/include \
-I/afs/rhic.bnl.gov/star/packages/SL22c \
-I/afs/rhic.bnl.gov/star/packages/SL22c/StRoot \
-I/afs/rhic.bnl.gov/star/packages/SL22c/.sl73_gcc485/include \
-I/afs/rhic.bnl.gov/star/ROOT/5.34.38/.sl73_gcc485/rootdeb/include" StRoot/D0JetAngularities/StEmcPosition2.h LinkDef.h
rootcint -f /tmp/testCint.cxx -c -DROOT_CINT -D__ROOT__ $INC StHIOverlayAngularities.h LinkDef.h
rootcint -f /tmp/testCint.cxx -c -DROOT_CINT -D__ROOT__ $INC StJet.h LinkDef.h
rootcint -f /tmp/testCint.cxx -c -DROOT_CINT -D__ROOT__ $INC StKaonPion.h LinkDef.h
rootcint -f /tmp/testCint.cxx -c -DROOT_CINT -D__ROOT__ $INC StPicoD0JetAnaMaker.h LinkDef.h

