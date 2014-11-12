REP=HigherPrecision
SOURCEFILE2=ancientSelectionApaBigQ
SOURCEFILE3=ancientSelectionApaSmallQ

all: wrap_lap.pyf wrap_lap.so
	make -C $(REP) all; cp $(REP)/$(SOURCEFILE2) . ; cp $(REP)/$(SOURCEFILE3) .   

wrap_lap.pyf:
	f2py -m wrap_lap -h wrap_lap.pyf spteqr.f;

wrap_lap.so: wrap_lap.pyf
	f2py  -c wrap_lap.pyf  spteqr.f -latlas -llapack -lblas;

clean:
	rm wrap_lap.so wrap_lap.pyf; make -C $(REP) clean


