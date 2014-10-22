all: transport.so

%.f90: %.F90 Makefile
	cpp -x f95-cpp-input -E -o $@ $<

%.so: %.f90 Makefile
	f2py -c --opt=-O3 --f90flags=-Wall --fcompiler=gnu95 $< -m $*


lxplus:
	-rm -r /afs/cern.ch/eng/sl/lintrack/Python-2.5/lib/python2.5/site-packages/pyviews/
	cp -a . /afs/cern.ch/eng/sl/lintrack/Python-2.5/lib/python2.5/site-packages/pyviews
	(make -C /afs/cern.ch/eng/sl/lintrack/Python-2.5/lib/python2.5/site-packages/pyviews transport-lxplus)

clean:
	-rm *.so

transport-lxplus:
	f2py -c --opt=-O3 --f90flags=-Wall --fcompiler=gnu95 transport.f90 -m transport

