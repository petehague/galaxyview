MYCPP := g++ -O3 -std=c++14 -g -D SERIAL
#MYCPP := icpc -O3 -std=c++14 -g -D SERIAL

AGENTFILES := RainfallMCMC/chain.cpp RainfallMCMC/pick.cpp RainfallMCMC/options.cpp
AGENTOPTS := -shared -ldl -fPIC -D MODULAR $(AGENTFILES)

all: binfolder castor pollux vfcomp stat abg.so gihist.so
	echo Make complete

tr: tring tring.so
	echo Tilted ring code complete

abg.so: source/abg.cpp
	$(MYCPP) $(AGENTOPTS) source/abg.cpp -oabg.so

gihist.so: source/gihist.cpp
		$(MYCPP) $(AGENTOPTS) source/gihist.cpp -ogihist.so

pollux: bin/pollux.o bin/hsl.o bin/geometry.o bin/stat.o bin/tasks.o
	$(MYCPP) bin/pollux.o bin/hsl.o bin/geometry.o bin/stat.o bin/tasks.o -opollux

stat: source/stat.cpp
	$(MYCPP) -D STANDALONE source/stat.cpp -ostat

vfcomp: source/vfcomp.cpp bin/projection.o bin/stat.o
	$(MYCPP) -D STANDALONE source/vfcomp.cpp bin/projection.o bin/stat.o -ovfcomp

castor: bin/castor.o bin/baryons.o bin/projection.o bin/geometry.o bin/observe.o bin/hsl.o
	$(MYCPP) bin/castor.o bin/baryons.o bin/projection.o bin/geometry.o bin/observe.o bin/hsl.o -ocastor

tring: bin/tring.o bin/projection.o bin/vfcomp.o bin/massmodels.o bin/hsl.o bin/stat.o
	$(MYCPP) bin/tring.o bin/projection.o bin/vfcomp.o bin/massmodels.o bin/hsl.o bin/stat.o -otring

tring.so: source/tring.cpp source/massmodels.cpp $(AGENTFILES)
		$(MYCPP) $(AGENTOPTS) source/tring.cpp source/massmodels.cpp source/vfcomp.cpp source/stat.cpp source/projection.cpp -otring.so

bin/stat.o: source/stat.cpp
	$(MYCPP) -c source/stat.cpp -obin/stat.o

bin/tasks.o: source/tasks.cpp
	$(MYCPP) -c source/tasks.cpp -obin/tasks.o

bin/hsl.o: source/hsl.cpp
	$(MYCPP) -c source/hsl.cpp -obin/hsl.o

bin/massmodels.o: source/massmodels.cpp
	$(MYCPP) -c source/massmodels.cpp -obin/massmodels.o

bin/observe.o: source/observe.cpp
	$(MYCPP) -c source/observe.cpp -obin/observe.o

bin/geometry.o: source/geometry.cpp
	$(MYCPP) -c source/geometry.cpp -obin/geometry.o

bin/castor.o: source/castor.cpp
	$(MYCPP) -c source/castor.cpp -obin/castor.o

bin/baryons.o: source/baryons.cpp
	$(MYCPP) -c source/baryons.cpp -obin/baryons.o

bin/projection.o: source/projection.cpp
	$(MYCPP) -c source/projection.cpp -obin/projection.o

bin/tring.o: source/tring.cpp
	$(MYCPP) -c source/tring.cpp -obin/tring.o

bin/vfcomp.o: source/vfcomp.cpp
	$(MYCPP) -c source/vfcomp.cpp -obin/vfcomp.o

bin/pollux.o: source/pollux.cpp
	$(MYCPP) -c source/pollux.cpp -obin/pollux.o

binfolder:
	mkdir -p bin

clean:
	rm -rf bin/*.o
	rm -f castor
	rm -f pollux
	rm -f tring
	rm -f tring.so
	rm -f vfcomp
	rm -rf *.dSYM*
	rm -f abg.so
