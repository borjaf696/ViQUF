# example GATB-core makefile
GATB=$(shell pwd)/lib/gatb-core/gatb-core/
GATB_LIB=$(GATB)/build/lib

#Gatb+sdsl

GATBFLAGS = -I$(GATB)/src  -I$(GATB)/build/include -I$(GATB)/thirdparty -I ~/include
CXXFLAGS = -std=c++0x -O3 -fpermissive
CXXFLAGS += ${GATBFLAGS}
LDFLAGS= -L$(GATB_LIB) -lgatbcore  -lpthread -lhdf5 -lz -std=c++0x -ldl -static
SDSLFLAGS = -L ~/lib -lsdsl -ldivsufsort -ldivsufsort64

#Graph
graphObj := ${patsubst %.cpp,%.o,${wildcard lib/graph/*.cpp}}
lib/graph/%.o: lib/graph/%.cpp lib/graph/*.h
	${CXX} -c ${CXXFLAGS} $< -o $@

#Main_code
m_main_obj := ${patsubst %.cpp,%.o,${wildcard src/*.cpp}}
src/%.o: src/%.cpp 
	${CXX} -c ${CXXFLAGS} ${SDSLFLAGS} $< -o $@

# Extra_code
extra_obj := ${patsubst %.cpp,%.o,${wildcard src/*.cpp}}
lib/extra/%.o: lib/extra/%.cpp lib/extra/%.h
	${CXX} -c ${CXXFLAGS} $< -o $@

mCodeRelease := src/test-meta.o
mCodeGatb := src/main.o

all: ${m_main_obj} ${graphObj}
	${CXX} ${graphObj} ${mCodeRelease} -o bin/output.out ${LDFLAGS} ${SDSLFLAGS}

meta: ${m_main_obj} ${graphObj}
	${CXX} ${graphObj} ${mCodeRelease} -o bin/meta.out ${LDFLAGS} ${SDSLFLAGS}

gatb :${mCodeGatb}
	${CXX} ${mCodeGatb} -o bin/basic.out ${LDFLAGS} ${SDSLFLAGS}


#Clean
clean:
	-rm ${graphObj}
	-rm ${m_main_obj}
	-rm ${extra_obj}
