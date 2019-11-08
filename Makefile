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
lib/graph/%.o: lib/graph/%.cpp graph/*.h lib/utils/*h
	${CXX} -c ${CXXFLAGS} $< -o $@

#Main_code
m_main_obj := ${patsubst %.cpp,%.o,${wildcard src/*.cpp}}
src/%.o: src/%.cpp 
	${CXX} -c ${CXXFLAGS} ${SDSLFLAGS} $< -o $@

#Read_Seqs
p_obj := ${patsubst %.cpp,%.o,${wildcard lib/utils/*.cpp}}
lib/utils/%.o: lib/utils/%.cpp lib/utils/*h
	${CXX} -c ${CXXFLAGS} ${LDFLAGS} $< -o $@

mCodeRelease := src/test-meta.o
mCodeGatb := src/main.o

all: ${m_main_obj}
	${CXX} ${mCodeRelease} -o bin/output.out ${LDFLAGS} ${SDSLFLAGS}

meta: ${m_main_obj} ${graphObj}
	${CXX} ${mCodeRelease} -o bin/meta.out ${LDFLAGS} ${SDSLFLAGS}

gatb :${mCodeGatb}
	${CXX} ${mCodeGatb} -o bin/basic.out ${LDFLAGS} ${SDSLFLAGS}


#Clean
clean:
	-rm ${graphObj}
	-rm ${m_main_obj}
