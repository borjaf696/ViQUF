# example GATB-core makefile
GATB=$(shell pwd)/lib/gatb-core/gatb-core/
GATB_LIB=$(GATB)/build/lib

#Gatb+sdsl

GATBFLAGS = -I$(GATB)/src  -I$(GATB)/build/include -I$(GATB)/thirdparty -I ~/include
CXXFLAGS = -std=c++0x -O3 -fpermissive
CXXFLAGS += ${GATBFLAGS}
LDFLAGS= -L$(GATB_LIB) -lgatbcore  -lpthread -lhdf5 -lz -std=c++0x -ldl -static
SDSLFLAGS = -L ~/lib -lsdsl -ldivsufsort -ldivsufsort64

#Main_code
m_main_obj := ${patsubst %.cpp,%.o,${wildcard src/*.cpp}}
src/%.o: src/%.cpp 
	${CXX} -c ${CXXFLAGS} ${SDSLFLAGS} $< -o $@

#Read_Seqs
p_obj := ${patsubst %.cpp,%.o,${wildcard Utils/*.cpp}}
Utils/%.o: Utils/%.cpp Utils/*h
	${CXX} -c ${CXXFLAGS} ${LDFLAGS} $< -o $@

m_code_release := src/test-meta.o

all: ${m_main_obj}
	${CXX} ${m_code_release} -o output.out ${LDFLAGS} ${SDSLFLAGS}

#Clean
clean:
	-rm ${m_main_obj}
