# example GATB-core makefile
GATB=$(shell pwd)/lib/gatb-core/gatb-core/
GATB_LIB=$(GATB)/build/lib
GATBFLAGS = -I$(GATB)/src  -I$(GATB)/build/include -I$(GATB)/thirdparty
CXXFLAGS =  -std=c++0x -O3 -fpermissive
CXXFLAGS += ${GATBFLAGS}
LDFLAGS=   -L$(GATB_LIB) -lgatbcore  -lpthread -lhdf5 -lz -std=c++0x -ldl -static

#Main_code
m_main_obj := ${patsubst %.cpp,%.o,${wildcard src/*.cpp}}
src/%.o: src/%.cpp 
	${CXX} -c ${CXXFLAGS} $< -o $@

#Read_Seqs
p_obj := ${patsubst %.cpp,%.o,${wildcard Utils/*.cpp}}
Utils/%.o: Utils/%.cpp Utils/*h
	${CXX} -c ${CXXFLAGS} ${LDFLAGS} $< -o $@

m_code_release := src/main.o

all: ${m_main_obj}
	${CXX} ${m_code_release} -o output.out ${LDFLAGS}

#Clean
clean:
	-rm ${m_main_obj}
