#tell make to create runme
all: runme


#compiler... note that icpc is smarter than g++ et al... something that compiles with this makefile using icpc won't necessarily work with g++
CPP = icpc
#CPP = c++
#CPP = g++
#CPP = gcc


#optimization
OPTFLAGS = -O3 -xhost -ipo
#OPTFLAGS = -O
#OPTFLAGS = -O3
#OPTFLAGS = -g


#list of types of OBJS which have different dependencies or are dependencies for other OBJS
ALLOBJS = \
	main.o \
	initialize.o \
	initial_conditions.o \
	dashpot.o \
	evolve_mts.o \
	evolve_kts.o \
	infinite_bc.o \
	update_system.o \
	print_mt_data.o \
	print_kt_data.o \
	restart.o

#objects to be compiled, but are dependencies of main.o
SUBOBJS = \
	initialize.o \
	initial_conditions.o \
	dashpot.o \
	evolve_mts.o \
	evolve_kts.o \
	infinite_bc.o \
	update_system.o \
	print_mt_data.o \
	print_kt_data.o \
	restart.o


#list dependencies
$(ALLOBJS): \
	common.h params.h protos1.h concentration.h concentration_methods.h kinetochore.h

main.o: \
	$(SUBOBJS) \
	protos2.h

initialize.o: \
	initial_conditions.o

evolve_mts.o: \
	protos2.h

restart.o: \
	initialize.o


#compile each .cpp file
.cpp.o:
	$(CPP) $(OPTFLAGS) -c -o $@ $<
#    $@  means "the target"
#    $<  means "whatever the dependencies are"


#create executable
#EXECUTABLE = runme
#$(EXECUTABLE): $(ALLOBJS)
runme: $(ALLOBJS)
	$(CPP) $(OPTFLAGS) -o $@ $(ALLOBJS)

