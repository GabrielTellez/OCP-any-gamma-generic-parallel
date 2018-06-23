OBJECTS=cresults.o partition.o textslice.o dispatch_partitions_binaryblock.o process_options.o ocp_functions.o
LIBRARIES=-ltbb -larprec -lboost_program_options -lgmpxx -lgmp -lmpfr
all:	$(OBJECTS) ocp
clean:
	rm -f $(OBJECTS) checksum-tbbparallel checksum-tbbparallel.o benchmarking benchmarking.o softdisk softdisk.o ocp.o softdisk_singledensity.o  softdisk_singledensity
checksum-tbbparallel: checksum-tbbparallel.o $(OBJECTS) 
	$(CXX) $(CXXFLAGS) checksum-tbbparallel.o $(OBJECTS) -o checksum-tbbparallel $(LIBRARIES) 
benchmarking: benchmarking.o $(OBJECTS) 
	$(CXX) $(CXXFLAGS) benchmarking.o $(OBJECTS) -o benchmarking $(LIBRARIES) 
debug: checksum-tbbparallel.o $(OBJECTS) softdisk_functions.h
	g++ checksum-tbbparallel.o $(OBJECTS) -o checksum-tbbparallel $(LIBRARIES)  -g 
checksum-tbbparallel-linux: checksum-tbbparallel.o $(OBJECTS) softdisk_functions.h
	g++ checksum-tbbparallel.o $(OBJECTS) -o checksum-tbbparallel-linux $(LIBRARIES) -lrt
softdisk: softdisk.o $(OBJECTS) 
	$(CXX) $(CXXFLAGS) $< $(OBJECTS) -o $@ $(LIBRARIES)
softdisk_singledensity: softdisk_singledensity.o $(OBJECTS) 
	$(CXX) $(CXXFLAGS) $< $(OBJECTS) -o $@ $(LIBRARIES)
ocp: ocp.o $(OBJECTS) 
	$(CXX) $(CXXFLAGS) $< $(OBJECTS) -o $@ $(LIBRARIES)
ocp2: ocp.o $(OBJECTS) 
	$(CXX) $(CXXFLAGS) $< $(OBJECTS) -o $@ $(LIBRARIES)
partition.o: partition.h partition.cpp partition_inst.hpp
dispatch_partitions_binaryblock.o: dispatch_partitions_binaryblock.h dispatch_partitions_binaryblock.cpp dispatch_partitions_binaryblock_inst.hpp
cresults.o: cresults.h cresults.cpp cresults_inst.hpp
checksum-tbbparallel.o: checksum-tbbparallel.cpp softdisk_functions.h
benchmarking.o: benchmarking.cpp softdisk_functions.h
softdisk.o: softdisk.cpp softdisk_functions.h
ocp.o: ocp_models.h ocp_models_base.h ocp_functions.h run_ocp_model.h ocp.cpp
ocp_functions.o: ocp_functions.h
