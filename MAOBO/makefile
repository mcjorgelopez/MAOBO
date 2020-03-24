CCOMP=g++
CFLAGS=-std=c++11 -Wall -fopenmp -O3
LFLAGS=-fopenmp -lm -O3
EJECUTABLE=./output/ejecutable

all: $(EJECUTABLE)
$(EJECUTABLE): ./obj/geometry/geometry.o ./obj/boundary/boundary.o ./obj/mesh/mesh.o ./obj/octree/octree_driver.o ./obj/comunicator/comunicator.o ./obj/main.o
	$(CCOMP) $(LFLAGS) ./obj/geometry/geometry.o ./obj/boundary/boundary.o ./obj/mesh/mesh.o ./obj/octree/octree_driver.o ./obj/comunicator/comunicator.o ./obj/main.o -o $(EJECUTABLE)
./obj/geometry/geometry.o: ./src/geometry/geometry.cc
	$(CCOMP) $(CFLAGS) -c ./src/geometry/geometry.cc -o ./obj/geometry/geometry.o
./obj/boundary/boundary.o: ./src/boundary/boundary.cc
	$(CCOMP) $(CFLAGS) -c ./src/boundary/boundary.cc -o ./obj/boundary/boundary.o
./obj/mesh/mesh.o: ./src/mesh/mesh.cc
	$(CCOMP) $(CFLAGS) -c ./src/mesh/mesh.cc -o ./obj/mesh/mesh.o
./obj/octree/octree_driver.o: ./src/octree/octree_driver.cc
	$(CCOMP) $(CFLAGS) -c ./src/octree/octree_driver.cc -o ./obj/octree/octree_driver.o
./obj/comunicator/comunicator.o: ./src/comunicator/comunicator.cc
	$(CCOMP) $(CFLAGS) -c ./src/comunicator/comunicator.cc -o ./obj/comunicator/comunicator.o
./obj/main.o: ./src/main.cc
	$(CCOMP) $(CFLAGS) -c ./src/main.cc -o ./obj/main.o
clean:
	rm -f ./*.msh
	rm -f ./*.txt
	rm -f ./output/ejecutable
	rm -f ./obj/main.o
	rm -f ./obj/geometry/geometry.o
	rm -f ./obj/boundary/boundary.o
	rm -f ./obj/mesh/mesh.o
	rm -f ./obj/octree/octree_driver.o
	rm -f ./obj/comunicator/comunicator.o
run:
	./output/ejecutable ./input/SOKET.msh SOKET.msh 8 4
