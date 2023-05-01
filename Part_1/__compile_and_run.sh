echo "Beginning compile"
echo "	---------------------------------"

g++ -o simulation simulation.cpp forces_and_integrators.cpp sysvec_utils.cpp vector_utils.cpp

echo "	---------------------------------"

echo "Compile complete. Running: \n"
echo "	---------------------------------"

./simulation 1 ./results/pt_1a.dat 100
./simulation 2 ./results/pt_1b.dat 100
./simulation 3 ./results/pt_1c.dat 100
./simulation 4 ./results/pt_1d.dat 2600 1 1 0.01 5

echo "	---------------------------------"
pause