echo off
echo "Beginning compile"
echo "	---------------------------------"

wsl g++ -o simulation simulation.cpp forces_and_integrators.cpp sysvec_utils.cpp vector_utils.cpp

echo "	---------------------------------"

echo "Compile complete. Running: \n"
echo "	---------------------------------"

wsl ./simulation 1 ./results/pt_1a.dat 100
wsl ./simulation 2 ./results/pt_1b.dat 100
wsl ./simulation 3 ./results/pt_1c.dat 100
wsl ./simulation 4 ./results/pt_1d.dat 2600 1 1 0.01 5

echo "	---------------------------------"
pause