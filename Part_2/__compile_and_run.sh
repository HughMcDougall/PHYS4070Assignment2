echo "Beginning compile"
echo "	---------------------------------"

g++ -o singlegrid singlegrid.cpp monte_carlo.cpp rand_utils.cpp vector_utils.cpp -fopenmp
g++ -o temp_sweep temp_sweep.cpp monte_carlo.cpp rand_utils.cpp vector_utils.cpp -fopenmp
g++ -o singlegrid_parallel singlegrid_parallel.cpp monte_carlo.cpp rand_utils.cpp vector_utils.cpp -fopenmp

echo "	---------------------------------"

echo "Compile complete. Running: \n"
echo "	---------------------------------"

./singlegrid

./temp_sweep  8 0.0 5.0 32 ./results/temp_sweep_8.dat
./temp_sweep 16 0.0 5.0 32  ./results/temp_sweep_16.dat
./temp_sweep  32 0.0 5.0 32 ./results/temp_sweep_32.dat
./temp_sweep  64 0.0 5.0 32 ./results/temp_sweep_64.dat

./temp_sweep  128 2.0 2.26918531421 32 ./results/temp_sweep_128_nearcrit.dat 100000 1000 1

./singlegrid_parallel


echo "	---------------------------------"
pause