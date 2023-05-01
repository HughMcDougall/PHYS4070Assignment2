echo off
echo "Beginning compile"
echo "	---------------------------------"

:: wsl g++ -o rand_test rand_test.cpp monte_carlo.cpp rand_utils.cpp vector_utils.cpp
:: wsl g++ -o singlegrid singlegrid.cpp monte_carlo.cpp rand_utils.cpp vector_utils.cpp -fopenmp
:: wsl g++ -o temp_sweep temp_sweep.cpp monte_carlo.cpp rand_utils.cpp vector_utils.cpp 
wsl g++ -o singlegrid_parallel singlegrid_parallel.cpp monte_carlo.cpp rand_utils.cpp vector_utils.cpp -fopenmp

echo "	---------------------------------"

echo "Compile complete. Running: \n"
echo "	---------------------------------"

:: wsl ./rand_test
:: wsl ./singlegrid
wsl ./singlegrid_parallel

:: wsl ./temp_sweep  8 0.0 5.0 32 ./results/temp_sweep_8.dat
:: wsl ./temp_sweep 16 0.0 5 32  ./results/temp_sweep_16.dat
:: wsl ./temp_sweep  32 0.0 5.0 32 ./results/temp_sweep_32.dat
:: wsl ./temp_sweep  64 0.0 5.0 32 ./results/temp_sweep_64.dat

:: wsl ./temp_sweep  128 2.0 2.26918531421 32 ./results/temp_sweep_128_nearcrit.dat 100000 1000 1

echo "	---------------------------------"
pause