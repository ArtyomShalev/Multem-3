# multem3 building ----------------
cd ../..
./clean.sh
cmake .
make
cp bin/multem3 examples/CompPhysComm2024_submission/.
# multem2 building ----------------
cd examples/CompPhysComm2024_submission/multem_different_versions/multem2
./clean.sh
cmake .
make
cp multem2 ../../. 
# multem2 lapack building ---------
cd ../multem2_lapack
./clean.sh
cmake .
make
cd multem2_lapack
cp multem2_with_lapack ../../../. 
# multem3_cerf
cd ../../multem3_cerf
./clean.sh
cmake .
make
cp bin/multem3_cerf ../../.
