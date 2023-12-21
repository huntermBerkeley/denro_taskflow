cd build
echo "1"
time mpirun -n 1 NLSigma/nlsmSolverNUTS ../NLSigma/par/nlsmCustom.par.json >> /dev/null
echo "2"
time mpirun -n 2 NLSigma/nlsmSolverNUTS ../NLSigma/par/nlsmCustom.par.json >> /dev/null
echo "4"
time mpirun -n 4 NLSigma/nlsmSolverNUTS ../NLSigma/par/nlsmCustom.par.json >> /dev/null
echo "8"
time mpirun -n 8 NLSigma/nlsmSolverNUTS ../NLSigma/par/nlsmCustom.par.json >> /dev/null
echo "16"
time mpirun -n 16 NLSigma/nlsmSolverNUTS ../NLSigma/par/nlsmCustom.par.json >> /dev/null
cd ..