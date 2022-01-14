# LU Decomposition

## C# Version usage

Generate NxN matrix and save it to a file 

```
LUDecomp.exe -g <size> <output file>
```

To run decomposition paralelly you need to have Microsoft MPI installed and then use the ```mpiexec``` tool

```
mpiexec -n <process amount> LUDecomp.exe <input file>
```

## C++ Version usage
TBD