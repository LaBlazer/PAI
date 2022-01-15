# LU Decomposition

## C# Version usage

Generate NxN matrix and save it to a file 

```
LUDecomp.exe -g <size> <output file>
```

To run decomposition parallelly you need to specify input matrix file and thread amount

```
LUDecomp.exe <input file> <thread amount>
```

To quickly test the program on matrix with size 2000x2000 run

```
LUDecomp.exe -g 2000 mtx.txt && LUDecomp.exe mtx.txt 4
```

## C++ Version usage

Basically the same as C# version

```
LUDecompCpp.exe <input file> <thread amount>
```