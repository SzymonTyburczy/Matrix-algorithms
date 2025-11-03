To compile benchmark_main.cpp use the following command in PowerShell:

```powershell
++ -std=c++17 -O2 -Wall -Wextra -I lab1\src `
>>   lab1\src\main_benchmark.cpp `
>>   lab1\src\helperFunctions.cpp `
>>   lab1\src\matrix_Strassen.cpp `
>>   lab1\src\matrix_Binet.cpp `
>>   lab1\src\matrix_ai.cpp `
>>   -lpsapi -o lab1\benchmark.exe;
```

To run the benchmark following tests:
iterative with argument `iterative`
Strassen with argument `strassen`
Binet with argument `binet`
AI recursive with argument `ai`
