@echo off
setlocal enableextensions

REM Initialize MSVC build environment (vcvars64)
call "C:\Program Files (x86)\Microsoft Visual Studio\2022\BuildTools\VC\Auxiliary\Build\vcvars64.bat" || goto :error

REM Configure and build with CMake (NMake Makefiles)
"C:\Program Files\CMake\bin\cmake.exe" -G "NMake Makefiles" -DCMAKE_BUILD_TYPE=Release .. || goto :error
"C:\Program Files\CMake\bin\cmake.exe" --build . || goto :error

REM Run sanity benchmarks
.\bench.exe --algo dijkstra --gen er --n 10000 --p 0.001 --threads 1 || goto :error
.\bench.exe --algo delta --gen er --n 10000 --p 0.001 --threads 8 --delta 4.0 || goto :error

echo Completed build and sanity runs successfully.
exit /b 0

:error
echo Build or run failed with error %ERRORLEVEL%.
exit /b %ERRORLEVEL%
