@echo off
setlocal EnableExtensions
cd /d "%~dp0\..\.."

set "MODE=standard"
if /I "%~1"=="paper" set "MODE=paper"
if not "%~1"=="" if /I not "%~1"=="paper" (
  echo Usage: scripts\windows\build_and_test.bat [paper]
  exit /b 2
)

set "BASH=C:\msys64\usr\bin\bash.exe"
if not exist "%BASH%" set "BASH=D:\msys64\usr\bin\bash.exe"
if not exist "%BASH%" (
  echo MSYS2 bash not found under C:\msys64 or D:\msys64.
  exit /b 2
)

if /I "%MODE%"=="paper" (
  set "BUILD_DIR=build-paper"
  set "PAPER_FLAG=ON"
) else (
  set "BUILD_DIR=build"
  set "PAPER_FLAG=OFF"
)

"%BASH%" -lc "export PATH=/ucrt64/bin:/mingw64/bin:/usr/bin:$PATH; cd \"$(cygpath -u '%CD%')\"; cmake -S . -B %BUILD_DIR% -G Ninja -DCMAKE_BUILD_TYPE=Release -DNANOKMC_PAPER_BUILD=%PAPER_FLAG% && cmake --build %BUILD_DIR% && ctest --test-dir %BUILD_DIR% --output-on-failure"
exit /b %ERRORLEVEL%
