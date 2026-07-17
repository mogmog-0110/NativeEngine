@echo off
REM Build the NativeEngine (no external dependencies -- pure C++17).
REM Usage:  build.bat        (from any shell; sets up MSVC then compiles)

set "VCVARS=C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvars64.bat"
if not exist "%VCVARS%" (
  echo ERROR: vcvars64.bat not found at "%VCVARS%"
  exit /b 1
)
call "%VCVARS%" >nul

if not exist "%~dp0build" mkdir "%~dp0build"
pushd "%~dp0build"

cl /nologo /std:c++17 /EHsc /O2 /W4 /WX ^
   /I "%~dp0src" ^
   "%~dp0src\world.cpp" "%~dp0src\physics_world.cpp" ^
   "%~dp0src\selftest.cpp" "%~dp0src\main.cpp" ^
   /Fe:NativeEngine.exe
set RC=%ERRORLEVEL%

popd
exit /b %RC%
