@echo off
set "VCVARS=C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvars64.bat"
call "%VCVARS%" >nul
if not exist "%~dp0build" mkdir "%~dp0build"
pushd "%~dp0build"
cl /nologo /std:c++17 /EHsc /O2 /W4 ^
   /I "%~dp0compat" /I "%~dp0..\PhysxRender" ^
   "%~dp0..\PhysxRender\bond_graph.cpp" "%~dp0..\PhysxRender\metrics.cpp" ^
   "%~dp0test\test_science_reuse.cpp" ^
   /Fe:test_science.exe
set RC=%ERRORLEVEL%
popd
exit /b %RC%
