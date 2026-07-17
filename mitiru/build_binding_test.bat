@echo off
call "C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvars64.bat" >nul
cd /d "%~dp0.."
if not exist build mkdir build
cl /nologo /std:c++17 /EHsc /O2 /W4 /I mitiru\compat /I mitiru /I src src\world.cpp src\physics_world.cpp mitiru\test_binding.cpp /Fe:build\test_binding.exe /Fo:build\
exit /b %ERRORLEVEL%
