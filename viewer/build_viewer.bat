@echo off
REM Build the NativeEngine viewer (modern OpenGL via freeglut).
REM Requires FREEGLUT_ROOT (or the default path below).
call "C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvars64.bat" >nul
cd /d "%~dp0.."
if not exist build mkdir build

set "FG=%FREEGLUT_ROOT%"
if "%FG%"=="" set "FG=E:\dev\freeglut\3.4_1.1"

cl /nologo /std:c++17 /EHsc /O2 /W4 /DFREEGLUT_LIB_PRAGMAS=0 ^
   /I viewer /I src /I "%FG%\include" ^
   viewer\viewer.cpp src\world.cpp ^
   /Fe:build\NativeViewer.exe /Fo:build\ ^
   /link /LIBPATH:"%FG%\lib\win64" freeglut.lib opengl32.lib
set RC=%ERRORLEVEL%

REM stage the freeglut DLL next to the viewer
copy /Y "%FG%\bin\win64\freeglut.dll" build\ >nul 2>&1
exit /b %RC%
