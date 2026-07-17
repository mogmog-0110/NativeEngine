@echo off
REM Build the NativeEngine Visual Debugger (modern OpenGL via freeglut + Dear ImGui).
REM Requires FREEGLUT_ROOT (or the default path below). ImGui is vendored under
REM viewer\third_party\imgui, so no other external dependency.
call "C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvars64.bat" >nul
cd /d "%~dp0.."
if not exist build mkdir build

set "FG=%FREEGLUT_ROOT%"
if "%FG%"=="" set "FG=E:\dev\freeglut\3.4_1.1"
set "IMGUI=viewer\third_party\imgui"

cl /nologo /std:c++17 /EHsc /O2 /W3 /D_CRT_SECURE_NO_WARNINGS /DFREEGLUT_LIB_PRAGMAS=0 ^
   /I viewer /I src /I "%IMGUI%" /I "%FG%\include" ^
   viewer\viewer.cpp src\world.cpp ^
   "%IMGUI%\imgui.cpp" "%IMGUI%\imgui_draw.cpp" "%IMGUI%\imgui_tables.cpp" "%IMGUI%\imgui_widgets.cpp" ^
   "%IMGUI%\backends\imgui_impl_glut.cpp" "%IMGUI%\backends\imgui_impl_opengl3.cpp" ^
   /Fe:build\NativeViewer.exe /Fo:build\ ^
   /link /LIBPATH:"%FG%\lib\win64" freeglut.lib opengl32.lib
set RC=%ERRORLEVEL%

REM stage the freeglut DLL next to the viewer
copy /Y "%FG%\bin\win64\freeglut.dll" build\ >nul 2>&1
exit /b %RC%
