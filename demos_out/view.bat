@echo off
REM View a NativeEngine demo in the modern OpenGL viewer.
REM Usage (from this folder):  view pbc_pair
REM   names: pile cylinders pendulum gas pbc_gas pbc_pair
if "%~1"=="" (
  echo usage: view ^<name^>   ^(pile cylinders pendulum gas pbc_gas pbc_pair^)
  exit /b 1
)
set "VIEWER=%~dp0..\build\NativeViewer.exe"
set "REC=%~dp0demo_%~1.pxrf"
if not exist "%REC%" ( echo not found: %REC%  ^(run: NativeEngine.exe demo %~1^) & exit /b 1 )
"%VIEWER%" "%REC%"
