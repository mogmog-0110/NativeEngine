@echo off
REM Play a NativeEngine demo recording in the OpenGL viewer.
REM Usage (from this folder):  play pbc_pair       (or: pile cylinders pendulum gas pbc_gas)
if "%~1"=="" (
  echo usage: play ^<name^>    where name is one of:
  echo   pile cylinders pendulum gas pbc_gas pbc_pair
  exit /b 1
)
set "PLAYER=%~dp0..\..\x64\ReleaseRender\PhysxRender.exe"
set "REC=%~dp0demo_%~1.pxrf"
if not exist "%REC%" (
  echo not found: %REC%   ^(run: NativeEngine.exe demo %~1^)
  exit /b 1
)
"%PLAYER%" "%REC%"
