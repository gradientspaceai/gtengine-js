@echo off
rem Builds oracle\build\oracle.exe from oracle\cpp with MSVC against the
rem upstream GTE headers.
rem
rem   oracle\build.cmd [case-file-glob]      e.g. oracle\build.cmd v19-*.cpp
rem
rem   GTE_PATH  upstream GTE folder (default D:\git\GeometricTools-upstream\GTE)
rem   VCVARS    path to vcvars64.bat (default: VS 2022 Community)
setlocal
if "%GTE_PATH%"=="" set "GTE_PATH=D:\git\GeometricTools-upstream\GTE"
if "%VCVARS%"=="" set "VCVARS=C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvars64.bat"
if not exist "%GTE_PATH%\Mathematics\Vector.h" (
    echo GTE headers not found at %GTE_PATH%
    exit /b 1
)
if not exist "%VCVARS%" (
    echo vcvars64.bat not found at %VCVARS%
    exit /b 1
)
set "CASES=%~1"
if "%CASES%"=="" set "CASES=*.cpp"
set "ROOT=%~dp0"
set "UPSTREAM=unknown"
for /f %%c in ('git -C "%GTE_PATH%" rev-parse HEAD') do set "UPSTREAM=%%c"
call "%VCVARS%" >nul 2>&1
if not exist "%ROOT%build\obj" mkdir "%ROOT%build\obj"
cl /nologo /std:c++17 /EHsc /O2 /fp:precise /bigobj /MP /W3 /wd4244 /wd4267 /wd4305 ^
   /DORACLE_UPSTREAM_COMMIT=\"%UPSTREAM%\" /I"%GTE_PATH%" /I"%ROOT%cpp" ^
   /Fo"%ROOT%build\obj\\" /Fe"%ROOT%build\oracle.exe" ^
   "%ROOT%cpp\main.cpp" "%ROOT%cpp\cases\%CASES%"
exit /b %ERRORLEVEL%
