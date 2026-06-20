@echo off
REM ============================================================
REM  Package influgen_2.exe with all runtime dependencies
REM  for deployment on Intel AND AMD servers.
REM  Run this AFTER building Release|x64 in Visual Studio.
REM ============================================================

set DEPLOY_DIR=influgen_deploy
set EXE_PATH=x64\Release\influgen_2.exe

if not exist "%EXE_PATH%" (
    echo ERROR: %EXE_PATH% not found. Build Release x64 first!
    exit /b 1
)

if not defined MKLROOT (
    echo ERROR: MKLROOT not set. Install Intel oneAPI MKL.
    exit /b 1
)

echo Creating deployment folder: %DEPLOY_DIR%
if exist %DEPLOY_DIR% rmdir /s /q %DEPLOY_DIR%
mkdir %DEPLOY_DIR%

REM --- Copy the exe ---
echo Copying influgen_2.exe ...
copy "%EXE_PATH%" "%DEPLOY_DIR%\"

REM --- Copy FSTI authorization DLL ---
echo Copying FSTI authorization ...
copy "lib\x64\fsti_netlic_server.dll" "%DEPLOY_DIR%\"

REM --- Copy MKL runtime DLLs ---
echo Copying MKL runtime DLLs ...
copy "%MKLROOT%\bin\mkl_core.2.dll"       "%DEPLOY_DIR%\"
copy "%MKLROOT%\bin\mkl_sequential.2.dll"  "%DEPLOY_DIR%\"
copy "%MKLROOT%\bin\mkl_def.2.dll"         "%DEPLOY_DIR%\"
copy "%MKLROOT%\bin\mkl_avx2.2.dll"        "%DEPLOY_DIR%\"
copy "%MKLROOT%\bin\mkl_mc3.2.dll"         "%DEPLOY_DIR%\"
copy "%MKLROOT%\bin\mkl_vml_def.2.dll"     "%DEPLOY_DIR%\"
copy "%MKLROOT%\bin\mkl_vml_avx2.2.dll"    "%DEPLOY_DIR%\"
copy "%MKLROOT%\bin\mkl_vml_cmpt.2.dll"    "%DEPLOY_DIR%\"

REM --- Copy MSVC runtime ---
echo Copying MSVC runtime ...
REM These are usually already installed via vc_redist.x64.exe on servers,
REM but including them ensures it works even without the redist installed.
set VCREDIST=C:\Program Files\Microsoft Visual Studio\18\Community\VC\Redist\MSVC\14.51.36231\x64\Microsoft.VC145.CRT
if exist "%VCREDIST%" (
    copy "%VCREDIST%\vcruntime140.dll"   "%DEPLOY_DIR%\"
    copy "%VCREDIST%\vcruntime140_1.dll" "%DEPLOY_DIR%\"
    copy "%VCREDIST%\msvcp140.dll"       "%DEPLOY_DIR%\"
    copy "%VCREDIST%\msvcp140_1.dll"     "%DEPLOY_DIR%\"
) else (
    echo WARNING: MSVC redist not found at expected path.
    echo          Install vc_redist.x64.exe on the target server instead.
)

echo.
echo ============================================================
echo  Deployment package ready: %DEPLOY_DIR%\
echo  Total files:
dir /b "%DEPLOY_DIR%" | find /c /v ""
echo.
echo  Copy the entire folder to the target server.
echo  Works on both Intel and AMD CPUs.
echo ============================================================
