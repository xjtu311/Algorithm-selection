@if "%DEBUG%" == "" @echo off
@rem ##########################################################################
@rem
@rem  portfolio-evaluator startup script for Windows
@rem
@rem ##########################################################################

@rem Set local scope for the variables with windows NT shell
if "%OS%"=="Windows_NT" setlocal

@rem Add default JVM options here. You can also use JAVA_OPTS and PORTFOLIO_EVALUATOR_OPTS to pass JVM options to this script.
set DEFAULT_JVM_OPTS="-Xmx2048m" "-Dsun.zip.disableMemoryMapping"

set DIRNAME=%~dp0
if "%DIRNAME%" == "" set DIRNAME=.
set APP_BASE_NAME=%~n0
set APP_HOME=%DIRNAME%..

@rem Find java.exe
if defined JAVA_HOME goto findJavaFromJavaHome

set JAVA_EXE=java.exe
%JAVA_EXE% -version >NUL 2>&1
if "%ERRORLEVEL%" == "0" goto init

echo.
echo ERROR: JAVA_HOME is not set and no 'java' command could be found in your PATH.
echo.
echo Please set the JAVA_HOME variable in your environment to match the
echo location of your Java installation.

goto fail

:findJavaFromJavaHome
set JAVA_HOME=%JAVA_HOME:"=%
set JAVA_EXE=%JAVA_HOME%/bin/java.exe

if exist "%JAVA_EXE%" goto init

echo.
echo ERROR: JAVA_HOME is set to an invalid directory: %JAVA_HOME%
echo.
echo Please set the JAVA_HOME variable in your environment to match the
echo location of your Java installation.

goto fail

:init
@rem Get command-line arguments, handling Windowz variants

if not "%OS%" == "Windows_NT" goto win9xME_args
if "%@eval[2+2]" == "4" goto 4NT_args

:win9xME_args
@rem Slurp the command line arguments.
set CMD_LINE_ARGS=
set _SKIP=2

:win9xME_args_slurp
if "x%~1" == "x" goto execute

set CMD_LINE_ARGS=%*
goto execute

:4NT_args
@rem Get arguments from the 4NT Shell from JP Software
set CMD_LINE_ARGS=%$

:execute
@rem Setup the command line

set CLASSPATH=%APP_HOME%\lib\hydra-1.1-development-cae8151.jar;%APP_HOME%\lib\fastrf-1.10.01-hydra.jar;%APP_HOME%\lib\mysqldbtae-2.08.01-hydra.jar;%APP_HOME%\lib\smac-smac-v2.10.02-hydra.jar;%APP_HOME%\lib\jDTs-jDTs.jar;%APP_HOME%\lib\jParFor-jParFor.jar;%APP_HOME%\lib\matlab-matlab.jar;%APP_HOME%\lib\Zilla-1.8-development-8fc9055.jar;%APP_HOME%\lib\commons-collections-3.2.1.jar;%APP_HOME%\lib\commons-lang3-3.3.1.jar;%APP_HOME%\lib\commons-collections4-4.0.jar;%APP_HOME%\lib\opencsv-2.4.jar;%APP_HOME%\lib\guava-18.0.jar;%APP_HOME%\lib\logback-access-1.1.2.jar;%APP_HOME%\lib\logback-classic-1.1.2.jar;%APP_HOME%\lib\logback-core-1.1.2.jar;%APP_HOME%\lib\jcommander-1.35-patched.jar;%APP_HOME%\lib\libsvm-3.17.jar;%APP_HOME%\lib\weka-stable-3.6.6.jar;%APP_HOME%\lib\spi-0.2.4.jar;%APP_HOME%\lib\commons-io-2.1.jar;%APP_HOME%\lib\commons-math-2.2.jar;%APP_HOME%\lib\exp4j-0.3.10.jar;%APP_HOME%\lib\jama-1.0.2.jar;%APP_HOME%\lib\jcip-annotations-1.0.jar;%APP_HOME%\lib\numerics4j-1.3.jar;%APP_HOME%\lib\mchange-commons-java-0.2.7.jar;%APP_HOME%\lib\mysql-connector-java-5.1.33.jar;%APP_HOME%\lib\c3p0-0.9.5-pre8.jar;%APP_HOME%\lib\c3p0-oracle-thin-extras-0.9.5-pre8.jar;%APP_HOME%\lib\metrics-core-3.1.0.jar;%APP_HOME%\lib\httpasyncclient-4.1.jar;%APP_HOME%\lib\metrics-jvm-3.1.0.jar;%APP_HOME%\lib\spring-web-4.1.5.RELEASE.jar;%APP_HOME%\lib\jna-4.0.0.jar;%APP_HOME%\lib\httpclient-4.4.1.jar;%APP_HOME%\lib\java-cup-0.11a.jar;%APP_HOME%\lib\httpcore-4.4.1.jar;%APP_HOME%\lib\httpcore-nio-4.4.1.jar;%APP_HOME%\lib\commons-logging-1.2.jar;%APP_HOME%\lib\spring-aop-4.1.5.RELEASE.jar;%APP_HOME%\lib\spring-beans-4.1.5.RELEASE.jar;%APP_HOME%\lib\spring-context-4.1.5.RELEASE.jar;%APP_HOME%\lib\spring-core-4.1.5.RELEASE.jar;%APP_HOME%\lib\aopalliance-1.0.jar;%APP_HOME%\lib\spring-expression-4.1.5.RELEASE.jar;%APP_HOME%\lib\commons-math3-3.3.jar;%APP_HOME%\lib\aeatk-2.08.01-hydra.jar;%APP_HOME%\lib\jackson-datatype-guava-2.5.1.jar;%APP_HOME%\lib\slf4j-api-1.7.7.jar;%APP_HOME%\lib\commons-codec-1.9.jar;%APP_HOME%\lib\jackson-databind-2.5.1.jar;%APP_HOME%\lib\jackson-core-2.5.1.jar;%APP_HOME%\lib\jackson-annotations-2.5.0.jar

@rem Execute portfolio-evaluator
"%JAVA_EXE%" %DEFAULT_JVM_OPTS% %JAVA_OPTS% %PORTFOLIO_EVALUATOR_OPTS%  -classpath "%CLASSPATH%" ca.ubc.cs.beta.hydra.executors.PortfolioExecutor %CMD_LINE_ARGS%

:end
@rem End local scope for the variables with windows NT shell
if "%ERRORLEVEL%"=="0" goto mainEnd

:fail
rem Set variable PORTFOLIO_EVALUATOR_EXIT_CONSOLE if you need the _script_ return code instead of
rem the _cmd.exe /c_ return code!
if  not "" == "%PORTFOLIO_EVALUATOR_EXIT_CONSOLE%" exit 1
exit /b 1

:mainEnd
if "%OS%"=="Windows_NT" endlocal

:omega
