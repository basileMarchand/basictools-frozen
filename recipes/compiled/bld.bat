if errorlevel 1 exit 1
%python% -m pip install --no-deps . -vv
rem if errorlevel 1 exit 1
rem move %SP_DIR%/libCppBasicTools%SHLIB_EXT% %STDLIB_DIR%/.
