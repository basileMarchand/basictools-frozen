if errorlevel 1 exit 1
%python% -m pip install --no-deps . -vv
if errorlevel 1 exit 1
move %SP_DIR%/libCppBasicTools%SHLIB_EXT% %STDLIB_DIR%/.
