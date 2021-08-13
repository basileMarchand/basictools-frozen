$PYTHON -m pip install --no-deps . -vv --use-feature=in-tree-build
echo $SP_DIR
echo $SHLIB_EXT
echo  $STDLIB_DIR
echo "symbolic link"
echo ln -s $SP_DIR/libCppBasicTools$SHLIB_EXT $STDLIB_DIR/../.
ln -s $SP_DIR/libCppBasicTools$SHLIB_EXT $STDLIB_DIR/../.
