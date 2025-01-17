#!/bin/bash
###################################################

# executable
var=xmeshfem3D

# configuration parameters
CONF_PARAM="--enable-vectorization --enable-openmp"

###################################################
testdir=`pwd`
me=`basename "$0"`

#checks if ROOT valid
if [ -z "${ROOT}" ]; then export ROOT=../../ ; fi

# sets source directory
cd $ROOT/
srcdir=`pwd`

cd $testdir/

# title
echo >> $testdir/results.log
echo "$me in: $testdir" >> $testdir/results.log
echo >> $testdir/results.log

#cleanup
rm -rf config.log config.status
rm -rf ./bin ./obj ./setup ./OUTPUT_FILES ./DATA

# configuration
# (out-of-source compilation)
echo "configuration: $srcdir/configure ${CONF_PARAM}" >> $testdir/results.log
$srcdir/configure ${CONF_PARAM} >> $testdir/results.log 2>&1

# checks exit code
if [[ $? -ne 0 ]]; then
  echo >> $testdir/results.log
  echo "configuration failed, please check..." >> $testdir/results.log
  exit 1
fi

# we need to output to console output, otherwise tests will fail by timeout in travis
sed -i "s:IMAIN .*:IMAIN = ISTANDARD_OUTPUT:" setup/constants.h >> $testdir/results.log

echo "" >> $testdir/results.log
echo "successful configuration" >> $testdir/results.log

# single compilation
echo "compilation: $var" >> $testdir/results.log
make clean >> $testdir/results.log 2>&1
make -j 4 $var >> $testdir/results.log 2>&1

# checks exit code
if [[ $? -ne 0 ]]; then
  echo >> $testdir/results.log
  echo "compilation failed, please check..." >> $testdir/results.log
  exit 1
fi

# checks binary
if [ ! -e bin/$var ]; then
  echo "compilation of $var failed, please check..." >> $testdir/results.log
  exit 1
else
  echo "binary exists: $var" >> $testdir/results.log
fi
echo "" >> $testdir/results.log

#cleanup
rm -rf ./bin/*

echo "successful compilation" >> $testdir/results.log


