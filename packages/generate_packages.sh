#!/bin/bash

mkdir -p dist

for p in omuse*
do
  cd $p || exit $?
  rm -rf dist || exit $?
  python setup.py sdist || exit $?
  # Issue with setuptools_scm root handling
  # https://github.com/pypa/setuptools_scm/issues/188
  # To be updated after new version of setuptools_scm using:
  # https://github.com/pypa/setuptools_scm/pull/870
  cd dist
  tarname=$(ls *.tar.gz)
  tar -xzf $tarname
  pname=$(basename $tarname .tar.gz)
  toml_file=$(ls $pname/*.toml)
  echo $toml_file
  sed -i '/root/d' $toml_file
  tar -czf $tarname $pname
  cd -
  cp dist/*.tar.gz ../dist/ || exit $?
  cd .. 
done
