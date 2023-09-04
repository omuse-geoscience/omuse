#!/bin/bash

mkdir -p dist

for p in omuse*
do
  cd $p || exit $?
  rm -rf dist || exit $?
  python setup.py sdist || exit $?
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
