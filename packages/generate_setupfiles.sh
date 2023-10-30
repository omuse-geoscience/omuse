#!/usr/bin/env bash

while read -r line; do
	pkgdir=$(echo "print('omuse-${line}'.lower())" | python);
	echo ${pkgdir}
	python setup_template.py ${line} > ${pkgdir}/setup.py;
	python pyproject_template.py ${line} > ${pkgdir}/pyproject.toml;
done < community_package_names
