#!/bin/sh

name="kalki"
new_name="legopipe"
mv bin/$name bin/$new_name
list_of_files=`find . -name "*" -type f | grep -v git | grep -v $0 | grep -v _temp`
for file in $list_of_files
do
	echo $file
	sed -e "s/$name/$new_name/g" $file > _temp
	mv _temp $file
done
