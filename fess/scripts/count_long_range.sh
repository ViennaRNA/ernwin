rm both.count; for file in `find -name 'temp.comp'`; do defines=`cat $file | grep define | wc -l`; longrange=`cat $file | grep longrange | wc -l`; echo $file $defines $longrange >> both.count; done;
