for fn in $1/*.cg; do
  basefn=$(basename $fn)
  pdb_id=${basefn%_*}
  file=$(echo ${pdb_id} | tr '[:upper:]' '[:lower:]')
  if [ ! -f ${file}.cif ]; then
     wget http://files.rcsb.org/download/${file}.cif.gz
  fi
done


