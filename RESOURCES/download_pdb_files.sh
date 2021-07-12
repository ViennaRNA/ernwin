# This script can be used to download the cif files from a list from the PDB database. As first argument it requires a csv file downloaded from the representative sets of rna structures: rna.bgsu.edu/rna3dhub/nrlist/

for file in  $(awk 'BEGIN{FS=","}{gsub(/"/, "", $2); print $2}' $1); do
  echo $file
  file=${file%%|*} 
  file_lc=${file,,}
  f2=${file_lc:1:2}
  #if [ ! -f ${file}.pdb ]; then
    #wget http://files.rcsb.org/download/${file}.pdb.gz 
    #wget ftp://ftp.wwpdb.org/pub/pdb/compatible/pdb_bundle/${f2}/${file_lc}/${file_lc}-pdb-bundle.tar.gz
  #fi
  if [ ! -f ${file}.cif ]; then
     wget http://files.rcsb.org/download/${file}.cif.gz
  fi
  #if [ ! -f ${file}.xml ]; then
  #  wget http://files.rcsb.org/download/${file}.xml.gz
  #fi
done
gunzip *.gz
