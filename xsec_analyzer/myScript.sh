
cat myFileList.txt | while read line 
do
#echo $line
file=$(basename "$line")
path=$(dirname "$line")
#echo $file
#echo $path
   
# Was CC1mu1p0pi while naming the output CC1mu1piXp_* -- it ran the wrong
# selection and produced files whose names claimed otherwise.
ProcessNTuples $line  numuMC CC1mu1piXp $path/ProcessNTuples/CC1mu1piXp_$file >> $path/ProcessNTuples/temp.${file}.out &

# do something with $line here
done


#ProcessNTuples  /data/uboone/new_numi_flux/Run1_fhc_new_numi_flux_fhc_pandora_ntuple.root  numuMC CC1mu1piXp Output/CC1mu1piXp_ProcessNTuple_Run1_fhc_new_numi_flux_fhc_pandora_ntuple.root > temp.out & 


