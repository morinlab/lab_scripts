
#samples=($(cut -f4 snp.txt | sort -u))
samples="$(echo snp | cut -f4 snp.txt | sort -u)"
echo ${#samples[@]}
testarr[0]=PT007
namet=${samples[1]}
echo $namet
#name=$(echo $namet | sed -e 's/\n//g')
#echo $name
name2=${testarr[0]}
echo ${samples[1]}=$name2
#name=PT007
target_directory=/Volumes/Shared/projects/dlbcl_relapse_exome/seg_files/TITAN
#directory=${directory/' '/' '}
[ -d $target_directory/${samples[1]} ] && echo "he he directory he"
for sample in "${samples[@]}"
do
    cd /
    directory="/Volumes/Shared/projects/dlbcl_relapse_exome/seg_files/TITAN/$sample"
    echo $directory
    if [ -d $directory ]
    then
        echo "directory exists..."
        cd $directory
        pwd
    fi
# cd /Volumes/Shared/projects/dlbcl_relapse_exome/seg_files/TITAN-all/PT${sample}/
#   cut -f 2-4,10 PT${sample}T*cluster_4_segs > /Users/jgrewal/Desktop/

done