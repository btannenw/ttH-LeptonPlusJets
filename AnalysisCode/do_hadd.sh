

dir=~/myeos/testout/

for sample in \
tt2b \
ttb  \
ttbb \
ttcc \
ttlf 
do

hadd -f ${dir}/${sample}.root ${dir}/${sample}_*.root 

done 


hadd -f ${dir}/data.root \
    ${dir}/tt2b.root \
    ${dir}/ttb.root  \
    ${dir}/ttbb.root \
    ${dir}/ttcc.root \
    ${dir}/ttlf.root 
