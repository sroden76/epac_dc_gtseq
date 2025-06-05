!/bin/bash

# Move to working directory
#cd /Users/aonoufriou/Documents/GitHub/Turtle_GTseq_SNPs/data-raw/sam.files/test
cd \Users\suzanne.roden\Documents\GitHub\Mnov.gtseq.analysis\data-raw\sam.files

#List of replicate IDs
replicates="$(cat ind.to.merge.txt)"


# Loop over replicate IDs
for i in $replicates; do
    ## Define temporary sam files
        sam1="./all/${i}_mapped2ref.sam"
        sam2="./all/${i}b_mapped2ref.sam"
        sam3="./all/${i}c_mapped2ref.sam"
        merged="./merged/${i}_merged.sam"

 # Check if input SAM files exist
    if [[ -f "$sam1" && -f "$sam2" ]]; then
     
     # Merge sam files
        samtools merge "$merged" "$sam1" "$sam2"

        echo "Merged: $merged"; 
    else
        echo "Warning: Missing files for replicate $i"
    fi
done

##Notes
# -f within the "if" statement is used to test if a file exists or not
# here it will only be true if both files are present (using the && operator)
# if one file is missing, the script will print the name of the missing file 
# and not run the merge command


