# Conditional statements and loops

## for
Using `for` loop, we can execute a set of commands for a finite number of times for every item in a list.
### Syntax
```
for VARIABLE in LIST
do
  command1
  command2
done
```
### Example
```
for query in /cluster/tufts/Yourlab/input/*.fasta
do
    echo "Running BLAST for $query"
    blastn -query $query -db /path/to/database -out ${query%.fasta}_blast_results.txt -outfmt 6
done
```

## while
The command next to `while` is evaluated. If it is successful or 0, then the commands between `do` and `done` are executed.
### Syntax
```
while CONDITION
do
    command1
    command2
    ...
done
```
### Example
```
#!/bin/bash

# File containing the list of sequence files
file_list="sequence_files.txt"
# Path to the BLAST database
blast_db="/path/to/blast_database"

# Read the list of files line by line
while IFS= read -r sequence_file
do
    # Check if the file exists
    if [ -f "$sequence_file" ]; then
        echo "Processing $sequence_file"
        # Run BLAST on the current sequence file
        blastn -query "$sequence_file" -db "$blast_db" -out "${sequence_file%.fasta}_blast_results.txt" -outfmt 6
        echo "Finished processing $sequence_file"
    else
        echo "File $sequence_file not found"
    fi
done < "$file_list"

echo "All files have been processed."
```

The content of `sequence_files.txt` is the name of all of the query fasta files: 
```
sequence1.fasta
sequence2.fasta
sequence3.fasta
```
## if: conditionals in a bash script
In bioinformatics, if statements are used to make conditional decisions in scripts, such as checking whether a file or result already exists before performing data analysis or processing tasks, thus optimizing workflows and avoiding redundant computations.
### Syntax
```
if [commands]
then
  [if-statements]
elif [commands]
  [elif-statements]
else
  [else-statements]
fi
```
### Example
```
#!/bin/bash

# Define the paths for input files and output file
input_bam="/path/to/genomic_data.bam"
output_vcf="/path/to/output_directory/variants.vcf"

# Check if the output VCF file already exists
if [ -f "$output_vcf" ]; then
    echo "Output file $output_vcf already exists. Skipping the variant calling step."
else
    echo "Output file $output_vcf does not exist. Running variant calling."

    # Run the variant calling tool (e.g., bcftools)
    bcftools mpileup -Ou -f /path/to/reference_genome.fasta "$input_bam" | \
    bcftools call -mv -Ov -o "$output_vcf"

    echo "Variant calling complete. Results saved to $output_vcf"
fi
```

This if statement ensures that the variant calling process is skipped if the output file already exists, preventing unnecessary computations and saving processing time. This approach is useful in bioinformatics workflows where redundant analyses can be avoided by checking for existing results.




[Previous: Commands](03_basictools.md)                                                                 

[Next: Advanced tools](05_advanced_tools.md)
