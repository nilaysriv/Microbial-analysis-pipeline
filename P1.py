import os
import subprocess

def main():
    """Main function to execute the pipeline."""

    # Update and install prerequisites
    os.system("sudo apt upgrade")
    os.system("sudo apt update")
    os.system("sudo apt install build-essential")
    os.system("sudo apt install R-base-core")

    # Prompt for input files
    fastq1 = input("Enter the first FASTQ file: ")
    fastq2 = input("Enter the second FASTQ file: ")

    # Perform profiling with Greengenes and SILVA databases
    for database in ("G", "S"):
        output_dir = f"Results_{'Greengenes' if database == 'G' else 'SILVA'}"
        os.system(f"PM-parallel-meta -R {fastq1} {fastq2} -o {output_dir} -l 150 -D {database}")

    # Create list.txt with paths to classification.txt files
    with open("list.txt", "w") as f:
        for database in ("Greengenes", "SILVA"):
            output_path = f"Results_{database}/classification.txt"
            f.write(f"{database}\t{output_path}\n")

    # Perform taxonomic classification and functional profiling
    os.system("PM-select-taxa -l list.txt -o taxa.txt -L 5")  # Genus level
    os.system("PM-select-func -l list.txt -o func.txt -L 2")  # KEGG pathway level 2

if __name__ == "__main__":
    main()
