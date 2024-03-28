import os
import subprocess
from concurrent.futures import ThreadPoolExecutor


def process_bed_file(bed_file, work_path, output_dir):
    # sort bed file
    subprocess.run(["bedtools", "sort", "-i", os.path.join(work_path, bed_file)])
    # Extract lines contain 'ID=gene'
    input_file_path = os.path.join(work_path, bed_file)
    output_file_path = os.path.join(output_dir, f"{os.path.splitext(bed_file)[0]}_gene.bed")

    with open(input_file_path, "r") as input_file, open(output_file_path, "w") as output_file:
        for line in input_file:
            if 'ID=gene' in line:
                output_file.write(line)


def main():
    work_path = "/bed_data"
    # work directory
    os.chdir(work_path)

    output_dir = "./gene_bed_pytest"
    os.makedirs(output_dir, exist_ok=True)

    bed_files = [file for file in os.listdir(work_path) if file.endswith(".bed")]

    with ThreadPoolExecutor() as executor:
        futures = []
        for bed_file in bed_files:
            futures.append(executor.submit(process_bed_file, bed_file, work_path, output_dir))

        for future in futures:
            future.result()


if __name__ == '__main__':
    main()
