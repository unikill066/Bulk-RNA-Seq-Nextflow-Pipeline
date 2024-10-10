import os, sys, csv

def parse_final_out_file(filepath):
    """
    Parse relevant information from a .final.out file.
    """
    data = dict()
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith("Number of input reads"):
                data["Number of input reads"] = line.split("|")[1].strip()
            elif line.startswith("Average input read length"):
                data["Average input read length"] = line.split("|")[1].strip()
            elif line.startswith("Uniquely mapped reads number"):
                data["Uniquely mapped reads number"] = line.split("|")[1].strip()
            elif line.startswith("Uniquely mapped reads %"):
                data["Uniquely mapped reads %"] = line.split("|")[1].strip()
            elif line.startswith("Average mapped length"):
                data["Average mapped length"] = line.split("|")[1].strip()
            elif line.startswith("Number of splices: Total"):
                data["Number of splices (Total)"] = line.split("|")[1].strip()
            elif line.startswith("Mismatch rate per base"):
                data["Mismatch rate per base (%)"] = line.split("|")[1].strip()
            elif line.startswith("Number of reads mapped to multiple loci"):
                data["Number of reads mapped to multiple loci"] = line.split("|")[1].strip()
            elif line.startswith("% of reads mapped to multiple loci"):
                data["% of reads mapped to multiple loci"] = line.split("|")[1].strip()
            elif line.startswith("Number of reads unmapped: too short"):
                data["Number of reads unmapped: too short"] = line.split("|")[1].strip()
            elif line.startswith("% of reads unmapped: too short"):
                data["% of reads unmapped: too short"] = line.split("|")[1].strip()
    
    return data

def parse_all_final_out_files(directory):
    """
    Process all .final.out files in the directory and save results to CSV.
    """
    final_data = list()
    header = None

    for filename in os.listdir(directory):
        if filename.endswith(".final.out"):
            filepath = os.path.join(directory, filename)
            parsed_data = parse_final_out_file(filepath)
            parsed_data["File"] = filename  # Add the filename at the end
            final_data.append(parsed_data)

            if header is None:
                header = list(parsed_data.keys())

    # Reorder header to move "File" to the front
    if header is not None:
        header = ["File"] + [col for col in header if col != "File"]

    return final_data, header

def write_to_csv(final_data, header, output_csv):
    with open(output_csv, 'w', newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=header)
        writer.writeheader()
        for row in final_data:
            writer.writerow(row)

if __name__ == "__main__":

    if len(sys.argv) != 3:
        print("Usage: python parse_final_out.py <config_directory> <directory_with_final_out_files>")
        sys.exit(1)

    config_directory = sys.argv[1]
    directory = sys.argv[2]

    if not os.path.isdir(directory):
        print(f"Error: {directory} is not a valid directory.")
        sys.exit(1)

    final_data, header = parse_all_final_out_files(directory)

    os.makedirs(os.path.join(config_directory, "2_1_map_metrics_output_qc"), exist_ok=True)
    output_csv = os.path.join(config_directory, "2_1_map_metrics_output_qc" , "final_out_summary.csv")
    write_to_csv(final_data, header, output_csv)
    print(f"Summary CSV file saved as {output_csv}")