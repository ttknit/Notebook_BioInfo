import pysam
import sys

def extract_info(bam_file, output_file):
    """
    Extracts read ID, ref ID, start, and end positions from a BAM file.
    
    Args:
        bam_file (str): Path to the input BAM file.
        output_file (str): Path to the output TSV file.
    """
    try:
        # Open BAM file for reading
        samfile = pysam.AlignmentFile(bam_file, "rb")
        
        # Open output file for writing
        with open(output_file, 'w') as out_f:
            # Iterate through each read in the BAM file
            for read in samfile.fetch():
                # Check if the read is mapped
                if not read.is_unmapped:
                    read_id = read.query_name
                    ref_id = samfile.get_reference_name(read.reference_id)
                    # Use pysam's built-in properties for start and end positions
                    start = read.reference_start + 1 # pysam is 0-based, SAM is 1-based
                    end = read.reference_end
                    
                    # Write to the output file
                    out_f.write(f"{read_id}\t{ref_id}\t{start}\t{end}\n")
    except FileNotFoundError:
        print(f"Error: BAM file '{bam_file}' not found.", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"An error occurred: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python extract_info.py <input.bam> <output.tsv>", file=sys.stderr)
        sys.exit(1)
    
    input_bam = sys.argv[1]
    output_tsv = sys.argv[2]
    
    extract_info(input_bam, output_tsv)
