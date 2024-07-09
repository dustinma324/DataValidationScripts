import pandas as pd
import argparse
import os

def convert_delimiters(input_filename):
    # Read the CSV file with comma delimiter
    data = pd.read_csv(input_filename, delimiter=',')
    
    # Generate the output file name
    base_name = os.path.splitext(input_filename)[0]
    output_filename = base_name + '_space_delimited.txt'
    
    # Save the data to the new file with space as delimiter
    data.to_csv(output_filename, sep=' ', index=False)
    
    return output_filename

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Convert CSV file delimiters from commas to spaces.')
    parser.add_argument('input_files', nargs='+', help='The input CSV files with comma delimiters')
    
    args = parser.parse_args()
    
    for input_file in args.input_files:
        output_file = convert_delimiters(input_file)
        print(f"File '{input_file}' has been converted and saved as '{output_file}'")

