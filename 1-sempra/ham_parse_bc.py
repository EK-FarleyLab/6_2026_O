import pickle
import pandas as pd
from Levenshtein import hamming

def process_fastq(result_path, r2_pos, bc_len, bc, bc_list, slide_window=1):

    r1_path = result_path + '_R1_001.fastq'
    r2_path = result_path + '_R2_001.fastq'
    cellBC_dict = {}
    no_map = {}

    total_sequences = 0
    discarded_sequences = 0
    mapped_no_sliding = 0
    mapped_sliding = 0
    unmapped_sequences = 0

    # Read EBC list
    bc_df = pd.read_csv(bc_list)
    bc_set = set(bc_df[bc])

    with open(r1_path, 'r') as r1, open(r2_path, 'r') as r2:
        line_num = 0
        while True:
            r1_line = r1.readline()
            r2_line = r2.readline()
            
            if not r1_line or not r2_line:
                break
            
            line_num += 1
            
            # Only process the sequence lines
            if line_num % 4 == 2:
                total_sequences += 1
                r1_seq = r1_line.strip()
                r2_seq = r2_line.strip()

                # Discard sequences containing 'N'
                if 'N' in r1_seq or 'N' in r2_seq:
                    discarded_sequences += 1
                    continue

                # Check length conditions
                if len(r1_seq) != 28 or len(r2_seq) != 90:
                    discarded_sequences += 1
                    continue

                # Parse r1 sequence into cellBC and UMI
                cellBC = r1_seq[:16]
                UMI = r1_seq[16:28]

                # Extract BC from r2 sequence
                # if cbc, 6
                # if ebc, 52
                bc_candidates = [r2_seq[r2_pos:r2_pos+bc_len]]
                if slide_window == 1:
                    bc_candidates.extend([r2_seq[r2_pos-1:r2_pos-1+bc_len], r2_seq[r2_pos+1:r2_pos+1+bc_len]])

                matched = False

                for idx, bc in enumerate(bc_candidates):
                    for ref_bc in bc_set:
                        if hamming(bc, ref_bc) <= 1:
                            matched = True
                            if idx == 0:
                                mapped_no_sliding += 1
                            else:
                                mapped_sliding += 1
                            bc = ref_bc  # Store the matched reference BC
                            break
                    if matched:
                        break

                if matched:
                    if cellBC not in cellBC_dict:
                        cellBC_dict[cellBC] = {}
                    if bc not in cellBC_dict[cellBC]:
                        cellBC_dict[cellBC][bc] = {}
                    if UMI not in cellBC_dict[cellBC][bc]:
                        cellBC_dict[cellBC][bc][UMI] = 1
                    else:
                        cellBC_dict[cellBC][bc][UMI] += 1
                else:
                    unmapped_sequences += 1
                    if cellBC not in no_map:
                        no_map[cellBC] = {}
                    if bc not in no_map[cellBC]:
                        no_map[cellBC][bc] = {}
                    if UMI not in no_map[cellBC][bc]:
                        no_map[cellBC][bc][UMI] = 1
                    else:
                        no_map[cellBC][bc][UMI] += 1

    # Save the dictionaries as pickle files
    with open(f"{result_path}/cellBC_dict.pkl", 'wb') as f:
        pickle.dump(cellBC_dict, f)
    with open(f"{result_path}/no_map.pkl", 'wb') as f:
        pickle.dump(no_map, f)

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description="Process two FASTQ files, extract barcode information, and save the result as a pickle.")
    parser.add_argument("result_path", help="Sample Name. Also path to save the resulting files.")
    parser.add_argument("r2_pos", type=int, help="The position on the R2 to extract barcode.")
    parser.add_argument("bc_len", type=int, help="Length of the barcode sequence to extract.")
    parser.add_argument("bc", help="Column on bc_list that contains all the barcode sequence. Should be either ebc or cbc.")
    parser.add_argument("bc_list", help="Path to the BC list CSV file. Should have column 'ebc' and 'cbc'.")
    parser.add_argument("--slide_window", type=int, default=1, help="Enable sliding window.")
    
    args = parser.parse_args()
    process_fastq(args.result_path, args.r2_pos, args.bc_len, args.bc, args.bc_list, args.slide_window)
