import argparse
import pandas as pd
import readfcs
import os
import flowkit as fk
import numpy as np


def main(fcs_path, wsp_path, output_path):
    indata = readfcs.read(fcs_path)
    df1 = pd.DataFrame(indata.X, index=indata.obs_names, columns=indata.var_names)
    
    # load workspace with fcs and wsp and create gating boolean matrix
    workspace = fk.Workspace(wsp_path, fcs_samples=fcs_path)
    workspace.analyze_samples(sample_id=os.path.basename(fcs_path))
    df_gate_membership = pd.DataFrame()
    # choose a sample to get the gates from, we assume all samples have the same gate tree
    for gate_name, gate_path in workspace.get_gate_ids(os.path.basename(fcs_path)):
        results = []
        result = workspace.get_gate_membership(
                os.path.basename(fcs_path),
                gate_name=gate_name,
                gate_path=gate_path
        )
        results.append(result)
        results = np.concatenate(results)
        df_gate_membership[':'.join(list(gate_path) + [gate_name])] = results
    # Output DataFrame as TSV file
    df = pd.concat([df1.reset_index(drop=True), df_gate_membership], axis=1)
    df.to_csv(output_path, sep='\t', index=False)
    print(f"DataFrame has been saved to {output_path}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert FCS and WSP file to DataFrame and save as TSV")
    parser.add_argument("fcs_path", help="Path to the input FCS file")
    parser.add_argument("wsp_path", help="Path to the input WSP file")
    parser.add_argument("output_path", help="Path to save the output TSV file")

    args = parser.parse_args()
    main(args.fcs_path, args.wsp_path, args.output_path)