#! /usr/bin/env python3
# CrossNN local script
import NN_model
import argparse
import MGMT_methylation_status as mgmt
import os


def load_model(model_path):
    model = NN_model.NN_classifier(model_path)
    return model

def predict(model, input_file_path, epic_annotation_path=None):
      #Check if input file is cov or bedMethyl
        if input_file_path.endswith('.cov') or input_file_path.endswith('.cov.gz'):
            assert epic_annotation_path is not None, "For .cov input files, please provide the --epic_annotation argument"
            return model.predict_from_cov(input_file_path, epic_annotation_path)
        else:
            return model.predict_from_bedMethyl(input_file_path)

def write_predictions(predictions, output_file_path):
    with open(output_file_path, 'w') as f:
        f.write(f"Number of features used: {predictions[2]}\n")
        for i in range(len(predictions[0])):
            f.write(f"{predictions[1][i]}\t{predictions[0][i]:.6f}\n")
   

def plot_predictions(predictions, plot_file_path):
    import matplotlib.pyplot as plt
    import numpy as np

    scores, class_labels, num_features = predictions

    # Select top 10
    top_n = 10
    top_scores = scores[:top_n]
    top_labels = class_labels[:top_n]

    plt.figure(figsize=(10, 6))
    y_pos = np.arange(len(top_labels))
    plt.barh(y_pos, top_scores, align='center', alpha=0.7, color='#00A389')
    plt.axvline(x=0.2, color='#F98300', linestyle='--')
    plt.yticks(y_pos, top_labels)
    plt.xlabel('Scores')
    plt.title('Top 10 Prediction Scores by Class')
    plt.text(max(top_scores)*0.7, 9, f'Number of features used: {num_features}', fontsize=10, verticalalignment='center')


    plt.savefig(plot_file_path)
    plt.close()

def get_mgmt_status(cov_file, output, reference="hg38"):
    report_file = output + ".mgmt_report.txt"
    coverage_DMR1, coverage_DMR2 = mgmt.load_coverage(cov_file, reference)
    status_DMR1, percent_methylation_DMR1, status_DMR2, percent_methylation_DMR2, used_positions, discarded_positions = mgmt.determine_mgmt_status(coverage_DMR1, coverage_DMR2)
    mgmt.write_report(report_file, status_DMR1, percent_methylation_DMR1, status_DMR2, percent_methylation_DMR2, used_positions, discarded_positions, cov_file)

def main():
    parser = argparse.ArgumentParser(description="Run CrossNN model")
    parser.add_argument("-m", "--model", required=True, help="Path to the model file, created with the training.py script")
    parser.add_argument("-i", "--input", required=True, help="Path to the sample bedMethyl file")
    parser.add_argument("-o", "--output", required=True, help="Path to the output file")
    parser.add_argument("--epic_annotation", required=True, help="Path to the EPIC annotation file")
    parser.add_argument("-r", "--reference", choices=["hg19", "hg38"], default="hg38", help="Reference genome version (default: hg38) for MGMT methylation status analysis")
    args = parser.parse_args()

    outfilename, _ = os.path.splitext(args.output)

    model = load_model(args.model)
    print("Model loaded. Making predictions...")
    predictions = predict(model, args.input, args.epic_annotation)
    print("Predictions made. Writing output...")
    write_predictions(predictions, outfilename + ".txt")
    plot_predictions(predictions, outfilename + ".png")
    if args.input.endswith('.cov') or args.input.endswith('.cov.gz'):
        print("Input file is .cov, running MGMT methylation status analysis...")
        get_mgmt_status(args.input, outfilename, args.reference)
    else:
        print("MGMT methylation status analysis is only available for .cov input files.")
    print("All done.")

if __name__ == "__main__":
    main()
