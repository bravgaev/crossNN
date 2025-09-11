#! /usr/bin/env python3
# CrossNN local script
import NN_model
import argparse


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
        for pred in predictions[0]:
            for class_label in predictions[1]:
                f.write(f"{class_label}\t{pred}\n")
        

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
    plt.barh(y_pos, top_scores, align='center', alpha=0.7)
    plt.yticks(y_pos, top_labels)
    plt.xlabel('Scores')
    plt.title('Top 10 Prediction Scores by Class (number of features used: {})'.format(num_features))


    plt.savefig(plot_file_path)
    plt.close()

def main():
    parser = argparse.ArgumentParser(description="Run CrossNN model")
    parser.add_argument("--model", required=True, help="Path to the model file, created with the training.py script")
    parser.add_argument("--input", required=True, help="Path to the sample bedMethyl file")
    parser.add_argument("--output", required=True, help="Path to the output file")
    parser.add_argument("--epic_annotation", required=True, help="Path to the EPIC annotation file")
    args = parser.parse_args()

    model = load_model(args.model)
    predictions = predict(model, args.input, args.epic_annotation)
    write_predictions(predictions, args.output)
    plot_predictions(predictions, args.output + ".png")


if __name__ == "__main__":
    main()
