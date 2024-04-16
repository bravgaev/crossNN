import argparse
import pickle

import h5py
import numpy as np
import pandas as pd
import torch
import torch.utils.data as data_utils
from sklearn.preprocessing import LabelEncoder
from torch import nn, optim
from torch.utils.data import DataLoader

class DNN(nn.Module):
    """Our proposed shallow neural network model."""

    def __init__(self,n_input, n_output):
        super(DNN, self).__init__()
        self.layer_out = nn.Linear( n_input,n_output ,bias=False)

    def forward(self, x):
        return self.layer_out(x)


def read_data(data_pwd: str) :
    #-> tuple(list, list, list, np.ndarray)
    """Read data from h5.file. h5 file can be generated from nanotrainingset repository."""
    fp = h5py.File(data_pwd)
    label = list(fp["Dx"][:])
    probeIDs = list(fp["probeIDs"][:])
    sampleIDs = list(fp["sampleIDs"][:])
    data = fp["betaValues"]
    probeIDs = [str(x).replace("b", "").replace("'", "") for x in probeIDs]
    label = [str(x).replace("b", "").replace("'", "") for x in label]
    # 

    return label, probeIDs, sampleIDs, data


def preprocessing(data: np.ndarray, probeIDs: list) -> pd.DataFrame:
    """Preprocess the data by binarizing the beta values and removing uniformed cpg features."""

    # load training data and preprocess
    x = data[:, :]

    # binarinize the beta values
    x = np.where(x > 0.6, 1, -1)
    x = pd.DataFrame(x.T, columns=probeIDs)

    # remove uniformed cpg features
    x = x.filter(x.columns[x.sum() != 2801])
    x = x.filter(x.columns[x.sum() != -2801])
    return x


def mask_input(val, mask_size):
    """Mask input data by randomly setting a subset of features to 0."""
    # masked values = 0

    masked_fea = np.random.choice(
        range(0, val.size()[1]),
        size=int((1 - mask_size) * val.size()[1]),
        replace=False,
        p=None,
    )
    val[:, np.tile(masked_fea, (val.size()[0], 1))[0]] = 0.0
    return val


def model_training(
    DM: nn.Module,
    train_All_loader: DataLoader,
    device: torch.device,
    epochs=1000,
    learning_rate=0.0001,
    weight_decay=0.0001,
    mask_size=0.0025,
):
    """Train the model."""
    criterion = nn.CrossEntropyLoss()

    optimizer = optim.Adam(DM.parameters(), lr=learning_rate, weight_decay=weight_decay)
    print("Begin training.")
    print(device)
    for _ in range(epochs):
        train_epoch_loss = 0
        DM.train()
        print(_)
        for X_train_batch, y_train_batch in train_All_loader:
            X_train_batch, y_train_batch = X_train_batch, y_train_batch.to(device)
            optimizer.zero_grad()
            y_train_pred = DM(mask_input(X_train_batch, mask_size=mask_size).to(device))

            train_loss = criterion(y_train_pred, y_train_batch)

            train_loss.backward()
            optimizer.step()

            train_epoch_loss += train_loss.item()
    return DM

def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Train a model on the Brain Tumor dataset."
    )
    parser.add_argument(
        "--data",
        type=str,
        default="data/Capper_et_al.h5",
        help="Path to the data file.",
    )
    parser.add_argument("--epochs", type=int, default=1000, help="Number of epochs.")
    parser.add_argument(
        "--learning_rate", type=float, default=0.0001, help="Learning rate."
    )
    parser.add_argument(
        "--weight_decay", type=float, default=0.0001, help="Weight decay."
    )
    parser.add_argument("--mask_size", type=float, default=0.0025, help="Mask size.")
    parser.add_argument("--batch_size", type=int, default=32, help="Batch size.")
    parser.add_argument(
        "--model_output",
        type=str,
        default="models/brain_tumor.pth",
        help="Path to the output file.",
    )
    parser.add_argument(
        "--pickle_output",
        type=str,
        default="models/trained_BrainTumor.pkl",
        help="Path to the output file.",
    )
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_arguments()
    
    # Load dataset
    label, probeIDs, sampleIDs, data = read_data(args.data)    
    
    # preprocess data
    X = preprocessing(data, probeIDs)
    print(X.shape)
    train_data_all = torch.tensor(X.to_numpy())
    train_data_all = train_data_all.float()
    
    
    
    # preprocess labels
    enc = LabelEncoder()
    enc.fit(label)
    enc_label = enc.transform(label)
    enc_label = torch.from_numpy(enc_label).type(torch.LongTensor)

    X = preprocessing(data,probeIDs)
    
    # Build training datasets
    train_all = data_utils.TensorDataset(train_data_all, enc_label)
    train_All_loader = DataLoader(
        train_all, batch_size=args.batch_size, drop_last=True, shuffle=True
    )
    # Build model
    DM = DNN(n_input = len(X.columns), n_output=len(set(label)))
    device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")    
    
    # Train the model

    model = model_training(
        DM,
        train_All_loader,
        device,
        args.epochs,
        args.learning_rate,
        args.weight_decay,
        args.mask_size,
        
    )
    
    # only save model
    checkpoint = model.state_dict()
    torch.save(checkpoint, args.model_output)

    # save pickle file containing model, enc(label encoder), example_bed(including features used in the model)
    example_bed = pd.DataFrame({"probe_id": X.columns.tolist()})
    with open(args.pickle_output, "wb") as file:
        pickle.dump([checkpoint, enc, example_bed], file)

        
        