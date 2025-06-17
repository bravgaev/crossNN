
## Description
crossNN is an explainable machine learning framework for cross-platform DNA methylation-based classification of cancer. 

## Model training

#### First, install dependencies
>`pip install -r requirements.txt`

crossNN has been developed and tested under standard Linux (CentOS, Ubuntu).

#### Run the training with custom parameters
> `python training.py --data your_data.h5 --epochs 1000 --learning_rate 0.0001 --weight_decay 0.00001 --mask_size 0.0025 --model_output ./model_output.pth --pickle_output ./model_pickle.p`

#### Parameters
> `--data`  Training data in .h5  
> `--epochs`    Number of training epochs  
> `--learning_rate` Learning rate of the model   
> `--weight_decay` Weight decay of the L2 norm of the weights  
> `--mask_size` Masking rate p  
> `--model_output` Trained Pytorch Model    
> `--pickle_output` Trained model with label encoder in Pickle file  

## Inference

We currently provide two pre-trained models:
- brain tumor model (based on Heidelberg v11b4 reference set)
- pan-cancer model (Berlin v5i assembly)

In addition to source code provided here, a sample implementation with a graphical user interface to crossNN models is available at [crossnn.charite.de](https://crossnn.charite.de).   
Inference using crossNN is also implemented in [nanoDx](https://gitlab.com/pesk/nanoDx), an end-to-end analysis pipeline for nanopore low-pass whole genome sequencing data.

## Citation
If you use the crossNN in your research, _please consider citing the following paper_ in your work:

Yuan, D., Jugas, R., Pokorna, P., Sterba, J., Slaby, O., Schmid, S., Siewert, C., Osberg, B., Capper, D., Halldorsson, S., Vik-Mo, E.O., Zeiner, P.S., Weber, K.J., Harter, P.N., Thomas, C., Albers, A., Rechsteiner, M., Reimann, R., Appelt, A., Schüller, U., Jabareen, N., Mackowiak, S., Ishaque, N., Eils, R., Lukassen, S., Euskirchen, P. (2025) _crossNN is an explainable framework for cross-platform DNA methylation-based classification of tumors_. Nat Cancer 1–12. https://doi.org/10.1038/s43018-025-00976-5


