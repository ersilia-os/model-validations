This folder contains dataset that would be used to check for bias in models currently in the hub. It was gotten from two databases, the <a href="https://coconut.naturalproducts.net/download">Coconut database</a> and <a href="https://github.com/ersilia-os/groverfeat">Grover</a>. 800 molecules was randomly sampled from Grover and 200 from Coconut. These molecules were standardized using the Standardiser package. They were processed to contain column of the SMILES, InchiKeys and Canonical SMILES

Validations.csv: This contains the final dataset that has been processed. It contains 1000 rows of data from both databases

Coconut_processed200: This contains the 200 randomly sampled Grover data

Grover.csv: This contains the 800 randomly sampled Grover data

Coconut_data.csv: The size of the original data gotten from Coconut was large > 1.7GB, in order to work with a smaller size, we created a new data from the original containing just the SMILES, InchiKey and molecular weight. If you're interested in obtaining the original data, it can be found here https://coconut.naturalproducts.net/download

reference_library.csv: This contains all the datasets gotten from the Grover database




