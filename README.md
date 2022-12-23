# Aoki_ICIrepertoire

Codes for reproducing results of ICI repertoire paper (Aoki et al., 2022).

Code for processing data is "221202_ICIrepertoire_preprocess_Fig3.ipynb."
Code for plotting figure is "221203_ICIrepertoire_plotting_Fig3.ipynb."

Raw data for single-cell RNAseq is stored in the directory "raw.data."
Raw data for single-cell TCRseq is stored in the directory "scTCR_rawdata."
Raw data for bulk TCRseq is stored in the directory "BulkTCR_rawdata."

Seurat object for creating figure is "Aoki_3rd_res0.6scTCRmerged.ABnt.AddSig.pseudotime.rda."

These codes are run in the Docker container haokimriid/singlecell
Please install the docker image by the following command.
```
docker pull haokimriid/singlecell
```

Launch Jupyter from the Docker image and work on Jupyter. 
For, Windows10 Pro, run the following command.
```
docker run --rm -p 8888:8888 -v /mnt/your-working -directory/:/opt/work haokimriid/singlecell jupyter notebook --allow-root
```
Then, Copy the http://127.0.0.1:8888/?token=~ part (URL) and launch it in your browser. 
