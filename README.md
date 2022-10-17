# Fitness seascapes are necessary for accurately predicting evolutionary response to drug therapy
Reproduce figures from paper (https://www.biorxiv.org/content/10.1101/2022.06.10.495696v1.full.pdf)
Steps to reproduce figures: 

1. clone this repository: git clone https://github.com/eshanking/seascapes_figures/
2. Install fears (pip install git+https://github.com/eshanking/fears.git)
3. Open and run desired experiment file (experiments/experiment_file.py)
4. Open figure_code and run desired figure function (make_fig) with experiment info path as the path to the relevant experiment info file (found in the results folder)
5. Figures are saved to the figures folder.

paper_figure_reproduction.py contains code to reproduce all of the experimental figures from the paper. This code will take several hours to run and will require ~2 GB of disk space.