CondaVer=3
unset PYTHONPATH
MINICONDA_DIR=./miniconda${CondaVer}  # change to wherever you like
if [[ "$OSTYPE" == "linux-gnu" ]]; then
    wget https://repo.continuum.io/miniconda/Miniconda${CondaVer}-latest-Linux-x86_64.sh
    bash Miniconda${CondaVer}-latest-Linux-x86_64.sh -b -p $MINICONDA_DIR
    rm Miniconda${CondaVer}-latest-Linux-x86_64.sh
elif [[ "$OSTYPE" == "darwin"* ]]; then
    wget https://repo.continuum.io/miniconda/Miniconda3-latest-MacOSX-x86_64.sh
    bash Miniconda${CondaVer}-latest-MacOSX-x86_64.sh -b -p $MINICONDA_DIR
    rm Miniconda${CondaVer}-latest-MacOSX-x86_64.sh
fi

cd $MINICONDA_DIR

source etc/profile.d/conda.sh
conda update -y -n base -c defaults conda

export PYTHONNOUSERSITE=true

cd ../

conda create -y -n mapping_env python=${CondaVer} pip
conda activate mapping_env

conda install -y -c conda-forge ROOT=6.20.2 #ROOT installation (use older version to avoid crash)
conda install -y -c conda-forge scikit-learn
conda install -y -c conda-forge pyyaml
conda install -y -c conda-forge pandas
conda install -y -c conda-forge matplotlib
conda install -y -c conda-forge root_numpy
