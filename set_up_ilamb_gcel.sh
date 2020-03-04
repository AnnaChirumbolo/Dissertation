#!/bin/bash -l

# to load the ilamb environment, run this script with
# . set_up_ilamb_gcel.sh

ILAMB_WORKING_DIR=$PWD

source /exports/csce/datastore/geos/groups/gcel/SOFTWARE_LIBS/conda/miniconda3_GCEL/etc/profile.d/conda.sh

# load lamb environment using conda
echo "Activating ilamb27"
conda activate ilamb27

# set the ilmab root directory to the current working directory
echo "Setting ILAMB root directory"
export ILAMB_ROOT=$ILAMB_WORKING_DIR

# make sure to load the School's MPI library - ILAMB uses this to process the results in parallel
echo "Loading mpi module"
module load mpi

echo "To check that ilamb is working, type: "which ilamb-run""

echo "To run ilamb, use the command"
echo "ilamb-run --config sample.cfg --model_root $ILAMB_ROOT/MODELS/ --regions global --clean"







