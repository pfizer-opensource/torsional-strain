# Calculation of Torsional Strain Energy

## Setting up the Environment
This step will get you setup with an evironment within which you can run the 
code in this repository using Orion.
1. Get access to OpenEye Orion and Magpie. If you are a licensed Orion user and 
   don't have access, please contact [OpenEye Support](support@eyesopen.com).
2. Activate a conda distribution (typically something like `source $CONDA/activate`)
3. Create an environment
        `conda env create -v -f environment.yaml -n torsion`
4. Activate the environment
        `conda activate torsion`
5. Install OpenEye packages
    a. Insert your magpie token into the first line of `requirements_dev.txt`:
        `--extra-index-url https://token:YOUR_TOKEN@magpie.eyesopen.com/simple/`
    b. Install packages with `pip install -r requirements_dev.txt`
6. Configure the orion CLI 
        `ocli --profile=default config profile`

You should now have an environment capable of running the orion code.

## Uploading the Floe to Orion
1. Activate the environment
2. Run `invoke clean package upload`

This will add a package called `torsional-strain 1.0.0` with a floe called `Torsional Strain Energy` to Orion. 
![Floe Screenshot](floe_screenshot.jpg?raw=true "Floe Screenshot")

The floe requires only one input: a dataset of molecules for which we want to calculate torsional strain.