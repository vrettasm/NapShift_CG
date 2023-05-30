# NapShift (version 0.0.1)

![Logo](./logos/NapShift_CG_logo.png)

This repository provides a "Python implementation" of the NapShift-CG (CG: Coarse Grain)
program to estimate the backbone atoms' chemical shift values from martinized PBD files.
It is based on the already published and tested NapShift program that works on full atomistic
proteins. For more information have a look at: https://github.com/vrettasm/NapShift.git.

M. Vrettas, PhD.

## Installation

There are two options to install the software.

1. The easiest way is to visit the GitHub web-page of the project and
[download the code](https://github.com/vrettasm/NapShift_CG/archive/master.zip)
in zip format. This option does not require a prior installation of git on the
computer.

2. Alternatively one can clone the project directly using git as follows:

    `$ git clone https://github.com/vrettasm/NapShift_CG.git`

## Required packages

The minimum version is **Python 3.7** (recommended >=3.8). The required packages
are given in the "requirements.txt": To simplify the installation of the packages just use:

    $ pip install -r requirements.txt

## Virtual environment (recommended)

It is highly advised to create a separate virtual environment to avoid
messing with the main Python installation. On Linux and macOS systems
type:

    $ python3 -m venv napshift_cg_venv

Note: "napshift_cg_venv" is an _optional_ name.

Once the virtual environment is created activate it with:

    $ source napshift_cg_venv/bin/activate

Make sure **pip** is updated:

    $ python3 -m pip install --upgrade pip

Then we can install all the requirements as above:

    $ pip install -r requirements.txt

or

    $ python3 -m pip install -r requirements.txt

N.B. For Windows systems follow the **equivalent** instructions.

## How to run

To execute the program (within the activated virtual environment), you can either
navigate to the main directory of the project (i.e. where the napshift_cg.py is located),
or locate it through the command line and then run the following command:

    $ ./napshift_cg.py -f path/to/filename.pdb

This is the simplest way to run NapShift. It will create a file named:
"prediction_filename_model_0_chain_A.tab" in the _current working directory_,
with the predicted chemical shift values for all backbone atoms (N, C, Ca, Cb, H, Ha).

   > **Hint**: If you want to run the program on multiple files (in the same directory) you
   > can use the '*' wildcard as follows:
   >  
   > $ ./napshift_cg.py -f path/to/*.pdb

This will run NapShift on all the files (in the directory) with the '.pdb' extension.

---

To explore all the options of NapShift, use:

    $ ./napshift_cg.py -h

You will see the following menu:

![Help](./logos/help_menu.png)

## References (and documentation)

The original work is described in detail at:

1. TBD.

### Contact

For any questions/comments (**_regarding this code_**) please contact me at:
vrettasm@gmail.com