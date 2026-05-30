# Setup guide

The hands-on sessions of the 2026 SSCHA School will be run on a virtual machine, with the SSCHA software pre-installed.
It is important that all the partecipants have the VM ready and working **before** the beginning of the school, 
as during the limited time available for the school we will not go through the software installation.

In the following we detail all the step required to be ready to run the tutorials and partecipate to the hands-on session

## The virtual machine

To run the tutorials, we recommend using a custom version of the Quantum Mobile virtual machine.
The download links for Intel or Apple Silicon chips are the following:

1. [Intel x86/AMD64](https://drive.google.com/file/d/12hFWtAFn83wqAwF9H6ClXbb9adOQfpeb/view?usp=sharing)
2. [Apple Silicon/ARM64](https://drive.google.com/file/d/1Xqogj1__GRIq2heB92s_2ff3QUF1nAZG/view?usp=sharing)

(A huge thanks to Marnik Bercx, Nicola Colonna and all the people involved from Quantum Mobile for building this version for the ICTP/MARVEL College: Materials simulations in the age of AI).

To start the virtual machine, you must import it inside Oracle VirtualBox. 

**NOTE: VirtualBox version 7.2 at least is required**. The virtual machine will not work on old VirtualBox versions.

Run until the machine starts. Note that running the VM for the first time requires a lot of free space on the disc (up to about 30 GB).

## Installing required force fields

Once inside the virtual machine, open a terminal.
Every tutorial will be run inside the SSCHA conda environment. To access it use

```bash
conda activate sscha
```

This needs to be done every time you open a new terminal. A `(sscha)` text should appear at the beggining of the prompt line.
Then install the force-fields required for running the tutorials:

```bash
pip install quippy-ase F3ToyModel
```

You are now ready for the hands-on sessions 1-4.


## EPIq installation 



