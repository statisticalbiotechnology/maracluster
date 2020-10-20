This repository holds the source code for clustering and merging of fragment spectra from shotgun proteomics experiments.

## MaRaCluster

### Installation

Installers for several operating systems can be found on the [Release page](https://github.com/statisticalbiotechnology/maracluster/releases). Both a command line as well as a graphical user interface are available for all platforms.

If you prefer to compile from source, or are running on a different operating system, [click here](#installation-from-source).

### Example

An example, including spectrum files and a shell script, can be downloaded from http://kaell.org/files/maracluster_sample.zip.

### Graphical User Interface (GUI)

To run the example, simply select the spectrum files after clicking the `Add ms2 files` button, select an output folder and click the `Run MaRaCluster` button.

A more elaborate guide for using the GUI can be found from the Help menu in the menu bar of the GUI.

### Command Line Interface (CLI)

The main functionality of MaRaCluster is provided by the `maracluster` command. This command has several sub-commands to execute different parts of the clustering and merging. The most important ones are `maracluster batch` and `maracluster consensus`. The first takes a list of ms2 spectra as input and outputs a list of clusters, the second takes one of these outputs and creates consensus spectra for each cluster.

To run `maracluster batch`, a flat text file with the absolute/relative path to each of the ms2 spectrum files (one per line) is needed. Such a file can easily be generated using a `ls -1` command, *e.g.* `ls -1 ms2/* > files.txt` on Linux or `dir <path> /b /s > files.txt` on Windows. Any ms2 spectrum format readable by ProteoWizard can be used as input. Use the following command to start clustering:
```
maracluster batch -b files.txt
```
This will create several files called `MaRaCluster.clusters_p<x>.tsv` in a subdirectory called `maracluster_output`, for a range of p-value thresholds `10e-<x>`. These output files contain one spectrum per line, with different clusters separated by an empty line. The spectrum is listed with the path to the spectrum file in the first column, the unmerged scannr (or scan index if no scannr is available) in the second column and the cluster index in the third column, all separated by tabs.

To run `maracluster consensus`, we take one of the cluster files as input, e.g.:
```
maracluster consensus -l maracluster_output/MaraCluster.clusters_p10.tsv
```

For more information and options run `maracluster -h` on the command line.

### Installation from source

To install MaRaCluster, you can use the provided installation script `./quickbuild.sh`, which will build the package in `./bin/build`, and install the executables in the `/usr/bin` folder (needs superuser rights). If you do not have superuser rights, or want to install the executable somewhere else, modify the script accordingly by setting the `-DCMAKE_INSTALL_PREFIX` flag to the desired location, and change the last line from `sudo make install` to `make install`.
