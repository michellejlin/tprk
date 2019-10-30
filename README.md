# tprK pipeline
This pipeline was designed to take Illumina and PacBio files straight off the sequencer to a final comparison table of all the different variable regions with their relative frequencies, as well as various pretty plots along the way.

## Table of Contents
* [Setup](#Setup)
* [Input Files](#Input-Files)
* [Usage](#Usage)
* [Arguments](#Arguments)
* [Running Parts Separately](#Running-Parts-Separately)

## Setup
**Installs**

Make sure you have a working installation of [Python](https://www.python.org/) and [Julia](https://julialang.org/downloads/). 

**Setting It Up**

This next section is a step-by-step breakdown of how to get this pipeline set up. If you know how to alias `tprk_pipeline.py`, you can skip this section entirely **except #9**. Otherwise, read on!

1. Download the tprK repository. If you have git installed on the command line you can open up a terminal window and type `git clone https://github.com/michellejlin/tprK.git`. Or if you're not comfortable doing that you can click the green button that says Clone or download and click on download zip. Then unzip the folder to wherever you want on your computer. 

Now we have to make sure the script can be run from anywhere. Navigate to your main tprK folder.

2. Open up a terminal window and change directories to the main tprK folder. If you downloaded the zip file to your downloads folder this would look something like: `cd /Users/username/Downloads/tprK-master/` 

3. Type in `pwd`. This will give you the current path. Copy this path.
4. Type in `nano ~/.bashrc`. For MAC OSX users, replace `bashrc` with `bash_profile`, keeping the rest of the punctuation.
5. This will open up an editor in the terminal window. Scroll down to get to the bottom of this file.
6. Type in ``alias tprk_pipeline.py="python INSERTPATHHERE/tprk_pipeline.py"`` into the terminal, where you replace INSERTPATHHERE with the path you just copied. An example might look like this: `alias tprk_pipeline.py="python /Users/uwvirongs/Documents/tprK-master/tprk_pipeline.py"`.
7. Hit Ctrl+X to quit out of the editor. Save when prompted (this may mean typing 'Y') and press Enter for any new prompts that show up.
8. Type in `source ~/.bashrc` (or replace `bashrc` with `bash_profile` for Mac OSX), to refresh the file.

Now we have to change the path of the Julia install in one of our files. 
9. **Do not skip this step, people who are skimming!** Open up `og_files_to_all_reads.R` in your favorite text editor. Change the path on line 3 to point to the path where your Julia install is. Here's what the beginning should look like:

```python
##This line must be set up for PacBio file processing!
##In a typical Mac installation, this path points to the Julia application in the Application folder.
julia <- julia_setup(JULIA_HOME = "/Applications/Julia-1.2.app/Contents/Resources/julia/bin/")
```


That's it! Once you've changed the path, you're ready to do some tprK analysis!

## Input Files
Put the following things in one folder:
- **All the sequence files to run analysis on**
    - PacBio Q20 files
    - Illumina files trimmed and run through Trimmomatic
    - By default the pipeline expects both PacBio and Illumina files for every sample. Running the pipeline with just PacBio or just Illumina files is possible with the `-pacbio` and `-illumina` flags respectively. However, some plots require both files to be generated and these plots will not be output.
- **Metadata file.** This should be a .csv with three columns: Sample Name, PacBio, Illumina, shown in the table below.
    - By default, the pipeline will look for a file named metadata.csv. You can rename and specify the name of the metadata file through the `-m` flag.
    - An example metadata file is provided. The general format of the metadata file should be three columns, separated by commas, as shown:

| Sample Name  | PacBio  | Illumina |
| ------------- | ------------- | ------------- |
| This will largely be the name used for generating tables and plots. | The PacBio file specified for the sample name. This must match exactly the name of the matching file in the folder. This should be a Q20 file.  | The Illumina file specified for the sample name. This must match exactly the name of the matching file in the folder. This should be a trimmed file run through Trimmomatic.|

## Usage
With setup done and input files ready, we can move on to actually using the pipeline!

This can be done in a few simple steps:
1. Open up terminal.
2. Navigate to the folder with your sequencer files and metadata.csv file using `cd`, i.e. `cd Users/uwvirongs/Documents/tprk-master/`.
3. Run the code! `tprk_pipeline.py metadata.csv`

Of course, that just runs the pipeline from start to finish, with completely default setup. To further specify arguments, refer to the [Arguments](#Arguments) section.

Additionally, you might not want to run the whole pipeline from start to finish every time. To see what the individual components of the pipeline do and how to run them separately, refer to [Running Parts Separately](#Running-Parts-Separately).

## Arguments
| Command | Description |
| --- | --- |
| `-m`,`--metadata_file` | Specify name of metadata.csv file containing sample name, Illumina, and PacBio files. By default, the pipeline will search for `metadata.csv`. |
| `-f`, `--relative_freq_filter` | Specify by what relative frequency an additional filtered final merged table and visualizations should be sorted at. By default, this is set to 0.2. |
| `-c`, `--count_filter` | Specify by what count an additional filtered final merged table and visualizations should be sorted at. By default, this is set to 5. |
| `-i`, `--illumina_filter` | Specify if PacBio reads should only include Illumina-supported reads that pass the filters given. By default, relative freq is set to 0.2 and count is set to 5. |
| `-pacbio` | Write this flag to specify that there are only PacBio files here. Comparison figures to Illumina will not be created. | 
| `-illumina` | Write this flag to specify that there are only Illumina files here. Comparison figures to PacBio will not be created. | 

## Running Parts Separately
You might not want to run the entire tprK pipeline from start to finish, every time. Here's a guide on what each component of this pipeline does and how to run them. 

In general, this information will be located near the top of each individual file in a commented out section. 
