# ImageCutouts
Fetch image cutouts from the legacy survey.

Command line code for fetching fits format image cutouts from the DESI legacy surveys sky viewer (legacysurvey.org).
Takes a csv file of targets containing position (in decimal degrees) and an object name, and a configurations file as inputs via command line arguments, run as:

*python3 cutout_fetching_lite.py target_list.csv config.txt*

Downloads requested cutouts and places them in specified output folder.



**Config File Description**

parameter | description
----------|------------
name_col | name of column used to ID the object in the target_file, propogates to downloaded file name
ra_col | name of the RA column in the target file (assumed to be in decimal degrees)
dec_col | name of the Dec column in the target file (assumed to be in decimal degrees)
survey | the key (see below) of the image survey required, e.g. 'sdss'
size_arcmin | cutout size in arcmin (i.e. the cutout will be size_arcmin * size_arcmin square)
outdir | directory to put cutouts in


**Available Surveys**

At present the survey provided in the config file needs to match the key used by the cutout urls, these are described in the table below.
Additionally some surveys download multi band images where the fits file obtained is a multi-band cube - e.g. an SDSS cutout contains 3 bands (g,r,i) in a single file.
If your input position is outside of the specified survey footprint an empty file is retrieved.


survey key | description
----------|------------
ls-dr8 | g/r/z cube from the 8th data release of the DESI legacy surveys (DECaLS, BASS, MzLS)
unwise-neo4 | w1/w2 cube from unWISE - NEO6
sdss | g/r/i cube from SDSS
hsc-dr2 | g/r/i cube from the 2nd data release of Hyper Supreme Cam
vlass1.2 | continuum cutout from the VLA sky survey epoch 1 quick look images - although the key is vlass1.2, all of epoch 1 is covered.
galex | nUV/fUV cube from GALEX
des-dr1 | g/r/z cube from DES DR1


The above list of surveys is representative, but incomplete - see legacysurvey.org/viewer for an up to date and complete list.
