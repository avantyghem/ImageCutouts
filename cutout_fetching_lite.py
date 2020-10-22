"""Script to fetch cutouts from the legacy survey"""

import os
import sys
import time
import numpy as np
import pandas as pd
import argparse
import wget
from typing import Any, Callable, Dict, Mapping, Optional, Union
from astropy.io import fits
from urllib.error import HTTPError, URLError


####################################################################################

def make_url(ra, dec, survey="dr8", s_arcmin=3, s_px=512, format="fits"):
    # Convert coords to string
    ra = str(np.round(ra, 5))
    dec = str(np.round(dec, 5))

    # Set pixscale
    s_arcsec = 60 * s_arcmin
    pxscale = s_arcsec / s_px

    # Convert image scales to string
    s_px, pxscale = str(s_px), str(np.round(pxscale, 4))

    url = (
        f"http://legacysurvey.org/viewer/cutout.{format}"
        f"?ra={ra}&dec={dec}"
        f"&layer={survey}&pixscale={pxscale}&size={s_px}"
    )
    return url


def make_filename(objname, prefix="", survey="DECaLS-DR8", format="fits"):
    # Take Julian coords of name to eliminate white space - eliminate prefix
    if len(objname.split(" "))>1:
        name = objname.split(" ")[1]
    else:
        name = objname
    filename = f"{prefix}{name}_{survey}.{format}"
    return filename


#####need to make arguments here configurable from config file
###target file is code arg
###name_col, ra_col, dec_col, survey, image_size(arcmin), output directory as config parms
###add in image band if only want e.g. r band from sdss
def grab_cutouts(
    target_file: Union[str, pd.DataFrame],
    config_file: str = '',
    imgsize_pix: int = 512,
    prefix: str = "",
    suffix: str = "",
    extra_processing: Optional[Callable] = None,
    extra_proc_kwds: Dict[Any, Any] = dict(),
) -> None:
    """Function to download image cutouts from any survey.

    Arguments:
        target_file {str, pd.DataFrame} -- Input file or DataFrame containing the list of target coordinates and names.
        config_file {str} -- contains config info such as column names in target_file, survey to get info from, band, cutout size etc.
    """
    config = pd.read_table(config_file,delim_whitespace=True).replace("'", "",regex=True)
    
    name_col = config['value'][np.where(config['parameter']=='name_col')[0][0]]
    ra_col = config['value'][np.where(config['parameter']=='ra_col')[0][0]]
    dec_col = config['value'][np.where(config['parameter']=='dec_col')[0][0]]
    survey = config['value'][np.where(config['parameter']=='survey')[0][0]]
    imgsize_arcmin = config['value'][np.where(config['parameter']=='size_arcmin')[0][0]]
    ###convert imgsize to float
    imgsize_arcmin = float(imgsize_arcmin)
    output_dir = config['value'][np.where(config['parameter']=='outdir')[0][0]]
    band = config['value'][np.where(config['parameter']=='band')[0][0]]
    
    if isinstance(target_file, str):
        targets = pd.read_csv(target_file)
    else:
        targets = target_file

    if suffix == "":
        suffix = survey

    for _, target in targets.iterrows():
        name = target[name_col]
        a = target[ra_col]
        d = target[dec_col]

        outfile = os.path.join(
            output_dir, make_filename(objname=name, prefix=prefix, survey=suffix)
        )
        grab_cutout(
            a,
            d,
            outfile,
            survey=survey,
            imgsize_arcmin=imgsize_arcmin,
            imgsize_pix=imgsize_pix,
            extra_processing=extra_processing,
        )


def grab_cutout(
    ra,
    dec,
    outfile,
    survey="vlass1.2",
    imgsize_arcmin=3.0,
    imgsize_pix=500,
    extra_processing=None,
    extra_proc_kwds=dict(),
):
    url = make_url(
        ra=ra, dec=dec, survey=survey, s_arcmin=imgsize_arcmin, s_px=imgsize_pix
    )
    if not os.path.exists(outfile):
        status = download_url(url, outfile)
        if status and (extra_processing is not None):
            extra_processing(outfile, **extra_proc_kwds)


def download_url(url: str, outfile: str, max_attempts: int = 5):
    # Often encounter the following error:
    # urllib.error.HTTPError: HTTP Error 504: Gateway Time-out
    # Repeat the download attempt for up to `max_attempts` tries
    # Return True if the download was successful
    for attempt in range(max_attempts):
        try:
            wget.download(url=url, out=outfile)
            return True
        except HTTPError as e:
            print(f"Failed attempt {attempt} to download {outfile} with an HTTPError")
        except URLError as e:
            print(f"Failed attempt {attempt} to download {outfile} with a URLError")
        time.sleep(1)

    print(f"Failed to download image {outfile}")
    return False


def parse_args():
    "parse input args, i.e. target and config file names"
    parser = argparse.ArgumentParser(description=
                                     "Download image cutouts from legacysurveys.org")
    parser.add_argument("target_file", help="list of positions [deg] to obtain cutouts for")
    parser.add_argument("--config", action='store', type=str, default='config.txt',
                        help="config file")
                        
    args = parser.parse_args()
    return args

####################################################################################

if __name__ == '__main__':
    args = parse_args()
    grab_cutouts(target_file=args.target_file, config_file=args.config)



