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
from tqdm import tqdm
from urllib.error import HTTPError, URLError


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
    name = objname.split(" ")[1]
    filename = f"{prefix}{name}_{survey}.{format}"
    return filename


def process_unwise(fname, band="w1"):
    # Remove axis from unWISE files that contains the band info
    banddict = {"w1": 0, "w2": 1}

    hdu = fits.open(fname)
    imdata = hdu[0].data[banddict[band]]

    # Fix header
    nhead = hdu[0].header.copy()
    nhead["NAXIS"] = 2
    nhead["BAND"] = band

    delkeys = ["NAXIS3", "BANDS", "BAND0", "BAND1"]
    for key in delkeys:
        del nhead[key]

    newhdu = fits.PrimaryHDU(imdata, header=nhead)
    nhl = fits.HDUList(newhdu)
    nhl.writeto(fname, overwrite=True)


def process_vlass_image(infile, outfile, ext=0, scale_unit=True, sfactor=1000):
    # Process the image to hdulist len==1 and 2D WCS header
    # Not needed if obtained from legacy survey
    hdu = fits.open(infile)
    data = hdu[ext].data.squeeze()
    header = hdu[ext].header

    if scale_unit:
        data = sfactor * data

    # Fix header to 2D
    hkeys = list(header.keys())

    crkeys = ["CTYPE", "CRVAL", "CDELT", "CRPIX", "CUNIT"]
    cr3 = [f"{c}3" for c in crkeys]
    cr4 = [f"{c}4" for c in crkeys]
    badkeys = cr3 + cr4 + ["NAXIS3", "NAXIS4"]

    for key in hkeys:
        if "PC3" in key or "PC4" in key or "_3" in key or "_4" in key:
            badkeys.append(key)
        if key in badkeys:
            del header[key]

    header["NAXIS"] = 2

    # Write fits file
    newhdu = fits.PrimaryHDU(data)
    newhdu.header = header
    nhlist = fits.HDUList(newhdu)
    nhlist.writeto(outfile)


def grab_vlass_unwise_cutouts(
    target_file, output_dir=None, vlass_dir="", unwise_dir="", **kwargs,
):
    if output_dir is not None:
        vlass_dir = output_dir
        unwise_dir = output_dir

    # Download VLASS cutouts
    grab_cutouts(
        target_file, output_dir=vlass_dir, survey="vlass1.2", suffix="VLASS", **kwargs
    )

    # Download unWISE cutouts
    grab_cutouts(
        target_file,
        output_dir=unwise_dir,
        survey="unwise-neo4",
        suffix="unWISE_NEO4",
        extra_processing=process_unwise,
        extra_proc_kwds={"band": "w1"},
        **kwargs,
    )


def grab_cutouts(
    target_file: Union[str, pd.DataFrame],
    name_col: str = "Component_name",
    ra_col: str = "RA",
    dec_col: str = "DEC",
    survey: str = "vlass1.2",
    output_dir: str = "",
    imgsize_arcmin: float = 3.0,
    imgsize_pix: int = 500,
    prefix: str = "",
    suffix: str = "",
    extra_processing: Optional[Callable] = None,
    extra_proc_kwds: Dict[Any, Any] = dict(),
) -> None:
    """Function to download image cutouts from any survey.

    Arguments:
        target_file {str, pd.DataFrame} -- Input file or DataFrame containing the list of target 
                                           coordinates and names.

    Keyword Arguments:
        name_col {str} -- The column name in target_file that contains the desired file name 
                         (default: {"Component_name"})
        ra_col {str} -- RA column name (default: {"RA"})
        dec_col {str} -- Dec column name (default: {"DEC"})
        survey {str} -- Survey name to pass to the legacy server (default: {"vlass1.2"})
        output_dir {str} -- Output path for the image cutouts (default: {""})
        prefix {str} -- Prefix for the output filename (default {""})
        suffix {str} -- Suffix for the output filename (default {survey})
        imgsize_arcmin {float} -- Image angular size in arcminutes (default: {3.0})
        imgsize_pix {int} -- Image size in pixels (default: {500})
    """
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

        url = make_url(
            ra=a, dec=d, survey=survey, s_arcmin=imgsize_arcmin, s_px=imgsize_pix
        )
        outfile = os.path.join(
            output_dir, make_filename(objname=name, prefix=prefix, survey=suffix)
        )
        if not os.path.exists(outfile):
            status = download_url(url, outfile)
            if status and (extra_processing is not None):
                extra_processing(outfile, **extra_proc_kwds)


"""
    # Download VLASS
    vurl = make_url(
        ra=a, dec=d, survey="vlass1.2", s_arcmin=imgsize_arcmin, s_px=imgsize_pix
    )
    vfile = os.path.join(vlass_dir, make_filename(objname=name, survey="VLASS"))
    if not os.path.exists(vfile):
        status = download_url(vurl, vfile)

    # Download unWISE
    uwurl = make_url(
        ra=a, dec=d, survey="unwise-neo4", s_arcmin=imgsize_arcmin, s_px=imgsize_pix
    )
    uwfile = os.path.join(
        unwise_dir, make_filename(objname=name, survey="unWISE-NEO4")
    )
    if not os.path.exists(uwfile):
        status = download_url(uwurl, uwfile)
        if status:
            process_unwise(fname=uwfile, band="w1")
"""


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
    """
    Parse input arguments
    """
    parser = argparse.ArgumentParser(description="Download VLASS and unWISE cutouts.")
    parser.add_argument("target_file", help="Path to the desired SOM to annotate")
    parser.add_argument(
        "-p",
        "--path",
        dest="path",
        help="Directory for output files",
        default=None,
        type=str,
    )
    parser.add_argument(
        "--vlass_path",
        dest="vlass_dir",
        help="Directory for output VLASS files",
        default="",
        type=str,
    )
    parser.add_argument(
        "--unwise_path",
        dest="unwise_dir",
        help="Directory for output unWISE files",
        default="",
        type=str,
    )
    parser.add_argument(
        "--img_size",
        dest="img_size",
        help="Image size in pixels",
        default=500,
        type=int,
    )
    parser.add_argument(
        "--ang_size",
        dest="ang_size",
        help="Image angular size in arcminutes",
        default=3.0,
        type=float,
    )
    parser.add_argument(
        "--name_col",
        dest="name_col",
        help="Name of the column containing the file name",
        default="Component_name",
        type=str,
    )
    parser.add_argument(
        "--ra_col",
        dest="ra",
        help="Name of the RA column in the input catalogue",
        default="RA",
        type=str,
    )
    parser.add_argument(
        "--dec_col",
        dest="dec",
        help="Name of the Dec column in the input catalogue",
        default="DEC",
        type=str,
    )

    args = parser.parse_args()
    return args


if __name__ == "__main__":

    args = parse_args()

    if args.path is not None:
        args.vlass_dir = args.path
        args.unwise_dir = args.path

    grab_vlass_unwise_cutouts(
        args.target_file,
        name_col=args.name_col,
        ra_col=args.ra,
        dec_col=args.dec,
        vlass_dir=args.vlass_dir,
        unwise_dir=args.unwise_dir,
        imgsize_arcmin=args.ang_size,
        imgsize_pix=args.img_size,
    )
