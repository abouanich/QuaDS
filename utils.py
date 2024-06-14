#!interpreter [optional-arg]
# -*- coding: utf-8 -*-

"""
Helper functions and useful other pieces of code

:project: DIVIS: Biological ontology Integration and Visualisation on RoseData
:author: Julie Bourbeillon
:author: Emmanuel Benoît
:copyright: Copyright 2020, INRAe http://www.inrae.fr
:credits: Julie Bourbeillon
:license: CeCILL \
See the LICENSE file in the project's top-level directory for details.
:version: 1.0
:maintainer: Julie Bourbeillon
:email: julie.bourbeillon@agrocampus-ouest.fr
:status: Beta
"""

################################################################################
# Import Libraries
################################################################################

# Data analysis and manipulation library
import pandas as pd

# Logging library
import logging

# System library to manipulate the file system
from os import path

################################################################################
# Helper functions and useful other pieces of code
################################################################################


def write_excel(
    filename: str, dataframe: pd.DataFrame, idx: bool = False
):
    """
    Function to write into excel file. Doesn't create new sheets when one of the
    same name already exists and replace it

    Args:
        filename: name of the Excel file
        dataframe: dataframe to write to the sheet
        idx: boolean to indicate if we have to write the line indexes or not
    """
    logging.info(f"-> Writing to Excel: {filename}")
    if not path.isfile(filename):
        logging.info("      -> File doesn't exists - creating")
        writer = pd.ExcelWriter(
            filename, engine="xlsxwriter", datetime_format="mmmm dd yyyy"
        )
        dataframe.to_excel(writer, na_rep="", index=idx)
        writer.save()
        return

def scale_df(df: pd.DataFrame, maxi: int) -> pd.DataFrame:
    """
    Function used to scale the values in a dataframe in the range provided as
    parameter
    Args:
        df: the dataframe to scale
        range: the range in which we want to scale the values
    Returns:
         the scaled dataframe
    """

    multiplier = maxi / df.max().max()
    new_df = df.multiply(multiplier)
    return new_df
