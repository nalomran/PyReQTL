from pathlib import Path, PurePath
from typing import List


def remove_vars(df,
                r_val: int,
                num_samples: int,
                cols: List[str]):
    """Removes variants that are homozygous variant or references in more that
    80% of the samples

    Parameters
    -----------
    df: data frame of the variants from readCounts

    r_val: number of reads corresponding to homozygous or homozygous reference

    num_samples: number of unique samples

    cols: list of the new columns names which will be created

    Returns
    -------
     a dataframe that meet the condition above

    """

    df[cols[0]] = (df['R'] == r_val).groupby(df['SNV']).transform('sum')
    df[cols[1]] = df[cols[0]] / num_samples

    df = df[df[cols[1]] < 0.8].reset_index()

    return df


def create_output_dir(dir_path: str) -> str:
    """Creates output directory if it doesn't exits

    Parameters
    -----------
    dir_path: the output file directory

    Return
    -------
    dir_path: the output file directory returned for later used


    """

    if not Path(dir_path).exists():
        print(f"creating output directory {dir_path}...\n")
        Path(dir_path).mkdir(parents=True, exist_ok=True)

    return dir_path


def extract_filename_tag(file_path: str) -> str:
    """extract file name tag alone without file extension

    Parameters
    ----------
    file_path: the path with the filename

    Return
    -------
    filename: str, the extracted filename


    """

    return Path(file_path).resolve().stem


def output_filename_generator(out_dir: str,
                              prefix: str,
                              file_suffix: str) -> str:
    """Concatenate output filename with file directory and file prefix pattern

    Parameters
    -----------
    out_dir: str, the directory of the file output

    prefix: the prefix of the file output

    file_suffix: the suffix of file output


    Return
    -------
    the concatenated output filename

    """

    return str(PurePath(out_dir, prefix + file_suffix))


def bool_conv_args(args: str) -> bool:
    """Convert string argument value to boolean True or False

    Parameters
    ----------
    args: argument value to represent True or False


    Return
    -------
    a converted string argument to a boolean value
    """

    # consider most possible truth values scenarios supplied by the user
    if args.lower() in ['yes', 'true', 't', 'y']:
        return True

    elif args.lower() in ['no', 'false', 'f', 'n']:
        return False

    else:
        raise argparse.ArgumentParser('Please make sure to enter a boolean '
                                      'value i.e. True or False.')