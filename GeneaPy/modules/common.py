""" Functions used across multiple scripts"""


def get_ensembl_release(hg_version):
    """ Change hg version to ensembl build release number"""
    if isinstance(hg_version, str):
        hg_ensembl_dict = {"hg19": 75, "hg38": 83}
        hg_version = correct_hg_version(hg_version)
        hg_ensembl = hg_ensembl_dict.get(hg_version.lower())
        return hg_ensembl
    elif isinstance(hg_version, int):
        return hg_version
    else:
        raise TypeError("hg_version must by str or int type")


def correct_hg_version(hg_version):
    """ Change genome build release from GRCH to HG"""
    hg_version_dict = {"grch37": "hg19", "grch38": "hg38"}
    if hg_version.lower().startswith("g"):
        hg_version = hg_version_dict.get(hg_version.lower())
    return hg_version
