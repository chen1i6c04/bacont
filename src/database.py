import os
import sys
import shutil
import tarfile
import hashlib
import subprocess
import requests
from loguru import logger
from alive_progress import alive_bar


def md5_hash(tarpath, buffer_size: int = 1024 * 1024):
    md5 = hashlib.md5()
    with open(tarpath, 'rb') as handle:
        data = handle.read(buffer_size)
        while data:
            md5.update(data)
            data = handle.read(buffer_size)
    return md5.hexdigest()


def uncompress(tarpath, outfile):
    with tarfile.open(tarpath, mode="r:gz") as archive, open(outfile, 'wb') as handle:
        member_name = 'data/ncbi_plasmid_full_seqs.fas'
        try:
            member = archive.getmember(member_name)
            handle.write(archive.extractfile(member).read())
        except ValueError:
            logger.error(f"Could not extract {member_name} to {outfile}")


@logger.catch(onerror=lambda _: sys.exit(1))
def check_database_installation(dbpath):
    """
    checks database is installed correctly
    :param dbpath: database directory
    """
    f = os.path.join(dbpath, 'ncbi_plasmids.mmi',)
    if not os.path.exists(f):
        logger.info(
            f"Database directory is missing file {f}. Database will be downloaded."
        )
        get_database(dbpath)


def building_database(db):
    try:
        subprocess.run(f"minimap2 -d {db}", shell=True, check=True)
        os.remove(db)
    except subprocess.CalledProcessError:
        logger.error("Could not building database.")


def get_database(dbpath):
    logger.info("Downloading MOB-suite database v. 3.1.8.")
    tarball = "MOB-suite_database_v_3.1.8.tar.gz"
    tarpath = os.path.join(dbpath, tarball)
    url = "https://zenodo.org/records/10304948/files/data.tar.gz"
    md5_value = '670e293e4ff1918db5b654278a60bea8'
    
    if os.path.exists(dbpath):
        shutil.rmtree(dbpath)
    os.mkdir(dbpath)

    try:
        with open(tarpath, 'wb') as handle, requests.get(url, stream=True) as r:
            total_length = r.headers.get("content-length")
            if total_length is not None:  # content length header is set
                total_length = int(total_length)
            with alive_bar(total=total_length, scale="SI") as bar:
                for data in r.iter_content(chunk_size=1024 * 1024):
                    handle.write(data)
                    bar(count=len(data))
    except IOError:
        logger.error(
            f"Could not download file from Zenodo! url={url}, path={tarpath}"
        )

    md5_sum = md5_hash(tarpath)
    if md5_sum == md5_value:
        logger.info(f"Database file download OK: {md5_sum}")
    else:
        logger.error(
            f"Corrupt database file! MD5 should be '{md5_value}' but is '{md5_sum}'"
        )
    database = os.path.join(dbpath, 'ncbi_plasmids')
    uncompress(tarpath, database)
    building_database(database)
    os.remove(tarpath)
