import os

from abseqPy.config import ABSEQROOT, EXTERNAL_DEP_DIR
from abseqPy.config import MEM_GB


# temporarily overrides PATH variable with EXTERNAL_DEP_DIR/bin, IGBLASTDB and IGDATA (if they exist)
class PriorityPath:
    def __init__(self):
        self.updated = False
        self.old_env = os.environ.copy()
        _env = os.environ.copy()

        # if the BIN directory exists, append it to the front of PATH variable
        override_path = os.path.abspath(os.path.join(ABSEQROOT, EXTERNAL_DEP_DIR, 'bin')) + os.path.sep
        if os.path.exists(override_path):
            _env['PATH'] = override_path + os.pathsep + _env['PATH']
            self.updated = True

        # if the igdata dir exists, override it irrespective of if there's already a IGDATA env
        override_igdata = os.path.abspath(os.path.join(ABSEQROOT, EXTERNAL_DEP_DIR, 'igdata')) + os.path.sep
        if os.path.exists(override_igdata):
            _env['IGDATA'] = override_igdata
            self.updated = True

        # if the igdb dir exists, override it irrespective of if there's already a IGBLASTDB env
        override_igdb = os.path.abspath(os.path.join(ABSEQROOT, EXTERNAL_DEP_DIR, 'databases')) + os.path.sep
        if os.path.exists(override_igdb):
            _env["IGBLASTDB"] = override_igdb
            self.updated = True

        if self.updated:
            os.environ.clear()
            os.environ.update(_env)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        if self.updated:
            os.environ.clear()
            os.environ.update(self.old_env)


def hasLargeMem(size=16):
    """
    tells if system has memory strictly larger than specified size in GB
    :param size: unit GB
    :return: true if virtual_memory > size
    """
    return MEM_GB > size
