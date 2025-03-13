"""
Contains utility methods for aeon modules.
"""

import os
import sys


class AeonUtil:
    """
    Utility static methods for aeon modules.
    """

    @staticmethod
    def get_aeon_dir():
        """
        Get the directory path for the aeon library
        """
        return os.path.dirname(sys.modules["aeon_ancestry"].__file__)

    @staticmethod
    def resolve_ref_filename(filename: str) -> str:
        """
        Resolves cli filename input. Checks if file exists in package first,
        then external. If both are unable to be found, raises ValueError.
        """
        ref_path = os.path.join(AeonUtil.get_aeon_dir(), filename)
        if os.path.exists(ref_path):
            return ref_path

        if os.path.exists(filename):
            return filename

        raise ValueError(f"Unable to find file: {filename}")
