"""Version information for :mod:`miadbviz`.

Run with ``python -m miadbviz.version``
"""

from __future__ import annotations

import logging
import os
from subprocess import CalledProcessError, check_output

__all__ = [
    "VERSION",
    "get_version",
    "get_git_hash",
]

VERSION = "0.0.1-dev"


def get_git_hash() -> str:
    """
    Get the :mod:`miadbviz` git hash.

    :returns: Git Hash.
    :rtype: str
    """
    with open(os.devnull, "w") as devnull:
        try:
            ret = check_output(
                ["git", "rev-parse", "HEAD"],
                cwd=os.path.dirname(__file__),
                stderr=devnull,
            )
        except CalledProcessError:
            return "UNHASHED"
        else:
            return ret.strip().decode("utf-8")[:8]


def get_version(with_git_hash: bool = False):
    """
    Get the :mod:`miadbviz` version string, including a git hash.

    :param with_git_hash: Flag to include git hash.
    :type with_git_hash: bool

    :returns: The version.
    :rtype: str
    """
    return f"{VERSION}-{get_git_hash()}" if with_git_hash else VERSION


if __name__ == "__main__":
    logging.debug(get_version(with_git_hash=True))
