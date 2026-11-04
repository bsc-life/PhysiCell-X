#!/usr/bin/env python3
# coding: utf-8

import os
import json
import platform
import urllib.request
import tarfile
import zipfile
import sys
import shutil
import subprocess

def get_os_arch():
    """Determine the operating system and architecture."""
    os_type = platform.system().lower()
    arch = platform.machine().lower()

    if os_type == 'linux':
        return 'linux-x64' if arch in ('x86_64', 'amd64') else 'linux-x86'
    elif os_type == 'darwin':
        if arch == 'arm64':
            return 'osx-arm64'
        elif arch == 'x86_64':
            return 'osx-x86_64'
        print(f"Unsupported architecture: {os_type} : {arch}")
        sys.exit(1)
    elif os_type.startswith('win'):
        return 'win64' if arch == 'amd64' else 'win32'
    else:
        print(f"Unsupported OS: {os_type}")
        sys.exit(1)

def ensure_directory_exists(path):
    """Create directory if it does not exist."""
    if not os.path.exists(path):
        os.makedirs(path)

def download_and_extract(url, dest_path, is_zip=False):
    """Download and extract a package from a URL."""
    # Ensure the destination directory exists
    ensure_directory_exists(dest_path)

    filename = url.split('/')[-1]
    file_path = os.path.join(dest_path, filename)

    # Download the file
    print(f"Downloading {filename} from {url}...")
    urllib.request.urlretrieve(url, file_path)

    # Extract the file
    print(f"Extracting {filename}...")
    if is_zip:
        with zipfile.ZipFile(file_path, 'r') as zip_ref:
            zip_ref.extractall(dest_path)
    elif filename.endswith(".deb"):
        result = subprocess.run(["dpkg-deb", "-x", file_path, dest_path])
        if result.returncode != 0:
            # Some shared filesystems reject timestamp updates on dangling
            # symlinks from split Debian development packages.  The dFBA build
            # only needs the headers and static archive, which dpkg-deb has
            # already extracted before reporting that metadata error.
            header = os.path.join(dest_path, "usr", "include", "sbml", "SBMLTypes.h")
            library = os.path.join(dest_path, "usr", "lib", "x86_64-linux-gnu",
                                   "libsbml-static.a")
            if not (os.path.isfile(header) and os.path.isfile(library)):
                result.check_returncode()
    else:
        with tarfile.open(file_path, 'r:gz') as tar_ref:
            tar_ref.extractall(dest_path)

    # Clean up
    os.remove(file_path)
    print(f"{filename} installed successfully.\n")

def normalize_libsbml_install(dest_path):
    """Expose the include/lib layout used by the dFBA Makefiles."""
    expected_header = os.path.join(dest_path, "include", "sbml", "SBMLTypes.h")
    expected_library = os.path.join(dest_path, "lib", "libsbml-static.a")
    if os.path.isfile(expected_header) and os.path.isfile(expected_library):
        return

    include_source = None
    library_source = None
    for root, _, files in os.walk(dest_path):
        if os.path.basename(root) == "sbml" and "SBMLTypes.h" in files:
            include_source = os.path.dirname(root)
        if "libsbml-static.a" in files:
            library_source = root

    if include_source and library_source:
        include_dest = os.path.join(dest_path, "include")
        library_dest = os.path.join(dest_path, "lib")
        if os.path.abspath(include_source) != os.path.abspath(include_dest):
            shutil.copytree(include_source, include_dest, dirs_exist_ok=True)
        if os.path.abspath(library_source) != os.path.abspath(library_dest):
            ensure_directory_exists(library_dest)
            # Only this archive is used by the dFBA Makefiles.  Copying the
            # complete Debian lib directory also copies libsbml.so, a dangling
            # link whose target belongs to the separate runtime package and
            # fails on some shared filesystems.
            shutil.copy2(os.path.join(library_source, "libsbml-static.a"),
                         expected_library)

def main():
    # Define paths
    current_folder = os.path.abspath(os.path.dirname(__file__))
    base_folder = os.path.abspath(os.path.join(current_folder, os.pardir))
    default_libs_path = os.path.join(base_folder, "addons", "dFBA", "ext")

    # Load package information
    json_packages = os.path.join(base_folder, 'beta', 'fba_packages.json')
    with open(json_packages) as fh:
        packages_dict = json.load(fh)

    # Determine OS and architecture
    arch = get_os_arch()

    # Install packages
    for pkg in ("coin-or", "libsbml"):
        pkg_path = os.path.join(default_libs_path, pkg)

        if pkg == "libsbml":
            normalize_libsbml_install(pkg_path)

        if pkg == "coin-or":
            installed = os.path.isfile(os.path.join(
                pkg_path, "include", "coin", "CoinPackedMatrix.hpp"))
        else:
            installed = (os.path.isfile(os.path.join(
                pkg_path, "include", "sbml", "SBMLTypes.h")) and
                os.path.isfile(os.path.join(pkg_path, "lib", "libsbml-static.a")))
        if installed:
            print(f"Package {pkg} is already installed. Skipping...")
            continue

        # Ensure the package directory exists
        ensure_directory_exists(pkg_path)

        pkg_dict = packages_dict.get(pkg, {}).get(arch)
        if not pkg_dict:
            print(f"No package information found for {pkg} on {arch}. Skipping...")
            continue

        url = pkg_dict['url']
        is_zip = url.endswith("zip")
        download_and_extract(url, pkg_path, is_zip)
        if pkg == "libsbml":
            normalize_libsbml_install(pkg_path)

if __name__ == "__main__":
    main()
