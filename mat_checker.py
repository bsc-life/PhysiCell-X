import argparse
import numpy as np
import h5py,glob
import scipy.io


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("data_folder", nargs='?', help="folder were the data is stored", default="./outputq")
    # temporarily works for a folder
    args = parser.parse_args()
    data_folder = args.data_folder

    for file in glob.glob(data_folder+"/*.mat"):
        print(file)
        try:
            f = h5py.File('somefile.mat','r')
            data = f.get('data/variable1')
        except:
            with open(data_folder+"/mat_files.txt",'a') as f:
                f.write("File "+file+ " not valid with h5py\n")
        try:
            mat = scipy.io.loadmat(file)

        except:
            with open(data_folder+"/mat_files.txt",'a') as f:
                f.write("File "+file+ " not valid with scipy\n")
        

