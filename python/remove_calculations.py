'''
Delete recursively folders and contents matching a list input.
'''

import shutil
import sys
import os

data_path = sys.argv[1]  # The path where data are located
remove_list = sys.argv[2:]  # A list of generic directories to be removed

print(remove_list)
# For all possible directories
for item in os.walk(data_path):

    path = item[0]  # The path
    for r in remove_list:  # Loop for folders to delete
        if r in path:  # Check if delete folder is in path
            shutil.rmtree(path)  # Remove folder and contents
            print('Removing: '+path)
