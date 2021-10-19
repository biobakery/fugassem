from template import tools
import os
import numpy as np

def transpose():
    inputFilePath = os.path.abspath('template/data/misc/TemplateData.txt')
    
    with open(inputFilePath) as file:
        lis = [x.replace('\n', '').split(' ') for x in file]

    x = np.array(lis)
    print("\nInput")
    print(x)
    print("\Output")
    print(np.transpose(x)) 