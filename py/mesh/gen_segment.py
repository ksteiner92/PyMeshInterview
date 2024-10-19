import numpy as np
from pymesh import *


class PolyFileReader:
    def __init__(self, path):
        nvert = 0
        nattrs = 0
        bndmarker = 0
        points = None
        with open(path, "r") as f:
            for il, line in enumerate(f.readlines()):
                tokens = line.split()
                if il == 0:
                    if len(tokens) != 4:
                        raise ValueError("First line does not have 4 parameters")
                    nvert = int(tokens[0])
                    points = np.empty((nvert, 2), dtype=np.float64)
                    nattrs = int(tokens[2])
                    bndmarker = int(tokens[3])
                    continue
                if il <= nvert:
                    points[int(tokens[0]), :] = np.array(tokens[1:3], dtype=np.float64)
