from cpp_modules._VoroCrustAlgorithm import dumpSeeds, determineSeedsInOut, enforceLipschitzness, randomSampleVolumeSeeds
from cpp_modules._PL_Complex import PL_Complex, VoroCrustVertex, VoroCrustEdge, VoroCrustFace
from cpp_modules._Vector3D import Vector3D
import cpp_modules._vorocrust_vtk as vorocrust_vtk
from cpp_modules._VoroCrust_kd_tree import VoroCrust_KD_Tree, VoroCrust_KD_Tree_Boundary, VoroCrust_KD_Tree_Ball
from cpp_modules._trees import Trees
from cpp_modules._RMPS import Seed, CornersRMPS, EdgesRMPS, FacesRMPS, SliverDriver

import logging

def get_logger(name, logfile, level=logging.INFO):
    logger = logging.getLogger(name=name)
    filehandler = logging.FileHandler(logfile, mode='w')
    formatter = logging.Formatter("%(name)s [%(levelname)s] : %(message)s")
    filehandler.setFormatter(formatter)
    
    logger.addHandler(filehandler)  
    logger.setLevel(level=level)
    return logger