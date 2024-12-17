import numpy as np
from matplotlib import pyplot as plt
from os import path
import os
import sys
import logging 

from __init__  import *

def run_voro_crust(*, input_dir, output_dir, theta, maxRadius, L_lip, alpha, max_number_of_iter, num_samples_of_edges, num_samples_of_faces):
    os.makedirs(output_dir, exist_ok=True)

    logger = logging.getLogger(name='run_voro_crust')
    filehandler = logging.FileHandler(path.join(output_dir, f'{path.basename(input_dir)}.log'), mode='w')
    formatter = logging.Formatter("%(name)s [%(levelname)s] : %(message)s")
    filehandler.setFormatter(formatter)
    
    logger.addHandler(filehandler)  
    logger.setLevel(level=logging.INFO)

    logger.info(f"Running VoroCrustRICH with the following parameters:")

    logger.info(f"Reading data from {input_dir}")

    vertices_file = path.join(input_dir, "vertices.txt")
    faces_file = path.join(input_dir, "faces.txt")

    vertices_data = np.loadtxt(vertices_file)
    faces_data = np.array(np.loadtxt(faces_file), dtype=int)

    vertices = [Vector3D(x=x, y=y, z=z) for [x, y, z] in vertices_data]

    logger.info("Building PL_Complex")
    plc = PL_Complex(vertices=vertices)

    logger.info(f"PLC bounding box = {plc.getBoundingBox()}, maxRadius = {maxRadius:g}")

    for indices in faces_data:
        plc.addFace(indices=indices)

    algorithm = VoroCrustAlgorithm( plc=plc, 
                                    sharpTheta=theta, 
                                    maxRadius=maxRadius, 
                                    L_Lipschitz=L_lip, 
                                    alpha=alpha, 
                                    maximal_num_iter=max_number_of_iter,
                                    num_of_samples_edges=num_samples_of_edges,
                                    num_of_samples_faces=num_samples_of_faces)

    logger.info("Running VoroCrust")
    algorithm.run()

    logger.info("Writing PL_Complex vtu file")
    vorocrust_vtk.write_vtu_PL_Complex(filename=path.join(output_dir, "plc.vtu"), plc=plc)
    logger.info("Writing Sharp Features Sampling to vtu file")
    vorocrust_vtk.write_vtu_trees(filename=path.join(output_dir, "trees_all.vtu"), trees=algorithm.trees)
    
    logging.info("Generating Algorithm Dump")
    algorithm.dump(dirname=output_dir)

    logger.info("Writing Boundary Ball Sampling")
    try:
        vorocrust_vtk.write_ballTree(filename=path.join(output_dir, "corner_sampling_all.vtp"), ball_tree=algorithm.trees.ball_kd_vertices)
        vorocrust_vtk.write_ballTree(filename=path.join(output_dir, "edges_sampling_all.vtp"), ball_tree=algorithm.trees.ball_kd_edges)
        vorocrust_vtk.write_ballTree(filename=path.join(output_dir, "faces_sampling_all.vtp"), ball_tree=algorithm.trees.ball_kd_faces)
    except RuntimeError as e:
        print(f"Error occurred during VoroCrust execution: {e}")
    except MemoryError as e:
        print(f"Not enough memory for writing trees: {e}")

    logger.info("Generating Seeds")
    seeds = algorithm.getSeeds()

    logger.info("Dumping Seeds")
    dumpSeeds(dirname=path.join(output_dir, "all_seeds.txt"), seeds=seeds)

    seeds_points = [seed.p for seed in seeds]
    try:
        logger.info("Writing Seeds vtu file")
        vorocrust_vtk.write_points(filename=path.join(output_dir, "seeds_all.vtu"), points=seeds_points)
        
    except RuntimeError as e:
        print(f"Error occurred during VoroCrust execution: {e}")
    except MemoryError as e:
        print(f"Not enough memory for writing seeds: {e}")
    
    logger.info("Determing Zone of Seeds")
    zone_plcs = [algorithm.plc]
    zone_seeds = determineZoneOfSeeds(seeds=seeds, zone_plcs=zone_plcs)

    logger.info("Dumping Zone Seeds")
    dumpSeeds(dirname=path.join(output_dir, "seeds_in"), seeds=zone_seeds[0])
    dumpSeeds(dirname=path.join(output_dir, "seeds_out"), seeds=zone_seeds[1])

    logger.info("Generate Volume Seeds")
    zone_volume_seeds = algorithm.randomSampleSeeds(zones_plcs=zone_plcs, zones_boundary_seeds=zone_seeds, maxSize=maxRadius)

    logger.info("Dumping Volume Seeds")
    dumpSeeds(dirname=path.join(output_dir, "zone_in_volume_seeds"), seeds=zone_volume_seeds[0])

    logger.info("Finished!")

if __name__ == "__main__":
    input_dir = "/home/itamarg/workspace/VoroCrustRICH/data/fox"
    output_dir = "/home/itamarg/workspace/VoroCrustRICH/output/fox_script"
    theta = np.pi * 45./180.
    maxRadius = 1.0
    L_lip = 0.3
    alpha = 0.13
    max_number_of_iter = 1
    num_samples_of_edges = int(1e5)
    num_samples_of_faces = int(1e6)

    run_voro_crust(
        input_dir=input_dir,
        output_dir=output_dir,
        theta=theta,
        maxRadius=maxRadius,
        L_lip=L_lip,
        alpha=alpha,
        max_number_of_iter=max_number_of_iter,
        num_samples_of_edges=num_samples_of_edges,
        num_samples_of_faces=num_samples_of_faces
    )
