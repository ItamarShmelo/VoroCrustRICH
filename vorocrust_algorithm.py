import numpy as np
import os
from os import path
import sys
import logging

from __init__ import *

logger = get_logger(name='VoroCrustAlgorithm', logfile=path.join('vorocrust_algorithm.log'), level=logging.INFO)

class VoroCrustAlgorithm:
    def __init__(self, *, 
                    input_dir=None, 
                    sharpTheta=None, 
                    maxRadius=None, 
                    L_Lipschitz=None, 
                    alpha=None,
                    maximal_num_iter=None,
                    num_of_samples_edges=None,
                    num_of_samples_faces=None):
        self.input_dir = input_dir
        self.sharpTheta = sharpTheta
        self.maxRadius = maxRadius
        self.L_Lipschitz = L_Lipschitz
        self.alpha = alpha
        self.maximal_num_iter = maximal_num_iter
        self.num_of_samples_edges = num_of_samples_edges
        self.num_of_samples_faces = num_of_samples_faces

        logger.info(f"Initializing VoroCrustAlgorithm with the following parameters:")
        for key, value in vars(self).items():
            logger.info(f"{key}: {value}")

        logger.info(f"reading data from {input_dir}")
        
        vertices_file = path.join(input_dir, "vertices.txt")
        faces_file = path.join(input_dir, "faces.txt")

        vertices_data = np.loadtxt(vertices_file)
        faces_data = np.array(np.loadtxt(faces_file), dtype=int)

        vertices = [Vector3D(x=x, y=y, z=z) for [x, y, z] in vertices_data]

        logger.info("Building PL_Complex")
        self.plc = PL_Complex(vertices=vertices)

        for face in faces_data:
            assert len(face) == 3, "Faces must have exactly 3 vertices"
            self.plc.addFace(indices=face)
        
        assert self.plc.checkAllVerticesAreUnique(), "Vertices are not unique"
        assert self.plc.checkAllVerticesAreOnFace(), "Not all vertices are part of a face not on a face"

        logger.info("Detecting Sharp Features")
        self.plc.detectSharpFeatures(sharpTheta=theta)

        self.trees = Trees()
        
        self.cornerDriver = CornersRMPS(maxRadius=self.maxRadius, 
                                        L_Lipschitz=self.L_Lipschitz, 
                                        sharpTheta=self.sharpTheta, 
                                        plc=self.plc)
        
        self.edgesDriver = EdgesRMPS(maxRadius=self.maxRadius, 
                                     L_Lipschitz=self.L_Lipschitz, 
                                     alpha=self.alpha, 
                                     sharpTheta=self.sharpTheta, 
                                     plc=self.plc)
        
        self.facesDriver = FacesRMPS(maxRadius=self.maxRadius, 
                                     L_Lipschitz=self.L_Lipschitz, 
                                     alpha=self.alpha, 
                                     sharpTheta=self.sharpTheta, 
                                     plc=self.plc)
        
        self.sliverDriver = SliverDriver(L_Lipschitz=self.L_Lipschitz)

    def generateBoundaryBallSamples(self):
        logger.info("Generating sharp features samples")
        self.trees.loadPLC(plc=self.plc, Nsample_edges=self.num_of_samples_edges, Nsample_faces=self.num_of_samples_faces)

        logger.info("Sharp Corner Ball Sampling")
        self.cornerDriver.loadCorners(sharp_corners=self.plc.sharp_corners)
        self.cornerDriver.doSampling(corner_ball_tree=self.trees.ball_kd_vertices, trees=self.trees)

        logger.info("Enforce Lipschitz Condition on Corner Balls")
        enforceLipschitzness(ball_tree=self.trees.ball_kd_vertices, L_Lipschitz=self.L_Lipschitz)

        for i in range(self.maximal_num_iter):
            logging.info(f"Iteration {i+1} of {self.maximal_num_iter}")
            
            logger.info("Doing Sharp Edge Ball Sampling")
            do_edge_sampling = True
            while do_edge_sampling:
            
                self.edgesDriver.loadEdges(sharp_edges=self.plc.sharp_edges)
                self.edgesDriver.doSampling(edges_ball_tree=self.trees.ball_kd_edges, trees=self.trees)
                self.trees.ball_kd_edges.remakeTree()
                logger.info("Enforce Lipschitz Condition on Edge Balls")
                do_edge_sampling = enforceLipschitzness(ball_tree=self.trees.ball_kd_edges, L_Lipschitz=L_lip)
            
            logger.info("Doing Face Ball Sampling")
            do_face_sampling = True
            while do_face_sampling:
                self.facesDriver.loadFaces(sharp_faces=self.plc.faces)
                self.facesDriver.doSampling(faces_ball_tree=self.trees.ball_kd_faces, trees=self.trees)
                self.trees.ball_kd_faces.remakeTree()
                logger.info("Enforce Lipschitz Condition on Face Balls")
                do_face_sampling = enforceLipschitzness(ball_tree=self.trees.ball_kd_faces, L_Lipschitz=L_lip)
            
            if i < max_number_of_iter - 1:
                # if no slivers were eliminated break the loop
                logger.info("Eliminating Slivers")
                if not self.sliverDriver.eliminateSlivers(trees=self.trees): break
                logger.info("Enforce Lipschitz Condition on all after sliver elimination")
                enforceLipschitzness(self.trees.ball_kd_vertices, L_Lipschitz=L_lip)
                enforceLipschitzness(self.trees.ball_kd_edges, L_Lipschitz=L_lip)
                enforceLipschitzness(self.trees.ball_kd_faces, L_Lipschitz=L_lip)

    def run(self):
        logger.info(f"Running VoroCrustAlgorithm")

        logger.info("Generating boundary balls")
        self.generateBoundaryBallSamples()

        logger.info(f"Generated {self.trees.ball_kd_vertices.size()} Vertex balls")
        logger.info(f"Generated {self.trees.ball_kd_edges.size()} Edge balls")
        logger.info(f"Generated {self.trees.ball_kd_faces.size()} Face balls")
        logger.info(f"Total number of balls: {self.trees.ball_kd_vertices.size() + self.trees.ball_kd_edges.size() + self.trees.ball_kd_faces.size()}")

        logger.info("Generating boundary seeds")
        self.boundary_seeds = self.sliverDriver.getSeeds(trees=self.trees)
        
        logger.info(f"Generated {len(self.boundary_seeds)} boundary seeds")

        self.inside_seeds, self.outside_seeds = determineZoneOfSeeds(seeds=self.boundary_seeds, zone_plcs=[self.plc])

        logger.info(f"Generated {len(self.inside_seeds)} inside seeds")
        logger.info(f"Generated {len(self.outside_seeds)} outside seeds")

        logger.info("Generating volume seeds")
        self.volume_seeds = randomSampleVolumeSeeds(zones_plcs=[self.plc], zones_boundary_seeds=[self.inside_seeds, self.outside_seeds], maxSize=self.maxRadius, trees=self.trees, L_Lipschitz=self.L_Lipschitz)

        logger.info(f"Generated {len(self.volume_seeds[0])} volume seeds")

        logger.info(
            "\n" +"-"*30 + 
            f"\nSummary:\n" +
            f"Total number of seeds: {len(self.inside_seeds) + len(self.outside_seeds) + len(self.volume_seeds[0])}\n" +
            f"Number of outside Seeds: {len(self.outside_seeds)}\n" +
            f"Number of inside boundary seeds: {len(self.inside_seeds)}\n" +
            f"Number of volume seeds: {len(self.volume_seeds[0])}\n" +
            "-"*30
            )
        
        logger.info(f"VoroCrust Finished!")

    def write_all_seeds_to_output_directory(self, *, output_dir):
        logger.info("Writing all seeds to output directory")
        os.makedirs(output_dir)

        assert self.inside_seeds is not None, "Inside seeds are not generated yet"
        write_seeds_to_dir(self.inside_seeds, output_dir=path.join(output_dir, "inside_seeds"))

        assert self.outside_seeds is not None, "Outside seeds are not generated yet"
        write_seeds_to_dir(self.outside_seeds, output_dir=path.join(output_dir, "outside_seeds"))

        assert self.volume_seeds is not None, "Volume seeds are not generated yet"
        write_seeds_to_dir(self.volume_seeds[0], output_dir=path.join(output_dir, "volume_seeds"))

def write_seeds_to_dir(seeds, output_dir):
    os.makedirs(output_dir)

    x = np.zeros(len(seeds), dtype=np.float64)
    y = np.zeros(len(seeds), dtype=np.float64)
    z = np.zeros(len(seeds), dtype=np.float64)
    r = np.zeros(len(seeds), dtype=np.float64)

    for i, seed in enumerate(seeds):
        x[i] = seed.p.x
        y[i] = seed.p.y
        z[i] = seed.p.z
        r[i] = seed.radius
    
    np.savetxt(path.join(output_dir, "x.txt"), x)
    np.savetxt(path.join(output_dir, "y.txt"), y)
    np.savetxt(path.join(output_dir, "z.txt"), z)
    np.savetxt(path.join(output_dir, "r.txt"), r)

if __name__ == "__main__":
    input_dir = sys.argv[1]
    dirname = path.basename(input_dir)
    output_dir = path.join("/home/itamarg/workspace/VoroCrustRICH/output/", sys.argv[2])
    theta = np.pi * 45./180.
    maxRadius = 1.0
    L_lip = 0.3
    alpha = 0.13
    max_number_of_iter = 1
    num_samples_of_edges = int(1e5)
    num_samples_of_faces = int(1e6)

    alg = VoroCrustAlgorithm(input_dir=input_dir,
                       sharpTheta=theta,
                       maxRadius=maxRadius,
                       L_Lipschitz=L_lip,
                       alpha=alpha,
                       maximal_num_iter=max_number_of_iter,
                       num_of_samples_edges=num_samples_of_edges,
                       num_of_samples_faces=num_samples_of_faces)

    alg.run()

    alg.write_all_seeds_to_output_directory(output_dir=path.join(output_dir, "seeds"))