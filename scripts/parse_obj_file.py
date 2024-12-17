import numpy as np
import os

def parse_obj_file(*, path_to_obj_file, output_path):
    os.makedirs(output_path, exist_ok=True)

    vertices = []
    faces = []
    with open(file=path_to_obj_file, mode='r') as obj_file:
        lines = obj_file.readlines()

        for line in lines:
            line = line.strip()
            words = line.split()
            if words[0] == 'v':
                coords = np.array(words[1:], dtype=float)
                vertices += [coords]
            
            if words[0] == 'f':
                face = np.array(words[1:], dtype=int) - 1
                faces += [face]


    vertices = np.array(vertices, dtype=float)
    faces = np.array(faces, dtype=int)

    np.savetxt(os.path.join(output_path, 'vertices.txt'), vertices)
    np.savetxt(os.path.join(output_path, 'faces.txt'), faces)

if __name__ =="__main__":
    path_to_obj_file = '/home/itamarg/workspace/VoroCrustRICH/obj_files/smiley_bad_nose.obj'
    output_path = '/home/itamarg/workspace/VoroCrustRICH/data/smiley_bad_nose/'

    parse_obj_file(path_to_obj_file=path_to_obj_file, output_path=output_path)