#!/usr/bin/env python3

from PyQt5 import QtGui  # Added to be able to import ovito

from ovito.modifiers import VoronoiAnalysisModifier
from ovito.io import import_file


def log_frames(name):
    '''
    Gather the number of frames from the log file.

    inputs:
        name = The name with the path of the LAMMPS log file
    outputs:
        frames = The total number of frames
    '''

    steps = []
    condition = False
    with open(name) as f:
        for line in f:
            line = line.strip().split(' ')

            if line[0] == 'Step':
                condition = True
                continue

            if line[0] == 'Loop':
                break

            if condition:
                steps.append(int(line[0]))

    frames = len(steps)

    return frames


def vp(name, frame, edge_threshold=0.1):
    '''
    Grap the Voronoi Polyhedra (VP) indexes.

    inputs:
        name = The name with the path of the trajectories
        edge_threshold = The threshold for the edges considered
        frame = The frame of the trajectories
    outputs:
        indexes = The VP indexes
    '''

    # Load input data and create an ObjectNode with a data pipeline.
    node = import_file(name, multiple_frames=True)

    voro = VoronoiAnalysisModifier(
                                   compute_indices=True,
                                   use_radii=False,
                                   edge_threshold=edge_threshold
                                   )

    node.modifiers.append(voro)

    out = node.compute(frame)
    indexes = out.particle_properties['Voronoi Index'].array

    return indexes
