# -*- coding: utf-8 -*-

# This code is part of Qiskit.
#
# (C) Copyright IBM 2020s.
#
# This code is licensed under the Apache License, Version 2.0. You may
# obtain a copy of this license in the LICENSE.txt file in the root directory
# of this source tree or at http://www.apache.org/licenses/LICENSE-2.0.
#
# Any modifications or derivative works of this code must retain this
# copyright notice, and modified files need to carry a notice indicating
# that they have been altered from the originals.


'''
The imports that follow are highly non-standard and require some explanation. 

This file is designed to run in both a modern, fully functioning Python
environment, with Python 3.x and the ability to use external libraries.
It is also designed to function using only the standard library (in
addition to MicroQiskit) in any Python from 2.7 onwards.

The deciding factor is whether Qiskit is available to be imported. If so,
the following external libraries are required dependencies:

qiskit
numpy
scipy
PIL

Otherwise, MicroQiskit will be used in place of Qiskit, and alternative
techniques using only the standard library will be used in place of the
other dependencies.

More information on Qiskit can be found at

https://qiskit.org

and information on MicroQiskit can be found at

https://github.com/qiskit-community/MicroQiskit
'''

import math

try:
    from qiskit import *
    simple_python = False
except:
    print('Unable to import Qiskit, so MicroQiskit will be used instead')
    from microqiskit import *
    simple_python = True

if not simple_python:
    import numpy as np
    from scipy.linalg import fractional_matrix_power
    from PIL.Image import new as newimage, Image
    
if simple_python:
    class Image():
        def __init__(self):
            self.mode = None
            self.size = None
            self.info = None
        def getpixel(self,xy):
            return self.info[xy]
        def putpixel(self, xy, value):
            self.info[xy] = value
        def show(self):
            for x in range(self.size[0]):
                for y in range(self.size[1]):
                    print('('+str(x)+','+str(y)+')'+': '+str(self.info[x,y]))

    def newimage(mode, size):
        img = Image()
        img.mode = mode
        img.size = size
        if mode=='L':
            blank = 0
        elif mode=='RGB':
            blank = (0,0,0)
        img.info = {(x,y):blank\
                    for x in range(size[0])\
                    for y in range(size[1])}
        return img
        

def _kron(vec0,vec1):
    '''
    Calculates the tensor product of two vectors.
    '''
    new_vec = []
    for amp0 in vec0:
        for amp1 in vec1:
            new_vec.append(amp0*amp1)
    return new_vec


def _combine_circuits(qc0,qc1):
    '''
    Combines a pair of circuits in parallel.
    
    For MicroQiskit, this only works if the circuits contain only
    initialization.
    '''
    
    if simple_python:
    
        warning = "Combined circuits should contain only initialization."
    
        # create a circuit with the combined number of qubits
        num_qubits = qc0.num_qubits + qc1.num_qubits
        combined_qc = QuantumCircuit(num_qubits)

        # extract statevectors for any initialization commands
        kets = [None,None]
        for j,qc in enumerate([qc0, qc1]):
            for gate in qc.data:
                assert gate[0]=='init', warning
                kets[j] = gate[1]

        # combine into a statevector for all the qubits
        ket = None
        if kets[0] and kets[1]:
            ket = _kron(kets[0], kets[1])
        elif kets[0]:
            ket = _kron(kets[0], [1]+[0]*(2**qc1.num_qubits-1))
        elif kets[1]:
            ket = _kron([1]+[0]*(2**qc0.num_qubits-1),kets[1])

        # use this to initialize
        if ket:
            combined_qc.initialize(ket)
                    
    else:

        circuits = [qc0,qc1]
        
        # make sure the quantum registers have distinct names
        qregs = [qc.qregs[0].name for qc in circuits]
        for j,qc in enumerate(circuits):
            qc.qregs[0].name = str(j)

        # combine the circuits in parallel
        combined_qc = circuits[0] + circuits[1]
        
        # restore the original names of the quantum registers
        for j,qc in enumerate(circuits):
            qc.qregs[0].name = qregs[j]
            
    return combined_qc


def _get_size(height):
    """
    Determines the size of the grid for the given height map.
    """
    Lx = 0
    Ly = 0
    for (x,y) in height:
        Lx = max(x+1,Lx)
        Ly = max(y+1,Ly)
    return Lx,Ly


def _circuit2probs(qc):
    """
    Runs the given circuit, and returns the resulting probabilities.
    """
    if simple_python:
        probs = simulate(qc,get='probabilities_dict')
    else:
        # separate circuit and initialization
        new_qc = qc.copy()
        new_qc.data = []
        initial_ket = [1]
        for gate in qc.data:
            if gate[0].name=='initialize':
                initial_ket = _kron(initial_ket,gate[0].params)
            else:
                new_qc.data.append(gate)
        # then run it
        ket = quantum_info.Statevector(initial_ket)
        ket = ket.evolve(new_qc)
        probs = ket.probabilities_dict()
    
    return probs
    

def _probs2height(qc, probs, log):
    """
    Determines the height map for a given circuit that encodes a height map
    and the a corresponding set of probabilities.
    """
    # get grid info
    (Lx,Ly) = eval(qc.name)
    grid,_ = make_grid(Lx,Ly)
    
    # set height to probs value, rescaled such that the maximum is 1
    max_h = max( probs.values() )   
    height = {(x,y):0 for x in range(Lx) for y in range(Ly)}
    for bitstring in probs:
        if bitstring in grid:
            height[grid[bitstring]] = float(probs[bitstring])/max_h
         
    # take logs if required
    if log:
        for pos in height:
            if height[pos]>0:
                height[pos] = max(math.log(log*height[pos])/math.log(log),0)
            else:
                height[pos] = 0
                        
    return height


def _image2heights(image):
    """
    Converts an rgb image into a list of three height dictionaries, one for
    each colour channgel.
    """
    Lx,Ly = image.size
    heights = [{} for j in range(3)]
    for x in range(Lx):
        for y in range(Ly):
            rgb = image.getpixel((x,y))
            for j in range(3):
                heights[j][x,y] = rgb[j]

    return heights


def _heights2image(heights):
    """
    Constructs an image from a set of three height dictionaries, one for each
    colour channel.
    """
    Lx,Ly = _get_size(heights[0])
    h_max = [max(height.values()) for height in heights]

    image = newimage('RGB',(Lx,Ly))
    for x in range(Lx):
        for y in range(Ly):
            rgb = []
            for j,height in enumerate(heights):
                if (x,y) in height:
                    h = float(height[x,y])/h_max[j]
                else:
                    h = 0
                rgb.append( int(255*h) )
            image.putpixel((x,y), tuple(rgb) )

    return image


def make_line ( length ):
    """
    Creates a list of bit strings of at least the given length, such
    that the bit strings are all unique and consecutive strings
    differ on only one bit.
    
    Args:
        length (int): Required length of output list.
    
    Returns:
        line (list): List of 2^n n-bit strings for n=⌊log_2(length)⌋
    """
    
    # number of bits required
    n = int(math.ceil(math.log(length)/math.log(2)))
    
    # iteratively build list
    line = ['0','1']
    for j in range(n-1):
        # first append a reverse-ordered version of the current list
        line = line + line[::-1]
        # then add a '0' onto the end of all bit strings in the first half
        for j in range(int(float(len(line)/2))):
            line[j] += '0'
        # and a '1' for the second half
        for j in range(int(float(len(line)/2)),int(len(line))):
            line[j] += '1'
            
    return line


def normalize(ket):
    """
    Normalizes the given statevector.
    
    Args:
        ket (list or array_like)
    
    Returns:
        ket (list or array_like)
    """
    N = 0
    for amp in ket:
        N += amp*amp.conjugate()
    for j,amp in enumerate(ket):
        ket[j] = float(amp)/math.sqrt(N)
    return ket


def make_grid(Lx,Ly=None):
    """
    Creates a dictionary that provides bit strings corresponding to
    points within an Lx by Ly grid.
    
    Args:
        Lx (int): Width of the lattice (also the height if no Ly is
            supplied).
        Ly (int): Height of the lattice if not Lx.
    
    Returns:
        grid (dict): Dictionary whose values are points on an
            Lx by Ly grid. The corresponding keys are unique bit
            strings such that neighbouring strings differ on only
            one bit.
        n (int): Length of the bit strings
        
    """
    # set Ly if not supplied
    if not Ly:
        Ly = Lx
    
    # make the lines
    line_x = make_line( Lx )
    line_y = make_line( Ly )
    
    # make the grid
    grid = {}
    for x in range(Lx):
        for y in range(Ly):
            grid[ line_x[x]+line_y[y] ] = (x,y)
            
    # determine length of the bit strings
    n = len(line_x[0]+line_y[0])
            
    return grid, n


def height2circuit(height, log=None):
    """
    Converts a dictionary of heights (or brightnesses) on a grid into
    a quantum circuit.
    
    Args:
        height (dict): A dictionary in which keys are coordinates
            for points on a grid, and the values are positive numbers of
            any type.
        log (int): If given, a logarithmic encoding is used with the
            given value as the base.
            
    Returns:
        qc (QuantumCircuit): A quantum circuit which encodes the
            given height dictionary.
    """
    # get bit strings for the grid
    Lx,Ly = _get_size(height)
    grid, n = make_grid(Lx,Ly)
    
    # create required state vector
    state = [0]*(2**n)
    max_h = max(height.values())
    for bitstring in grid:
        (x,y) = grid[bitstring]
        if (x,y) in height:
            h = height[x,y]
            if log:
                state[ int(bitstring,2) ] = math.sqrt(log**(float(h)/max_h))
            else:
                state[ int(bitstring,2) ] = math.sqrt( h )
    state = normalize(state)
        
    # define and initialize quantum circuit            
    qc = QuantumCircuit(n,n)
    if simple_python:
        # microqiskit style
        qc.initialize(state)
    else:
        qc.initialize(state,range(n))
    qc.name = '('+str(Lx)+','+str(Ly)+')'

    return qc

    
def circuit2height(qc, log=None):
    """
    Extracts a dictionary of heights (or brightnesses) on a grid from
    the quantum circuit into which it has been encoded.
    
    Args:
        qc (QuantumCircuit): A quantum circuit which encodes a height
            dictionary.
        log (int): If given, a logarithmic decoding is used with the
            given value as the base.
            
    Returns:
        height (dict): A dictionary in which keys are coordinates
            for points on a grid, and the values are floats in the
            range 0 to 1.
    """
    
    probs = _circuit2probs(qc)
    return _probs2height(qc, probs, log)

        
def swap_heights(height0, height1, fraction, log=None):
    """
    Given a pair of height maps for the same sized grid, a set of partial
    swaps is applied between corresponding qubits in each circuit.
    
    Args:
        height0, height1 (dict): Dictionaries in which keys are coordinates
            for points on a grid, and the values are floats in the range 0
            to 1.
        fraction (float): Fraction of swap gates to apply.
        log (int): If given, a logarithmic decoding is used with the
            given value as the base.
            
    Returns:
        new_height0, new_height1 (dict): As with the height inputs.
    """
    assert _get_size(height0)==_get_size(height1), \
    "Objects to be swapped are not the same size"
    
    # convert heights to circuits
    circuits = [height2circuit(height,log=log)\
                for height in [height0,height1]]
    
    # get number of qubits in each circuit
    num_qubits = circuits[0].num_qubits
    
    # combine the two circuits into one (in parallel)
    combined_qc = _combine_circuits(circuits[0], circuits[1])
    
    # apply the partial swap
    if not simple_python:
        U = np.array([
        [1, 0, 0, 0],
        [0, 0, 1, 0],
        [0, 1, 0, 0],
        [0, 0, 0, 1]
        ])
        U = fractional_matrix_power(U,fraction)
    for q in range(num_qubits):
        if not simple_python:
            combined_qc.unitary(U, \
                                [combined_qc.qregs[0][q],\
                                 combined_qc.qregs[1][q]],\
                                 label='partial_swap')
        else:
            q0 = q
            q1 = num_qubits + q
            combined_qc.cx(q1,q0)
            combined_qc.h(q1)
            combined_qc.cx(q0,q1)
            combined_qc.rz(-math.pi*fraction/2,q1)
            combined_qc.cx(q0,q1)
            combined_qc.rz(math.pi*fraction/2,q1)
            combined_qc.h(q1)
            combined_qc.cx(q1,q0)
    
    # get the probability distributions for each register
    p = _circuit2probs(combined_qc)
    marginals = [{},{}]
    for string in p:
        substrings = [string[0:num_qubits], string[num_qubits::]]
        for j,substring in enumerate(substrings):
            if substring in marginals[j]:
                marginals[j][substring] += p[string]
            else:
                marginals[j][substring] = p[string]
    
    # convert the prob dists to heights
    new_heights = []
    for j,marginal in enumerate(marginals):
        new_heights.append( _probs2height(circuits[j],marginal,log) )
        
    return new_heights[0], new_heights[1]


def height2image(height):
    """
    Converts a dictionary of heights (or brightnesses) on a grid to
    an image.

    Args:
        height (dict): A dictionary in which keys are coordinates
                for points on a grid, and the values are positive
                numbers of any type.

    Returns:
        image (Image): Monochrome image for which the given height
            dictionary determines the brightness of each pixel. The
            maximum value in the height dictionary is always white.
    """
    Lx,Ly = _get_size(height)
    h_max = max(height.values())

    image = newimage('L',(Lx,Ly))
    for x in range(Lx):
        for y in range(Ly):
            if (x,y) in height:
                h = float(height[x,y])/h_max
            else:
                h = 0
            image.putpixel((x,y), int(255*h) )

    return image


def swap_images(image0, image1, fraction, log=None):
    """
    Given a pair of same sized grid images, a set of partial swaps is applied
    between corresponding qubits in each circuit.
    
    Args:
        image0, image1 (Image): RGB encoded images.
        fraction (float): Fraction of swap gates to apply.
        log (int): If given, a logarithmic decoding is used with the
            given value as the base.
            
    Returns:
        new_image0, new_image1 (Image): RGB encoded images.
    """
    heights0 = _image2heights(image0)
    heights1 = _image2heights(image1)

    new_heights0 = []
    new_heights1 = []
    for j in range(3):
        nh0, nh1 = swap_heights(heights0[j], heights1[j], fraction, log=log)
        new_heights0.append(nh0)
        new_heights1.append(nh1)

    new_image0 = _heights2image(new_heights0)
    new_image1 = _heights2image(new_heights1)

    return new_image0, new_image1

def image2circuits(image, log=None):
    """
    Converts an image to a set of three circuits, with one corresponding to
    each RGB colour channel.

    Args:
        image (Image): An RGB encoded image.
        log (int): If given, a logarithmic encoding is used with the
            given value as the base.

    Returns:
        circuits (list): A list of quantum circuits encoding the image.
    """

    heights = _image2heights(image)

    circuits = []
    for height in heights:
        circuits.append( height2circuit(height, log=log) )

    return circuits


def circuits2image(circuits, log=None):
    """
    Extracts an image from list of circuits encoding the RGB channels.

    Args:
        circuits (list): A list of quantum circuits encoding the image.
        log (int): If given, a logarithmic decoding is used with the
            given value as the base.

    Returns:
        image (Image): An RGB encoded image.
    """

    heights = []
    for qc in circuits:
        heights.append( circuit2height(qc, log=log) )

    return _heights2image(heights)


def row_swap_images(image0, image1, fraction, log=None):
    """
    A variant of `swap_images` in which the swap process is done on each line
    of the images individually, rather than with the images as a whole. This
    makes it much faster.
    
    Args:
        image0, image1 (Image): RGB encoded images.
        fraction (float): Fraction of swap gates to apply.
        log (int): If given, a logarithmic decoding is used with the
            given value as the base.
            
    Returns:
        new_image0, new_image1 (Image): RGB encoded images.
    """
    images = [image0, image1]

    Lx,Ly = images[0].size

    # create separate images for each row
    rows = [[],[]]
    for j in range(2):
        for y in range(Ly):   
            rows[j].append(newimage('RGB',(Lx,1)))
            for x in range(Lx):
                rows[j][y].putpixel((x,0),images[j].getpixel((x,y)))


    # do the swap on the row images
    for y in range(Ly):
        rows[0][y], rows[1][y] = swap_images(rows[0][y], rows[1][y], fraction, log=log)

    # reconstruct the full images
    new_images = [newimage('RGB',(Lx,Ly)) for _ in range(2)]
    for j in range(2):
        for y in range(Ly):
            for x in range(Lx):
                new_images[j].putpixel((x,y),rows[j][y].getpixel((x,0)))

    return new_images[0], new_images[1]