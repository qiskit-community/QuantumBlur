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

import numpy as np
from scipy.linalg import fractional_matrix_power

from qiskit import *

from PIL import Image

def _get_size(height):
    Lx = 0
    Ly = 0
    for (x,y) in height:
        Lx = max(x+1,Lx)
        Ly = max(y+1,Ly)
    return Lx,Ly

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
    n = int(np.ceil(np.log(length)/np.log(2)))
    
    # iteratively build list
    line = ['0','1']
    for j in range(n-1):
        # first append a reverse-ordered version of the current list
        line = line + line[::-1]
        # then add a '0' onto the end of all bit strings in the first half
        for j in range(int(len(line)/2)):
            line[j] += '0'
        # and a '1' for the second half
        for j in range(int(len(line)/2),int(len(line))):
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
        N += amp*np.conj(amp)
    for j,amp in enumerate(ket):
        ket[j] = amp/np.sqrt(N)
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
                state[ int(bitstring,2) ] = np.sqrt( log**(h/max_h) )
            else:
                state[ int(bitstring,2) ] = np.sqrt( h )
    state = normalize(state)
        
    # define and initialize quantum circuit            
    qc = QuantumCircuit(n,n)
    qc.initialize(state,range(n))
    qc.name = '('+str(Lx)+','+str(Ly)+')'

    return qc


def _circuit2probs(qc):
    
    # blank copy of circuit
    new_qc = qc.copy()
    new_qc.data = []
    
    initial_ket = np.array([1])
    for gate in qc.data:
        if gate[0].name=='initialize':
            initial_ket = np.kron(initial_ket,gate[0].params)
        else:
            new_qc.data.append(gate)
    
    ket = quantum_info.Statevector(initial_ket)
    ket = ket.evolve(new_qc)
    
    return ket.probabilities_dict()
    

def _probs2height(qc, probs, log):
    
    # get grid info
    (Lx,Ly) = eval(qc.name)
    grid,_ = make_grid(Lx,Ly)
    
    # set height to probs value, rescaled such that the maximum is 1
    max_h = max( probs.values() )   
    height = {}
    for bitstring in probs:
        if bitstring in grid:
            height[grid[bitstring]] = probs[bitstring]/max_h
         
    # take logs if required
    if log:
        for pos in height:
            if height[pos]>0:
                height[pos] = max( np.log(log*height[pos])/np.log(log), 0)
            else:
                height[pos] = 0
                        
    return height
    
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
    
    image = Image.new('L',(Lx,Ly))
    for x in range(Lx):
        for y in range(Ly):
            if (x,y) in height:
                h = height[x,y]/h_max
            else:
                h = 0
            image.putpixel((x,y), int(255*h) )
            
    return image


def image2circuits(image, log=None):
    """
    Converts an image to a set of three circuits, with one corresponding to each RGB colour channel.
    
    Args:
        image (Image): An RGB encoded image.
        log (int): If given, a logarithmic encoding is used with the
            given value as the base.
        
    Returns:
        circuits (list): A list of quantum circuits encoding the image.
    """
        
    # turn rgb into heights dictionaries
    Lx,Ly = image.size
    heights = [{} for j in range(3)]
    for x in range(Lx):
        for y in range(Ly):
            rgb = image.getpixel((x,y))
            for j in range(3):
                heights[j][x,y] = rgb[j]
                
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
       
    Lx,Ly = _get_size(heights[0])
    h_max = [max(height.values()) for height in heights]
    
    image = Image.new('RGB',(Lx,Ly))
    for x in range(Lx):
        for y in range(Ly):
            rgb = []
            for j,height in enumerate(heights):
                if (x,y) in height:
                    h = height[x,y]/h_max[j]
                else:
                    h = 0
                rgb.append( int(255*h) )
            image.putpixel((x,y), tuple(rgb) )
            
    return image


def swap_heights(height0, height1, fraction, log=None):
    
    # convert heights to circuits
    circuits = [height2circuit(height,log=log) for height in [height0,height1]]
    
    # make sure the quantum registers have distinct names
    qregs = [qc.qregs[0].name for qc in circuits]
    for j,qc in enumerate(circuits):
        qc.qregs[0].name = str(j)
    
    # combine them (in parallel)
    double_qc = circuits[0] + circuits[1]
    
    # apply the partial swap
    U = np.array([
        [1, 0, 0, 0],
        [0, 0, 1, 0],
        [0, 1, 0, 0],
        [0, 0, 0, 1]
    ])
    U = fractional_matrix_power(U,fraction)
    for q in range(double_qc.qregs[0].size):
        double_qc.unitary(U, [double_qc.qregs[0][q],double_qc.qregs[1][q]], label='partial_swap')
    
    # get the probability distributions for each register
    p = _circuit2probs(double_qc)
    marginals = [{},{}]
    nq0 = circuits[0].qregs[0].size
    for string in p:
        substrings = [string[0:nq0], string[nq0::]]
        for j,substring in enumerate(substrings):
            if substring in marginals[j]:
                marginals[j][substring] += p[string]
            else:
                marginals[j][substring] = p[string]
                
    # restore the original names of the quantum registers
    for j,qc in enumerate(circuits):
        qc.qregs[0].name = qregs[j]
    
    # convert the prob dists to heights
    new_heights = []
    for j,marginal in enumerate(marginals):
        new_heights.append( _probs2height(circuits[j],marginal,log) )
        
    return new_heights[0], new_heights[1]