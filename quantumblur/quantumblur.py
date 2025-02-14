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


"""
The imports that follow are highly non-standard and require some explanation. 

This file is designed to run in both a modern, fully functioning Python
environment, with Python 3.x and the ability to use external libraries.
It is also designed to function using only the standard library (in
addition to MicroMoth) in any Python from 2.7 onwards.

The deciding factor is whether Qiskit is available to be imported. If so,
the following external libraries are required dependencies:

qiskit
numpy
scipy
PIL

Otherwise, MicroMoth will be used in place of Qiskit, and alternative
techniques using only the standard library will be used in place of the
other dependencies.

More information on Qiskit can be found at

https://qiskit.org

and information on MicroMoth can be found at

https://github.com/moth-quantum/MicroMoth
"""

import math
import random

# determine whether qiskit can be used, or whether to default to
# MicrMoth and the standard library
try:
    from qiskit import QuantumCircuit, quantum_info
    from qiskit_aer import AerSimulator
    from qiskit_aer.library import SaveStatevectorDict
    simple_python = False
except:
    print('Unable to import Qiskit, so MicroMoth will be used instead')
    from micromoth import QuantumCircuit, simulate
    simple_python = True

    
# this is overwritten by the PIL class if available
class Image():
    """
    A minimal reimplementation of the the PIL Image.Image class, to allow all
    image based tools to function even when only the standard library is
    available.
    
    To initialize an Image oject, use the `newimage` function.
    
    Attributes:
        mode (str): If L, pixel values are a single integer. If 'RGB', they
            are a tuple of three integers.
        size (tuple): Specifies width and height.
    """
    def __init__(self):
        self.mode = None
        self.size = None
        self._image_dict = None
    def getpixel(self,xy):
        """
        Returns pixel value at the given coordinate.
        """
        return self._image_dict[xy]
    def putpixel(self, xy, value):
        """
        Sets the pixel value at the given coordinate.
        """
        self._image_dict[xy] = value
    def todict(self):
        """
        Returns dictionary of pixel values with coordinates as keys.
        Not present in PIL version.
        """
        return self._image_dict
    def show(self):
        """
        If the PIL version of this class is used, this function creates a PNG
        image and displays it. This version instead simply prints all
        coordinates and pixel values.
        """
        for x in range(self.size[0]):
            for y in range(self.size[1]):
                print('('+str(x)+','+str(y)+')'+': '+str(self._image_dict[x,y]))
    def resize(self, new_size, method):
        print("This functionality has not been implemented.")

# this is overwritten by the PIL function if available               
def newimage(mode, size):
    """
    A minimal reimplementation of the the PIL Image.new function.
    Creates an Image object for the given mode and size.
    """
    img = Image()
    img.mode = mode
    img.size = size
    if mode=='L':
        blank = 0
    elif mode=='RGB':
        blank = (0,0,0)
    img._image_dict = {(x,y):blank\
                for x in range(size[0])\
                for y in range(size[1])}
    return img

# if external libraries can be used, import the ones we need
if not simple_python:
    import numpy as np
    from scipy.linalg import fractional_matrix_power
    from PIL.Image import new as newimage, Image


def _kron(vec0,vec1):
    """
    Calculates the tensor product of two vectors.
    """
    new_vec = []
    for amp0 in vec0:
        for amp1 in vec1:
            new_vec.append(amp0*amp1)
    return new_vec


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


def circuit2probs(qc):
    """
    Runs the given circuit, and returns the resulting probabilities.
    """
    if simple_python:
        probs = simulate(qc,get='probabilities_dict')
    else:
        qc_run = qc.copy()
        qc_run.append(SaveStatevectorDict(qc.num_qubits),qc.qregs[0])
        rawamps = AerSimulator().run(qc_run,shots=1).result().data()['statevector_dict']
        probs = {}
        for string, amp in rawamps.items():
            probs[str(bin(int(string,16))[2::].zfill(qc.num_qubits))] = np.abs(amp)**2
    
    return probs


def _image2heights(image):
    """
    Converts an rgb image into a list of three height dictionaries, one for
    each colour channel.
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
        for j in range(int(float(len(line))/2)):
            line[j] += '0'
        # and a '1' for the second half
        for j in range(int(float(len(line))/2),int(len(line))):
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


def make_strip(offset, Lx, rotate=False):
    """
    Creates a dictionary that provides bit strings corresponding to
    points within a grid along a strip of width Lx and height Ly.
    At each height y, each row of the strip is offset by an amount
    offset[y].
    
    Args:
        offset (list): List of horizontal offsets.
        Lx (int): Width of the lattice.
        rotate (bool): If True, the above will refer to vertical
        offsets and height of the lattice, respectively.
    
    Returns:
        grid (dict): Dictionary whose values are points on grid
        whose size is defined by the input. The corresponding keys
        are unique bit strings such that neighbouring strings differ
        on only one bit.
        n (int): Length of the bit strings
        
    """
    # get the required input info
    Ly = len(offset)

    # make the lines
    line_x = make_line( Lx )
    line_y = make_line( Ly )

    # make the grid
    grid = {}
    for y in range(Ly):
        for x in range(offset[y], offset[y]+Lx):
            string = line_x[x%Lx]+line_y[y]
            if not rotate:
                grid[string] = (x,y)
            else:
                grid[string] = (y,x)

    # determine length of the bit strings
    n = len(line_x[0]+line_y[0])
            
    return grid, n


def height2circuit(height, log=False, eps=1e-2, grid=None):
    """
    Converts a dictionary of heights (or brightnesses) on a grid into
    a quantum circuit.
    
    Args:
        height (dict): A dictionary in which keys are coordinates
            for points on a grid, and the values are positive numbers of
            any type.
        log (bool): If given, a logarithmic encoding is used.
        grid (dict): By defeault, a grid generated by `make_grid` is
            used. However, a custom grid can be provided.
            
    Returns:
        qc (QuantumCircuit): A quantum circuit which encodes the
            given height dictionary.
    """
    # get bit strings for the grid
    Lx,Ly = _get_size(height)
    if grid == None:
        grid, n = make_grid(Lx,Ly)
    else:
        n = len(list(grid.keys())[0])
    
    # create required state vector
    state = [0]*(2**n)
    if log:
        # normalize heights
        max_h = max(height.values())
        height = {pos:float(height[pos])/max_h for pos in height}
        # find minimum (not too small) normalized height
        min_h = min([height[pos] for pos in height if height[pos] > eps])
        # this minimum value defines the base
        base = 1.0/min_h
    for bitstring in grid:
        (x,y) = grid[bitstring]
        if (x,y) in height:
            h = height[x,y]
            if log:
                state[ int(bitstring,2) ] = math.sqrt(base**(float(h)/min_h))
            else:
                state[ int(bitstring,2) ] = math.sqrt( h )
    state = normalize(state)
        
    # define and initialize quantum circuit            
    qc = QuantumCircuit(n)
    if simple_python:
        # micromoth style
        qc.initialize(state)
    else:
        qc.initialize(state, range(n))
    qc.name = '('+str(Lx)+','+str(Ly)+')'

    return qc


def probs2height(probs, size=None, log=False, grid=None, min_h=None):
    """
    Extracts a dictionary of heights (or brightnesses) on a grid from
    a set of probabilities for the output of a quantum circuit into
    which the height map has been encoded.
    
    Args:
        probs (dict): A dictionary with results from running the circuit.
            With bit strings as keys and either probabilities or counts as
            values.
        size (tuple): Size of the height map to be created. If not given,
            the size is deduced from the number of qubits (assuming a
            square image).
        log (bool): If given, a logarithmic decoding is used.
        grid (dict): By defeault, a grid generated by `make_grid` is
            used. However, a custom grid can be provided.
            
    Returns:
        height (dict): A dictionary in which keys are coordinates
            for points on a grid, and the values are floats in the
            range 0 to 1.
    """
    
    # get grid info
    if size:
        (Lx,Ly) = size
    else:
        Lx = int(2**(len(list(probs.keys())[0])/2))
        Ly = Lx
    if grid == None:
        grid, n = make_grid(Lx,Ly)
    
    # set height to probs value, rescaled such that the maximum is 1
    max_h = max( probs.values() )   
    height = {(x,y):0.0 for x in range(Lx) for y in range(Ly)}
    for bitstring in probs:
        if bitstring in grid:
            height[grid[bitstring]] = float(probs[bitstring])/max_h
         
    # take logs if required
    if log:
        min_h = min([height[pos] for pos in height if height[pos] > 1e-100])
        base = 1/min_h
        for pos in height:
            if height[pos]>1e-100:
                height[pos] = max(math.log(height[pos]/min_h)/math.log(base),0)
            else:
                height[pos] = 0.0
                        
    return height

    
def circuit2height(qc, log=False, grid=None):
    """
    Extracts a dictionary of heights (or brightnesses) on a grid from
    the quantum circuit into which it has been encoded.
    
    Args:
        qc (QuantumCircuit): A quantum circuit which encodes a height
            dictionary. The name attribute should hold the size of
            the image to be created (as a tuple cast to a string).
        log (bool): If given, a logarithmic decoding is used.
        grid (dict): By defeault, a grid generated by `make_grid` is
            used. However, a custom grid can be provided.
            
    Returns:
        height (dict): A dictionary in which keys are coordinates
            for points on a grid, and the values are floats in the
            range 0 to 1.
    """
    
    probs = circuit2probs(qc)
    try:
        # get size from circuit
        size = eval(qc.name)
    except:
        # if not in circuit name, infer it from qubit number
        L = int(2**(qc.num_qubits/2))
        size = (L,L)
    return probs2height(probs, size=size, log=log, grid=grid)


def combine_circuits(qc0,qc1):
    """
    Combines a pair of initialization circuits in parallel
    Creates a single register circuit with the combined number of qubits,
    initialized with the tensor product state.
    """

    warning = "Combined circuits should contain only initialization."

    # create a circuit with the combined number of qubits
    num_qubits = qc0.num_qubits + qc1.num_qubits
    combined_qc = QuantumCircuit(num_qubits)

    # extract statevectors for any initialization commands
    kets = [None,None]
    for j,qc in enumerate([qc0, qc1]):
        for gate in qc.data:
            if simple_python:
                assert gate[0]=='init', warning
                kets[j] = gate[1]
            else:
                assert gate[0].name=='initialize', warning
                kets[j] = gate[0].params

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
        if simple_python:
            combined_qc.initialize(ket)
        else:
            combined_qc.initialize(ket,range(num_qubits))
    
    # prevent circuit name from being used for size determination
    combined_qc.name = 'None'
            
    return combined_qc


def partialswap(combined_qc, fraction):
    """
    Apply a partial swap to a given combined circuit (made up of two equal
    sized circuits combined in parallel) by the given fraction.
    """
    num_qubits = int(combined_qc.num_qubits/2)
    
    if not simple_python:
        U = np.array([
        [1, 0, 0, 0],
        [0, 0, 1, 0],
        [0, 1, 0, 0],
        [0, 0, 0, 1]
        ])
        U = fractional_matrix_power(U,fraction)
    for q in range(num_qubits):
        q0 = q
        q1 = num_qubits + q
        if not simple_python:
            combined_qc.unitary(U, [q0,q1],\
                                 label='partial_swap')
        else:
            combined_qc.cx(q1,q0)
            combined_qc.crx(math.pi*fraction,q0,q1)
            combined_qc.cx(q1,q0)  

            
def probs2marginals(combined_qc, probs):
    """
    Given a probability distribution corresponding to a given combined
    circuit (made up of two equal sized circuits combined in parallel),
    this function returns the two marginals for each subcircuit.
    """
    num_qubits = int(combined_qc.num_qubits/2)
    
    marginals = [{},{}]
    for string in probs:
        substrings = [string[0:num_qubits], string[num_qubits::]]
        for j,substring in enumerate(substrings):
            if substring in marginals[j]:
                marginals[j][substring] += probs[string]
            else:
                marginals[j][substring] = probs[string]
    
    return marginals


def swap_heights(height0, height1, fraction, log=False, ):
    """
    Given a pair of height maps for the same sized grid, a set of partial
    swaps is applied between corresponding qubits in each circuit.
    
    Args:
        height0, height1 (dict): Dictionaries in which keys are coordinates
            for points on a grid, and the values are floats in the range 0
            to 1.
        fraction (float): Fraction of swap gates to apply.
        log (bool): If given, a logarithmic decoding is used.
            
    Returns:
        new_height0, new_height1 (dict): As with the height inputs.
    """

    assert _get_size(height0)==_get_size(height1), \
    "Objects to be swapped are not the same size"   
    
    # set up the circuit to be run
    circuits = [height2circuit(height) for height in [height0,height1]]
    combined_qc = combine_circuits(circuits[0], circuits[1])
    partialswap(combined_qc, fraction)
    
    # run it an get the marginals for each original qubit register
    p = circuit2probs(combined_qc)           
    marginals = probs2marginals(combined_qc, p)     
    
    # convert the marginals to heights
    new_heights = []
    for j,marginal in enumerate(marginals):
        new_heights.append( probs2height(marginal,size=eval(circuits[j].name),log=log) )
        
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


def swap_images(image0, image1, fraction, log=False):
    """
    Given a pair of same sized grid images, a set of partial swaps is applied
    between corresponding qubits in each circuit.
    
    Args:
        image0, image1 (Image): RGB encoded images.
        fraction (float): Fraction of swap gates to apply.
        log (bool): If given, a logarithmic decoding is used.
            
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

def image2circuits(image, log=False, grid=None):
    """
    Converts an image to a set of three circuits, with one corresponding to
    each RGB colour channel.

    Args:
        image (Image): An RGB encoded image.
        log (bool): If given, a logarithmic encoding is used.
        grid (dict): By defeault, a grid generated by `make_grid` is
            used. However, a custom grid can be provided.

    Returns:
        circuits (list): A list of quantum circuits encoding the image.
    """

    heights = _image2heights(image)

    circuits = []
    for height in heights:
        circuits.append( height2circuit(height, log=log, grid=grid) )

    return circuits


def circuits2image(circuits, log=False):
    """
    Extracts an image from list of circuits encoding the RGB channels.

    Args:
        circuits (list): A list of quantum circuits encoding the image.
        log (bool): If given, a logarithmic decoding is used.

    Returns:
        image (Image): An RGB encoded image.
    """

    heights = []
    for qc in circuits:
        heights.append( circuit2height(qc, log=log) )

    return _heights2image(heights)


def row_swap_images(image0, image1, fraction, log=False):
    """
    A variant of `swap_images` in which the swap process is done on each line
    of the images individually, rather than with the images as a whole. This
    makes it much faster.
    
    Args:
        image0, image1 (Image): RGB encoded images.
        fraction (float): Fraction of swap gates to apply.
        log (bool): If given, a logarithmic decoding is used.
            
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


def blur_height(height, xi, axis='x', circuit=None, log=False, grid=None):
    """
    Applies a predetermined blur effect designed for a smooth blur.
    
    Args:
        height (dict): A dictionary in which keys are coordinates
            for points on a grid, and the values are positive numbers of
            any type.
        xi (float): Fraction of pi rotation to apply on the qubit for
            which the largest rotation is aplied
        axis (string): `rx` rotations are used when this is `'x'`, and
            `ry` rotations are used otherwise.
        circuit (QuantumCircuit): Rotations are applied to the given circuit
            if supplied. Otherwise one is made from `height`.
        log (bool): If True, a logarithmic encoding is used.
        grid (dict): By defeault, a grid generated by `make_grid` is
            used. However, a custom grid can be provided.
            
    Returns:
        circuit (QuantumCircuit): Circuit on which the blur effect has been
            added.
    """
    
    # get size and bit strings for the grid
    Lx,Ly = _get_size(height)
    if grid == None:
        grid, n = make_grid(Lx,Ly)
    else:
        n = len(list(grid.keys())[0])
    # invert grid dict to have coords as keys
    coord_grid = {grid[string]:string for string in grid}
    
    rates = [0]*n
    for x in range(Lx):
        for y in range(Ly):
            # for this point, go through all neighbours
            # and find all bits on which address differs
            string  = coord_grid[x,y]
            axes = []
            for (dx,dy) in [(0,1),(0,-1),(1,0),(-1,0)]:
                if (x+dx,y+dy) in coord_grid:
                    nstring  = coord_grid[x+dx,y+dy]
                    for j,b in enumerate(nstring):
                        if b!=string[j]:
                            axes.append(n-j-1)
            # add the height at this point to the rates for each of these
            for j in axes:
                if (x,y) in height:
                    rates[j] += height[x,y]

    # normalize the rates
    max_rate = max(rates)
    for j in range(n):
        rates[j] /= max_rate
        
    # make the circuit the rotation
    qc_rot = QuantumCircuit(n)
    for j in range(n):
        theta = np.pi*rates[j]*np.pi*xi
        if axis=='x':
            qc_rot.rx(theta,j)
        else:
            qc_rot.ry(theta,j)
            
    # add to initial circuit
    if circuit:
        if simple_python:
            circuit = circuit + qc_rot
        else:
            circuit = circuit.compose(qc_rot)
    else:
        if simple_python:
            circuit = height2circuit(height,log=log) + qc_rot
        else:
            circuit = circuit = height2circuit(height,log=log).compose(qc_rot)

    circuit.name = '('+str(Lx)+','+str(Ly)+')'
        
    return circuit


def blur_image(image, xi, circuits=None, axis='x',log=False):
    """
    Applies a predetermined blur effect designed for a smooth blur.
    
    Args:
        image (Image): An RGB encoded image.
        xi (float): Fraction of pi rotation to apply on the qubit for
            which the largest rotation is aplied
        axis (string): `rx` rotations are used when this is `'x'`, and
            `ry` rotations are used otherwise.
        circuits (list): Rotations are applied to the given circuits
            if supplied. Otherwise they are made from `image`.
        log (bool): If True, a logarithmic encoding is used.
            
    Returns:
        circuits (list): Circuits on which the blur effect has been
    """       

    heights = _image2heights(image)
    
    if circuits==None:
        circuits=[None,None,None]
    
    for j,height in enumerate(heights):
        circuits[j] = blur_height(height, xi, axis=axis, circuit=circuits[j], log=log)
        

    return circuits


def dotdot(L,diamond=0,delta=0,depth=0):
    """
    Creates a circuit for a heightmap depicting the letteer 'I'.
    The circuit is designed to be low-depth and compatible with near-term devices.
    
    Args:
        L (int): Linear size of the image. Must be at least 32.
        diamond (float): Degree to which filling is done with a diamond effect (maximum at 1).
        delta (float): Strength of random variations to all rotations made (no randomness when 0).
        depth (int): Depth of a final random circuit using rotations, using rotations of strength
            delta and entangling gates.
            
    Returns:
        qc (QuantumCircuit): Quantum circuit encoding the height map
        line (list): A list of qubits, ordered such that all entangling gates are between neighbours.
    """
    
    def dt():
        if delta!=0:
            return delta*(0.5-random.random())*math.pi
        else:
            return 0
    
    # make sure L is a compatible value (a power of 2 and >= 32)
    assert L>=32, 'Size must be at least 32x32 pixels.'
    L = int(2**math.ceil(math.log(L)/math.log(2)))

    # qubits per register for I shape
    ni = 5

    # qubits per register for fill
    nf = int(math.log(L/2**ni)/math.log(2))

    # total qubits per register
    nr = ni + nf
    # and total
    n = 2*nr

    # entangling gates on the following pairs are required
    pair_extend = (1,nr+2)
    pairs_fill = [(ni+j,ni+nr+j) for j in range(nf)]

    # we assume qubits are on a line such that these pairs are neigbours
    line = list(pair_extend)
    for pair in pairs_fill:
        line += list(pair)
    for j in range(n):
        if j not in line:
            line.append(j)

    # two sets of disjoint pairs cover this line
    r1 = [(line[j-1],line[j]) for j in range(1,n,2)]
    r2 = [(line[j-1],line[j]) for j in range(2,n,2)]

    qc = QuantumCircuit(n)

    # fill
    fill_qubits = []
    for c,t in pairs_fill:
        qc.ry(math.pi/2+dt(),c)
        # for full fill, use theta=pi/2 here
        theta = (1-diamond)*math.pi/2
        qc.rx(theta,t)
        # otherwise, this cx gives a diamond effect
        if diamond!=0:
            qc.cx(c,t)

    # horizontal lines
    qc.rx(math.pi+dt(),ni-1)
    for j in [0,2,3]:
        qc.rx(math.pi/2+dt(),j)
    qc.rx(0.39*math.pi+dt(),1)

    # central bar
    qc.rx(math.pi+dt(),nr+1)
    qc.ry(math.pi/2+dt(),nr+4)
    qc.ry(math.pi/2+dt(),nr)
    qc.ry(math.pi/2+dt(),nr+3)

    # extend at top
    qc.rx(math.pi+dt(),1)
    if simple_python:
        qc.crx(math.pi/2,1,nr+2)
    else:
        qc.ch(1,nr+2)
    qc.x(1)

    # cover with cz gates to entangle everything
    for c,t in r1+r2:
        qc.cz(c,t)
        
    if delta!=0:
        for _ in range(depth):
            for j in range(n):
                qc.rx(dt(),j)
            for c,t in r1+r2:
                if simple_python:
                    qc.h(t)
                    qc.cx(c,t)
                    qc.h(t)
                else:
                    qc.cz(c,t)

    return qc, line
