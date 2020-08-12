from quantumblur import *

def height_test():
    
    eps = 0.01

    # heights that we'll test
    h0 = {(0,0):1,(0,1):0,(1,0):0,(1,1):1}
    h1 = {(0,0):0,(0,1):1,(1,0):1,(1,1):0}

    # correct answers for the rotation we'll do
    rotated = {}
    rotated[0] = {(0, 0): 1.0, (0, 1): 0.2984464104095249, (1, 0): 0.2984464104095249, (1, 1): 1.0}
    rotated[1] = {(0, 0): 0.2984464104095249, (0, 1): 1.0, (1, 0): 1.0, (1, 1): 0.2984464104095249}

    for j,h in enumerate([h0,h1]):

        qc = height2circuit(h)

        # check that the height comes out of the circuit unchanged
        new_h = circuit2height(qc)
        for pos in h:
            assert abs(h[pos]-new_h[pos])<eps

        # check that the rotation works correctly
        for q in range(qc.num_qubits):
            qc.rx(0.5,q)  
        new_h = circuit2height(qc)
        for pos in h:
            assert abs(rotated[j][pos]-new_h[pos])<eps

    # correct answers for the swap we'll do
    swapped = {}
    swapped[0] = {(0, 0): 1.0, (0, 1): 0.33333333333333315, (1, 0): 0.33333333333333315, (1, 1): 1.0}
    swapped[1] = {(0, 0): 0.33333333333333315, (0, 1): 1.0, (1, 0): 1.0, (1, 1): 0.33333333333333315}

    # check that swap works correctly
    new_h,_ = swap_heights(h0, h1, 1.0/3)
    for pos in h0:
        assert abs(swapped[0][pos]-new_h[pos])<eps
    
    print(':)')


height_test()