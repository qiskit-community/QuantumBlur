# QuantumBlur

A tool for doing quantum things with height maps and images. All explained in [this blog](https://medium.com/qiskit/introducing-procedural-generation-using-quantum-computation-956e67603d95) and the corresponding [paper](https://arxiv.org/abs/2007.11510).

In this respository is a Python version of the source code, but there is also a [Unity implementation](https://github.com/TigrisCallidus/QuantumBlurUnity/blob/master/README.md).

## Requirements

Either of the following:
* Python 3.7 and beyond with Qiskit, NumPy, SciPy and PIL.
* Any Python from 2.7 and beyond, with only the standard library and [MicroMoth](https://github.com/moth-quantum/MicroMoth) (A minimal quantum SDK, previously known as [MicroQiskit](https://github.com/qiskit-community/MicroQiskit));

The former is recommended, but the latter is more flexible to running everything in strange places.

## Installation

You just need the file [quantumblur.py](quantumblur/quantumblur.py) somewhere importable. You can do this with a simple copy/paste, but you can also pip install this repository with the following command.
```
pip install git+https://github.com/qiskit-community/QuantumBlur.git
```

## How to use

See the [quick start guide](QuickStart.ipynb).
