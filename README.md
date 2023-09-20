# DAQC_simulator
Classical simulator for Digital-Analog Quantum Computing circuits. Writen in Matlab.

This repository provides 3 simulators. A digital quantum simulator, a banged-DAQC simulator and a stepwise-DAQC.

The programs here are a direct implementation of the methods described in:
Paula García-Molina, Ana Martin, Mikel Garcia de Andoin and Mikel Sanz, "Noise in digital and digital-analog quantum computation", arXiv::2107.12969 (2022)
https://arxiv.org/abs/2107.12969

The functions for the simulators are given in the x_simulator.m files.
The inputs can be handcrafted DAQC schedules.
Another way of simulating DAQC circuits is to generate a DAQC schedule from a digital circuit, by using QASM_parser.m.
QASM parser generates a digital circuit description that can be simulated with D_simulator.m.
Conversions between different types of descriptions can be made using the xxx2yyy.m functions.

For any question, contact mikelgda@gmail.com

Licensed under CC BY 4.0
