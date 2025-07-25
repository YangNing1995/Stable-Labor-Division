# Stable Division of Labor

This repository contains the code for our paper:
**"Synthetic Conscription: Achieving Stable Labor Division with Rational Design of Gene Circuitry"**
*Zhi Sun, Ning Yang, Baiyi Jiang, Kang Xia, Lei Zhang, Weijie Li, Shanwen Chen, Chunbo Lou, Chao Tang, Xiaojing Yang*

## Abstract
Division of labor within populations is a game-changer for complex functions in synthetic biology, dramatically bolstering system efficiency, stability, and robustness. However, different functions inevitably introduce different burden, resulting different fitness of different functional subgroups, which destabilizes the system. Inspired by the real-world “conscription” systems, we propose a “synthetic conscription” strategy: by resetting the population to a homogeneous state for labor re-division after stimulus response, the whole population would collectively share the burden that division of labor might bring to subgroups, thereby ensuring long-term stability. To achieve this goal, we first formalized the steady-state and kinetic requirements for “synthetic conscription”. We then exhaustively enumerated all feasible two-node transcriptional regulatory topologies satisfying these criteria. Leveraging parameter configurations enabling essential dynamics identified through this systematic search, we engineered a gene circuit implementing the optimal topology. The resulting circuit enables robust, high-fidelity division of labor across repeated stimulus cycles. Stochastic simulations coupled with global sensitivity analysis identified critical parameters governing the functional details of the conscription strategy, with experimental validation. Finally, integrating this circuit into probiotic Escherichia coli Nissle 1917 enabled targeted drug delivery and effective mitigation of inflammatory bowel disease (IBD). This proof-of-concept demonstration underscores the potential of synthetic conscription for chronic disease intervention.

## Repository Structure

### Topology Enumeration
This folder contains code for the systematic enumeration of gene circuit topologies. The enumeration algorithms explore all possible network architectures to identify circuits capable of achieving stable division of labor in synthetic biological systems.

### Stochastic Simulation
This folder contains the code for stochastic simulations of gene circuit dynamics, including the fitting process of experimental data, bifurcation diagram analysis, and parameter swap simulations.

## Software Requirements
This repository uses both Python and MATLAB for circuit design, simulation, and analysis.

### MATLAB
- MATLAB R2024a 

### Python Packages
The Python components rely on the following libraries:

- `intersect==1.2`
- `jupyterlab==4.4.4`
- `matplotlib==3.10.3`
- `numpy==2.2.6`
- `scipy==1.15.3`
- `seaborn==0.13.2`

You can install all Python dependencies using:

```bash
pip install intersect==1.2 jupyterlab==4.4.4 matplotlib==3.10.3 numpy==2.2.6 scipy==1.15.3 seaborn==0.13.2
