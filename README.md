# `phq`: Temperature Dependent Phonon Quasiparticle Dispersions of Complex Crystals from First Principles

[TOC]

This repo is created by [Zhen Zhang](mailto:zhenzhang305@hotmail.com). The associated paper is [here](https://www.sciencedirect.com/science/article/pii/S0010465519301523). The program files DOI is: [http://dx.doi.org/10.17632/sk4jsjc6p9.1](https://dx.doi.org/10.17632/sk4jsjc6p9.1).

## Installation Instructions

The `phq` program can be installed on Linux, macOS and Windows systems by gfortran, `ifort` or other fortran compilers. ifort is recommended for faster running speed and smaller memory requirement. The source codes are in the `src` folder and the makefiles are in the system folder. Simply copy one of the suitable makefiles from system folder and paste it in the src folder and rename it as `makefile`. Run:

```shell
$ make
```

The phq executable file is compiled.

The source codes work as following:

1. `configure.f90` contains 3 auxiliary modules, which configure the physical quantities, parameters and text readin settings.

2. `main.f90` the main module, which performs the normal mode projection and anharmonic phonon extraction and so forth.

3. `phq.f90`  runs the program and records the running time.

## Input Files

There are 4 input files needed, which can be prepared according to the example provided:

1. `input` general controlling and settings. Parameters:

    `dt` MD time step in atomic unit.  `step_md_use` number of MD steps needed for the calculation of mode-projected VAFs. `correlation_time` desired decay time for VAFs in unit of `dt` which should be less than the `step_md_use`. `pole` parameter used in the maximum entropy method to filter high frequency components which should be less than the ``step_md_use``. `supercell` size of the supercell. `temperature` MD simulation temperature. `method` choose renormalized frequencies by which method when building the effective harmonic dynamical matrices: `0` represents the fitting approach. `1` represents the FT. `2` represents the MEM. If `method` is not specified or integers other than `1` and `2` are entered, the fitting approach will be used, which is the default method.

2. `scf.out` structure information of the primitive cell from self-consistent calculation. Parameters:

   `ntype` number of elements.  `natom` number of atoms. `mass` atomic mass of each element following the symbol of the element. `lattice_parameter` scale of lattice vectors in unit of Bohr radius. `cell_parameters` lattice vectors in cartesian cooradinates in unit of `lattice_parameter`. `atomic_positions` atomic positions of atoms in reduced coordinates following the symbol of element.

3. `dyn.out` harmonic phonon results from Quantum ESPRESSO ph.x output. Parameters:

    `q` **q**-point in cartesian coordinates in unit of 2$\pi$/`lattice_parameter`.They are required as program input. `freq` harmonic phonon frequencies. They are required as program input as well as the following six columns of eigenvectors. Eigenvector of each atom has x, y and z components and each component has a real part and an imaginary part. The order of atoms should be the same as that of the atoms entered in `atomic_positions` in `scf.out`. `Dynamical Matrix` dynamical matrices, which are not required as program input if other program is used. `Dielectric Tensor` whether this is needed depends on whether LO-TO splitting needs to be considered in the system. `Effective Charges` whether they are needed depends on whether LO-TO splitting needs to be considered in the system.

4. `md.out` MD information. Parameters:

   `total_step` total actual MD steps. In practice, recorded MD steps should be configurations after reaching thermal equilibrium and less than this number.`atomic_positions` initial atomic positions in reduced coordinates of the supercell lattice vectors. `md_step` recorded MD steps, followed by the instantaneous `atomic_md_positions` in the MD simulation. `atomic_md_positions` atomic positions in the MD simulation in reduced coordinates of the supercell lattice vectors. Atoms in each one of the primitive cell should be together instead of atoms of the same element being together, and the order of atoms in each of the primitive cell should be in the same order as provided in `scf.out`.

## Executing the Program

After putting the phq executable in the same folder with the input files, run:

```shell
$ ./phq < input
```

## Output Files

The main output files of `phq` are as following:

1. `corr.vaf` velocity autocorrelation functions of each normal mode .
2. `corr_fit.vaf` fitted curves of correlation functions according to phonon quasiparticle and anharmonic perturbation theory.
3. `corr_fourier.vaf` Fourier transformation of correlation functions.
4. `frequency.freq` harmonic, fitted renormalized, Fourier transformed renormalized and maximum entropy method renormalized phonon frequencies.
5. `tau_fit.tau` fitted phonon quasiparticles' lifetime.
6. `vector_q.out` eigenvectors of the primitive cell.
7. `vector.out` eigenvectors of the supercell.
8. `harmonic_matrix.mat` harmonic force constants.
9. `gamma_matrix.mat` effective harmonic force constant matrix.
10. `dynamical_matrix_md.mat` effective harmonic dynamical matrices.
11. `dynmatmd` renormalized phonon information with the same format of `ph.x` executable's output in the Quantum ESPRESSO suite. All post processing can be started from here.

## Post Processing

The output `dynmatmd` files from phq can be read in by the `q2r.x` executable in the [Quantum ESPRESSO](https://www.quantum-espresso.org/) suite. Effective harmonic force constants, renormalized phonon dispersions, anharmonic entropy and free energy can be obtained from `q2r.x` output. If other Fourier interpolation program is used, either use the effective harmonic dynamical matrices file `dynamical_matrix_md.mat` or use the effective harmonic force constant matrix file `gamma_matrix.mat` to fit the format.

## Example

There is an example of diamond silicon in the example folder. Postprocessing files are also in the silicon example folder.

## License

[GNU General Public License v3](./LICENSE.txt)

