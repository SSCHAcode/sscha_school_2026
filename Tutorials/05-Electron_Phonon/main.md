# Hands-on Session 8: EPIq - Anharmonicity in Electron-Phonon Coupling Related Properties

## Introduction

In this hands-on session, we learn how to include anharmonic effects calculated within the SSCHA in the calculation of electron-phonon coupling-related properties using [EPIq](https://the-epiq-team.gitlab.io/epiq-site/).

![SSCHA EPIq Logo](figures_08/SSCHA_epiq_Logo.png)

In some systems, the first-principles calculation of electron-phonon coupling matrix elements can be demanding. EPIq (Electron-Phonon Wannier Interpolation over k- and q-points) is an open-source software that allows speeding up the calculation of electron-phonon coupling-related properties using the Wannier interpolation technique. Details on the interpolation scheme can be found [here](https://the-epiq-team.gitlab.io/epiq-site/docs/th_found). Within EPIq, it is possible to include anharmonic corrections to the dynamical matrices as calculated within the SSCHA.

## Requirements

In the interest of time, the following starting data are at your disposal for this hands-on session:

1. Electron-phonon matrix elements $g^{\nu}_{m,n}(\mathbf{k},\mathbf{q})$ computed from first principles.
2. Wannier interpolation files encoding the transformation to the optimally smooth subspace, $U_{mn}$: $\ket{\mathbf{R}n} = \frac{1}{\sqrt{N_k^w}} \sum_{\bf k=1}^{N_k^w}\sum_{m=1}^{N_{\rm w}}e^{-i\mathbf{k}\cdot \mathbf{R}}U_{mn}(\mathbf{k})|\psi_{{\bf k}m}\rangle$
3. Anharmonic dynamical matrices $D^{SCHA}_{\mu,\nu}(\mathbf{k},\mathbf{q})$.
4. Harmonic dynamical matrices (as a reference) $D^{HARM}_{\mu,\nu}(\mathbf{k},\mathbf{q})$.

> **Attention:**
> A folder prepared for you with these data for the present tutorial can be downloaded from the `08_EPIq` folder in the shared cloud. Right-click on `tutorial_data` on the navigation bar and download the whole folder. It contains:
>
> 1. The electron-phonon coupling matrix elements can be found in the `mat_elem` folder. Each file corresponds to a different $\mathbf{q}$-point in the first Brillouin zone.
> 2. The `Wannier` folder contains the files (`.eig`, `.chk`).
> 3. The dynamical matrices are stored in the `dyn_mat` directory. `dynq*` files are harmonic while `MoS2.Hessian.dyn*` are anharmonic dynamical matrices computed with the SSCHA code.
>
> Place all the downloaded material where you intend to run the tutorial. A suggested structure is, for example:
>
> ```console
> handson_8/
>     |
>     + ----- epiq/
>     |   |
>     |   + ----- bin/
>     |   |
>     |   + ----- src/
>     |
>     |
>     + ----- MoS2/
>         |
>         + ----- MoS2.eig
>         |
>         + ----- MoS2.chk
>         |
>         + ----- mat_elem/
>         |   |
>         |   + ----- MoS2_elph.mat.1_q*
>         |
>         + ----- dyn_mat/
>             |
>             + ----- dynq*
>             |
>             + ----- MoS2.Hessian.dyn*
> ```

## About EPIq

![EPIq site and paper QR codes](figures_08/qr-code-epiq-both.png)

The Electron-Phonon Interpolation over q package exploits Wannier interpolation to obtain many properties in solids. Bloch's theorem allows describing an infinite system in real space with a continuum of Hamiltonians for different k-points in the Brillouin zone.

The properties of a material refer to a certain observable averaged over the sample. Thanks to Bloch's theorem, it is often convenient to perform the average in reciprocal k-space, in other words, as a sum over the whole Brillouin zone of the quantities defined at each k-point.

$$
\langle O \rangle = \frac{1}{N_k}\sum_{\mathbf{k}} F_O(\mathbf{k})
$$

The quality of the averaging depends on the fineness of the sampling $N_k$ and the smoothness of the integrand function $F_O(\mathbf{k})$ (a constant function is fully described with just one k-point).

![Sampling illustration](figures_08/sampling.png)

Wannier interpolation is an efficient way to refine the Brillouin zone sampling.

![Wannier interpolation illustration](figures_08/wannier_interpolation.png)

### Calculations Available in EPIq

1. Adiabatic (static) and non-adiabatic (dynamic) force constant matrices.
2. **Electron-phonon contribution to the phonon linewidth** and related quantities.
3. Isotropic and **anisotropic Eliashberg equations**.
4. Double Resonant Raman scattering.
5. Electron lifetime and relaxation time.

### EPIq Workflow

![EPIq workflow](figures_08/workflow.png)

The core steps of any calculation employing EPIq are performed using the main executable `epiq.x` and consist of two main stages:

1. A preliminary step where electron-phonon coupling matrix elements and the Hamiltonian are Fourier-transformed to real space and written to file.
2. The electron-phonon coupling matrix elements and the Hamiltonian are transformed back to reciprocal space to compute the property of interest at arbitrary k- and q-values.

### EPIq Input File

The input file is divided into three namelists:

1. `&control`, specifying what calculation the code will perform
2. `&electrons`, specifying the electronic parameters of the property to be computed
3. `&phonons`, specifying the phonon parameters for the property to be computed

Finally, the last lines of the input indicate the electron momentum (k-) and phonon momentum (q-) meshes on which the matrix elements are interpolated.

## Let's Practice: Calculation of Electron-Phonon Coupling Related Properties for Doped Monolayer ${\rm MoS}_2$

![Monolayer MoS2](figures_08/mos2.png)

In this tutorial, we calculate electron-phonon coupling-related properties for doped monolayer ${\rm MoS}_2$ and evaluate the effect of anharmonicity on them.

### Wannier Interpolation

As explained in the previous section, any calculation within EPIq starts with a preliminary step called `dump`. During this step, the Hamiltonian and the electron-phonon coupling matrix elements in real space are computed and written to file, to be used in any subsequent calculation. This step is performed only once. The input file is the following:

```fortran
&control
    dump_gR=.true.,
    prefix='MoS2',
    elphmat_dir="./mat_elem/"
  /
  &electrons
  /
  &phonons
    nq1=8,
    nq2=8,
    nq3=1,
  /
```

Notice, in particular, the following parameters:

- `dump_gR` in the namelist `control` tells the program to save the auxiliary file containing the Hamiltonian and the electron-phonon coupling matrix elements in real space. Since we are not calculating any property of the system, the namelists `electrons` and `phonons` are essentially empty. Only `nq1,nq2,nq3` must be supplied to specify the q-points mesh where the input electron-phonon coupling matrix elements were computed (in this case, an 8x8x1 q-grid).

> **Note:** To keep everything tidy, you can keep the el-ph matrix elements in a separate folder and use the variable `elphmat_dir="<path-to-folder>"` in `&control`.

Run this preliminary calculation using:

```console
mpirun -n <NPROC> {$path_to_epiq}/bin/epiq.x -inp dump.in > dump.out
```

After the dump has finished, some output files have been produced.

- The output file (binary, named `G_and_H.bin`) contains the Hamiltonian and the electron-phonon coupling matrix elements in the Wannier representation. If `ascii_G_and_H=.true.` is added in the `&control` namelist, then the produced output file is readable and named `G_and_H.asc`.
- `.dat` files containing the real space localization of the Hamiltonian and the electron-phonon coupling matrix elements.

### Check Real Space Localization

If the Wannier transformation is well converged, the matrix elements are optimally localized in real space. Always check their localization.

> **Exercise:**
>
> Using the two-column files:
>
> - `MoS2_gw_R_ph.mu*.dat.pe_1`: $|R| \qquad \sum_{m,n} | g^{\nu}_{m,n}(0,\mathbf{R})|^2$
> - `MoS2_gw_R_el.mu*.dat.pe_1`: $|r| \qquad \sum_{m,n} | g^{\nu}_{m,n}(\mathbf{r},0)|^2$
>
> Plot the averaged modulus of the electron-phonon matrix elements as a function of the distance in real space $|R|$. Are they localized?

### Phonon Linewidth Calculation

We will now focus on one of the system properties that can be calculated with EPIq: the phonon linewidth $\gamma_{\mathbf{q},\nu}$. We will consider the Allen phonon linewidth, defined by the following equation:

$$
\gamma_{\mathbf{q},\nu} = \frac{4 \pi \omega_{\mathbf{q},\nu}}{N_k}\sum_{m,n}\sum_{\mathbf{k}}|g^{\nu}_{m,n}(\mathbf{k},\mathbf{q})|^2\delta(\epsilon_{\mathbf{k}+\mathbf{q},m}-\epsilon_{F})\delta(\epsilon_{\mathbf{k},n}-\epsilon_{F})
$$

where the electron-phonon coupling is defined from the deformation potential as:

$$
g^{\nu}_{m,n}(\mathbf{k},\mathbf{q}) = \sum_s\mathbf{e}^s_{\mathbf{q},\nu}\cdot\mathbf{d}^s_{m,n}(\mathbf{k},\mathbf{q})/\sqrt{2M_s\omega_{\mathbf{q},\nu}}
$$

#### Example of Input File for Linewidth Calculation of Monolayer ${\rm MoS}_2$

We first calculate the mode-resolved $\gamma$ at the M-point of the Brillouin zone. In the following example, we perform the calculation for a phonon of momentum $\mathbf{q}=\mathbf{M}$ for two values of the electronic smearing. The input parameters are explained in detail in the [EPIq manual](https://the-epiq-team.gitlab.io/epiq-site/docs/manual/). Here is the input file:

```fortran
&control
  prefix='MoS2',
  calculation='ph_linewidth',
  read_dumped_gr=.true.,
  dump_gR=.false.,
  elphmat_dir="./mat_elem/"
  out2json=.true.
 /

 &electrons
      ngauss=0,
      sigma_min=0.01,
      sigma_max=0.05
      nsigma=5,
  /

  &phonons
      use_alternative_dyn=.true.
      prefix_alt_dyn='./dyn_mat/dynq'
      Fourier_interp_dyn=.true.,
      nq1=8,nq2=8,nq3=1,
  /

  k-points
  automatic
  4 4 1 0 0 0

  q-points
  crystal
  1
  0.5 0 0 1 ! M
```

The linewidth calculation is then started by:

```console
mpirun -n <NPROC> {$path_to_epiq}/bin/epiq.x -inp lw.in > lw.out
```

Note that this calculation is done within the harmonic approximation, as we are using the dynamical matrices indicated by the variable `prefix_alt_dyn='./dyn_mat/dynq'`.

#### Parameters

- The linewidth calculation is selected by setting the `calculation` parameter equal to `'ph_linewidth'` in the `&control` namelist.
- In the namelist `control`, `read_dumped_gr` is set to `.true.`, indicating that electron-phonon coupling matrix elements can be read from `G_and_H.bin`.
- In the namelist `electrons`, `sigma_min`, `sigma_max`, `ngauss`, and `nsigma` specify the maximum, minimum, type, and number of electronic smearing values to use. `efermi` and `ef_from_input` specify the initial guess for the Fermi level (the Fermi level calculated by Quantum ESPRESSO is usually a good guess) and whether the Fermi level should be recalculated by EPIq (`ef_from_input` equal to `.false.`) or set from input (`ef_from_input` equal to `.true.`). The variable `thr_compute_k` is used to restrict the calculation only to k-points possessing at least one eigenvalue in the specified energy region near the Fermi level.
- In the namelist `phonons`, `Fourier_interp_dyn` equal to `.true.` asks to interpolate the dynamical matrices for the q-points that do not belong to the Wannier grid. Alternatively, EPIq gives the opportunity to read eigenvalues and eigenvectors produced by `matdyn.x` of the Quantum ESPRESSO package (`matdyn.eig` file), by setting `Fourier_interp_dyn=.false.` and `read_modes=.true.` in input, or even to directly read dynamical matrices from a `matdyn.dyn` file, by setting `Fourier_interp_dyn=.false.` and `read_modes=.false.`.

#### Output Files

If the `out2json` flag is set to `.true.`, the file `MoS2_lambda.json` will be produced. It can be automatically parsed using Python as in the following lines, where the variable `q_pts` is a dictionary whose entries are the results of the calculation for each q-point.

```python
import json
ff = './MoS2_lambda.json'
with open(ff,'r') as f:
    q_pts = json.load(f)

first_q = q_pts["1"]

print("Fraction coordinate of the q point:", first_q["xq_frac"])
print("Electronic temperatures used:", first_q["T"])
print("Results for the first mode:", first_q["1"])
print("Frequency of the second mode in meV:", first_q["2"]["freq"])
```

Here is a simple Python script to plot the linewidth estimated with the Allen formula:

```python
import matplotlib.pyplot as plt
import numpy as np
for q in list(q_pts.values())[:]:
    for mod in range(1,10):
        plt.plot( q["T"],q[f"{mod}"]["gamma_allen"],label=r"$\omega_"+f"{mod}$={q[f'{mod}']['freq']:.0f} (meV)",linestyle='--' )
    plt.xlabel('T  (eV)')
    plt.ylabel(r'$\gamma(T)$  (meV)')
    plt.legend(title=f"q={ np.round(np.array(q['xq_frac']),2) }")
    plt.show()
```

The other output file, `MoS2_lambda.d`, contains all the properties calculated by EPIq. The way this file is formatted is specified in the output file `lw.out`. Notice in particular that the first column contains the electronic smearing, while the second column contains the Allen linewidth.

> **Exercise:**
>
> Perform a convergence study of the linewidth at $\mathbf{M}$ as a function of the smearing and the k-mesh density. What is the minimum temperature at which the linewidth of the eighth mode is considered converged with a k-mesh of 8x8x1? And with 12x12x1?

> **Exercise:**
>
> Which mode shows a larger smearing dependence?

### Inclusion of Anharmonicity

Now, we want to observe what changes with the inclusion of **anharmonic effects**. To this aim, we need to correctly specify the prefix of the free-energy Hessian matrices calculated within the SSCHA using `prefix_alt_dyn='./dyn_mat/MoS2.Hessian.dyn'`.

> **Exercise:**
>
> Compute the linewidth at $\mathbf{q}=\mathbf{M}$ using the anharmonic dynamical matrix computed with the SSCHA code. Why does it seem like the first and the second mode are exchanged with respect to the harmonic dynamical matrices? Which modes present a larger anharmonic correction?

### Dispersion Along a Line

Now we want to perform the calculation of phonon linewidth along a certain crystalline direction and produce a plot like this one:

![Phonon linewidth dispersion](figures_08/lw.png)

where the thickness of the lines is proportional to the calculated phonon linewidth. To do this, we repeat the phonon linewidth calculation, this time considering the whole $\Gamma$--M path:

```fortran
&control
   dump_gR=.false.,
   read_dumped_gr=.true.,
   prefix='MoS2',
   calculation='ph_linewidth',
   elphmat_dir="./mat_elem/"
&end
&electrons
   ngauss=0,
   sigma_min=0.01,
   sigma_max=0.02
   nsigma=2,
&end
&phonons
   use_alternative_dyn=.true.
   prefix_alt_dyn='./dyn_mat/MoS2.Hessian.dyn'
   Fourier_interp_dyn=.true.,
   nq1=8,nq2=8,nq3=1,
&end
k-points
automatic
8 8 1  0  0  0
q-points
crystal
11
0.001 0 0 1
0.05 0 0 1
0.10 0 0 1
...
...
0.50 0 0 1
```

Once the linewidth calculation has finished, we can obtain a plottable file using the `linewidth_path.x` post-processing tool. Here is an example of input file:

```fortran
&input_lambda
  prefix='MoS2'
  lkp_sequential=.true.
  sigma_min=0.01,
  sigma_max=0.02,
  nsigma=2,
  chosen_sigma=0.01
&end
 3.159998   0.000000   0.000000
-1.579999   2.736639   0.000000
 0.000000   0.000000  19.141895
 crystal
 11
 0.001 0 0 1
 0.05 0 0 1
 0.10 0 0 1
 ...
 ...
 0.50 0 0 1
```

Notice the parameter `chosen_sigma`, which specifies what smearing will be used to produce the plottable file, and the lattice parameters at the end of the namelist. The post-processing is executed as follows:

```console
{$path_to_epiq}/bin/linewidth_path.x < path.in > path.out
```

Finally, use the following gnuplot script to plot the q-resolved linewidth for the acoustic modes:

```gnuplot
set ylabel '{/Symbol w}(meV)'
set xlabel 'Gamma - M'
set style fill transparent solid 0.25
pl for [i=0:9] 'MoS2_lw_path.d' every 9::i u 1:2 w l lt rgb 'black' title ''
repl  for [i=0:9] 'MoS2_lw_path.d' every 9::i u ($1):($2-$3/2):($2+$3/2)\
w filledc lt rgb 'red' title ''
```

> **Exercise:**
>
> Try to produce two plots of the whole phonon spectrum, comparing the harmonic and the anharmonic results. Do you observe any differences?

## Advanced Tutorial: Migdal-Eliashberg Calculation Using SSCHA Hessian Matrices

EPIq also allows solving the anisotropic Eliashberg equations on the imaginary axis in order to calculate the $\mathbf{k}$-resolved superconducting gap. We give an example for the superconducting gap of doped monolayer $\textrm{MoS}_2$ at T = 1 K.

```fortran
&control
    dump_gR=.false.,
    read_dumped_gr=.true.,
    prefix='MoS2',
    calculation='migdal_eliashberg',
    elphmat_dir="./mat_elem/"
    &end
    &electrons
    theta=1,
    efermi=-2.0,
    ef_from_input=.false.,
&end
&phonons
    use_alternative_dyn=.false.
    prefix_alt_dyn='./dyn_mat/MoS2.Hessian.dyn'
    read_modes=.false.,
    nq1=8,
    nq2=8,
    nq3=1,
&end
&input_migdal
    initialize='step',
    gap_threshold=0.025,
    sigma_me=0.1,
    ME_Fermi_thickness=0.4,
    mustar=0.1,
    nmatsu=64,
    nitermax=100,
    alpha_mix_me=0.5,
    gap_init=4.0,
&end
nkfs
64 64 1
40
```

Run the EPIq calculation as follows:

```console
mpirun -n 4 {$path_to_epiq}/bin/epiq.x <input_ME > out_ME
```

Note that for the `migdal_eliashberg` calculation, the number of $\mathbf{k}$-points (40 in this case) must be a multiple of the MPI processes (4 in this case).

Notice the following parameters in the input file:

- In the `&electrons` namelist, `theta=1.0` specifies the temperature of the calculation in Kelvin.
- In the `&input_migdal` namelist, `sigma_me` specifies the broadening (in eV) to be used in the calculation and `ME_Fermi_thickness` specifies the energy range around the Fermi level where electron eigenvalues are considered in the calculation. `nmatsu` indicates the number of Matsubara frequencies to be employed in the sum (see [here](https://the-epiq-team.gitlab.io/epiq-site/docs/calc/eliash/) for further details).
- The code solves the equation by generating a certain number of random $\mathbf{k}$-points, defined by the `64 64 1` grid in input and having an eigenvalue on the Fermi surface.

EPIq also calculates the Fermi surface using the grid specified after `nkfs` (64 64 1 here). After EPIq is done, two output files are produced:

- `MgB2.bxsf` contains the Fermi surface in the XCrysDen .bxsf format.
- `MgB2_ME.d` contains the Fermi-surface-resolved Eliashberg gap.

You can produce a plot of the Fermi-surface-resolved superconducting gap using the `plot_ME_fs.x` post-processing tool. Execute it as:

```console
{$path_to_epiq}/bin/plot_ME_fs.x
```

This is, for example, what you get for ${\rm MgB}_2$, with its famous double gap:

![Fermi surface of MgB2](figures_08/fs_mgb2.png)

Try to obtain the same kind of plot for ${\rm MoS}_2$ using [fermisurfer](https://mitsuaki1987.github.io/fermisurfer/) ([GitHub repository](https://github.com/mitsuaki1987/fermisurfer.git/)) by typing:

```console
fermisurfer MoS2.frmsf
```

> **Exercise:**
>
> Perform a convergence study of the average superconducting gap as a function of the number of Matsubara frequencies and k-points. What is the required value of `nmatsu` and k-points to have a precision better than 0.1 meV?

> **Note:**
> To perform a complete calculation of anharmonic electron-phonon coupling-related properties within SSCHA + EPIq, the following steps are required:
>
> 1. The [SSCHA code](http://sscha.eu/) $\rightarrow$ First, compute the free energy Hessian within the stochastic self-consistent harmonic approximation (SSCHA).
> 2. The [Quantum ESPRESSO code](https://www.quantum-espresso.org/) $\rightarrow$ Then, compute electron-phonon coupling matrix elements following the instructions reported on the [EPIq site](https://the-epiq-team.gitlab.io/epiq-site/docs/tutorials/step1/).
> 3. [Wannier90 code](http://www.wannier.org/) $\rightarrow$ Identify the unitary transformation connecting the Bloch and the maximally-smooth gauge (required to interpolate the electron-phonon coupling matrix elements in the Wannier basis).
> 4. [EPIq code](https://the-epiq-team.gitlab.io/epiq-site/) $\rightarrow$ Perform the electron-phonon coupling interpolation in the Wannier basis and calculate physical properties including the anharmonic correction from SSCHA.
>
> All these open-source software can be downloaded following the instructions on each website.
>
> If you want to know more about the full procedure, please take a look at the tutorials on the [EPIq site](https://the-epiq-team.gitlab.io/epiq-site/docs/tutorials/tutorials/).
