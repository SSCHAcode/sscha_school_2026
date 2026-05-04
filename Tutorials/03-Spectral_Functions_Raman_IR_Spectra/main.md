# Hands-on Session 5 - Raman and Infrared Spectra with the Time-Dependent Self-Consistent Harmonic Approximation

In the previous tutorial, you learned how to compute the spectral function by integrating the bubble in Fourier space, using the dynamical ansatz formulated by [Bianco et al., Physical Review B, 96, 014111, 2017](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.96.014111). Here, we will instead employ the Lanczos algorithm within the Time-Dependent Self-Consistent Harmonic Approximation (TD-SCHA) [Monacelli, Mauri, Physical Review B 103, 104305, 2021](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.103.104305).

For this, we need the package `tdscha` (it is recommended to configure it with the Julia speedup to run faster; see the installation guide).

## Computing the IR Signal in Ice

We use an already computed ensemble of phase XI of ice (low-temperature ice at ambient pressure, a prototype of standard cubic ice) to obtain the IR spectrum.

Inside the `data` directory, we find an already calculated ensemble of ice XI at 0 K with the corresponding original dynamical matrix `start_dyn_ice1` used to generate the ensemble and the dynamical matrix `final_dyn_ice1` after the SSCHA minimization.

### An Introduction

The infrared spectrum is related to the dipole-dipole response function:

$$
\chi_{MM}(\omega) = \int_{-\infty}^{\infty}dt e^{-i\omega t}\left<M(t) M(0)\right>
$$

where the average $\left<M(t)M(0)\right>$ is the quantum average at finite temperature.

Exploiting the TD-SCHA formalism introduced in the previous lecture, this response function can be written as:

$$
\chi_{MM}(\omega) = \boldsymbol{r}(M) \boldsymbol{G}(\omega) \boldsymbol{q}(M)
$$

where $\boldsymbol{G}(\omega)$ is the TD-SCHA Green's function, while $\boldsymbol{r}$ and $\boldsymbol{q}$ are vectors that quantify the perturbation and response, respectively.

In particular, if we neglect two-phonon effects (nonlinear coupling with light), we obtain:

$$
\chi_{MM}(\omega) = \sum_{ab}\frac{\mathcal Z_{\alpha a} {\mathcal Z}_{\alpha b}}{\sqrt{m_am_b}} G_{ab}(\omega)
$$

where ${\mathcal Z}_{\alpha a}$ is the Born effective charge of atom $a$ with polarization $\alpha$, and $G_{ab}(\omega)$ is the one-phonon Green's function (its imaginary part is precisely the spectral function).

Indeed, we need to compute the effective charges. This can be done directly by Quantum ESPRESSO using linear response theory (ph.x).

> **Exercise:**
>
> Use your knowledge of CellConstructor to extract a structure file from the final dynamical matrix to submit the calculation of the dielectric tensor, effective charges, and Raman tensor in Quantum ESPRESSO.
>
> Hint: The structure is the attribute `structure` of the Phonons object. The structure in the SCF file can be saved with the `save_scf` method of the Structure object.
>
> You can then attach the structure to the header of the ESPRESSO `ir_raman_header.pwi`.
>
> Notice that we are using norm-conserving pseudopotentials and the LDA exchange-correlation functional, as the Raman tensor in Quantum ESPRESSO is implemented only with them. However, this is usually an excellent approximation.
>
> You must run the pw.x code and the ph.x code (`ir_raman_complete.phi`), which performs the phonon calculation.
>
> We provide the final output file in `ir_raman_complete.pho`.

### Prepare the Infrared Response

We need to attach the Raman tensor and effective charges computed inside `ir_raman_complete.pho` to the final dynamical matrix. We will use this to initialize the response function calculation, as in the equation above.

To attach the content of an ESPRESSO ph calculation (only dielectric tensor, Raman tensor, and Born effective charges) to a specific dynamical matrix, use:

```python
dyn.ReadInfoFromESPRESSO("ir_raman_complete.pho")
```

If you save the dynamical matrix in Quantum ESPRESSO format, before the frequencies and the diagonalization, there will be the dielectric tensor:

```
 Dielectric Tensor:

      1.890128098000           0.000000000000           0.000000000000
      0.000000000000           1.912811137000           0.000000000000
      0.000000000000           0.000000000000           1.916728724000
```

Followed by the effective charges and the Raman tensor.

### Submitting the IR Calculation

With the following script, we submit a TD-SCHA calculation for the IR.

```python
import numpy as np
import cellconstructor as CC, cellconstructor.Phonons
import sscha, sscha.Ensemble
import tdscha, tdscha.DynamicalLanczos as DL

# Load the starting dynamical matrix
dyn_start = CC.Phonons.Phonons("start_dyn_ice")

# Load the ensemble
temperature = 0 # K
population = 2
n_configs = 10000

ensemble = sscha.Ensemble.Ensemble(dyn_start, temperature)
ensemble.load("data", population, n_configs)

# Load the final dynamical matrix
final_dyn = CC.Phonons.Phonons("final_dyn_ice")
final_dyn.ReadInfoFromESPRESSO("ir_raman_complete.pho")

# Update the ensemble weights for the final dynamical matrix
ensemble.update_weights(final_dyn, temperature)

# Setup the TD-SCHA calculation with the Lanczos algorithm
lanczos = DL.Lanczos(ensemble)
lanczos.ignore_v3 = True
lanczos.ignore_v4 = True

# If you have julia-enabled tdscha installed uncomment
# lanczos.mode = DL.MODE_FAST_JULIA
# for a x10-x15 speedup.

lanczos.init()


# Setup the IR response
polarization = np.array([1,0,0])  # Polarization of light
lanczos.prepare_ir(pol_vec = polarization)


# Run the algorithm
n_iterations = 1000
lanczos.run_FT(n_iterations)
lanczos.save_status("ir_xpol")
```

**Congratulations!** You ran your first TD-SCHA calculation. You can plot the results by using:

```console
tdscha-plot.py ir_xpol.npz
```

The script `tdscha-plot.py` is automatically installed with the tdscha package.

![IR spectrum with both *include_v3* and *include_v4* set to False.](05_raman_ir/IR_v2.png)

Additionally, `tdscha-plot.py` takes three more parameters: the range of frequencies to display and the smearing.

### Deep Dive into the Calculation

Let us dive a bit into the calculation. The beginning of the script should be almost self-explanatory, as we are just loading dynamical matrices, dielectric tensors, and effective charges.

The line

```python
ensemble.update_weights(final_dyn, temperature)
```

deserves special attention. Here, we are changing the weights of the configurations inside the ensemble to simulate the specified dynamical matrix and temperature, even if they differ from those used to generate the ensemble. This is useful for computing the spectrum at several temperatures without extracting and calculating a new ensemble each time.

```python
# Setup the TD-SCHA calculation with the Lanczos algorithm
lanczos = DL.Lanczos(ensemble)
lanczos.ignore_v3 = True
lanczos.ignore_v4 = True
lanczos.init()
```

Then we initialize the Lanczos algorithm for the TD-SCHA, passing the ensemble.

The `ignore_v3` and `ignore_v4` are flags that, if set to True, make the 3-phonon and 4-phonon scattering be ignored during the calculation.

As you can see from the output, our IR signal had very sharp peaks because we ignored any phonon-phonon scattering process that could give rise to a finite lifetime.

By setting only `ignore_v4` to True, we reproduce the behavior of the bubble approximation. Notably, while the four-phonon scattering is exceptionally computationally and memory demanding in free energy Hessian calculations, within the Lanczos algorithm, accounting for the four-phonon scattering is only a factor of two more expensive than using just the third order, without requiring any additional memory.

```python
# Setup the IR response
polarization = np.array([1,0,0])  # Polarization of light
lanczos.prepare_ir(pol_vec = polarization)
```

Here we tell the Lanczos which kind of calculation we want to do. In other words, we set the $\boldsymbol{r}$ and $\boldsymbol{q}$ vectors for the Lanczos calculation.

- `prepare_ir`
- `prepare_raman`
- `prepare_mode`
- `prepare_perturbation`

The names are intuitive; besides the Raman and IR, `prepare_mode` allows you to study the response function of a specific phonon mode, and `prepare_perturbation` enables defining a custom perturbation function.

```python
# Run the algorithm
n_iterations = 1000
lanczos.run_FT(n_iterations)
lanczos.save_status("ir_xpol")
```

Here we start the calculation of the response function. The number of iterations indicates how many Lanczos steps are required. Each step adds a new pole to the Green's function. Therefore, many steps are necessary to converge broad spectrum features, while much fewer are needed if the peaks are sharp. We save the status so that we can restore it later.

Finally, the commented line:

```python
lanczos.mode = DL.MODE_FAST_JULIA
```

This line only works if Julia and PyCall are correctly set up on your machine; in that case, run the script with `python-jl` instead of python. It will provide a massive speedup of a factor between 10x and 15x. The calculation can also be run in parallel using `mpirun` before calling the Python executable (or `python-jl`). In this case, to work correctly, you should have mpi4py installed and working.

> **Exercise:**
>
> Compute the Lanczos with the bubble approximation and without any approximation, and check the differences.

![IR signal accounting for the three-phonon scattering](05_raman_ir/IR_v3.png)

![IR signal accounting for all anharmonic scattering. The peaks appearing slightly below 2500 cm-1 is a combination mode known to be present in ice. See [Cherubini et al., J Chem Phys 155, 184502, 2021](https://pubs.aip.org/aip/jcp/article-abstract/155/18/184502/199619/The-microscopic-origin-of-the-anomalous-isotopic?redirectedFrom=fulltext)](05_raman_ir/IR_v4.png)

> **Exercise:**
>
> Try to see how different polarizations of the light affect the result.

## Analyze the Output

In the last part, we used the script `tdscha-plot.py` to display the simulation result. This is a quick way to show the results of a calculation.

Here, we will dive deeper into the calculation output file to extract the response function and obtain the results.

The Lanczos algorithm provides a set of coefficients $a_i$, $b_i$, and $c_i$ through which the Green's function is evaluated using a continued fraction:

$$
G(\omega) = \frac{1}{a_1 - (\omega + i\eta)^2 +  \frac{b_1c_1}{a_2 - (\omega+i\eta)^2 + \frac{c_2b_2}{a_3 - \cdots}}}
$$

Each iteration of the algorithm adds a new set of coefficients written to the standard output. Thanks to this expression, we only need the series of coefficients to compute the dynamical Green's function at any frequency and with any smearing. The Green's function can be computed with:

```python
green_function = lanczos.get_green_function_continued_fraction(frequencies, smearing=smearing)
```

Here, `frequencies` is an array in Rydberg. The response function is the negative of the imaginary part of the Green's function; thus, to reproduce the plot, we have:

```python
import tdscha, tdscha.DynamicalLanczos
import cellconstructor as CC, cellconstructor.Units
import numpy as np
import matplotlib.pyplot as plt

# Load the result of the previous calculation
lanczos = tdscha.DynamicalLanczos.Lanczos()
lanczos.load_status("ir_xpol_v4")

# Get the green function
W_START = 0
W_END = 3700
N_W = 10000
SMEARING = 10

frequencies = np.linspace(W_START, W_END, N_W)

# Convert in RY
frequencies_ry = frequencies / CC.Units.RY_TO_CM
smearing_ry = SMEARING / CC.Units.RY_TO_CM

# Compute the green function
green_function = lanczos.get_green_function_continued_fraction(frequencies_ry,
        smearing=smearing_ry)

# Get the response function
ir_response_function = - np.imag(green_function)

# Plot the data
plt.plot(frequencies, ir_response_function)
plt.show()
```

The previous script plots the data, precisely like `tdscha-plot.py`; however, now you have full access to the response function, both its imaginary and real parts.

> **Exercise:**
>
> Plot the IR data at various smearings and as a function of the number of steps (50, 100, 200, 300, and 1000). How does the signal change with smearing and the number of steps? When is it converged?

## Raman Response

The Raman response is very similar to the IR. Raman probes the fluctuations of the polarizability instead of those of the polarization, and it occurs when the sample interacts with two light sources: the incoming electromagnetic radiation and the outgoing one. The outgoing radiation has a frequency that is shifted with respect to the incoming one by the energy of the scattering phonons. The signal on the red side of the pump is called Stokes, while the signal on the blue side is the Anti-Stokes. Since the outgoing radiation has higher energy than the incoming one in the Anti-Stokes, it is generated only by existing (thermally excited) phonons inside the sample, and therefore it has a lower intensity than the Stokes.

On the Stokes side, the intensity of the scattered light with a frequency redshift of $\omega$ is:

$$
I(\omega) \propto \left<\alpha_{xy}(\omega)\alpha_{xy}(0)\right> (n(\omega) + 1)
$$

where $\alpha$ is the polarizability along the $xy$ axis. We can do a linear expansion around the equilibrium position of the polarizability, and obtain:

$$
\alpha_{xy}(\omega) = \sum_{a = 1}^{3N}\frac{\partial \alpha_{xy}}{\partial R_a(\omega)} (R_a(\omega) - \mathcal R_a)
$$

$$
\alpha_{xy}(\omega) = \sum_{a = 1}^{3N}\Xi_{xya} (R_a(\omega) - \mathcal R_a)
$$

If we insert this into the expression for the intensity, the average between the positions is the atomic Green's function divided by the square root of the masses, giving:

$$
I(\omega) \propto \sum_{ab} \frac{\Xi_{xy a} \Xi_{xy b}}{\sqrt{m_a m_b}} G_{ab}(\omega)(n(\omega) + 1)
$$

where $G_{ab}(\omega)$ is the atomic Green's function on atoms $a$ and $b$, while $\Xi_{xy a}$ is the Raman tensor along the electric fields directed in $x$ and $y$ on atom $a$.

The multiplication factor $n(\omega) + 1$ comes from the observation of the Stokes nonresonant Raman (it would be just $n(\omega)$ for the Anti-Stokes).

As we did for the IR signal, we can prepare the calculation of the Raman scattering by computing the polarizability-polarizability response function.

```python
# Setup the polarized Raman response
polarization_in = np.array([1,0,0])
polarization_out = np.array([1,0,0])
lanczos.prepare_raman(pol_vec_in=polarization_in,
        pol_vec_out=polarization_out)
```

Note that here we have to specify two polarizations of the light: the incoming radiation and the outgoing radiation.

> **Exercise:**
>
> Compute and plot the intensity of the Raman in the Stokes and Anti-Stokes configurations. Try with different polarizations and even orthogonal polarizations; what changes?
>
> The Bose-Einstein factor $n(\omega)$ can be computed with the following function:

```python
# n(w) Bose-Einstein occupation number:
# w is in Ry, T is in K
n_w = tdscha.DynamicalLanczos.bose_occupation(w, T)
```

## Unpolarized Raman and IR

In the previous section, we saw how to compute Raman and IR with specific polarizations of the incoming and outgoing radiation on oriented crystals (single crystals). However, the most common situation is a powder sample probed with unpolarized light.

In this case, we need to look at the Raman and IR response for unpolarized samples. While this is just the average of the x, y, and z directions for the IR signal, the Raman is more complex. In particular, the unpolarized Raman signal can be computed from the so-called *invariants*, where the perturbations in the polarizations are:

$$
I_A = \frac{1}{3}( {xx} + {yy} + {zz})^2/9
$$

$$
I_{B_1} = ({xx} - {yy})^2 / 2
$$

$$
I_{B_2} = ({xx} - {zz})^2/2
$$

$$
I_{B_3} = ({yy} - {zz})^2/2
$$

$$
I_{B_4} = 3({xy})^2
$$

$$
I_{B_5} = 3({yz})^2
$$

$$
I_{B_6} = 3({xz})^2
$$

The total intensity of unpolarized Raman is:

$$
I_{unpol}(\omega) = 45 \cdot I_a(\omega) + 7 \cdot \sum_{i=1}^6 I_{B_i}(\omega)
$$

The tdscha code implements a way to compute each perturbation separately. For example, the Raman response related to $I_A$ is calculated with:

```python
lanczos.prepare_raman(unpolarized=0)
```

While the $I_{B_i}$ components are computed using index $i$. For example, to compute $I_{B_5}$:

```python
# To compute I_B5 we do
lanczos.prepare_raman(unpolarized=5)
```

To obtain the total spectrum, you need to add the scattering factor $n(\omega) + 1$ and sum all these perturbations with the correct prefactor (45 for $I_A$ and 7 for the sum of all $I_B$).

To reset a calculation and start a new one, you can use:

```python
lanczos.reset()
```

which may be called before preparing the perturbation.

> **Exercise:**
>
> Compute the unpolarized Raman spectrum of ice and plot the results.

![Unpolarized Raman of ice, Stokes.](05_raman_ir/raman_unpolarized.png)

You should use a supercell size sufficiently large to properly converge the simulation. In this case, the 1x1x1 supercell is too small to converge the calculation and obtain meaningful results.
