# Hands-on Session 5 - Raman and Infrared Spectra with the Time-Dependent Self-Consistent Harmonic Approximation


## Theoretical background

This tutorial focuses on the computation of dynamical properties of materials.
The key quantity which is measured by experiments is the *dynamical response function* $\chi(\omega)$.
The response function probes how the material responds to a time-dependent external perturbation.
We can model the result of any experiment as: the material is in equilibrium for $t < t_0$. Then, at $t_0$ a new perturbation is turned on,
and we measure a property $A$ of the material at a certain time $t$. This response will be the convolution of the perturbation at all times between $t_0$ and $t$, weighted by how the perturbation propagates in time, which is precisely the response function $\chi(t-t')$:
$$
A(t) = A_0 + \int_{t_0}^t dt' \chi(t - t') F(t')
$$
where $F(t')$ is the external time-dependent perturbation, and $A_0$ is the value of the property $A$ in equilibrium.
Alternatively, by passing in Fourier space, the convolution becomes a simple product:
$$
A(\omega) = A_0 + \chi(\omega) F(\omega)
o$$

This is very general, and it applies to any experiment. In this tutorial, we focus on IR and Raman spectroscopy.
For the case of IR, the perturbation is the electric field of the light, which interacts with the material by generating an oscillating dipole. The response of the material is the *emitted* (or absorbed) electric field, which is related to how the dipole moment of the material oscillates in time.
So the observable $A(t)$ is the dipole moment of the material, and the perturbation $F(t)$ is the amplitude of the electric field of the light.

The Kubo equation for the response function is
$$
\chi(t) = \frac{i}{\hbar}\theta(t)\left<M(t)M(0)\right>
$$
where $\hat M(t)$ is the dipole moment operator in the Heisenberg picture, and $\left<\cdot\right>$ is the quantum average at finite temperature.
The SSCHA can compute the response function $\chi(t)$ by using the TD-SCHA formalism. Here we will use the Lanczos algorithm to compute the response function in Fourier space, which is more efficient than computing it in time and then Fourier transforming it.

The first step is to convert the diopole operator in phonons creation and annihilation operators. To this aim, we can approximate the dipole moment as a linear function of the atomic displacements $\hat u$:
$$
\hat M(t) = \left.\frac{\partial M}{\partial u}\right|_{u = 0} \hat u(t) + O(u^2)
$$
The derivative of the dipole moment with respect to the atomic displacements is precisely the Born effective charge:
$$
Z^i_{\alpha\beta} = \frac{\partial M_\alpha}{\partial u^i_\beta} = \frac{\partial^2 \mathcal E}{\partial u^i_\beta \partial E_\alpha}
$$
where $E_\alpha$ is the electric field in the $\alpha$ Cartesian direction, $u^i_\beta$ is the displacement of atom $i$ in the $\beta$ Cartesian direction, and $\mathcal E$ is the total energy of the system (Born-Oppenheimer).
Using this in the expression of the suscieptibility, we obtain:
$$
\chi_{\text{IR}\alpha}(t) = \sum_{ij}\sum_{\beta\gamma} \frac{Z^i_{\alpha\beta} Z^j_{\alpha\gamma}}_{\sqrt{m_i m_j}}\sqrt{m_im_j} \left<\hat u^i_\beta(t) \hat u^j_\gamma(0)\right>
$$
The last term is the so-called phonon Green's function, which represent the propagation of a phonon created at time $t=0$ and destroyed at time $t$.
$$
G^{ij}_{\alpha\beta}(t) = \sqrt{m_i m_j}\left<\hat u_\alpha^i(t) \hat u^j_\beta(0)\right>
$$

Analogously, the Raman spectrum is related to the response of the material to two electric fields, which interact with the material by generating an oscillating polarizability. The observable $A(t)$ is the polarizability of the material, and the perturbation $F(t)$ is the amplitude of the two overlapping electric fields which interferee within the material. The Raman susceptibility can be written as:
$$
\chi_\text{Raman}(t) = \left<\hat\alpha_{\alpha\beta}(t) \hat\alpha_{\alpha\beta}(0)\right>
$$
where $\hat \alpha_{\alpha\beta}$ is the polarizability operator, which can be approximated as a linear function of the atomic displacements:
$$
\hat \alpha_{\alpha\beta}(t) = \left.\frac{\partial \alpha_{\alpha\beta}}{\partial u^i_\gamma}\right|_{u = 0} \hat u^i_\gamma(t) + O(u^2)
$$
$$
\Xi^i_{\alpha\beta\gamma} = \frac{\partial \alpha_{\alpha\beta}}{\partial u^i_\gamma} = \frac{\partial^3 \mathcal E}{\partial u^i_\gamma \partial E_\alpha \partial E_\beta}
$$
where $\Xi^i_{\alpha\beta\gamma}$ is the Raman tensor, which is the third derivative of the total energy (Born-Oppenheimer) with respect to the atomic displacements and the two electric fields (incoming-outgoing).

Notably, the Raman requires two electric fields because it is a scattering, where incoming and outcoming radiation are different. The IR instead is an absorption/emission, where incoming and outcoming radiation are the same.

The Raman suscieptibility is then:
$$
\chi_\text{Raman}(t) =  \left<\hat\alpha_{\alpha\beta}(t) \hat\alpha_{\alpha\beta}(0)\right> = \sum_{ij}\sum_{\gamma\delta} \frac{\Xi^i_{\alpha\beta\gamma} \Xi^j_{\alpha\beta\delta}}{\sqrt{m_i m_j}}\sqrt{m_im_j} \left<\hat u^i_\gamma(t) \hat u^j_\delta(0)\right>
$$

Also here, the last term is the phonon Green's function, which represent the propagation of a phonon created at time $t=0$ and destroyed at time $t$.
Therefore, to compute the IR and Raman spectra, we need to compute the phonon Green's function.
Going in Fourier space, we get
$$
\chi_{\text{IR}\,\alpha}(\omega) = \sum_{ij}\sum_{\beta\gamma} \frac{Z^i_{\alpha\beta} Z^j_{\alpha\gamma}}{\sqrt{m_i m_j}}\sqrt{m_im_j} G^{ij}_{\beta\gamma}(\omega)
$$
$$
\chi_{\text{Raman}\,\alpha\beta}(\omega) = \sum_{ij}\sum_{\gamma\delta} \frac{\Xi^i_{\alpha\beta\gamma} \Xi^j_{\alpha\beta\delta}}{\sqrt{m_i m_j}}\sqrt{m_im_j} G^{ij}_{\gamma\delta}(\omega)
$$

### The phonon Green's function

The phonon Green's function can be computed within the TD-SCHA. In general, we can write the Green's function as a non-interacting green function plus a self-energy:
$$
G^{-1}(\omega) = G_0^{-1}(\omega) - \Pi(\omega)
$$
where $G_0(\omega)$ is a non-interacting Green's function for *harmonic*-like phonons, while $\Pi(\omega)$ is the phonon self-energy, which depends on the frequency.
In the TD-SCHA theory, the non-interacting Green's function is defined as the Green's function of the auxiliary SSCHA harmonic Hamiltonian.

$$
G_0^{-1}(\omega) = I\omega^2 - D_\text{SSCHA}
$$
where $I$ is the identity matrix and $D_\text{SSCHA}$ is the dynamical matrix of the SSCHA Hamiltonian.

The self-energy $\Pi(\omega)$ correspond to the same expression derived for the Free energy Hessian; see [Bianco et al., Physical Review B, 96, 014111, 2017](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.96.014111).
$$
\Pi_{ab}(\omega) = -\frac 12\sum_{cd\mu\nu} \overset{(3)}{\Phi}_{acd} \Lambda^{cd}_{\mu\nu}(\omega)\left[ 1 +\frac 12 \overset{(4)}{\Phi}\Lambda(\omega)\right]^{-1} \overset{(3)}{\Phi}_{\nu b}
$$
where $\Lambda^{cd}_{\mu\nu}(\omega)$ is the Fourier transform of the two-phonon propagator, which can be computed from the equilibrium non interacting Green's function $G_0(\omega)$ as:
$$
\Lambda^{ijlm}_{\alpha\beta\gamma\delta}(t) = \sqrt{m_i m_j m_l m_m}\left<\hat u^i_\alpha(t) \hat u^j_\beta(t) \hat u^l_\gamma(0) \hat u^m_\delta(0)\right>_0
$$

Computing the full inversion of the self-energy for every value of the frequency is computationally extremely expensive.
It is possible, with an efficient algorithm (Lanczos), to prove that the Green function can be obtained by inverting a special matrix, see [Monacelli, Mauri, Physical Review B 103, 104305, 2021](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.103.104305).

In this tutorial, we will use this Lanczos algorithm to compute the phonon Green's function and, consequently, the IR spectrum of CsPbI3.

For this, we need the package `tdscha` (it is recommended to configure it with the Julia speedup to run faster; see the installation guide).

## The Phonon Green's function

We need to first relax a complete SSCHA calculation, exactly as for the free energy hessian.
You find the script `sscha_relax.py`, which obtain the CsPbI3 structure and final dynamical matrix at 500 K.

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
