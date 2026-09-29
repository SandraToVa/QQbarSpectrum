# QQbarSpectrum

MATLAB code for the charmonium spectrum of quarkonium and hybrids
(spin 0 and spin 1), with optional quarkonium–hybrid mixing and hyperfine
splitting of the spin-1 hybrids. It also fits the hyperfine parameters
A and B to the lattice hybrid spectrum.

Authors: Sandra Tomàs (spectrum codes, fit), Ruben Oncala (original
subroutines and potentials).

## Quick start

```matlab
% Full spectrum, wave-function plots and chi^2 at a fixed (A, B)
QQbarSpectrum

% Best (A, B), chi^2/dof and errors
QQbarSpectrum_Optimize
```

A single family can also be called directly:

```matlab
[E, W, x, info] = GQQbarS1Jcal1(1.496, A, B);
% E    : energies of the 7 lowest states (GeV)
% W    : W(c,:,k) = radial wave function u_c(r) of channel c of state k
% x    : radial points of W (GeV^-1)
% info : mesh, eigenvectors, channel weights (info.weights), potential handle
```

## Files

| File | What it does |
|---|---|
| `QQbarS0J{0,1,2}.m`, `QQbarS1J{0,1,2}.m` | Quarkonium, spin 0 / spin 1, for each J |
| `GQQbarS0Jcal{0,1,2}.m` | Spin-0 hybrids, for each Jcal |
| `GQQbarS1Jcal{0,1,2,3}.m` | Spin-1 hybrids with hyperfine terms (depend on A, B) |
| `quarkoniumConfig.m` | Flags for the quarkonium files: `Eg`, `Eval`, `mix` |
| `hybridConfig.m` | Flags for the hybrid files: `Eg`, `Eval`, `mix`, `hf`, `glamb1`, `glamb3`, `sig`, `pm`, `r0` |
| `hybridAssignments.m` | Lattice data `t`, `e` and which computed state each data point is compared with (one table for mix=false, one for mix=true) |
| `QQbarSpectrum.m` | Computes everything and plots the wave functions |
| `QQbarSpectrumTable.m` | Spectrum tables as in mixing.pdf: n, S, J, j, L, E, J^PC, dominant component, %Q, %H of every state |
| `QQbarSpectrum_Optimize.m` | Fit of A and B |
| `spinAverages.m` | Spin averages and splittings of each multiplet vs the data at given (A, B), with a χ² breakdown |
| `proposeMixAssignment.m` | Suggests the mix=true assignment table from the mix=false one (wave-function overlaps), as a cross-check |
| `solveRadialSchrodinger.m` | Solver: energies and wave functions of a coupled-channel radial problem |
| `radialMesh.m`, `radialHamiltonian.m` | Mesh and Hamiltonian matrix used by the solver |

Each physics file contains its own potential matrix (`potentialMatrix`) and
potentials. The number of channels is read from the potential matrix, so it
follows the `mix` flag automatically.

## Configuration

All flags live in the two config files. You can edit the defaults there, or
set them at the top of a script:

```matlab
hybridConfig('reset');
hybridConfig('mix', true);          % hybrid files
quarkoniumConfig('mix', true);      % quarkonium files
cfg = hybridConfig;                 % read the current values
```

The settings persist in memory until `hybridConfig('reset')` or
`clear functions` / `clear all`. `QQbarSpectrum_Optimize` resets and sets them
explicitly each time it runs.

Conventions:

- Units are GeV and GeV^-1.
- The potential matrices already include the centrifugal terms and the factor m.
  The solver solves `-u'' + V u = m E u` and returns `E` in GeV.
- With `Eg = true` (and `m = 1.496`), `2m - Eval` is added to the potentials,
  so `E` is the mass of the state.
- The hyperfine terms are next-to-leading order in 1/m, which is why they
  carry no extra factor m. `Vsa` and `Vsb` are evaluated at `r0`, as the
  interpolation model requires.

## Numerical method

`solveRadialSchrodinger` uses a regularized Lagrange–Laguerre mesh
(D. Baye, *Phys. Rep.* 565 (2015) 1). The radial functions behave as r at the
origin and decay exponentially at large r. The whole coupled-channel problem
becomes one symmetric matrix, and a single `eig` gives every level and wave
function at once. This means:

- no shooting, bracketing or matching, and no overflow or NaNs;
- degenerate levels come out correctly;
- the energies depend smoothly on the parameters, which the fit needs.

Default mesh: `N = 100` points per channel, last point at `rmax = 25` GeV^-1
(the fit uses `N = 60`). The energies converge to about 1e-6 GeV, far below
the lattice errors. The fit repeats this check at the minimum with a larger
mesh. On problems with exact solutions (harmonic oscillator, Coulomb, linear)
the solver agrees to about 1e-13.

Options (`opts` struct, last argument of any physics file):

| Field | Meaning | Default |
|---|---|---|
| `N` | mesh points per channel | 100 |
| `rmax` | last mesh point (GeV^-1) | 25 |
| `x` | points where `W` is evaluated | `linspace(0,22,1000)` |

`GQQbar...(m, A, B, 'potential')` returns only the potential handle in
`info.V`, without solving. The fit uses this.

This replaces the earlier CP shooting solver (`computeEigenvalues_vS`,
`shootForEigenvalue_vS`, ...), which you can still find in the git history.

## The fit (`QQbarSpectrum_Optimize.m`)

1. **Model.** The spin-1 hybrid Hamiltonian is exactly linear in A and B,
   `H(A,B) = H0 + A*HA + B*HB`. H0, HA and HB are built once per Jcal; the
   script checks the linearity. Each (A, B) evaluation is then only a
   diagonalization. Spin-0 hybrid energies do not depend on A or B and are
   computed once.
2. **Grid scan** over `[lb, ub]` (default 25×25) for the starting point and
   the chi^2 map.
3. **Local fit** with `lsqnonlin`, using the exact Jacobian from
   Hellmann–Feynman: `dE/dA = <psi|HA|psi>/m`.
4. **Errors.** `cov = inv(J'J)` at the minimum (Δχ² = 1). Also printed:
   the errors scaled by `sqrt(chi^2/dof)` when chi^2/dof > 1 (PDG
   convention), the correlation, the p-value, and the pull of each data
   point. The plots show the fit against the data and the chi^2 map with the
   68% and 95% contours (Δχ² = 2.30, 6.18).

A warning is printed if the minimum lies on the boundary of `[lb, ub]`.

### State assignment

Each row of the table in `hybridAssignments.m` is `{spin, Jcal, index}`:
the data point `t(k)` is compared with the `index`-th state, in energy
order, of `GQQbarS<spin>Jcal<Jcal>.m`, **at the best fit**.

When levels of different decoupled channel blocks cross as A and B change,
pure energy order would swap states during the fit and put kinks in chi^2.
So the fit follows each level as "k-th state of its channel block". At the
minimum it checks the energy order again, and if it changed, it refits from
there until it is stable. The printed assignment table shows, for every data
point, the block (channels) and the energy at the best fit, so you can check
that each point is matched to the intended state.

With `mix = true` quarkonium states enter the spectra, so the indices
change. That case uses its own table. `proposeMixAssignment.m` can suggest
one from wave-function overlaps with the mix=false states, but the final
table should be checked against the wave functions.

## Known limitations

- `QQbarS1J{0,1,2}.m` implement only `mix = true`. With `mix = false` they
  stop with an error.
- `GQQbarS1Jcal*.m` with `mix = false, hf = false` is not implemented (it
  would be the spin-0 case). The solver stops with an error if a potential
  matrix is identically zero.
- `cF` and the `Eg` shift are only switched on when `m == 1.496` exactly.

## Status (Eg = true, hf = true)

- mix = false: A = −0.0447 ± 0.0140, B = 0.0014 ± 0.0078 (errors scaled by
  √(χ²/dof)), ρ = 0.57, χ²/dof = 33.0/12 = 2.75. The four spin-0 points, which
  do not depend on A or B, contribute 7.5 to the χ².
- mix = true: A = −0.0920 ± 0.0150, B = −0.0397 ± 0.0071, ρ = 0.01,
  χ²/dof = 43.8/12 = 3.65. The mix = true assignment table should be checked at
  this minimum.

The fitted energies get the error σ_E = √(σ_fit² + 0.03²) GeV, where σ_fit is
the (A, B) covariance propagated to each energy (checked on the 1σ ellipse).
