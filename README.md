# Replication details for "Inducing Spectral Gaps for the Cohomological Laplacians of $\text{Sp}_{2n}(\mathbb{Z})$"

This is the replication code for providing lower bounds for the spectral gaps of the cohomological Laplacians $\Delta_1$ related to the symplectic groups $\text{Sp}\_\{2n\}(\mathbb{Z})$. More precisely, we provide lower bounds for $\lambda>0$ such that $\Delta_1-\lambda I$ is a sum of squares for $Sp\_\{2n\}(\mathbb{Z})$ for $n\in\{2,3\}$ and a specific quotient of $\text{Sp}\_\{2n\}(\mathbb{Z})$ for $n\geq 2$. The latter is the quotient by the normal subgroup generated by the commutators $[Z_i,Z_i^T]$, where $Z_i=I_{2n}+E_{i,i+n}$ with $E_{i,j}$ denoting the $2n\times 2n$ matrix filled with zeros except for the entry $(i,j)$. The detailed description of the presentation (i.e. the *Steinberg presentation*) taken into account in order to define $\Delta_1$ can be found in our article in section 2.

For the computations we used julia in version `1.11.2`.

## Getting Julia 1.11.2

To get the correct Julia version, first install [Juliaup](https://github.com/JuliaLang/juliaup), a cross-platform installer for Julia.

After installing Juliaup, install Julia in version `1.11.2` by running

```bash
juliaup add 1.11.2
```

## Obtaining code
To obtain the code for the replication, you can either download it directly from [Zenodo](https://zenodo.org/records/15225937), or use git for this. In the latter case, first clone this repository via
```bash
git clone https://github.com/jtszy/SP_2N_Cohomology.git
```

## Setting up the environment
First, run julia in the terminal in `SP_2N_Cohomology` folder
```bash
julia --project=.
```
Next, to set up the proper environment for the replication run in Julia REPL
```julia
julia> using Pkg; Pkg.instantiate()
```
This command installs and precompiles, if needed, all the necessary dependences,
so it may take a while.
Note that this step needs to be executed only once per installation.

## Running actual replication

Our scripts perform the necessary optimizations to find such sums of squares decomposition. In all the execution commands below, there is a possibility to calculate the desired sum of squares decomposition from the already precomputed solution which shall substantially reduce the execution time (we skip in these cases solving the semi-definite optimization problem and focus on the *certification* procedure only, see Appendix B in our paper). In the case you choose using the precomputed solution, it is necessary to download them from [Zenodo](https://zenodo.org/records/15225937) since we have not uploaded them on GitHub due to its memory constraints. 

In order to compute everything from scratch, set the `precomputed` flag to `false`. In the case you wish to use the precomputed solutions, set `precomputed` to `true`.

### $\Delta_1-\lambda I$ is a sum of squares for $\text{Sp}_4(\mathbb{Z})$ and $\text{Sp}_6(\mathbb{Z})$
We wish to prove that for the Steinberg presentations of $\text{Sp}\_4(\mathbb{Z})$ and $\text{Sp}\_6(\mathbb{Z})$ on $8$ and $18$ generators respectively (as defined in section 2 of our paper)
$\Delta_1-\lambda I_{8}$ and $\Delta_1-\lambda I_{18}$ is a sum of squares for $\lambda=0.82$ and $\lambda=0.99$ respectively.

The following command needs to be executed in the terminal in `SP_2N_Cohomology` folder:
```bash
julia --project=. ./scripts/Sp_2n_Steinberg.jl n precomputed delta_1
```

The running time of the script will be approximately `25` minutes and `115` minutes on a standard laptop computer for the cases $n=2$ and $n=3$ respectively. Therefore, in the latter case, we encourage to use the precomputed solution, focusing on the part providing the rigorous proof only - the running time in such a case will be approximately `70` minutes on a standard laptop computer.

### $\Delta_1-\lambda I$ is a sum of squares for $\text{Sp}_{2n}(\mathbb{Z})/\langle [Z_i,Z_i^T]\rangle$
This is the part responsible for the proof for the lower bound of the spectral gap of $\Delta_1$ for $\text{Sp}\_\{2n\}(\mathbb{Z})/\langle [Z_i,Z_i^T]\rangle$. The key idea is to use the induction technique and reduce the proof to two concrete calculations for $\text{Sp}\_6(\mathbb{Z})/\langle [Z_i,Z_i^T]\rangle$. More precisely, one distinguishes two particular summands of $\Delta_1^-$, denoted by $\text{Sq}^-$ and $\text{Adj}^-$. The two calculations mentioned before are then providing lower bounds for $\lambda$ in expressions $\text{Sq}^-+\Delta_1^+-\lambda I$ and $\text{Adj}^-+\Delta_1^+-\lambda I$, the lower bounds being $0.99$ and $0.24$ respectively. The latter is much more important since it is responsible for the proof for all $n$ starting from a particular one, which can be set to $3$ by the first calculation. The details can be found in our article. 

In order to replicate the computations, run in the terminal in `SP_2N_Cohomology` folder: 
```bash
julia --project=. ./scripts/Sp_2n_Steinberg.jl 3 precomputed sq_adj
```
The `sq_adj` flag corresponds to $\text{Sq}^-+\Delta_1^+-\lambda I$ and $\text{Adj}^-+\Delta_1^+-\lambda I$ respectively. Two corresponding options for this flag are: `"sq"` and `"adj"`. If one launches the replication from precomputed solution, the running times shall be approximately `2` hours. In the case one wish to run the whole semi-definite optimization, it will be about `3` hours.
