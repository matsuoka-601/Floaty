# Floaty
Soft body & fluid simulation based on Position Based Dynamics (PBD). 

![demo](img/demo.gif)

[Try Demo Here!](https://floaty-fluid.netlify.app)

- **[Position Based Dynamics (PBD)](https://matthias-research.github.io/pages/publications/posBasedDyn.pdf) and [Position Based Fluids (PBF)](https://mmacklin.com/pbf_sig_preprint.pdf)** is used for simulating soft bodies and fluid.
- **For two-way coupling of soft body and fluid simulation**, a method described in [Unified Particle Physics for Real-Time Applications](https://mmacklin.com/uppfrta_preprint.pdf) is used. 
    - **Buoyancy is naturally simulated** thanks to this paper.
- The simulation is **accelerated with multithreading using [wasm-bindgen-rayon](https://github.com/RReverser/wasm-bindgen-rayon)**.

## How to run
```
npm install
npm run build
npm run dev
```
If you have trouble running the repo, feel free to open an issue.
## Demo Video
[![alt設定](http://img.youtube.com/vi/lzzW41bLtIo/0.jpg)](https://www.youtube.com/watch?v=lzzW41bLtIo)
## Overview of The Simulation
### Soft Bodies
Soft bodies are simulated using [distance constraint](https://github.com/matsuoka-601/PBDRust/blob/f6eacb716e3d9488d0906437ac7badc4df9f14a1/solver.ts#L9) and [area constraint](https://github.com/matsuoka-601/PBDRust/blob/f6eacb716e3d9488d0906437ac7badc4df9f14a1/solver.ts#L46). 

Distance constraints limit adjacent particles of each soft body from moving too far apart or too close together. And area constraints limit the movement of particles so that soft bodies do not collapse.

For more detail, reading section 3.3 (Constraint Projection) and 4.4 (Cloth Balloons) in the original PBD paper would help.
### Fluid
The simulation of fluid is based on [Position Based Fluids](https://mmacklin.com/pbf_sig_preprint.pdf). That is, $i$ th fluid particle is constrained to satisfy the following density constraint（ $\rho_0$ is the target density and $\rho_i$ is the density of $i$ th particle）.

$$
    C_i(\boldsymbol{p}_1,\dots, \boldsymbol{p}_n) = \frac{\rho_i}{\rho_0}-1=0
$$

My initial implementation is based on [the one by Matthias Müller](https://github.com/matthias-research/pages/blob/master/challenges/fluid2d.html). 

To simulate surface tension, [cohesion and curvature term is used](https://github.com/matsuoka-601/PBDRust/blob/ee13d82bf866fcea8be5faf9e5fdbc9db9085f41/src/lib.rs#L303), which is described in [Versatile Surface Tension and Adhesion for SPH Fluids](https://cg.informatik.uni-freiburg.de/publications/2013_SIGGRAPHASIA_surfaceTensionAdhesion.pdf) by Akinci et al.
### Two-Way Coupling of Fluid and Soft Bodies
Following the framework described in [Unified Particle Physics for Real-Time Applications](https://mmacklin.com/uppfrta_preprint.pdf), two-way coupling of fluid and soft bodies can be naturally simulated.

In PBF calculations, fluid particles and soft body particles are treated in almost the same way. Buoyancy is automatically achieved using the mass weighted version of PBD (Eq. (4) and (5) in the paper).
## Some Tricks to Stabilize the Simulation
### Prevent Soft Bodies from Get Entangling
The softbodies can easily get entangled when they collide. To prevent this problem, [repulsive forces are applied between overlapping soft bodies](https://github.com/matsuoka-601/PBDRust/blob/13f187b69d1259da6d52b0a4078263530c859ddd/PBD.ts#L369). They still sometimes get entangled even with this trick, but it's less frequent than it would be without it.

In addition, [soft body particles within other soft bodies do not interact with those soft bodies](https://github.com/matsuoka-601/PBDRust/blob/main/PBD.ts#L294). This also prevents soft bodies from becoming entangled with each other.
### Prevent Fluid Particles from Entering Soft Bodies
To prevent fluid particles from remaining inside soft bodies, [fluid particles inside soft bodies do not interact with soft body particles](https://github.com/matsuoka-601/PBDRust/blob/13f187b69d1259da6d52b0a4078263530c859ddd/src/lib.rs#L209). Fluid particles naturally go out of the soft bodies thanks to this trick.
## References
- [Position Based Dynamics](https://matthias-research.github.io/pages/publications/posBasedDyn.pdf) by Müller et al.
- [Position Based Fluids](https://mmacklin.com/pbf_sig_preprint.pdf) by Macklin et al.
- [Unified Particle Physics for Real-Time Applications](https://mmacklin.com/uppfrta_preprint.pdf) by Macklin et al.
- [Small Steps in Physics Simulation](https://mmacklin.com/smallsteps.pdf) by Macklin et al.
