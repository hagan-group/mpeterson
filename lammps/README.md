A collection of extensions, tools, and scripts for working with
[LAMMPS](https://github.com/lammps/lammps).

# Installing LAMMPS

LAMMPS can be configured and installed using the `install-lammps.py` script. This will clone the
LAMMPS repository, install the extensions in this project's `extensions` directory, configure CMake,
and build and install LAMMPS. This script requires `python` (version `3.5+`) to run. See the
script's `--help` flag for options. Note that you might need to build LAMMPS manually for more
complex installations.


## Available Extensions

The following extensions are included when installing LAMMPS using this script:

- `compute active/atom`: Computes tangential active forces on particles. The direction of the
  active force is based on the bonds that the particle is a part of; that is, the active force
  points from particles with a lower ID to those with higher IDs by default. This is designed
  to work well for polymer-like molecules that are linear chains of particles. Note that this
  requires `fix addforce` to actually apply the computed active force to each particle.
- `fix bond/flip`: Dynamically flips bonds and rearranged dihedrals based on their current energy.
  This is designed to allow for fluidized membranes.
