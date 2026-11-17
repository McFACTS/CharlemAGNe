# CharlemAGNe: C++ Hydrostatic Annulus Radially Localized Equations Module for Active Galactic NuclEi

This module solves analytic disk models including the
  [Sirko-Goodman](https://academic.oup.com/mnras/article/341/2/501/1105444) and
  [Thompson et al.](https://scixplorer.org/abs/2005ApJ...630..167T/abstract) 
  disks as a function of radius from 
  the supermassive black hole (SMBH)
  of an AGN disk.
There are two reasons for writing this in C++:
- The intention is that this module will be fast enough to generate a disk
  in a small fraction of a second. 
  Prior to this, 
  [McFACTS](https://github.com/mcfacts/mcfacts)
  took a comparable amount of time to generate a disk using existing Python
  implementations and simulate one galaxy for one Myr.
  A faster disk model will allow for Monte-Carlo sampling of disk properties
  (such as SMBH mass and eddington rate).
- I would like to use this module to set up the initial conditions for
  hydrodyanamics simulations, and those simulations are written in C++
  (E.g. [KHARMA](https://github.com/parthenon-hpc-lab/kharma/)).
  In order to set up those initial conditions with the same opacity model
  which is used in-situ in the simulation, CharlemAGNe itself
  must use the same opacity module.
  The Parthenon-based hydro codes I am interested in working with all use
  [singularity-opac](https://github.com/lanl/singularity-opac).
  Writing this module in C++ allows me to use singularity-opac, and to include
  the C++ library in a problem generator for the initial conditions
  of a KHARMA simulation.

## Authors and acknowledgment

Authors:
- Vera Delfavero
- Cora Prather

This module is a spiritual successor to pAGN:
 [github](https://gitlab.com/aatrani/pagn),
 [paper](https://academic.oup.com/mnras/article/530/4/3689/7659378?login=true),
 [arxiv](https://arxiv.org/abs/2403.00060).

## A note on the name

With codes like "pAGN" and "AGNostic" in the environment now, it seems
  like words containing "AGN" with religious tones are a motif the
  community is sticking with.
CharlemAGNe fits that motif nicely, and I find a great deal of joy in the pun.
However, I must include a disclaimer that we do not endorse the actions
  of the historic figure Charlemagne, nor do we endorse the actions
  of the Holy Roman Empire.

## Project Milestones

Reproduce SG03 and TQM05 in a faster C++ code compatible with singularity-opac opacity tables: :white_check_mark:

Python interface and tests: :white_check_mark:

Upstream support for $R = (\rho / (T / 10^6 K))$ coordinates in singularity-opac: 
 in-progress [(see pull request).](https://github.com/lanl/singularity-opac/pull/88) 

Support for alternative opacity tables: not started

Direct singularity-opac Python surface: not started

Cylindrical coordinates and 2D disk profiles: not started

## LLM Policy

Large Language Models (LLMs) may be used for learning, ideas,
 to find bugs, and to check work.
They may not be used to generate code or pull requests without supervision.
They may not be used to generate documentation.
Significant design choices should be made by humans.
LLMs may assist with getting implementation details right 
 (I have found them very helpful in getting meson-python 
 and CMake to play nicely with one another).

## Contributing

Please reach out to me at xevra86@gmail.com before opening a pull request.


## License

You may distribute the C++ module without the pAGN opacity table files under
  the [BSD-3-CLAUSE license](LICENSE-BSD-3-CLAUSE).

The Python wheel (for the time being) includes an opacity table
  which is derived from pAGN's opacity files.
The Python module is therefore distributed under the
  [GNU Public License v3](LICENSE-GPL-3).
