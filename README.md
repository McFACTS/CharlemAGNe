# CharlemAGNe: C++ Hydrostatic Annulus Radially Localized Equations Module for Active Galactic NuclEi


## Authors and acknowledgment

Authors:
- Vera Delfavero
- Cora Prather

This module is a spiritual successor to pAGN:
 [github](https://gitlab.com/aatrani/pagn),
 [paper](https://academic.oup.com/mnras/article/530/4/3689/7659378?login=true),
 [arxiv](https://arxiv.org/abs/2403.00060).

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
