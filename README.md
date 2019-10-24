# LArContent
[![Build Status](https://travis-ci.org/PandoraPFA/LArContent.svg?branch=master)](https://travis-ci.org/PandoraPFA/LArContent)
[![Coverity Scan Build Status](https://scan.coverity.com/projects/13057/badge.svg)](https://scan.coverity.com/projects/pandorapfa-larcontent)
[![codecov](https://codecov.io/gh/PandoraPFA/LArContent/branch/master/graph/badge.svg)](https://codecov.io/gh/PandoraPFA/LArContent)

Pandora algorithms and tools for LAr TPC event reconstruction

LArContent is distributed under the [GPLv3 License](http://www.gnu.org/licenses/gpl-3.0.en.html)

[![License](https://www.gnu.org/graphics/gplv3-127x51.png)](https://www.gnu.org/licenses/gpl-3.0.en.html)

## License and Copyright
Copyright (C), LArContent Authors

LArContent is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

# Non-Standard Features
## ProtoDUNEDataAnalysis Algorithm

The ProtoDUNE data analysis algorithm is designed to extract parameters from the reconstruction that may be of use when evaluating ProtoDUNE
data.  This is needed as the current event validation algorithms rely on the presence of MCParticles in the event.  This algorithm is
designed to run after the reconstruction is complete (parallel with the current event validation algorithms) and takes in the PFO and MC
particle list names as input.  The output is a root file containing the varaibles of interest.

This algorithm can run on either pure data (no mc), modified data (has a single MC particle representing the beam trigger information) or
pure simulation (complete MC).  The algorithm can run for neutrino experiments also, which may be of use to examine the properties of
reconstructed cosmic rays, however, the beam related varaibles in the output will be null.  A visual display of the triggered beam
information is available for debug purposes.


