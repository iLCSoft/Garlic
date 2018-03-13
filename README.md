# Garlic
[![Build Status](https://travis-ci.org/iLCSoft/Garlic.svg?branch=master)](https://travis-ci.org/iLCSoft/Garlic)

Garlic is a Marlin Processor to identify photons and electrons.

Garlic is distributed under the [GPLv3 License](http://www.gnu.org/licenses/gpl-3.0.en.html)

[![License](https://www.gnu.org/graphics/gplv3-127x51.png)](https://www.gnu.org/licenses/gpl-3.0.en.html)


## Instructions
to build Garlic:
```
mkdir build
cd build
cmake -C /some_path/ILCSoft.cmake ..
make
make install
```

example sttering file is given in `example_garlic_steer.xml`

Suggested cluster properties file given in `data/garlicCuts_v1.txt`.
This data is obtained from the 30-layer silicon ecal in model `ILD_o1_v06`


## License and Copyright
Copyright (C), Garlic Authors

Garlic is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License long with this program.  If not, see <http://www.gnu.org/licenses/>.
