# Bites

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://h-mip.github.io/Bites.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://h-mip.github.io/Bites.jl/dev)
[![Build Status](https://github.com/h-mip/Bites.jl/workflows/CI/badge.svg)](https://github.com/h-mip/Bites.jl/actions)
[![Coverage](https://codecov.io/gh/h-mip/Bites.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/h-mip/Bites.jl)

Bites is a Julia package for analyzing mosquito-host biting networks. In addition to bipartite mosquito-host models, Bites now supports quadripartite simulations with birds, mosquitoes, humans, and horses, suitable for West Nile Virus (WNV) transmission modeling. The package includes features for simulating incubation periods and customizing host-specitic bite probabilities.


## Installation
```julia
(@v1.6) pkg> add https://github.com/h-mip/Bites.jl
```

## About

Bites.jl was developed as part of the [Host-Mosquito Interaction Project (H-MIP)](https://h-mip.github.io/), funded by the European Research Council (ERC) under the European Union’s Horizon 2020 research and innovation programme (Grant agreement No. 853271). The quadripartite network components used for West Nile Virus modeling were added as part of the [E4Warning project](https://www.e4warning.eu) (Eco-Epidemiological Intelligence for early Warning and response to mosquito-borne disease risk in Endemic and Emergence settings) funded by the European Union’s Horizon Europe programme under Grant Agreement 101086640.

Copyright 2021-2025 John R.B. Palmer

Bites.jl is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Bites.jl is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see http://www.gnu.org/licenses.