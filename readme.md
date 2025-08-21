In your home directory,
* `mkdir QCFD`
and download `CLBM_ZhixingSong.tar.gz` and untar it.
You may rename this folder to whatever you want but make sure it is consistent with `QCFD_HOME`
* `vi ~/.bashrc` and add the following,
```bash
export QCFD_HOME=$HOME/QCFD/CLBM_ZhixinSong/
export QCFD_SRC=$QCFD_HOME/src/
```
then
`cd QCFD_SRC/CLBM`

`Julia`

```julia
julia> using Pkg
julia> Pkg.activate("../..")  # Activate the main project environment from src/CLBM
julia> include("clbm_run.jl")       
```
