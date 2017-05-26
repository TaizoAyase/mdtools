mdtools
==========

# About
scripts for analysis of MD trajectory

using python, ruby, Gromacs-built-in-tools(>=5.0.5) etc...

## Required python package

* [mdtraj](http://mdtraj.org/) (>= 1.8.0)
* numpy (>= 1.10.4)
* scipy (>= 0.15.1)
* seaborn (>= 0.7.0)
* matplotlib (>= 1.5.0)
* tqdm (>= 4.10.0)

Some script is only compatible with MDAnalysis module.
For using scripts, **python 2.x is required**
Please see the import statement in the scripts.

* [MDAnalysis](http://www.mdanalysis.org/) (>=0.13.0)

Need to set your PYTHONPATH:

```
export PYTHONPATH=/path/to/mdtools/directory:$PYTHONPATH
```
(for bash)

## Required ruby package

* bioruby(>= 1.5.0)

## License

As the license of MDAnalysis module indicates, all source code in this repository is available under the GNU General Public License, version 2 (or an open source licence compatible with GPLv2).
