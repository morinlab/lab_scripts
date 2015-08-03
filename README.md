Morin Lab Scripts
=================

This repository serves to consolidate scripts for shared use in the Morin lab. Each script is housed in a separate subdirectory, which contains all related files, if applicable. 

Tags
----

Tags are used to track specific releases/versions of this repository. This allows for reproducibility when using these scripts for projects. This repository adopts the format `X.Y` (_e.g._ `1.0`), where X is the major version and Y is the minor version. When updates are made to existing scripts, increment the minor version number (_e.g._ `1.0 -> 1.1`). When new scripts are added to the repository, increment the major version (_e.g._ `1.0 -> 2.0`).  

Here's how to tag commits and push them to GitHub. 

```bash
git tag 1.1
git push --tags
```

And then on another machine, you can checkout that tag as follows. 

```bash
git clone git@github.com:morinlab/lab_scripts.git
cd lab_scripts
git checkout -b tags/1.1 tags/1.1
```
