# EBAME7 tutorial

Welcome to the EBAME7 hands-on for exploring the functional uncharacterised fraction of genomes and metagenomes. Go to the [**wiki**](https://github.com/genomewalker/ebame7/wiki/Unknowns) to get to the anvi'o part of the tutorial

___

Code for the AGNOSTOS tutorial

To recreate the analyses follow the instructions below:

  > This has been tested on R 4.1.1

Clone the repo with:

  ```bash
git clone https://github.com/genomewalker/ebame7.git
cd ebame7
```

Then let's install the packages we will need for the tutorial. First start R to get renv installed:

```
R
```

> If you open the project file `ebame7.Rproj` in Rstudio it will perform the same steps.

If everything went well, [renv](https://rstudio.github.io/renv/articles/renv.html) will be installed and you will get a message like:

```
* Installing renv 0.14.0 ... Done!
Successfully installed and loaded renv 0.14.0.
* ~/repos/ebame7. [renv 0.14.0]
```

And restore the environment:

```r
renv::restore()
q()
```
