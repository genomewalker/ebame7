# EBAME 6 tutorial

Welcome to the EBAME4 hands-on for exploring the functional uncharacterised fraction of genomes and metagenomes. Go to the [**wiki**](https://github.com/genomewalker/ebame6/wiki) to get to the tutorial

___

Code for the AGNOSTOS tutorial

To recreate the analyses follow the instructions below:

  > This has been tested on R 4.1.1

Clone the repo with:

  ```bash
git clone https://github.com/genomewalker/ebame6.git
cd ebame6
```

Then let's install the packages we will need for the tutorial. First start R to get renv installed:

```
R
```

> If you open the project file `ebame6.Rproj` in Rstudio it will perform the same steps.

If everything went well, [renv](https://rstudio.github.io/renv/articles/renv.html) will be installed and you will get a message like:

```
* Installing renv 0.14.0 ... Done!
Successfully installed and loaded renv 0.14.0.
* ~/repos/ebame6. [renv 0.14.0]
```

And restore the environment:

```r
renv::restore()
q()
```
