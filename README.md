For questions please email rsmeht4@emory.edu.

# DeCoTUR
## Detecting Coevolving Traits Using Relatives
This package computes the "coevolution score" presented in [this preprint](https://www.biorxiv.org/content/10.1101/2022.03.14.484367v1).

## Installation
This package is implemented in R. If you don't have `devtools`, then run
`install.packages('devtools')`.

To install DeCoTUR, run

`library(devtools)`

`install_github("rohansmehta/decotur")`

## Using DeCoTUR
The primary function you would want to use is

`get_score`

Everything else is a helper function. The function `get_score` has many arguments. The two most important are the presence-absence matrix and the distance matrix.

For the presence-absence matrix, each row is a trait, and each column is a sample. For example

| Trait | Sample1 | Sample2 | ... | SampleN |
| ----- | ------- | ------- | --- | -------- |
| Trait1 | 1 | 0 | ... | 1 |
| Trait2 | 1 | 1 | ... | 0 |
| Trait3 | 1 | 1 | ... | 1 |
| ... | ... | ... | ... | ... |
| TraitM | 0 | 0 | ... | 1 |

For the distance matrix, each row and column is a sample. For example

| | Sample1 | Sample2 | ... | SampleN |
| ----- | ------- | ------- | --- | -------- |
| Sample1 | 0 | 0.001 | ... | 1.01 |
| Sample2 | 2.03 | 0 | ... | 0.023 |
| Sample3 | 6.07 | 1.97 | ... | 1.45 |
| ... | ... | ... | ... | ... |
| SampleN | 0.001 | 0.0015 | ... | 0 |

The other inputs to this function are: A) a method for determining the close-pair cutoff. This can be 
1. Auto. I don't recommend this. It looks at the pairwise distance distribution and tries to find a cutoff that encompasses the first peak. Very complicated, not worth it, only in here as a curiosity and if someone wants to be "super rigorous" about determining the cutoff.
2. Distance. You give it a distance cutoff, it uses it.
3. Fraction. You give it a fraction (i.e. I want the 5% closest pairs, so my fraction is 0.05).
4. Fixed Number. I.e. I want 500 close pairs.

Each method requires particular parameters detailed in the function descriptions. Put those parameters in a list, and that's another input into the function. For instance, if you want to use the Fixed Number method, that requires 1) a number of close pairs (let's say 500), a random seed (let's say 1), and whether or not you want to see the histogram of the distance distribution (which I would recommend), so it's TRUE.

So for the closepair_params part of the argument, this example would need 
`closepair_params = list(500, 1, TRUE)`

B) A block size for computation. If you divide up the traits you want to evaluate into blocks, the computation goes much faster. Provided that you had a lot of traits to being with. If not, set this number to the number of traits.

C) Whether or not you want to downweight bushes. I downweight them and say it's a good idea. You might not want to, for whatever reason.

D) Whether or not you want p values with your results. I would recommend keeping this at the default of TRUE; unless it really slows down your analysis, it won't hurt.

E) Whether or not you want text output for progress. I would recommend keeping this at the default of TRUE.

F) Whether or not you want to run the relatively "speed" optimized version or the relatively "memory" optimized version. In my experience, there is no reason to choose "memory" over "speed" so long as you have a small enough "blocksize", but the option is there in case you run into memory issues.
