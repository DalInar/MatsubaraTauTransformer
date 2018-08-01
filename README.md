# MatsubaraTauTransformer
This code transforms greens functions and self-energies between matsubara frequency and tau space.

The meat of the code is the Python script `MatsubaraTauTransformer.py`.  This code will read in a file containing a Greens function or Self-Energy in either Matsubara frequency or Tau space, and transform it to the other space.  The transform from tau to matsubara is straightforward, but the transform from tau to matsubara requires careful consideration of high-frequency tails.

The repo also contains some test data, as well as a Jupyter notebook that introduces the reader to the theory of transforming with high-frequency tails and reproduces much of the same code.  

## Using the Code
The code requires 5 input command line arguments: the input file, the output file, the type of function (greens or self-energy), and dataspace of the input data (matsubara or tau), and the value of beta.  You can get a usage information by running `python MatsubaraTauTransformer.py -h`.

For example, suppose you want to transform a self-energy in matsubara space into tau space. The test data included in the repo is from DMFT on the single band hubbard model with U=3, mu=-1, and beta=5, and the matsubara self-energy is in `selfenergy_9`.  If we want the result to be saved to `SE_tau_out`, then we would run:

`python MatsubaraTauTransformer.py selfenergy_9  SE_tau_out selfenergy matsubara 5`

If instead we wanted to transform the Greens function in tau space into matsubara space, we would run:

`python MatsubaraTauTransformer.py G_tau_9 G_omega_out greens tau 5`

In general, for tau data the input file should be formatted in 2 columns (tau, f(tau)).  For matsubara data, the input file should be formatted in 3 columns (mat_freq, Re(f), Im(f))
