# greensolverfd
 
    Given an arbitrary potential this code can solve the corresponding 
one-dimensional Schrodinger equation using a finite difference method plus 
either a shooting Green's function method for fixed chemical potential or 
a quadratic eigenvalue problem for fixed number of electrons.

    Also given a one-dimensional density this code can find the 
corresponding one-dimensional potential using either of the two methods 
mentioned above.

    Run the setup.m script to add subfolders to the path and ensure that 
the resp_mex function has been compiled.  Examples are available in
Examples folder and runtests.m should complete in under 10 seconds.
