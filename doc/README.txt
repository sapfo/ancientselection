*"to install"
0) requirements:
- python
- numpy
- scipy
- matplotlib (if you want to use the plotting function)
Finally, for large gamma, also install:
- mpack
http://mplapack.sourceforge.net/

1) unzip ancientselection-0.2.tar.gz
download and copy ancientselection-0.2.tar.gz under whereyoudownloadedthecode
the e.g.
cd whereyoudownloadedthecode;
tar -zxvf ancientselection-0.2.tar.gz;

2) create wrapper for the lapack function spteqr etc.
cd  whereyoudownloadedthecode/ancientselection-0.2/ancientselection
make

3) create a symbolic link to where you want to run the code (there are two main scripts)
cd yourfavoritepath;
ln -s whereyoudownloadedthecode/ancientselection-0.2/ancientselection/ancientselection.py .;
ln -s whereyoudownloadedthecode/ancientselection-0.2/ancientselection/plotting.py .;

4) run
Examples of data are given under the directory "tests"

0) to get the help screen
python ancientselection.py;
or 
python ancientselection.py -h;
(idem for plotting.py)

An example dataset is given under tests:
(cd tests; then run one of the followings)

There are two main ways to compute the likelihood:
A) Computing it over fixed parameter values, essentially a cube.
In this case the likelihoods are saved for all input paramaters and can be plotted as a contour plot for example.
e.g. (under tests)

2 dimensions (fixing Ne=1000):
python ancientselection.py  -i ASIP_data_fixedNe.py --run exhaust_dim2_run1 --exhaust cube --T0dim 20 --Gammadim 20 --codedir ../ancientselection;
(takes around 9m)
then plot:
python plotting.py -i Out/Lik_array_ASIP_data_fixedNe_exhaust_dim2_run1.npy;

1 dimension (fixing gamma=0, Ne=1000)
python ancientselection.py  -i ASIP_data_sel0.py --run exhaust_dim1_run1 --exhaust cube --T0dim 100 --codedir ../ancientselection;
(takes around 2min)
then plot:
python plotting.py -i Out/Lik_array_ASIP_data_sel0_exhaust_dim1_run1.npy;

B) With an optimization algorithm (nelder mead). 
In this case the optima are saved printed to the screen.
python ancientselection.py  -i ASIP_data_fixedNe.py --run optim_run1 --codedir ../ancientselection;
(takes around: 6min)

The format for the saved output is npy format (https://github.com/numpy/numpy/blob/master/doc/neps/npy-format.txt).



