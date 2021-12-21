To add new materials on a local installation in a conda environment.

1. Activate your desired conda environment; if you don't already have
a dedicated conda environment for wptherml, create one like this:
`conda create --name wptherml_materials python=3.9`

2. Activate your conda environment
`conda activate wptherml_materials`

3. If this is a new environment, make sure you have dependencies installed
`conda install numpy scipy matplotlib`
`conda install -c conda-forge notebook`

4. Change directories to the location of your wptherml repository folder.
If you don't have a local repository of wptherml, you can clone one
from github:

`git clone https://github.com/FoleyLab/wptherml.git`

5. Change directories to the wptherml/wptherml/datalib folder.  This
is where the data files for different materials are, and also where the python file is that reads these data files and stores them to refractive index arrays for the materials you have specified in your main wptherml script.

6. For each new material you would like to add, put a text file here with a
descriptive name that contains the refractive index data vs wavelength in the following format:
column 1: wavelength in meters
column 2: real part of refractive index vs wavelength (aligned with column 1)
column 3: imaginary part of refractive index vs wavelength (aligned with column 1)

7. Once all data files have been added, open datalib.py in an editor and make the following changes:

a. Add the chemical names of your new material to the list of supported_materials near the top of the file.  For concreteness, I will pretend I am adding the data file for silicon nitride, so I will add 'SiN' to the list of supported materials.

b. Find the comment that says "#All materials with data files need a if-statement here!".  In any one of those if statements (or in a new one), you need
to add an if statement that evaluates to true if the user has specified 
'SiN' as one or more of the materials in their main wptherml script.
A before and after of these lines of code follows:

# Before addition
    elif (arg=='W' or arg=='Re' or arg=='Rh' or arg=='Ru' or arg=='Al'):
        n = Read_RI_from_File(lam, arg)

# After addition
    elif (arg=='SiN' or arg=='W' or arg=='Re' or arg=='Rh' or arg=='Ru' or arg=='Al'):
        n = Read_RI_from_File(lam, arg)

8. Find the function 'def Read_RI_from_File(lam, matname): and add a line that will
catch if the matname matches your new material 'SiN'... Again a Before and After example

# Before addition
def Read_RI_from_File(lam, matname):
    if (matname=='W'):
        file_path = path + 'W_Palik_RI_f.txt'
        a = np.loadtxt(file_path)
    elif (matname=='TiO2'):
        file_path = path + 'TiO2_Siefke.txt'

# After addition
def Read_RI_from_File(lam, matname):
    if (matname=='W'):
        file_path = path + 'W_Palik_RI_f.txt'
        a = np.loadtxt(file_path)
    elif (matname=='SiN'):
        file_path = path + 'SiN_RI.txt'
        a = np.loadtxt(file_path)
    elif (matname=='TiO2'):
        file_path = path + 'TiO2_Siefke.txt'

8.  Once added, you need to build wptherml... from the main wptherml directory, simply type

`python setup.py install`
