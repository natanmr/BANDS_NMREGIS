####################################
## Bands Strucuture from VASP     ##
## Author: Natan Moreira Regis    ##
## E-mail: n.m.regis@df.ufscar.br ##
## Licence: GPL-V3                ##                               
## Last update: August 6, 2024    ##
####################################

########################
## External Libraries ##
########################
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd

###################
## IO ultilities ##
###################
def open_file(filename):
    """
    Open a file for reading.

    Arguments:
        filename (str): The name of the file to open.

    Returns:
        file object: The opened file object.

    Raises:
        SystemExit: If the file is not found.
    """

    # Attempt to open the file in read mode
    try:
        file = open(filename, 'r')
        return file
    # Handle the case where the file is not found    
    except FileNotFoundError:
        print(f"File {filename} not found. Exiting...")
        exit()

def read_poscar(path):
    """
    Read the POSCAR file from a given directory.

    Args:
        path (str): The directory path where the POSCAR file is located.

    Returns:
        list: Lines of the POSCAR file.
    """
    # Open and read the POSCAR file
    return open_file(f"{path}/POSCAR").readlines()


def write_bands_files(bands_obj) -> None:
    """
    Write bands.dat and other informations in gap.txt file

    Arguments:
        bands_obj (object): bands object:

    Returns: 
        None        
    """

    


class ReadEIGENVAL:
    '''
    This class recives the file argument
    and return the Content of Eigenval
    optionally, If you were performed 
    hybrid calculations, you can exclude 
    the fake bands at the beginning (exclude>0)
    of at the end (exclude<0) of the EIGENVAL
    file.
    '''

    def __init__(self, path="./", file="EIGENVAL", exclude=0):
        """
        Initialize an instance of ReadEIGENVAL.  

        Args:
            path (str): Path to the EIGENVAL file. Defaults to "./"
            file (str): File with EIGENVAL data. Defaults to "EIGENVAL"
            exclude (int): Points to exclude of bands
        """

        # Open EIGENVAL file
        EIGENVALFile = self.OpenEIGENVAL(path=path, file=file)

        # Read header informations
        self.ParseEIGENVALHeader(EIGENVALFile)

        # Read bands data
        self.ParseEIGENVALValues(EIGENVALFile)

        # Exclude specific k-points:
        self.eingenvals, self.mesh = exclude_k_points(exclude, self.NKPOINTS, self.eingenvals)

    def OpenEIGENVAL(self, path, file):
        """
        Open EIGENVAL file 
        """

        # Read lines of EIGENVAL file
        return open_file(f"{path}/{file}").readlines()

    def ParseEIGENVALHeader(self, EIGENVALFile):
        """
        Read Header of EIGENVAL file
        """

        # Spin polarization of the system
        self.NATOMS, _, _, self.ISPIN = map(int, EIGENVALFile[0].split())

        # Volume and lattice parameters
        self.v0, self.a0, self.b0, self.c0, _ = map(float, EIGENVALFile[1].split())

        # Name of the system
        self.name = EIGENVALFile[4].split()[0]

        # Number of electrons, k-points and bands
        self.NELECT, self.NKPOINTS, self.NBANDS = map(int, EIGENVALFile[5].split())

    def ParseEIGENVALValues(self, EIGENVALFile):
        """ 
        Read eingenvalues in each k-point 
        """
        # Store the k-points and eigenvals
        self.eingenvals = []

        # Run over all k-points
        for i in range(7, self.NKPOINTS * (self.NBANDS + 2), self.NBANDS + 2):

            # Read the k-point and the weight
            k1, k2, k3, _ = map(float, EIGENVALFile[i].split())

            # Temporary variable to store the eigenvalues 
            tmp = []

            # Run over the nbands
            for j in range(i + 1, i + 1 + self.NBANDS):
                # append each band's engenvals
                tmp.append(EIGENVALFile[j].split()[:])

                # Save the k-point and the correspondent bands (eingenvals)
            self.eingenvals.append([[k1, k2, k3], tmp])

class ReadPROCAR:
    '''
    Class to read PROCAR file 
    This class suport only lorbit = 10
    and SOC calculations
    '''

    def __init__(self, NATOM, NBANDS, NKPOINTS, file="PROCAR", exclude=0, path="./"):
        """
        Initialize an instance of ReadPROCAR.  

        Args:
            path (str): Path to the EIGENVAL file. Defaults to "./"
            file (str): File PROCAR data. Defaults to "PROCAR"
            exclude (int): Points to exclude of bands
        """

        # Open PROCAR
        PROCAR = self.OpenPROCAR(path=path, file=file)

        self.ReadOrbitalProj(PROCAR=PROCAR, NATOMS=NATOM, NBANDS=NBANDS, NKPOINTS=NKPOINTS)

    def OpenPROCAR(self, path, file):
        """
        Open PROCAR file 
            Returns
                list: lines of PROCAR
        """
        return open_file(f"{path}/{file}").readlines()

    def ReadOrbitalProj(self, PROCAR, NKPOINTS, NBANDS, NATOMS):
        """
        Reads and processes the PROCAR file to extract k-point, tot, sx, sy, and sz values.

        Args:
            PROCAR (list[str]): Procar file lines
            NKPOINTS (int): Number of k-points
            NBANDS (int): Number of bands
            NATOMS (int): Number of atoms

        Returns:
            list[dict]: A list of dictionaries, each containing the k-point and associated tot, sx, sy, and sz values.
        """

        # Save the projection's values
        Projections = []

        # K-points and bands counter
        nk = -1
        nb = 0

        for line in PROCAR:
            # Split the line into non-empty elements
            elements = [ele for ele in line.strip().split() if ele]

            # Ignore blank lines
            if not elements:
                continue

            # Check k-point line, where begins new data block
            if elements[0] == "k-point":

                # Update number of k-points
                nk += 1

                # Reset number of bands
                nb = 0

                # Save k-point number
                k = nk

                # Reset projections saving variable
                tot = [];
                sx = [];
                sy = [];
                sz = []

                # Aux variable to save results only once
                aux = 0

            elif elements[0] == "band":

                # Update band counter
                nb += 1

                # Update projection counter
                counter = 0

            elif elements[0] == "tot":
                # Update the current dictionary with tot, sx, sy, sz values
                if counter == 0:
                    tot.append(elements[5])
                elif counter == 1:
                    sx.append(elements[5])
                elif counter == 2:
                    sy.append(elements[5])
                elif counter == 3:
                    sz.append(elements[5])
                counter += 1

            # Save the bands values
            if nb == NBANDS and aux < 1:
                Projections.append([k, sz])
                aux += 1

        self.Projections = Projections

    def Projection(self):
        return self.Projections

###########################
## Structural ultilities ##
###########################
def lv2_rv(poscar):
    """
    Calculate the reciprocal lattice vectors given a POSCAR file content.

    Arguments:
        poscar (list): List of strings representing the lines of a POSCAR file.

    Returns:
        tuple: Reciprocal lattice vectors b1, b2, and b3 as numpy arrays.
    """
    
    # Get lattice parameters
    a1, a2, a3 = get_lattice_vectors(poscar=poscar)

    # Calculate the volume of the unit cell
    volume = np.dot(a1, np.cross(a2, a3))

    # Calculate the reciprocal lattice vectors
    b1 = 2 * np.pi * np.cross(a2, a3) / volume
    b2 = 2 * np.pi * np.cross(a3, a1) / volume
    b3 = 2 * np.pi * np.cross(a1, a2) / volume

    return b1, b2, b3

def get_lattice_vectors(poscar):
    """
    Get lattice vectors from VASP's POSCAR file

    Arguments: 
        poscar (list): List of strings representing the lines of a POSCAR file.
    
    Returns: 
        tuple: lattice vectors.

    """
    # Parse the lattice vectors from the POSCAR file
    a1 = np.array([float(element) for element in poscar[2].split()])
    a2 = np.array([float(element) for element in poscar[3].split()])
    a3 = np.array([float(element) for element in poscar[4].split()])

    return a1, a2, a3

################
## Ultilities ##
################
def exclude_k_points(exclude, nkpoints, values):
    '''
    Exclude some k-points from bands structure
    This is userfull for DFT hybrid calculations

    Args:
        exclude (int): The k-points for exclude.
        use positive values for exclude begining points 
        or negative for ending points. 

        nkpoints (int): Number of k-points.

        values (list): List of values to exclude.

    Returns:
        values (list): a list with updated k-points
        mesh (list): The removed k-points

    Raises:
        Garantees that at least one k-point is returning
    '''

    # List with excluded points
    mesh = []

    # Check if the argument exclude is different from zero and if the value is lowest than the number of k-points
    if exclude != 0 and abs(exclude) < nkpoints:

        # Case the k-points to exclude is in beginning of EIGENVAL file
        if exclude > 0:
            mesh = values[:exclude]
            values = values[exclude:]

        # Case the k-points to exclude is in end of EIGENVAL file
        if exclude < 0:
            mesh = values[exclude:]
            values = values[:exclude]

    # Case without exclude k-points
    else:
        pass

    return values, mesh

######################
## Bands ultilities ##
######################

class Bands:
    """ 
    Class for get band structure and other informations
        Atributes
            - path (str): path for reading data
            - npath (int): Number of k-points. Set the number of k-points in case of hybrid calculations.
            - file (str, optional): EIGENVAL file. Defaults to "EIGENVAL". 
            - normalize (Bool, optional): Normalise the k-path from 0 to 1. Defaults to "True"
            - exclude (int, optional): Exclude bands. Important for hybrid calculations. Set positive for excluding bands in begining of EIGENVAL, or negative in the end. Defauls to 0, i.e. do not exclude any bands

    ## Methods
        -  
    """

    def __init__(self, path, file="EIGENVAL", normalize=True, exclude=0, **kwargs):

        self.path = path

        self.read_eigenval = ReadEIGENVAL(path=path, file=file)  # Read eingenval file

        self.name = self.read_eigenval.name  # name of system

        self.NBANDS = self.read_eigenval.NBANDS  # Number of bands

        self.NATOMS = self.read_eigenval.NATOMS

        self.mesh = self.read_eigenval.mesh  # Mesh points

        self.NKPOINTS = self.read_eigenval.NKPOINTS  # Number of k-poins

        self.eingenvals = self.read_eigenval.eingenvals  # Read eigenval

        self.normalize = normalize  # Normalize k-path

        self.exclude = exclude  # Bands to exclude

        self.ISPIN = self.read_eigenval.ISPIN  # System spin

        # Read POSCAR and get lattice parameters
        try:
            with open(f"{path}/POSCAR", 'r') as poscar:
                self.poscar = poscar.readlines()
                [self.b1, self.b2, self.b3] = lv2_rv(self.poscar)
        except FileNotFoundError:
            print("POSCAR not found")
            exit()

        self.EcEv()  # Read VBM and CBM in each k-point

    def EcEv(self):
        """
        Calculate the states in the valence and conduction band in each k-point.

        This function determines the conduction band minimum (CBM) and valence band maximum (VBM)
        for each k-point based on the eigenvalues and occupancy values in the self.eingenvals and
        self.mesh attributes.

        If ISPIN is not equal to 1, the function will print a message and exit.
        For hybrid calculations (when self.exclude is not 0), it calculates the band indices
        of the transition points using mesh points. For non-hybrid calculations, it directly 
        determines the CBM and VBM for each k-point.

        Results are stored in the self.CB_VB attribute as a list of lists, where each inner list 
        contains the k-point identifier, VBM, and CBM.

        Attributes:
            self.CB_VB (list): List to store the results of CBM and VBM for each k-point.
        """

        self.CB_VB = []  # Store the results

        # Check if ISPIN is supported
        if self.ISPIN != 1:
            print("Not implemented yet...")
            exit()

        if self.exclude != 0:
            Ncb, Nvb = 0, 0

            # Calculate the band of transition in the mesh points
            for i in range(self.exclude):
                for j in range(self.NBANDS):
                    if float(self.mesh[i][1][j][2]) == 0.0:
                        if Ncb != 0 and Nvb != 0:
                            if Ncb != int(self.mesh[i][1][j][0]) and Nvb != int(self.mesh[i][1][j - 1][0]):
                                print("Band of the gap change!")
                        Ncb = int(self.mesh[i][1][j][0])
                        Nvb = int(self.mesh[i][1][j - 1][0])
                        break

            # Calculate the valence and conduction states in each k-point given the band index
            for k in range(self.NKPOINTS - abs(self.exclude)):
                cb = float(self.eingenvals[k][1][Ncb - 1][1])
                vb = float(self.eingenvals[k][1][Nvb - 1][1])
                self.CB_VB.append([self.eingenvals[k][0], vb, cb])
        else:
            # In case of non-hybrid calculations, calculate the transition
            for i in range(self.NKPOINTS):
                for j in range(1, self.NBANDS):
                    if float(self.eingenvals[i][1][j][2]) == 0.0:
                        cb = float(self.eingenvals[i][1][j][1])
                        vb = float(self.eingenvals[i][1][j - 1][1])
                        self.CB_VB.append([self.eingenvals[i][0], vb, cb])
                        break

    def VBM(self):
        """
        Calculate the VBM

        Returns:
            Float: VBM value
        """
        valid_points = range(self.NKPOINTS - abs(self.exclude))
        self.k_vbm, self.vbm = max(((k, self.CB_VB[k][1]) for k in valid_points), key=lambda x: x[1])
        return self.vbm

    def kvbm(self):
        self.VBM()
        return self.eingenvals[self.k_vbm][0]
    
    def kcbm(self):
        self.CBM()
        return self.eingenvals[self.k_cbm][0]

    def CBM(self):
        """
        Calculate the CBM

        Returns:
            Float: CBM value
        """
        valid_points = range(self.NKPOINTS - abs(self.exclude))
        self.k_cbm, self.cbm = min(((k, self.CB_VB[k][2]) for k in valid_points), key=lambda x: x[1])
        return self.cbm

    # Calculate the averange of the states in the VB and CB
    def av_VB_CB(self):
        mean_v = 0
        mean_c = 0
        for k in range(self.NKPOINTS - abs(self.exclude)):
            mean_v += self.CB_VB[k][1]
            mean_c += self.CB_VB[k][2]

        return mean_c / (k + 1) - mean_v / (k + 1)

    def gap(self):
        """
        Calculate the fundamental gap

        Returns:
            Float: Fundamental gap value
        """
        return self.CBM() - self.VBM()

    def gap_dir(self):
        """
        Calculate the direct gap

        Returns:
            Float: Direct fundamental gap value
        """
        # Initialize the minimum direct gap with a large value
        gap_d = 30

        # Iterate through valid k-points
        valid_points = range(self.NKPOINTS - abs(self.exclude))
        for k in valid_points:
            v = self.CB_VB[k][1]  # Valence band maximum at k-point
            c = self.CB_VB[k][2]  # Conduction band minimum at k-point

            # Update the minimum direct gap if a smaller value is found
            if c - v < gap_d:
                gap_d = c - v
                self.k_dir = k  # Store the k-point with the smallest direct gap

        # Return the minimum direct gap
        return gap_d

    def k_gap_dir(self):
        self.gap_dir()
        return self.k_dir

    def VBM_dir(self):
        """
        Get the VBM at the k-point with the smallest direct gap

        Returns:
            Float: VBM value at the k-point with the smallest direct gap
        """
        # Ensure the k-point with the smallest direct gap is calculated
        self.gap_dir()

        # Return the VBM value at the k-point with the smallest direct gap
        return self.CB_VB[self.k_dir][1]

    def CBM_dir(self):
        """
        Get the CBM at the k-point with the smallest direct gap

        Returns:
            Float: CBM value at the k-point with the smallest direct gap
        """
        # Ensure the k-point with the smallest direct gap is calculated
        self.gap_dir()

        # Return the CBM value at the k-point with the smallest direct gap
        return self.CB_VB[self.k_dir][2]

    # Calculate the band structure
    def bands(self, orbital="none"):

        self.bandas = []  # Save the results

        if orbital == 'sz':

            self.cols = ["KPOINT", "Energy", "sigma z"]

            PROCAR = ReadPROCAR(path=self.path, NBANDS=self.NBANDS, NKPOINTS=self.NKPOINTS,
                                NATOM=self.NATOMS).Projection()

            ks = self.K2Kpath()  # Get parametrized k-points
            for j in range(self.NBANDS):
                for k in range(self.NKPOINTS - abs(self.exclude)):
                    self.bandas.append([ks[k], float(self.eingenvals[k][1][j][1]), float(PROCAR[k][1][j])])

        else:
            self.cols = ["KPOINT", "Energy"]

            ks = self.K2Kpath()  # Get parametrized k-points
            for j in range(self.NBANDS):
                for k in range(self.NKPOINTS - abs(self.exclude)):
                    self.bandas.append([ks[k], float(self.eingenvals[k][1][j][1])])

        return self.bandas, self.cols

    def k2BZ(self, type="VBM", special_points = []):

        if type == "VBM":
            k = self.kvbm()
        elif type == "CBM":
            k = self.kcbm()
        elif type == "dir":
            k = self.k_gap_dir()
        else:
            print("Provite the type of gap that you want. Exiting...")
            exit()

        print("Type k-gap: ", type)

        # Get unitcell type and special points
        a, b, c = get_lattice_vectors(self.poscar)
        
        ks = self.K2Kpath()

        return ""

    def K2Kpath(self):
        ks = [0]
        dk = [[0., 0., 0.]]
        for k in range(self.NKPOINTS - abs(self.exclude)):
            if k == 0:
                kmod = 0
            else:
                dk.append([self.b1[0] * self.eingenvals[k][0][0] + self.b2[0] * self.eingenvals[k][0][1] + self.b3[0] *
                           self.eingenvals[k][0][2],
                           self.b1[1] * self.eingenvals[k][0][0] + self.b2[1] * self.eingenvals[k][0][1] + self.b3[1] *
                           self.eingenvals[k][0][2],
                           self.b1[2] * self.eingenvals[k][0][0] + self.b2[2] * self.eingenvals[k][0][1] + self.b3[2] *
                           self.eingenvals[k][0][2]
                           ])
                kmod += ((dk[k - 1][0] - dk[k][0]) ** 2 + (dk[k - 1][1] - dk[k][1]) ** 2 + (
                        dk[k - 1][2] - dk[k][2]) ** 2) ** (1. / 2.)
                ks.append(kmod)

        if self.normalize == True:
            for k in range(self.NKPOINTS - abs(self.exclude)):
                ks[k] = ks[k] / ks[-1]

        return ks
    
################
## Plot utils ##
################
def plot_band_axis(ax, bands_obj, shift_fermi = True, xmin = 0, xmax = 1, ymin = -4, ymax = 4, k_labels = None):

    data = bands_obj.bands()
    data = pd.DataFrame(data = data[0], columns=data[1])

    # Get VBM:
    vbm_tbm = bands_obj.VBM()

    special_points_num = []

    kold = 0
    for k, kp in enumerate(data["KPOINT"]):
            
        if k > 1 and data["KPOINT"][k] == data["KPOINT"][k-1]:
            if data["KPOINT"][k] not in special_points_num:
                        special_points_num.append(kp)
        if kp == 1.0:
            ax.plot(data["KPOINT"][kold:k+1], data["Energy"][kold:k+1]-vbm_tbm, color = "black", linestyle = "-", label = "Bands")
            kold = k+1

    # set lines on special points
    for sp in special_points_num:
        ax.axvline(x = sp, color="black", linewidth=1.0, linestyle="dashed")

    # Line on fermi level
    ax.axhline(y = 0, color="black", linewidth=1.0, linestyle="dotted", label = r"$E_{Fermi}$")

    # set xticks
    special_points_num.append(data["KPOINT"][k])
    special_points_num.insert(0, data["KPOINT"][0])

    print(k_labels)
    if k_labels is None:
        ax.set_xticks(ticks=special_points_num)
        ax.set_xticklabels([])
    else:
        k_labels_list = [rf"${k.strip()}$" for k in k_labels.split(",")]
        if len(k_labels_list) == len(special_points_num):
            ax.set_xticks(ticks=special_points_num)
            ax.set_xticklabels(k_labels_list)
        else:
            print("Wrong number of k-labels... Using ticks with no labels.")
            ax.set_xticks(ticks=special_points_num)
            ax.set_xticklabels([])
            
    # axis limits
    ax.set_xlim(xmin, xmax); ax.set_ylim(ymin, ymax)

    return ax

#####################
## Menu ultilities ##
#####################
def menu() -> None:

    print("#########################")
    print("## VASP bands analysis ##")
    print("## Author: NMRegis     ##")
    print("## Licence: GPL-V3     ##")
    print("#########################")
    
    print("Which analysis do you wanna to do?")
    print("1. BANDS  ")
    print("2. Bands with Hybrid Functional ")
    print("3. BANDS with SOC ")

    type_bands = input("Enter you answer (Defaults to 1): ") or "1"
    
    path_calc = input("Enter the path for bands calculations (Default: ./): ") or "./"

    normalize = input("Do you want to normalize the k-path (Defaults to True)? ") or True

    k_labels = input("Write the k-labels of special points in BZ (Ex.: \Gamma, K, M, \Gamma) : ") or None

    if type_bands == "1":
        BANDS_data = Bands(path=path_calc, file = "EIGENVAL", exclude=0, normalize=normalize)
        fig, ax, = plt.subplots()
        ax = plot_band_axis(ax=ax, bands_obj=BANDS_data, k_labels=k_labels)
        plt.savefig("bands.png", dpi = 100)

    elif type_bands == "2":
        print("Comming soon...")
        exit()
    elif type_bands == "3":
        print("Comming soon...")
        exit()

    write_files = input("Do you want to write the bands files (y/n. Defaults to y)? ") or "y"
    if write_files == "y":
        print("Comming soon...")
        exit()
        write_bands_files(bands_obj = BANDS_data)
    
if __name__ == '__main__':
    menu()

