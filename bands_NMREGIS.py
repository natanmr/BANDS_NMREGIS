####################################
## Bands Strucuture from VASP     ##
## Author: Natan Moreira Regis    ##
## E-mail: n.m.regis@df.ufscar.br ##
## Licence: GPL-V3                ##
## Last update: May 31th, 2025    ##
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


import pandas as pd

def write_bands_files(path, bands_obj, shift_vbm=True, **kwargs) -> None:
    """
    Write bands.dat and other informations in gap.txt file.

    Arguments:
        bands_obj (object): bands object.
        shift_vbm (bool): Shift the VBM level to 0. Defaults to True.
    Returns: 
        None        
    """

    vbm = bands_obj.VBM() if shift_vbm else 0

    if bands_obj.ISPIN == 1:
        data = bands_obj.bands()
        df = pd.DataFrame(data=data[0], columns=data[1])
        if shift_vbm:
            df["Energy"] -= vbm

        with open(f"{path}bands.dat", "w") as f:
            prev_k = df["KPOINT"].iloc[0]
            for i, row in df.iterrows():
                curr_k = row["KPOINT"]
                if curr_k < prev_k:
                    f.write("\n")
                f.write(f"{curr_k:.6f} {row['Energy']:.6f}\n")
                prev_k = curr_k
    elif bands_obj.ISPIN == 2:
        bands_data, cols = bands_obj.bands()
        df_up = pd.DataFrame(data=bands_data[0], columns=["KPOINT", "Energy_up"])
        df_down = pd.DataFrame(data=bands_data[1], columns=["KPOINT", "Energy_down"])

        if shift_vbm:
            df_up["Energy_up"] -= vbm
            df_down["Energy_down"] -= vbm

        # Combine spin-up and spin-down data into a single DataFrame
        # Assuming KPOINTs are aligned for both spin channels
        combined_data = []
        for i in range(len(df_up)):
            combined_data.append([
                df_up["KPOINT"].iloc[i],
                df_up["Energy_up"].iloc[i],
                df_down["Energy_down"].iloc[i]
            ])
        
        df_combined = pd.DataFrame(combined_data, columns=["KPOINT", "Energy_up", "Energy_down"])

        # Write combined data to a single file
        with open(f"{path}bands.dat", "w") as f:
            # Write header
            f.write(f"{'KPOINT':<10} {'Energy_up':<15} {'Energy_down':<15}\n")
            
            prev_k = df_combined["KPOINT"].iloc[0]
            for i, row in df_combined.iterrows():
                curr_k = row["KPOINT"]
                if curr_k < prev_k:
                    f.write("\n") # Add a blank line for segment breaks
                f.write(f"{curr_k:.6f} {row['Energy_up']:.6f} {row['Energy_down']:.6f}\n")
                prev_k = curr_k
    else:
        print(f"Warning: ISPIN = {bands_obj.ISPIN} not supported for writing files.")
        
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
        self.eingenvals = []
        
        bands_per_spin = self.NBANDS // self.ISPIN if self.ISPIN == 2 else self.NBANDS


        for i in range(7, self.NKPOINTS * (self.NBANDS + 2), self.NBANDS + 2):
            k1, k2, k3, _ = map(float, EIGENVALFile[i].split())
            
            if self.ISPIN == 2:
                tmp_spin_up = []
                tmp_spin_down = []

                for j in range(i + 1, i + 1 + self.NBANDS):
                    band_data = EIGENVALFile[j].split()
                    band_index = int(band_data[0])
                    energy = float(band_data[1])
                    occupancy = float(band_data[2])

                    if band_index <= bands_per_spin:
                        tmp_spin_up.append([band_index, energy, occupancy])
                    else:
                        tmp_spin_down.append([band_index - bands_per_spin, energy, occupancy])
                self.eingenvals.append([[k1, k2, k3], [tmp_spin_up, tmp_spin_down]])
            else: # ISPIN = 1 or other unsupported ISPIN values will be treated as non-spin-polarized
                tmp = []
                for j in range(i + 1, i + 1 + self.NBANDS):
                    tmp.append(EIGENVALFile[j].split()[:])
                self.eingenvals.append([[k1, k2, k3], tmp])

class ReadPROCAR:
    '''
    Class to read PROCAR file 
    This class suport only lorbit = 10
    and SOC calculations
    '''

    def __init__(self, NATOM, NBANDS, NKPOINTS, file="PROCAR", exclude=0, path="./", ISPIN=1):
        """
        Initialize an instance of ReadPROCAR.  

        Args:
            path (str): Path to the EIGENVAL file. Defaults to "./"
            file (str): File PROCAR data. Defaults to "PROCAR"
            exclude (int): Points to exclude of bands
            ISPIN (int): Spin polarization from EIGENVAL, for PROCAR interpretation.
        """

        # Open PROCAR
        PROCAR = self.OpenPROCAR(path=path, file=file)

        self.ReadOrbitalProj(PROCAR=PROCAR, NATOMS=NATOM, NBANDS=NBANDS, NKPOINTS=NKPOINTS, ISPIN=ISPIN)

    def OpenPROCAR(self, path, file):
        """
        Open PROCAR file 
            Returns
                list: lines of PROCAR
        """
        return open_file(f"{path}/{file}").readlines()

    def ReadOrbitalProj(self, PROCAR, NKPOINTS, NBANDS, NATOMS, ISPIN):
        """
        Reads and processes the PROCAR file to extract k-point, tot, sx, sy, and sz values.

        Args:
            PROCAR (list[str]): Procar file lines
            NKPOINTS (int): Number of k-points
            NBANDS (int): Number of bands
            NATOMS (int): Number of atoms
            ISPIN (int): Spin polarization from EIGENVAL, for PROCAR interpretation.

        Returns:
            list[dict]: A list of dictionaries, each containing the k-point and associated tot, sx, sy, and sz values.
                        For ISPIN=2, this will be structured to hold spin-up and spin-down projections.
        """

        Projections = []

        nk = -1
        nb_counter = 0 # This will count total bands parsed within a k-point (including spin channels)
        
        # For ISPIN=2, NBANDS here is the total (spin-up + spin-down), so bands_per_spin is more accurate
        bands_per_spin_channel = NBANDS // ISPIN if ISPIN == 2 else NBANDS

        current_k_projections = [] # To store projections for the current k-point
        
        # Aux variables for ISPIN=2
        current_spin_up_bands = []
        current_spin_down_bands = []

        for line in PROCAR:
            elements = [ele for ele in line.strip().split() if ele]

            if not elements:
                continue

            if elements[0] == "k-point":
                if nk >= 0: # If it's not the very first k-point
                    # Before moving to the next k-point, process and save the data
                    if ISPIN == 2:
                        Projections.append([k, [current_spin_up_bands, current_spin_down_bands]])
                    else:
                        Projections.append([k, current_k_projections])

                nk += 1
                k = nk
                nb_counter = 0
                current_k_projections = []
                current_spin_up_bands = []
                current_spin_down_bands = []


            elif elements[0] == "band":
                if nb_counter > 0: # If it's not the very first band for the k-point
                    # Save the accumulated projections for the previous band
                    if ISPIN == 2:
                        # Determine if the previous band was spin-up or spin-down based on nb_counter
                        if nb_counter <= bands_per_spin_channel: # If it was a spin-up band
                             current_spin_up_bands.append({'band': (nb_counter - 1), 'tot': tot, 'sx': sx, 'sy': sy, 'sz': sz})
                        else: # If it was a spin-down band
                             current_spin_down_bands.append({'band': (nb_counter - 1 - bands_per_spin_channel), 'tot': tot, 'sx': sx, 'sy': sy, 'sz': sz})
                    else:
                        current_k_projections.append({'band': (nb_counter - 1), 'tot': tot, 'sx': sx, 'sy': sy, 'sz': sz})

                nb_counter += 1
                tot = []; sx = []; sy = []; sz = [] # Reset for the new band
                

            elif elements[0] == "tot":
                # Assuming the order is tot, sx, sy, sz
                if len(elements) > 5: # Make sure there are enough elements to parse
                    tot.append(float(elements[5]))
                    sx.append(float(elements[6]))
                    sy.append(float(elements[7]))
                    sz.append(float(elements[8]))
        
        # After the loop, add the projections for the very last band of the last k-point
        if nk >= 0:
            if ISPIN == 2:
                if nb_counter <= bands_per_spin_channel: # If it was a spin-up band
                    current_spin_up_bands.append({'band': (nb_counter - 1), 'tot': tot, 'sx': sx, 'sy': sy, 'sz': sz})
                else: # If it was a spin-down band
                    current_spin_down_bands.append({'band': (nb_counter - 1 - bands_per_spin_channel), 'tot': tot, 'sx': sx, 'sy': sy, 'sz': sz})
                Projections.append([k, [current_spin_up_bands, current_spin_down_bands]])
            else:
                current_k_projections.append({'band': (nb_counter - 1), 'tot': tot, 'sx': sx, 'sy': sy, 'sz': sz})
                Projections.append([k, current_k_projections])

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
        for each k-point based on the eigenvalues and occupancy values.

        For ISPIN=2, it calculates VBM/CBM for both spin-up and spin-down channels.
        Results are stored in self.CB_VB as a list of lists.
        If ISPIN=1: [[k_point, VBM, CBM], ...]
        If ISPIN=2: [[k_point, [VBM_up, CBM_up], [VBM_down, CBM_down]], ...]
        """

        self.CB_VB = []

        if self.ISPIN == 1:
            if self.exclude != 0:
                Ncb, Nvb = 0, 0
                for i in range(abs(self.exclude)):
                    for j in range(self.NBANDS):
                        if float(self.mesh[i][1][j][2]) == 0.0:
                            if Ncb != 0 and Nvb != 0:
                                if Ncb != int(self.mesh[i][1][j][0]) and Nvb != int(self.mesh[i][1][j - 1][0]):
                                    print("Band of the gap change!")
                            Ncb = int(self.mesh[i][1][j][0])
                            Nvb = int(self.mesh[i][1][j - 1][0])
                            break
                for k in range(self.NKPOINTS - abs(self.exclude)):
                    cb = float(self.eingenvals[k][1][Ncb - 1][1])
                    vb = float(self.eingenvals[k][1][Nvb - 1][1])
                    self.CB_VB.append([self.eingenvals[k][0], vb, cb])
            else:
                for i in range(self.NKPOINTS):
                    for j in range(1, self.NBANDS):
                        if float(self.eingenvals[i][1][j][2]) == 0.0:
                            cb = float(self.eingenvals[i][1][j][1])
                            vb = float(self.eingenvals[i][1][j - 1][1])
                            self.CB_VB.append([self.eingenvals[i][0], vb, cb])
                            break
        elif self.ISPIN == 2:
            bands_per_spin = self.NBANDS // self.ISPIN
            if self.exclude != 0:
                print("Warning: Hybrid calculations with ISPIN=2 are not fully implemented in EcEv yet. Proceeding without 'exclude' logic.")
                pass 

            for i in range(self.NKPOINTS):
                # Process spin-up bands
                vbm_up = -np.inf
                cbm_up = np.inf
                for j in range(bands_per_spin):
                    energy = float(self.eingenvals[i][1][0][j][1])
                    occupancy = float(self.eingenvals[i][1][0][j][2])
                    if occupancy > 0.0 and energy > vbm_up:
                        vbm_up = energy
                    if occupancy == 0.0 and energy < cbm_up:
                        cbm_up = energy
                
                # Process spin-down bands
                vbm_down = -np.inf
                cbm_down = np.inf
                for j in range(bands_per_spin):
                    energy = float(self.eingenvals[i][1][1][j][1])
                    occupancy = float(self.eingenvals[i][1][1][j][2])
                    if occupancy > 0.0 and energy > vbm_down:
                        vbm_down = energy
                    if occupancy == 0.0 and energy < cbm_down:
                        cbm_down = energy

                self.CB_VB.append([self.eingenvals[i][0], [vbm_up, cbm_up], [vbm_down, cbm_down]])
        else:
            print(f"ISPIN = {self.ISPIN} not supported yet.")
            exit()

    def VBM(self, spin_channel='both'):
        """
        Calculate the VBM.

        Args:
            spin_channel (str): 'up', 'down', or 'both'. Defaults to 'both'.

        Returns:
            Float: VBM value. If 'both', returns the maximum VBM across both spin channels.
        """
        valid_points = range(self.NKPOINTS - abs(self.exclude))
        
        if self.ISPIN == 1:
            self.k_vbm, self.vbm = max(((k, self.CB_VB[k][1]) for k in valid_points), key=lambda x: x[1])
            return self.vbm
        elif self.ISPIN == 2:
            vbm_values = []
            for k in valid_points:
                if spin_channel == 'up' or spin_channel == 'both':
                    vbm_values.append(self.CB_VB[k][1][0])
                if spin_channel == 'down' or spin_channel == 'both':
                    vbm_values.append(self.CB_VB[k][2][0])
            
            if not vbm_values:
                return None 

            self.vbm = max(vbm_values)
            
            if spin_channel == 'both':
                for k in valid_points:
                    if self.CB_VB[k][1][0] == self.vbm:
                        self.k_vbm = k
                        self.vbm_spin = 'up'
                        break
                    if self.CB_VB[k][2][0] == self.vbm:
                        self.k_vbm = k
                        self.vbm_spin = 'down'
                        break
            elif spin_channel == 'up':
                 for k in valid_points:
                    if self.CB_VB[k][1][0] == self.vbm:
                        self.k_vbm = k
                        self.vbm_spin = 'up'
                        break
            elif spin_channel == 'down':
                for k in valid_points:
                    if self.CB_VB[k][2][0] == self.vbm:
                        self.k_vbm = k
                        self.vbm_spin = 'down'
                        break
            return self.vbm
        return None # Should not reach here if ISPIN is handled
        
    def kvbm(self):
        self.VBM()
        return self.eingenvals[self.k_vbm][0]
    
    def kcbm(self):
        self.CBM()
        return self.eingenvals[self.k_cbm][0]

    def CBM(self, spin_channel='both'):
        """
        Calculate the CBM.

        Args:
            spin_channel (str): 'up', 'down', or 'both'. Defaults to 'both'.

        Returns:
            Float: CBM value. If 'both', returns the minimum CBM across both spin channels.
        """
        valid_points = range(self.NKPOINTS - abs(self.exclude))

        if self.ISPIN == 1:
            self.k_cbm, self.cbm = min(((k, self.CB_VB[k][2]) for k in valid_points), key=lambda x: x[1])
            return self.cbm
        elif self.ISPIN == 2:
            cbm_values = []
            for k in valid_points:
                if spin_channel == 'up' or spin_channel == 'both':
                    cbm_values.append(self.CB_VB[k][1][1])
                if spin_channel == 'down' or spin_channel == 'both':
                    cbm_values.append(self.CB_VB[k][2][1])

            if not cbm_values:
                return None 
            
            self.cbm = min(cbm_values)
            
            if spin_channel == 'both':
                for k in valid_points:
                    if self.CB_VB[k][1][1] == self.cbm:
                        self.k_cbm = k
                        self.cbm_spin = 'up'
                        break
                    if self.CB_VB[k][2][1] == self.cbm:
                        self.k_cbm = k
                        self.cbm_spin = 'down'
                        break
            elif spin_channel == 'up':
                for k in valid_points:
                    if self.CB_VB[k][1][1] == self.cbm:
                        self.k_cbm = k
                        self.cbm_spin = 'up'
                        break
            elif spin_channel == 'down':
                for k in valid_points:
                    if self.CB_VB[k][2][1] == self.cbm:
                        self.k_cbm = k
                        self.cbm_spin = 'down'
                        break
            return self.cbm
        return None # Should not reach here if ISPIN is handled

    # Calculate the averange of the states in the VB and CB
    def av_VB_CB(self):
        mean_v = 0
        mean_c = 0
        
        # This method needs to be adapted for ISPIN=2 if you want averages per spin channel
        # For now, it will only work for ISPIN=1 or return 0 for ISPIN=2 as it's not adapted
        if self.ISPIN == 1:
            for k in range(self.NKPOINTS - abs(self.exclude)):
                mean_v += self.CB_VB[k][1]
                mean_c += self.CB_VB[k][2]
            return mean_c / (k + 1) - mean_v / (k + 1)
        else:
            print("Warning: av_VB_CB is not adapted for ISPIN=2. Returning 0.")
            return 0


    def gap(self, spin_channel='both'):
        """
        Calculate the fundamental gap.

        Args:
            spin_channel (str): 'up', 'down', or 'both'. Defaults to 'both'.

        Returns:
            Float: Fundamental gap value. If 'both', returns the minimum gap across both spin channels.
        """
        if self.ISPIN == 1:
            return self.CBM() - self.VBM()
        elif self.ISPIN == 2:
            gap_up = self.CBM(spin_channel='up') - self.VBM(spin_channel='up')
            gap_down = self.CBM(spin_channel='down') - self.VBM(spin_channel='down')
            
            if spin_channel == 'up':
                return gap_up
            elif spin_channel == 'down':
                return gap_down
            elif spin_channel == 'both':
                return min(gap_up, gap_down)
        return None # Or raise an error
        

    def gap_dir(self, spin_channel='both'):
        """
        Calculate the direct gap.

        Args:
            spin_channel (str): 'up', 'down', or 'both'. Defaults to 'both'.

        Returns:
            Float: Direct fundamental gap value. If 'both', returns the minimum direct gap across both spin channels.
        """
        gap_d = 30.0 # Initialize with a large value

        valid_points = range(self.NKPOINTS - abs(self.exclude))
        
        if self.ISPIN == 1:
            for k in valid_points:
                v = self.CB_VB[k][1]
                c = self.CB_VB[k][2]
                if c - v < gap_d:
                    gap_d = c - v
                    self.k_dir = k
            return gap_d
        elif self.ISPIN == 2:
            gap_values = []
            for k in valid_points:
                if spin_channel == 'up' or spin_channel == 'both':
                    v_up = self.CB_VB[k][1][0]
                    c_up = self.CB_VB[k][1][1]
                    gap_values.append(c_up - v_up)
                
                if spin_channel == 'down' or spin_channel == 'both':
                    v_down = self.CB_VB[k][2][0]
                    c_down = self.CB_VB[k][2][1]
                    gap_values.append(c_down - v_down)
            
            if not gap_values:
                return None
            
            min_gap_overall = min(gap_values)
            self.k_dir = None
            self.dir_gap_spin = None

            for k in valid_points:
                if spin_channel == 'up' or spin_channel == 'both':
                    v_up = self.CB_VB[k][1][0]
                    c_up = self.CB_VB[k][1][1]
                    if (c_up - v_up) == min_gap_overall:
                        self.k_dir = k
                        self.dir_gap_spin = 'up'
                        break 
                
                if spin_channel == 'down' or spin_channel == 'both':
                    v_down = self.CB_VB[k][2][0]
                    c_down = self.CB_VB[k][2][1]
                    if (c_down - v_down) == min_gap_overall:
                        self.k_dir = k
                        self.dir_gap_spin = 'down'
                        break
            return min_gap_overall
        return None # Or raise an error

    def k_gap_dir(self):
        self.gap_dir()
        return self.eingenvals[self.k_dir][0]

    def VBM_dir(self, spin_channel='both'):
        """
        Get the VBM at the k-point with the smallest direct gap

        Returns:
            Float: VBM value at the k-point with the smallest direct gap
        """
        self.gap_dir(spin_channel=spin_channel) # Ensure direct gap is calculated and k_dir is set

        if self.ISPIN == 1:
            return self.CB_VB[self.k_dir][1]
        elif self.ISPIN == 2:
            if self.dir_gap_spin == 'up':
                return self.CB_VB[self.k_dir][1][0]
            elif self.dir_gap_spin == 'down':
                return self.CB_VB[self.k_dir][2][0]
            elif spin_channel == 'both': # Return the VBM corresponding to the overall direct gap
                if self.dir_gap_spin == 'up': return self.CB_VB[self.k_dir][1][0]
                if self.dir_gap_spin == 'down': return self.CB_VB[self.k_dir][2][0]
        return None

    def CBM_dir(self, spin_channel='both'):
        """
        Get the CBM at the k-point with the smallest direct gap

        Returns:
            Float: CBM value at the k-point with the smallest direct gap
        """
        self.gap_dir(spin_channel=spin_channel) # Ensure direct gap is calculated and k_dir is set

        if self.ISPIN == 1:
            return self.CB_VB[self.k_dir][2]
        elif self.ISPIN == 2:
            if self.dir_gap_spin == 'up':
                return self.CB_VB[self.k_dir][1][1]
            elif self.dir_gap_spin == 'down':
                return self.CB_VB[self.k_dir][2][1]
            elif spin_channel == 'both': # Return the CBM corresponding to the overall direct gap
                if self.dir_gap_spin == 'up': return self.CB_VB[self.k_dir][1][1]
                if self.dir_gap_spin == 'down': return self.CB_VB[self.k_dir][2][1]
        return None

    # Calculate the band structure
    def bands(self, orbital="none"):

        ks = self.K2Kpath()  # Get parametrized k-points

        if self.ISPIN == 1:
            self.bandas = []
            self.cols = ["KPOINT", "Energy"]
            if orbital == 'sz': # This case is typically for SOC in non-spin-polarized calculations
                PROCAR = ReadPROCAR(path=self.path, NBANDS=self.NBANDS, NKPOINTS=self.NKPOINTS,
                                     NATOM=self.NATOMS, ISPIN=self.ISPIN).Projection()
                self.cols = ["KPOINT", "Energy", "sigma z"]
                for j in range(self.NBANDS):
                    for k in range(self.NKPOINTS - abs(self.exclude)):
                        self.bandas.append([ks[k], float(self.eingenvals[k][1][j][1]), PROCAR[k][1][j]['sz'][0]]) # Assuming 'sz' for band 'j'
            else:
                for j in range(self.NBANDS):
                    for k in range(self.NKPOINTS - abs(self.exclude)):
                        self.bandas.append([ks[k], float(self.eingenvals[k][1][j][1])])
            return self.bandas, self.cols
        
        elif self.ISPIN == 2:
            self.bandas_up = []
            self.bandas_down = []
            # For ISPIN=2, the 'bands' method will return two separate lists of bands for plotting.
            # The columns will be implicit in how they are used by the plotting function.
            self.cols = ["KPOINT", "Energy"] # This is just a placeholder for the structure of each sub-list
            bands_per_spin = self.NBANDS // self.ISPIN

            if orbital != 'none':
                print("Warning: Orbital projection with ISPIN=2 is not fully implemented yet for this method.")
                print("Returning energies only.")
            
            for j in range(bands_per_spin): # Iterate through bands per spin channel
                for k in range(self.NKPOINTS - abs(self.exclude)):
                    energy_up = float(self.eingenvals[k][1][0][j][1])
                    energy_down = float(self.eingenvals[k][1][1][j][1])
                    
                    self.bandas_up.append([ks[k], energy_up])
                    self.bandas_down.append([ks[k], energy_down])

            # Return a tuple of (spin_up_bands, spin_down_bands)
            return (self.bandas_up, self.bandas_down), self.cols
        else:
            print(f"Warning: ISPIN = {self.ISPIN} not supported for bands plotting. Returning empty data.")
            return [], []


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
                # Need to be careful here if self.eingenvals[k][0] structure changes for ISPIN=2
                # It currently retrieves k-point coordinates from the first element of the inner list
                # which should be consistent regardless of ISPIN
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
            max_k = ks[-1] if ks else 1.0 # Avoid division by zero if ks is empty
            for k_idx in range(len(ks)):
                ks[k_idx] = ks[k_idx] / max_k

        return ks
    
################
## Plot utils ##
################
def plot_band_axis(ax, bands_obj, shift_vbm = True, xmin = 0, xmax = 1, ymin = -4, ymax = 4, k_labels = None):

    # Get VBM (for shifting if needed)
    vbm_to_shift = bands_obj.VBM() if shift_vbm else 0

    special_points_num = []

    if bands_obj.ISPIN == 1:
        data, cols = bands_obj.bands()
        df = pd.DataFrame(data=data, columns=cols)
        
        kold = 0
        for k_idx, kp in enumerate(df["KPOINT"]):
            if k_idx > 0 and df["KPOINT"][k_idx] < df["KPOINT"][k_idx-1]: # Check for segment break
                special_points_num.append(df["KPOINT"][k_idx-1])
                ax.plot(df["KPOINT"][kold:k_idx], df["Energy"][kold:k_idx]-vbm_to_shift, color="black", linestyle="-", label="Bands")
                kold = k_idx
        # Plot the last segment
        ax.plot(df["KPOINT"][kold:], df["Energy"][kold:]-vbm_to_shift, color="black", linestyle="-", label="Bands")
        # Add the last k-point if it's not already there and the DataFrame is not empty
        if not df.empty and df["KPOINT"].iloc[-1] not in special_points_num:
            special_points_num.append(df["KPOINT"].iloc[-1])

    elif bands_obj.ISPIN == 2:
        bands_data, cols = bands_obj.bands()
        df_up = pd.DataFrame(data=bands_data[0], columns=["KPOINT", "Energy"])
        df_down = pd.DataFrame(data=bands_data[1], columns=["KPOINT", "Energy"])

        # Flags to ensure labels appear only once in the legend
        first_spin_up_plot = True
        first_spin_down_plot = True

        # Plot spin-up bands
        kold_up = 0
        for k_idx, kp in enumerate(df_up["KPOINT"]):
            if k_idx > 0 and df_up["KPOINT"][k_idx] < df_up["KPOINT"][k_idx-1]: # Check for segment break
                special_points_num.append(df_up["KPOINT"][k_idx-1])
                ax.plot(df_up["KPOINT"][kold_up:k_idx], df_up["Energy"][kold_up:k_idx]-vbm_to_shift, 
                        color="blue", linestyle="-", label="Spin Up" if first_spin_up_plot else None)
                first_spin_up_plot = False
                kold_up = k_idx
        if not df_up.empty:
            ax.plot(df_up["KPOINT"][kold_up:], df_up["Energy"][kold_up:]-vbm_to_shift, 
                    color="blue", linestyle="-", label="Spin Up" if first_spin_up_plot else None)
            first_spin_up_plot = False


        # Plot spin-down bands
        kold_down = 0
        for k_idx, kp in enumerate(df_down["KPOINT"]):
            if k_idx > 0 and df_down["KPOINT"][k_idx] < df_down["KPOINT"][k_idx-1]: # Check for segment break
                if df_down["KPOINT"][k_idx-1] not in special_points_num: # Avoid duplicate special points
                    special_points_num.append(df_down["KPOINT"][k_idx-1])
                ax.plot(df_down["KPOINT"][kold_down:k_idx], df_down["Energy"][kold_down:k_idx]-vbm_to_shift, 
                        color="red", linestyle="--", label="Spin Down" if first_spin_down_plot else None)
                first_spin_down_plot = False
                kold_down = k_idx
        if not df_down.empty:
            ax.plot(df_down["KPOINT"][kold_down:], df_down["Energy"][kold_down:]-vbm_to_shift, 
                    color="red", linestyle="--", label="Spin Down" if first_spin_down_plot else None)
            first_spin_down_plot = False

        # Add the last k-point if it's not already there and the DataFrame is not empty
        if not df_up.empty and df_up["KPOINT"].iloc[-1] not in special_points_num:
            special_points_num.append(df_up["KPOINT"].iloc[-1])
        
        # Sort and unique special points
        special_points_num = sorted(list(set(special_points_num)))

    # set lines on special points
    for sp in special_points_num:
        ax.axvline(x=sp, color="black", linewidth=1.0, linestyle="dashed")

    # Line on fermi level
    ax.axhline(y=0, color="black", linewidth=1.0, linestyle="dotted", label=r"$E_{Fermi}$")

    # set xticks
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
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.legend() # Add a legend to differentiate spin up/down

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
    print("1. BANDS (Detects ISPIN automatically from EIGENVAL) ")
    print("2. Bands with Hybrid Functional (Not fully implemented for ISPIN=2 yet) ")
    print("3. BANDS with SOC (Not fully implemented for ISPIN=2 yet) ")

    type_analysis = input("Enter your answer (Defaults to 1): ") or "1" # Renamed from type_bands

    path_calc = input("Enter the path for bands calculations (Default: ./): ") or "./"

    normalize_input = input("Do you want to normalize the k-path (Defaults to True)? (y/n): ") or "y"
    normalize = normalize_input.lower() == 'y'

    k_labels = input("Write the k-labels of special points in BZ (Ex.: \Gamma, K, M, \Gamma) : ") or None

    shift_vbm_input = input("Shift the VBM level to zero (Defaults to True)? (y/n): ") or "y"
    shift_vbm = shift_vbm_input.lower() == 'y'

    write_files_input = input("Do you want to write the bands files (y/n. Defaults to y)? ") or "y"
    write_files = write_files_input.lower() == 'y'

    if type_analysis == "1":
        BANDS_data = Bands(path=path_calc, file="EIGENVAL", exclude=0, normalize=normalize)
        
        figure_filename = "bands.png"
        if BANDS_data.ISPIN == 2:
            figure_filename = "bands_spin_polarized.png"

        fig, ax = plt.subplots()
        ax = plot_band_axis(ax=ax, bands_obj=BANDS_data, k_labels=k_labels, shift_vbm=shift_vbm)
        plt.savefig(figure_filename, dpi=100)
        
        if write_files:
            write_bands_files(path=path_calc, bands_obj=BANDS_data, shift_vbm=shift_vbm)
            if BANDS_data.ISPIN == 2:
                print(f"Bands data written to {path_calc}bands.dat (combined spin channels).")
            else:
                print(f"Bands data written to {path_calc}bands.dat")

    elif type_analysis in ["2", "3"]:
        print("Coming soon for these analysis types...")
        # If you implement these, remember to consider ISPIN automatically.
        exit()
    else:
        print("Invalid option. Exiting.")
        exit()

if __name__ == '__main__':
    menu()