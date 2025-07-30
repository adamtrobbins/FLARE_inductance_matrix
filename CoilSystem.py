import numpy as np

class CoilSystem: 

    def __init__(self, data_dir):

        self.L = np.loadtxt(f'{data_dir}/L.csv', delimiter=',')
        self.R = np.loadtxt(f'{data_dir}/R.csv')
        self.coil_names = np.genfromtxt(f'{data_dir}/coil_names.csv', dtype=str)
        self.coil_coords = np.loadtxt(f'{data_dir}/coil_coords.csv', delimiter=',')

        if (len(set([len(self.L), len(self.R), len(self.coil_names), len(self.coil_coords)])) != 1):
            raise ValueError("Inconsistent # of coils in data files")
        
        self.NUM_COILS = len(self.R)

        self.enabled_mask = np.loadtxt(f'{data_dir}/enabled_coils.csv', dtype=bool)
        self.__disable_coils()

    def prune_coil_array(self, arr):
        if len(arr) != len(self.enabled_mask):
            raise ValueError(f"Array length does not match number of coils: {len(arr)} != {len(self.enabled_mask)}")
        
        return arr[self.enabled_mask]
    
    def get_cindex(self):
        return {name: idx for idx, name in enumerate(self.coil_names)}

    def __disable_coils(self):
        enabled_mask = self.enabled_mask
        if len(enabled_mask) != self.NUM_COILS:
            raise ValueError("Enabled mask length does not match number of coils")
        
        disabled_coils_names = self.coil_names[~enabled_mask]

        self.L = self.L[enabled_mask][:, enabled_mask]
        self.R = self.R[enabled_mask]
        self.coil_names = self.coil_names[enabled_mask]
        self.coil_coords = self.coil_coords[enabled_mask]

        self.NUM_COILS = len(self.R)
