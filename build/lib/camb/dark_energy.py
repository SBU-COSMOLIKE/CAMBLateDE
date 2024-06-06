from .baseconfig import F2003Class, fortran_class, numpy_1d, CAMBError, np, \
    AllocatableArrayDouble, f_pointer
from ctypes import c_int, c_double, byref, POINTER, c_bool


class DarkEnergyModel(F2003Class):
    """
    Abstract base class for dark energy model implementations.
    """
    _fields_ = [
        ("__is_cosmological_constant", c_bool),
        ("__num_perturb_equations", c_int),
        ("w_lam", c_double),
        ("wa",c_double),
        ("cs2_lam", c_double),
        ("no_perturbations",c_bool)]

    def validate_params(self):
        return True

@fortran_class
class LateDE(DarkEnergyModel):
    """
    DHFS: Abstract base class for models using w and wa parameterization with use w(a) = w + (1-a)*wa parameterization,
    or bin w.
    """
    
    _fortran_class_module_ = 'LateDE'
    _fortran_class_name_ = 'TLateDE'

    _fields_ = [
        ("DEmodel", c_int, "select one model among five: (1) w=cte, (2) CPL, (3) 3 bins w, (4) 5 bins w (5) 10 bins w"),
        ("max_num_of_bins", c_int, "Maximum number of bins"),
        ("z_knot", AllocatableArrayDouble, "Array of redshift bins"),
        ("w_knot", AllocatableArrayDouble, "Array of w in each knot"),        
        # Equation of State
        ("w0", c_double, "Bin w parameter: EoS for the 1st bin"),
        ("w1", c_double, "Bin w parameter: EoS for the 2nd bin"),
        ("w2", c_double, "Bin w parameter: EoS for the 3rd bin"),
        ("w3", c_double, "Bin w parameter: EoS for the 4rd bin"),
        ("w4", c_double, "Bin w parameter: EoS for the 5th bin"),
        ("w5", c_double, "Bin w parameter: EoS for the 6th bin"),
        ("w6", c_double, "Bin w parameter: EoS for the 7th bin"),
        ("w7", c_double, "Bin w parameter: EoS for the 8th bin"),
        ("w8", c_double, "Bin w parameter: EoS for the 9th bin"),
        ("w9", c_double, "Bin w parameter: EoS for the 10th bin"),
        # Redshift 
        ("z1", c_double, "Bin w parameter: redshift for the 1st bin"),
        ("z2", c_double, "Bin w parameter: redshift for the 2nd bin"),
        ("z3", c_double, "Bin w parameter: redshift for the 3rd bin"),
        ("z4", c_double, "Bin w parameter: redshift for the 4th bin"),
        ("z5", c_double, "Bin w parameter: redshift for the 5th bin"),
        ("z6", c_double, "Bin w parameter: redshift for the 6th bin"),
        ("z7", c_double, "Bin w parameter: redshift for the 7th bin"),
        ("z8", c_double, "Bin w parameter: redshift for the 8th bin"),
        ("z9", c_double, "Bin w parameter: redshift for the 9th bin"),
        ("z10", c_double, " Bin w parameter: redshift for the 10th bin"),
        ("sigma", c_double, "Transition width"),        
    ]

    def set_params(self, DEmodel, max_num_of_bins, z_knot, w_knot,
                     w0=-1, w1=-1, w2=-1, w3=-1, w4=-1, w5=-1, w6=-1, w7=-1, w8=-1, w9=-1,
                     z1=0.7, z2=1.4, z3=2.1, z4=2.8, z5=3.5, z6=4.2, z7=4.9, z8=5.6, z9=6.3, z10=7.0,
                     sigma=0.1):

        self.DEmodel=DEmodel
        self.max_num_of_bins=max_num_of_bins 
        self.z_knot=z_knot
        self.w_knot=w_knot        
        self.w0=w0 
        self.w1=w1 
        self.w2=w2
        self.w3=w3
        self.w4=w4
        self.w5=w5
        self.w6=w6
        self.w7=w7
        self.w8=w8
        self.w9=w9
        self.z1=z1
        self.z2=z2 
        self.z3=z3
        self.z4=z4 
        self.z5=z5
        self.z6=z6 
        self.z7=z7
        self.z8=z8
        self.z9=z9
        self.z10=z10
        self.sigma=sigma
@fortran_class
class DarkEnergyPPF(LateDE):
    """
    VM: CLASS IMPLEMENTS w, w0wa and binw

    """
    # cannot declare c_Gamma_ppf directly here as have not defined all fields in DarkEnergyEqnOfState (TCubicSpline)
    _fortran_class_module_ = 'DarkEnergyPPF'
    _fortran_class_name_ = 'TDarkEnergyPPF'


# short names for models that support w/wa
F2003Class._class_names.update({'ppf': DarkEnergyPPF})
