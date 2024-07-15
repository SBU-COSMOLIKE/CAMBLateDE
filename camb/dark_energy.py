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

# JVR Modification: Implementing the PC model
# PC is just a mapping between the alphas and the w_i
# So just use the same model binw and set alpha_i instead of w_i
PCs = np.array([
    [0.97597922, 0.21318807, 0.04285546, 0.01218697, 0.00550414],
    [-0.21443583,  0.90363374,  0.35153986,  0.10769025,  0.04781089],
    [ 0.0380459 , -0.36577559,  0.81290758,  0.39122488,  0.22557723],
    [-0.00584294,  0.06483369, -0.45591547,  0.67855926,  0.57224196],
    [-4.55250983e-04,  1.31278706e-03,  7.68471741e-02, -6.12172124e-01, 7.86980223e-01]
])
PCs_10 = np.array([
    [0.69386562, 0.59814274, 0.33365475, 0.18008562, 0.10753549, 0.05628191,
    0.03782315, 0.02204437, 0.01732748, 0.00843283],
    [-0.60436776,  0.20618268,  0.49883458,  0.4430748,   0.29462262,  0.18247331,
    0.12416477,  0.07917829,  0.06338614,  0.03569602],
    [-0.34508494,  0.55774976,  0.14261381, -0.34185917, -0.38398621, -0.35246186,
    -0.28775366, -0.20919111, -0.16139707, -0.09131441],
    [ 0.14131458, -0.35304698,  0.23454367,  0.39952541,  0.13953395, -0.35507667,
    -0.40626604, -0.41568027, -0.33053273, -0.2193479 ],
    [-0.10131702,  0.32668632, -0.53135655,  0.00557896,  0.45977156,  0.23657403,
    0.12708682, -0.40878027, -0.31738624, -0.22202804],
    [ 0.05687853, -0.20357395,  0.40383485, -0.23399916, -0.32981993,  0.4101677,
    0.44553287, -0.40646979, -0.23764009, -0.2006876 ],
    [-0.00390202 , 0.00949467, -0.02064439,  0.0625282,  -0.0933946,  -0.00548024,
    0.06822735,  0.65336832, -0.5772576,  -0.47120548],
    [ 0.02645227, -0.1207529,   0.3307611,  -0.63354939,  0.56320007,  0.13487063,
    -0.3128708,   0.13016026, -0.13090139,  0.08112874],
    [-1.75930258e-06, -1.01668057e-02,  4.09255186e-02, -1.00118467e-01,
    9.91660724e-02, -3.19610064e-02, -4.47020220e-02, -1.54769989e-02,
    5.93318124e-01, -7.89353159e-01],
    [ 0.00462894, -0.03234056,  0.0869642,  -0.15703697,  0.28070739, -0.68540542,
    0.645063,   -0.01742014, -0.02804128,  0.03059421]
])

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
        ("C_0", c_double, "First Chebyshev term"),        
        ("C_1", c_double, "Second Chebyshev term"),        
        ("C_2", c_double, "Third Chebyshev term"),        
        ("C_3", c_double, "Fourth Chebyshev term"),        
    ]

    def set_params(self, DEmodel, max_num_of_bins=0, z_knot=[], w_knot=[],
                     w0=-1, w1=-1, w2=-1, w3=-1, w4=-1, w5=-1, w6=-1, w7=-1, w8=-1, w9=-1,
                     z1=0.7, z2=1.4, z3=2.1, z4=2.8, z5=3.5, z6=4.2, z7=4.9, z8=5.6, z9=6.3, z10=7.0,
                     sigma=0.1, C_0=0, C_1=0, C_2=0, C_3=0, alpha1=None, alpha2=None, alpha3=None, alpha4=None, alpha5=None,
                     alpha6=None, alpha7=None, alpha8=None, alpha9=None, alpha10=None):

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
        self.C_0=C_0
        self.C_1=C_1
        self.C_2=C_2
        self.C_3=C_3
        if DEmodel == 5 and alpha1 is not None:
            # Override with PC w
            alphas = np.array([alpha1, alpha2, alpha3, alpha4, alpha5])
            ws = -1*np.ones(5) + np.dot(PCs.T, alphas)
            self.w0 = ws[0]
            self.w1 = ws[1]
            self.w2 = ws[2]
            self.w3 = ws[3]
            self.w4 = ws[4]
        elif DEmodel == 6 and alpha1 is not None:
            # Override with PC w
            alphas = np.array([alpha1, alpha2, alpha3, alpha4, alpha5, alpha6, alpha7, alpha8, alpha9, alpha10])
            ws = -1*np.ones(10) + np.dot(PCs.T, alphas)
            self.w0 = ws[0]
            self.w1 = ws[1]
            self.w2 = ws[2]
            self.w3 = ws[3]
            self.w4 = ws[4]
            self.w5 = ws[5]
            self.w6 = ws[6]
            self.w7 = ws[7]
            self.w8 = ws[8]
            self.w9 = ws[9]

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
