import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from sklearn.experimental import enable_iterative_imputer
from sklearn.impute import IterativeImputer
from sklearn.linear_model import BayesianRidge
from sklearn.gaussian_process import GaussianProcessRegressor, kernels
from sklearn.metrics import r2_score, mean_squared_error
from sklearn.model_selection import train_test_split

class polymer_data():
    def __init__(self, df=None):
        self.df = df

    def log_cal(self):
        gases = ['He', 'H2', 'N2', 'O2', 'CO2', 'CH4']
        for gas in gases:
            self.df['log10_'+gas] = self.df[gas].apply(lambda x: np.log10(x) if not pd.isna(x) else np.nan)
        return self.df

    def imputator(self, max_iter=200, random_state=0):
        train_data = self.df
        log_gases = ['log10_He', 'log10_H2', 'log10_N2', 'log10_O2', 'log10_CO2', 'log10_CH4']
        model =  IterativeImputer(max_iter=max_iter, random_state=random_state, sample_posterior=True)
        log_perms = train_data[log_gases].to_numpy()
        log_full = model.fit_transform(log_perms)
        new_cols = [log_gas+'_'+'Bayes' for log_gas in log_gases]
        new_df = pd.DataFrame(data=log_full, columns=new_cols)
        r_df = pd.concat([train_data, new_df], axis=1)
        
        return r_df

    def get_FP(self):

        self.df['FP'] = self.df.SMILES.apply(lambda m: AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(m), radius=2, nBits=2048))
        fp = np.array(self.df['FP'].tolist())
        
        return fp
    

def GPR_model(X_train, y_train, length_scale, noise_level, a, b):

    kernel = a*kernels.RBF(length_scale=length_scale
                                  ) + b*kernels.WhiteKernel(
                                      noise_level=noise_level)
    model= GaussianProcessRegressor(kernel=kernel, n_restarts_optimizer=0, normalize_y=True).fit(X_train, y_train)
    return model




# def log_cal(df):
#     gases = ['He', 'H2', 'N2', 'O2', 'CO2', 'CH4']
#     for gas in gases:
#         df['log10_'+gas] = df[gas].apply(lambda x: np.log10(x) if not pd.isna(x) else np.nan)
#     return df

# def imputator(df, max_iter=200, random_state=0):
#     train_data = df
#     log_gases = ['log10_He', 'log10_H2', 'log10_N2', 'log10_O2', 'log10_CO2', 'log10_CH4']
#     model =  IterativeImputer(max_iter=max_iter, random_state=random_state, sample_posterior=True)
#     log_perms = train_data[log_gases].to_numpy()
#     log_full = model.fit_transform(log_perms)
#     new_cols = [log_gas+'_'+'Bayes' for log_gas in log_gases]
#     new_df = pd.DataFrame(data=log_full, columns=new_cols)
#     r_df = pd.concat([train_data, new_df], axis=1)
    
#     return r_df

# def get_FP(data):

#     data['ECFP'] = data.SMILES.apply(lambda m: AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(m), radius=2, nBits=2048))
#     fp = np.array(data['ECFP'].tolist())
    
#     return fp