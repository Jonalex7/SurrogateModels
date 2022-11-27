import numpy as np
from datetime import datetime
from os import path, makedirs
from sklearn.gaussian_process import GaussianProcessRegressor # is it needed?
from sklearn.gaussian_process.kernels import RBF # is it needed?

class PC_Kriging():
    def __init__(self, config=None):
        if config is None:
            config = {"pol_type": ['hermite', 'legendre'],
                      "order": [2, 1]}
        assert "pol_type" in config and \
               "order" in config and \
            "Missing config"

        self.date_record = datetime.now().strftime("%Y_%m_%d_%H%M%S")

    def fit(self, experiments):
        pass

    def predict(self, points):
        pass