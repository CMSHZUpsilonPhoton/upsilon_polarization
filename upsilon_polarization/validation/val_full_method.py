import numpy as np
import matplotlib.pyplot as plt

from upsilon_polarization import full_method

if __name__ == "__main__":
    cms_data = full_method.load_CMS_data()

    print(f"Weight: {full_method.get_weight(1, 1, 1, 1, 1, 1, 1, 1, 1, cms_data, 1, 1, 0)}")

