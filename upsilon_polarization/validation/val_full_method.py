import numpy as np

from upsilon_polarization import full_method

if __name__ == "__main__":
    cms_data = full_method.load_CMS_data()

    print(
        f"Dummy Weight (Scalars - NOMINAL): {full_method.get_weight(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)}"
    )

    a = np.arange(6)
    print(
        f"Dummy Weight (NOMINAL): {full_method.get_weight(a, a, a, a, a, a, a, a, a, a, a)}"
    )

    print(
        f"Dummy Weight (PLUS): {full_method.get_weight_plus(a, a, a, a, a, a, a, a, a, a, a)}"
    )

    print(
        f"Dummy Weight (MINUS): {full_method.get_weight_minus(a, a, a, a, a, a, a, a, a, a, a)}"
    )

    N_SAMPLES = 6000
    wgts = full_method.get_weight(
                                    np.random.uniform(1,1000,N_SAMPLES), 
                                    np.random.uniform(-2.5,2.5,N_SAMPLES), 
                                    np.random.uniform(-3,3,N_SAMPLES), 
                                    np.random.uniform(1,1000,N_SAMPLES), 
                                    np.random.uniform(-2.5,2.5,N_SAMPLES), 
                                    np.random.uniform(-3,3,N_SAMPLES),
                                    np.random.uniform(1,1000,N_SAMPLES), 
                                    np.random.uniform(-2.5,2.5,N_SAMPLES), 
                                    np.random.uniform(-3,3,N_SAMPLES), 
                                    np.random.uniform(-8,8,N_SAMPLES), 
                                    np.random.uniform(1,500,N_SAMPLES), 
        )
    print(
        f"""Dummy Weight (Large sample - NOMINAL) min and max: {wgts.min()}, {wgts.max()}"""
    )
