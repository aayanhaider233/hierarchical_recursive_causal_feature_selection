import pandas as pd
import statsmodels.api as sm

CLOCKS = {
    "horvath2013": "horvath2013_EAA",
    "hannum": "hannum_EAA",
    "grimage2": "grimage2_EAA",
}

def compute_eaa(train_metadata, test_metadata):
    train_metadata = train_metadata.copy()
    test_metadata = test_metadata.copy()
    for clock, eaa_col in CLOCKS.items():
        X_train = sm.add_constant(train_metadata["age"])
        y_train = train_metadata[clock]

        model = sm.OLS(y_train, X_train).fit()

        train_metadata[eaa_col] = model.resid
        train_metadata.drop(columns=[clock], inplace=True)
        X_test = sm.add_constant(test_metadata["age"])
        test_metadata[eaa_col] = test_metadata[clock] - model.predict(X_test)
        test_metadata.drop(columns=[clock], inplace=True)
    return train_metadata, test_metadata