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

    train_eaa = pd.DataFrame()
    test_eaa = pd.DataFrame()

    train_eaa["sample_id"] = train_metadata["sample_id"]
    test_eaa["sample_id"] = test_metadata["sample_id"]

    for clock, eaa_col in CLOCKS.items():
        X_train = sm.add_constant(train_metadata["age"])
        y_train = train_metadata[clock]

        model = sm.OLS(y_train, X_train).fit()

        train_eaa[eaa_col] = model.resid

        X_test = sm.add_constant(test_metadata["age"])
        test_eaa[eaa_col] = test_metadata[clock] - model.predict(X_test)

    train_metadata = train_metadata.drop(columns=list(CLOCKS.keys()))
    test_metadata = test_metadata.drop(columns=list(CLOCKS.keys()))

    return train_metadata, test_metadata, train_eaa, test_eaa