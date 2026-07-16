import numpy as np
import pandas as pd

def beta_to_m(concatenated_df, sample_id_col="sample_id", epsilon=1e-6):

    sample_ids = concatenated_df[[sample_id_col]]
    beta_values = concatenated_df.drop(columns=[sample_id_col]).astype(float)

    beta_values = beta_values.clip(epsilon, 1 - epsilon)

    m_values = np.log2(beta_values / (1 - beta_values))

    m_df = pd.concat([sample_ids, m_values], axis=1)

    return m_df