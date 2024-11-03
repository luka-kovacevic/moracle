import pandas as pd

stage_success_probabilities = {
    "Preclinical": 0.05,
    "-1.0": 0,
    "1.0": 0.08,
    "2.0": 0.10,
    "3.0": 0.17,
    "4.0": 0.50
}

def compute_prob_clin_success(df: pd.DataFrame, stage_probs: dict[str, float] = stage_success_probabilities):
    # df.loc[:, "Trial Phase"] = df.loc[:, "Trial Phase"].map({'nan': "Preclinical", 1.0: "Stage 1", 2.0: "Stage 2", 3.0: "Stage 3", 4.0: "Approved"})
    

    weighted_probs = []
    total_weight = 0

    if df.shape[0] > 0:
        for i in range(df.shape[0]):
            stage = df.iloc[i]["Trial Phase"]
            weight = df.iloc[i]["Sim. Score"]  # Default to 1 if no score
            prob = stage_probs[str(stage)]
            weighted_probs.append(prob * weight)
            total_weight += weight

    return sum(weighted_probs) / total_weight if total_weight > 0 else 0