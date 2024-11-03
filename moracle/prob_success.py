import pandas as pd

stage_success_probabilities = {
    "-1.0": 0,
    "0.5": 0.02,
    "1.0": 0.08,
    "2.0": 0.10,
    "3.0": 0.17,
    "4.0": 0.50
}

def compute_prob_clin_success(df: pd.DataFrame, stage_probs: dict[str, float] = stage_success_probabilities):
    weighted_probs = []
    total_weight = 0

    for i in range(df.shape[0]):
        stage = df.iloc[i]["Trial Phase"]
        weight = df.iloc[i]["Sim. Score"]  # Default to 1 if no score
        prob = stage_probs[str(stage)]
        weighted_probs.append(prob * weight)
        total_weight += weight

    return sum(weighted_probs) / total_weight if total_weight > 0 else 0