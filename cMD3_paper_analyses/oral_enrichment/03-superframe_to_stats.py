#!/usr/bin/env python

import pandas as pd
from scipy import stats
from sklearn.metrics import roc_auc_score
import argparse

def compute_wilcoxon_auc(df, score_col, dataset_col, col_group, positive_class="positive", negative_class="negative"):
    results = []
    datasets = df[dataset_col].unique()
    for ds in datasets:
        ds_df = df[df[dataset_col] == ds]
        try:
            positives = ds_df[ds_df[col_group] == positive_class][score_col]
            negatives = ds_df[ds_df[col_group] == negative_class][score_col]

            if len(positives) == 0 or len(negatives) == 0:
                continue

            # Wilcoxon rank sum test (Mann-Whitney)
            wilco_stat, wilco_p = stats.ranksums(positives, negatives)

            # AUC score - binary labels: positive=1, negative=0
            labels = ds_df[col_group].apply(lambda x: 1 if x == positive_class else 0)
            scores = ds_df[score_col]
            auc = roc_auc_score(labels, scores)

            results.append({
                dataset_col: ds,
                'wilcoxon_statistic': wilco_stat,
                'wilcoxon_pvalue': wilco_p,
                'auc_score': auc
            })
        except Exception as e:
            # Optionally print error for debugging
            # print(f"Error processing dataset {ds}: {e}")
            continue

    return pd.DataFrame(results)

def main():
    parser = argparse.ArgumentParser(description="Compute Wilcoxon and AUC per dataset")
    parser.add_argument('input_file', help='Input TSV file with data')
    parser.add_argument('output_file', help='Output TSV file for results')
    parser.add_argument('--score_col', default='Oral_Richness', help='Column name for score (default: Oral_Richness)')
    parser.add_argument('--dataset_col', default='DATASET', help='Column name for dataset (default: DATASET)')
    parser.add_argument('--group_col', default='CLASS', help='Column name for group labels (default: CLASS)')
    parser.add_argument('--positive_class', default='positive', help='Positive class label (default: positive)')
    parser.add_argument('--negative_class', default='negative', help='Negative class label (default: negative)')

    args = parser.parse_args()

    df = pd.read_csv(args.input_file, sep='\t', header=0, index_col=0, low_memory=False)

    results = compute_wilcoxon_auc(
        df,
        score_col=args.score_col,
        dataset_col=args.dataset_col,
        col_group=args.group_col,
        positive_class=args.positive_class,
        negative_class=args.negative_class
    )

    results.to_csv(args.output_file, sep='\t', index=False)

if __name__ == "__main__":
    main()
