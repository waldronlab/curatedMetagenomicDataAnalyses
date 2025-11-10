#!/usr/bin/env python

import pandas as pd
import numpy as np
from scipy import stats as sts
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
import seaborn as sns

sns.set_style("whitegrid")
sns.set_context(
    "paper",
    rc={
        "font.size": 12,
        "axes.titlesize": 12,
        "axes.labelsize": 12,
        "axes.ticklabelsize": 12,
        "xtick.labelsize": 12,
        "ytick.labelsize": 12,
    },
)

if __name__ == '__main__':
    keys = [
        'sex clr', 'sex arcsin', 'age clr', 'age arcsin', 'bmi clr', 'bmi arcsin'
    ]
    tabs = [
        'meta_analyses/clr_sex_species.tsv',
        'meta_analyses/arcsin_sex_species.tsv',
        'meta_analyses/clr_age_species.tsv',
        'meta_analyses/arcsin_age_species.tsv',
        'meta_analyses/clr_bmi_species.tsv',
        'meta_analyses/arcsin_bmi_species.tsv'
    ]

    cols_cat = ['RE_Effect']
    cols_reg = ['RE_Correlation']

    fig, axes = plt.subplots(1, 3, figsize=(16, 5))  # 1 row, 3 columns

    for idx, (col, col_set, ax) in enumerate(
        zip(['sex', 'age', 'bmi'], [cols_cat, cols_reg, cols_reg], axes)
    ):
        frame = pd.read_csv(
            f'meta_analyses/clr_{col}_species.tsv',
            sep='\t',
            header=0,
            index_col=0
        )[[col_set[0]]]
        if col == "age":
            frame.drop(
                'k__Bacteria|p__Actinobacteria|c__Actinobacteria|o__Actinomycetales|f__Actinomycetaceae|g__Actinomyces|s__Actinomyces_massiliensis',
                inplace=True,
                errors='ignore'
            )
        frame.dropna(subset=[col_set[0]], inplace=True)
        frame = frame.rename(columns={col_set[0]: f'clr_{col_set[0]}'})

        frame_arc = pd.read_csv(
            f'meta_analyses/arcsin_{col}_species.tsv',
            sep='\t',
            header=0,
            index_col=0
        )[[col_set[0]]]
        frame_arc.dropna(subset=[col_set[0]], inplace=True)
        frame_arc = frame_arc.rename(columns={col_set[0]: f'arc_{col_set[0]}'})
        frame_arc = frame_arc.loc[frame.index]

        frame = pd.concat([frame, frame_arc], axis=1)

        sns.regplot(
            data=frame,
            x=f'clr_{col_set[0]}',
            y=f'arc_{col_set[0]}',
            ax=ax,
            fit_reg=True,
            scatter=True
        )

        rho, rhop = sts.spearmanr(frame[f'clr_{col_set[0]}'], frame[f'arc_{col_set[0]}'])
        rho_x, rhop_x = sts.pearsonr(frame[f'clr_{col_set[0]}'], frame[f'arc_{col_set[0]}'])

        ax.set_title(f"{col} r={rho:.2f}({rhop:.2e}) / {rho_x:.2f}({rhop_x:.2e})")
        ax.set_xlabel("CLR-derived effect size")
        ax.set_ylabel("Arcsin sqrt-derived effect size")

    plt.tight_layout()
    [plt.savefig(f'images/scatters_of_effSizes_clr_vs_asin.{fmt}', dpi=200) for fmt in ["svg", "png"]]
    plt.close()
