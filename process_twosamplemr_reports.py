import pandas as pd
from scipy.stats import false_discovery_control as fdc
from pathlib import Path
from io import StringIO


def parse_twosamplemr_reports(reports):

    for report in reports:
        with open(report) as rep_file:
            metrics_table = rep_file.readlines()[15:21]

        table = (
            pd.read_csv(
                StringIO("".join(metrics_table)),
                sep="|",
                header=0,
                index_col=1,
                skipinitialspace=True,
            )
            .dropna(axis=1, how="all")
            .iloc[1:]
        )

        table["padj"] = fdc(table["pval"].astype(float), method="bh")

        report_dir = report.parent.name

        table.to_csv(f"results/tables/{report_dir}.csv")


if __name__ == "__main__":
    reports = Path("./results/reports/").glob("**/*md")

    parse_twosamplemr_reports(reports)
