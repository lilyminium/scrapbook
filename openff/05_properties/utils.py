

from openff.evaluator.datasets import PhysicalProperty


def dataset_from_csv(csv: str = "01_mnsol_data.csv"):
    import pandas as pd
    from openff.evaluator.datasets import PhysicalPropertyDataSet

    df = pd.read_csv(csv)
    dataset = PhysicalPropertyDataSet.from_pandas(df)
    return dataset


def property_from_csv(csv: str = "01_mnsol_data.csv",
                      index: int = 0):
    dataset = dataset_from_csv(csv)
    return dataset.properties[index]
