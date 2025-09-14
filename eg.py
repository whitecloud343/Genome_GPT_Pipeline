import pandas as pd
import ast
import matplotlib.pyplot as plt

# Load harmonized dataset
df = pd.read_parquet("genome_phoenix.parquet")


df.to_csv("genome_harmonized1.csv", index=False)
