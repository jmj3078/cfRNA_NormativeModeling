"""Long-format result tables shared by every evaluation. A (method, scenario) pair is
replaced as a whole on rewrite, so reruns never duplicate rows."""
import pandas as pd


def append(path, df):
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        old = pd.read_csv(path)
        keys = set(zip(df["method"], df["scenario"]))
        old = old[[k not in keys for k in zip(old["method"], old["scenario"])]]
        df = pd.concat([old, df], ignore_index=True)
    df.to_csv(path, index=False)
