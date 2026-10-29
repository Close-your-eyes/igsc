# python -m pip install umap-learn pandas
import argparse
import pandas as pd
from umap import UMAP

parser = argparse.ArgumentParser()
parser.add_argument("input", help="Input CSV with numeric feature columns")
parser.add_argument("output", help="Output CSV")
parser.add_argument("--neighbors", type=int, default=15)
parser.add_argument("--min-dist", type=float, default=0.1)
parser.add_argument("--metric", default="euclidean")
args = parser.parse_args()

data = pd.read_csv(args.input)

embedding = UMAP(
    n_components=2,
    n_neighbors=args.neighbors,
    min_dist=args.min_dist,
    metric=args.metric,
    random_state=42,
).fit_transform(data)

pd.DataFrame(embedding, columns=["UMAP1", "UMAP2"]).to_csv(
    args.output, index=False
)
print(f"Saved {len(embedding)} rows to {args.output}")