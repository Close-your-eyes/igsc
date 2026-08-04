import pandas as pd
from tcrdist.repertoire import TCRrep


def calculate_tcr_distances(
    data,
    chain="TRB",
    species="human",
    db_file="alphabeta_gammadelta_db.tsv",
):
    chain = chain.upper()

    configurations = {
        "TRA": {
            "chain": "alpha",
            "cdr3": "cdr3_a_aa",
            "v": "v_a_gene",
            "j": "j_a_gene",
            "full_attribute": "pw_alpha",
            "cdr3_attribute": "pw_cdr3_a_aa",
        },
        "TRB": {
            "chain": "beta",
            "cdr3": "cdr3_b_aa",
            "v": "v_b_gene",
            "j": "j_b_gene",
            "full_attribute": "pw_beta",
            "cdr3_attribute": "pw_cdr3_b_aa",
        },
    }

    if chain not in configurations:
        raise ValueError("chain must be 'TRA' or 'TRB'")

    if species not in {"human", "mouse"}:
        raise ValueError("species must be 'human' or 'mouse'")

    config = configurations[chain]
    df = pd.DataFrame(data).copy()

    required = [config["cdr3"], config["v"], config["j"]]
    missing = [column for column in required if column not in df.columns]

    if missing:
        raise ValueError(
            f"{chain} requires these missing columns: {', '.join(missing)}"
        )

    df = df.dropna(subset=required).copy()

    df[config["cdr3"]] = (
        df[config["cdr3"]]
        .astype(str)
        .str.strip()
        .str.upper()
    )

    df[config["v"]] = df[config["v"]].astype(str).str.strip()
    df[config["j"]] = df[config["j"]].astype(str).str.strip()

    if "count" not in df.columns:
        df["count"] = 1

    if df.empty:
        raise ValueError(
            "No complete TCR records remain after removing missing values"
        )

    repertoire = TCRrep(
        cell_df=df,
        organism=species,
        chains=[config["chain"]],
        db_file=db_file,
    )

    labels = (
        repertoire.clone_df[config["cdr3"]]
        .astype(str)
        .tolist()
    )

    return {
        "full_distance": getattr(
            repertoire,
            config["full_attribute"],
        ),
        "cdr3_distance": getattr(
            repertoire,
            config["cdr3_attribute"],
        ),
        "labels": labels,
        "clone_data": repertoire.clone_df,
    }
    
    
    
    
import pandas as pd
from tcrdist.repertoire import TCRrep


def calculate_tcr_distances2(
    data,
    chain="TRB",
    species="human",
    db_file="alphabeta_gammadelta_db.tsv",
    cpus=1,
):
    """
    Calculate dense TCR distance matrices.

    This explicitly calls TCRrep.compute_distances(), including for
    repertoires containing more than 10,000 clones.

    Parameters
    ----------
    data : pandas.DataFrame
        TCR repertoire data. An R data.frame supplied through
        reticulate is automatically converted to a pandas DataFrame.
    chain : str
        "TRA" or "TRB".
    species : str
        "human" or "mouse".
    db_file : str
        tcrdist3 database filename or path.
    cpus : int
        Number of CPUs used by tcrdist3.

    Returns
    -------
    dict
        Dense full and CDR3-only distance matrices, labels, and the
        cleaned clone data.
    """
    chain = str(chain).upper()
    species = str(species).lower()
    cpus = int(cpus)

    configurations = {
        "TRA": {
            "tcrdist_chain": "alpha",
            "cdr3_column": "cdr3_a_aa",
            "v_column": "v_a_gene",
            "j_column": "j_a_gene",
            "full_attribute": "pw_alpha",
            "cdr3_attribute": "pw_cdr3_a_aa",
        },
        "TRB": {
            "tcrdist_chain": "beta",
            "cdr3_column": "cdr3_b_aa",
            "v_column": "v_b_gene",
            "j_column": "j_b_gene",
            "full_attribute": "pw_beta",
            "cdr3_attribute": "pw_cdr3_b_aa",
        },
    }

    if chain not in configurations:
        raise ValueError("chain must be 'TRA' or 'TRB'")

    if species not in {"human", "mouse"}:
        raise ValueError("species must be 'human' or 'mouse'")

    if cpus < 1:
        raise ValueError("cpus must be at least 1")

    config = configurations[chain]

    # Defensive copy, including when data originates from R.
    df = pd.DataFrame(data).copy()

    required = [
        config["cdr3_column"],
        config["v_column"],
        config["j_column"],
    ]

    missing = [
        column
        for column in required
        if column not in df.columns
    ]

    if missing:
        raise ValueError(
            f"{chain} requires these missing columns: "
            + ", ".join(missing)
        )

    # Remove incomplete records.
    df = df.dropna(subset=required).copy()

    if df.empty:
        raise ValueError(
            "No complete TCR records remain after removing missing values"
        )

    # Clean sequence and gene columns.
    df[config["cdr3_column"]] = (
        df[config["cdr3_column"]]
        .astype(str)
        .str.strip()
        .str.upper()
    )

    df[config["v_column"]] = (
        df[config["v_column"]]
        .astype(str)
        .str.strip()
    )

    df[config["j_column"]] = (
        df[config["j_column"]]
        .astype(str)
        .str.strip()
    )

    if "count" not in df.columns:
        df["count"] = 1

    # Construct without automatically calculating distances.
    repertoire = TCRrep(
        cell_df=df,
        organism=species,
        chains=[config["tcrdist_chain"]],
        db_file=db_file,
        compute_distances=False,
    )

    repertoire.cpus = cpus

    # Explicit dense calculation. This bypasses the automatic >10,000
    # clone safeguard, but still requires enough system memory.
    repertoire.compute_distances()

    labels = (
        repertoire.clone_df[config["cdr3_column"]]
        .astype(str)
        .tolist()
    )

    full_distance = getattr(
        repertoire,
        config["full_attribute"],
    )

    cdr3_distance = getattr(
        repertoire,
        config["cdr3_attribute"],
    )

    return {
        "full_distance": full_distance,
        "cdr3_distance": cdr3_distance,
        "labels": labels,
        "clone_data": repertoire.clone_df,
        "chain": chain,
        "species": species,
        "cpus": cpus,
    }