configfile: "configs/config.yaml"


def expand_datasets(dataset_configs):
    """
    Expand dataset configurations into individual dataset instances.
    Handles both single instance and range specifications.

    Returns a dict mapping dataset_id to dataset parameters.
    """
    expanded = {}

    for ds in dataset_configs:
        endo_genome = ds.get("endogenous_genome")
        cont_genome = ds.get("contaminant_genome")
        bact_db = ds.get("bacterial_db")
        proportions = ds.get("proportions")
        if not (len(proportions) == 3):
            raise ValueError(f"{proportions} should be 3 float values with sum 1")

        endo_vars = ds.get("endogenous_variation")
        if not isinstance(endo_vars, list):
            endo_vars = [endo_vars]

        read_lengths = ds.get("read_length")
        if not isinstance(read_lengths, list):
            read_lengths = [read_lengths]

        for endo_var, read_len in product(endo_vars, read_lengths):
            dataset_id = f"{endo_genome}_{cont_genome}_{bact_db}_{endo_var}_len{read_len}_0p7"

            expanded[dataset_id] = {
                "endogenous_genome": endo_genome,
                "contaminant_genome": cont_genome, 
                "endogenous_variation": endo_var,
                "bacterial_db": bact_db,
                "read_length": read_len,
                "proportions": proportions,
                "endogenous_prop": proportions[0],
                "contaminant_prop": proportions[1],
                "bacterial_prop": proportions[2],
            }

    return expanded


def get_dataset_param(wildcards, param_name):
    """Retrieve a parameter value for a given dataset_id."""
    dataset_id = wildcards.dataset_id
    if dataset_id in DATASET_CONFIGS:
        return DATASET_CONFIGS[dataset_id][param_name]
    raise ValueError(f"Dataset {dataset_id} not found in DATASET_CONFIGS")


# Expand all datasets from config
DATASET_CONFIGS = expand_datasets(config.get("datasets", []))
DATASET_IDS = list(DATASET_CONFIGS.keys())