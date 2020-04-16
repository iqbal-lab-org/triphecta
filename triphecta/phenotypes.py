import csv
import json
import logging

from triphecta import utils

data_lookup = {
    "NA": None,
    ".": None,
    "": None,
    "N/A": None,
    "NULL": None,
    "NONE": None,
    "T": True,
    "TRUE": True,
    "F": False,
    "FALSE": False,
    "R": True,
    "S": False,
    "RESISTANT": True,
    "SUSCEPTIBLE": False,
}


class Phenotypes:
    def __init__(self, input_phenos_tsv):
        self.input_phenos_tsv = input_phenos_tsv
        self.phenos, self.pheno_types = Phenotypes._load_phenotypes_tsv_file(
            self.input_phenos_tsv
        )
        self.pheno_types, self.bool_pheno_types = Phenotypes._get_pheno_types(
            self.pheno_types
        )

    @classmethod
    def convert_one_variable_string(cls, s):
        try:
            return data_lookup[s.upper().strip()]
        except KeyError:
            return float(s)
        except ValueError:
            raise TypeError(
                f"Error inferring type of '{s}' when loading phenotype file"
            )

    @classmethod
    def _load_phenotypes_tsv_file(cls, infile):
        phenos = {}

        with utils.open_file(infile) as f:
            reader = csv.DictReader(f, delimiter="\t")
            if "sample" not in reader.fieldnames:
                raise RuntimeError(
                    f"Must have a 'sample' column in phenotypes file. Not found in file {infile}"
                )
            pheno_types = {x: set() for x in reader.fieldnames if x != "sample"}

            for row in reader:
                if row["sample"] in phenos:
                    raise RuntimeError(
                        f"Duplicate sample name '{row['sample']}' in phenotypes file {infile}"
                    )

                phenos[row["sample"]] = {
                    x: Phenotypes.convert_one_variable_string(row[x])
                    for x in row
                    if x != "sample"
                }
                for p in pheno_types:
                    pheno_types[p].add(type(phenos[row["sample"]][p]))

        return phenos, pheno_types

    @classmethod
    def _get_pheno_types(cls, types):
        errors = []
        infer_type = {}
        bools = set()

        for pheno, types_set in types.items():
            if float in types_set and bool in types_set:
                errors.append(f"Mixed types bool and float in phenotype '{pheno}'")
            elif types_set == {type(None)}:
                logging.warning(f"WARNING: all values are null for phenotype '{pheno}'")
                infer_type[pheno] = type(None)
            elif float in types_set:
                infer_type[pheno] = float
            else:
                infer_type[pheno] = bool
                bools.add(pheno)

        if len(errors):
            raise RuntimeError("Problems with phenotype data:\n" + "\n".join(errors))

        assert len(infer_type) == len(types)
        return infer_type, bools

    def write_template_constraints_json(self, outfile):
        constraints = {}
        for pheno, pheno_type in self.pheno_types.items():
            constraints[pheno] = {"must_be_same": True, "params": {}}
            if pheno_type == bool:
                constraints[pheno]["method"] = "equal"
            elif pheno_type == float:
                constraints[pheno]["method"] = "range"
                constraints[pheno]["params"] = {"low": 0, "high": 1}
            else:
                raise TypeError

        with utils.open_file(outfile, "w") as f:
            json.dump(constraints, f, sort_keys=True, indent=2)

    def __getitem__(self, sample):
        return self.phenos[sample]
