class PhenotypeCompare:
    def __init__(self, constraints, count_unknown_as_diff=True):
        self.compare_functions = {
            "equal": PhenotypeCompare._compare_method_equal,
            "range": PhenotypeCompare._compare_method_range,
            "abs_distance": PhenotypeCompare._compare_method_abs_distance,
        }
        self.constraints = constraints
        errors = self._sanity_check_constraints()
        if len(errors):
            raise RuntimeError("Errors in constraints:\n" + "\n".join(errors))
        self.count_unknown_as_diff = count_unknown_as_diff
        self.required_diff_keys = {
            k for k in self.constraints if not self.constraints[k]["must_be_same"]
        }

    def _sanity_check_constraints(self):
        found_must_be_same = False
        errors = []

        for d in self.constraints.values():
            if d["must_be_same"]:
                found_must_be_same = True

            if d["method"] not in self.compare_functions:
                errors.append(f"Unknown method {d}")
                continue

            if d["method"] == "equal":
                if "params" not in d:
                    d["params"] = {}
                elif len(d["params"]) > 0:
                    errors.append(f"method is 'equal', params supplied: {d}")
            if d["method"] == "range" and (
                "low" not in d["params"] or "high" not in d["params"]
            ):
                errors.append(f"method is 'range', low and high not supplied: {d}")
            if d["method"] == "abs_distance" and ("max_dist" not in d["params"]):
                errors.append(f"method is 'abs_distance', max_dist not supplies: {d}")

        if not found_must_be_same:
            errors.append(
                f"Must have at least one constraint where 'must_be_same' is True"
            )

        return errors

    def satisfy_required_differences(self, pheno1, pheno2):
        for key in self.required_diff_keys:
            if PhenotypeCompare._phenos_equal_account_for_none(
                pheno1[key],
                pheno2[key],
                self.compare_functions[self.constraints[key]["method"]],
                False,
                **self.constraints[key]["params"],
            ):
                return False
        return True

    @staticmethod
    def _compare_method_equal(p1, p2):
        return p1 == p2

    @staticmethod
    def _compare_method_range(p1, p2, low=None, high=None):
        return (low <= p1 <= high) == (low <= p2 <= high)

    @staticmethod
    def _compare_method_abs_distance(p1, p2, max_dist=None):
        return abs(p1 - p2) <= max_dist

    @classmethod
    def _phenos_equal_account_for_none(
        cls, p1, p2, compare_function, count_unknown_as_diff, **kwargs
    ):
        if p1 is None or p2 is None:
            return not count_unknown_as_diff
        else:
            return compare_function(p1, p2, **kwargs)

    def phenos_agree_on_one_feature(self, pheno1, pheno2, key):
        return PhenotypeCompare._phenos_equal_account_for_none(
            pheno1[key],
            pheno2[key],
            self.compare_functions[self.constraints[key]["method"]],
            self.count_unknown_as_diff,
            **self.constraints[key]["params"],
        )

    def phenos_agree_on_features(self, pheno1, pheno2, keys):
        for key in keys:
            if not self.phenos_agree_on_one_feature(pheno1, pheno2, key):
                return False
        return True

    def differences(self, pheno1, pheno2):
        """Returns number of differences between the two phenotypes.
        Assumes that satisfy_required_differences(pheno1, pheno2) is True.
        (Or at least doesn't care if it's True or False.)
        Counts the differences only from he constraints where
        'must_be_same' is True"""
        differences = 0
        for key, constraint in self.constraints.items():
            if constraint["must_be_same"] is False:
                continue
            elif not self.phenos_agree_on_one_feature(pheno1, pheno2, key):
                differences += 1

        return differences
