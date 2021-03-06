import collections
import os

from triphecta import distances, utils, variant_counts


class Genotypes:
    def __init__(
        self,
        file_of_vcf_filenames=None,
        distance_matrix_file=None,
        variant_counts_file=None,
        check_vcf_files_exist=True,
        testing=False,
    ):
        self.distance_matrix_file = (
            None
            if distance_matrix_file is None
            else os.path.abspath(distance_matrix_file)
        )
        self.file_of_vcf_filenames = (
            None
            if file_of_vcf_filenames is None
            else os.path.abspath(file_of_vcf_filenames)
        )
        self.variant_counts_file = (
            None
            if variant_counts_file is None
            else os.path.abspath(variant_counts_file)
        )
        self.check_vcf_files_exist = check_vcf_files_exist

        if testing:
            self.sample_names_list = []
            self.distances = {}
            self.vcf_variant_counts = []
            self.vcf_files = {}
        else:
            self.load_all_data()

        self._make_sample_name_to_index()
        self.excluded_samples = dict()

    def _make_sample_name_to_index(self):
        self.sample_name_to_index = {
            name: i for i, name in enumerate(self.sample_names_list)
        }

    def load_all_data(self):
        if self.file_of_vcf_filenames is None:
            raise RuntimeError("Must provide file_of_vcf_filenames")
        else:
            self.vcf_files = utils.load_file_of_vcf_filenames(
                self.file_of_vcf_filenames,
                check_vcf_files_exist=self.check_vcf_files_exist,
            )

        if self.distance_matrix_file is None:
            raise RuntimeError("Must provide distance matrix file")
        else:
            (
                self.sample_names_list,
                self.distances,
            ) = distances.load_distance_matrix_file(self.distance_matrix_file)

        if self.variant_counts_file is not None:
            self.vcf_variant_counts = variant_counts.load_variant_count_list_from_tsv(
                self.variant_counts_file
            )

    def distance(self, sample1, sample2):
        key = tuple(
            sorted(
                [self.sample_name_to_index[sample1], self.sample_name_to_index[sample2]]
            )
        )
        return self.distances[key]

    def sample_names(self):
        for sample in self.sample_names_list:
            yield sample

    def distance_dict(self, sample, top_n=None):
        all_distances = {
            other: self.distance(sample, other)
            for other in self.sample_names()
            if sample != other and other not in self.excluded_samples
        }
        if top_n is None:
            return all_distances
        else:
            value_counts = collections.Counter(all_distances.values())
            total = 0
            for max_value, count in sorted(value_counts.items()):
                total += count
                if total >= top_n:
                    break

            return {k: v for k, v in all_distances.items() if v <= max_value}

    def update_excluded_samples_using_variant_counts(
        self, minimum_percent_hom_calls=90.0, count_het_to_hom_as_hom=True
    ):
        for i, counts in enumerate(self.vcf_variant_counts):
            total_calls = sum(counts)
            hom_calls = counts.hom + (
                counts.het_to_hom if count_het_to_hom_as_hom else 0
            )
            if 100 * hom_calls / total_calls < minimum_percent_hom_calls:
                if i not in self.excluded_samples:
                    self.excluded_samples[i] = set()
                self.excluded_samples[i].add("Too few hom calls")
