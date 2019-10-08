"""
From https://stackoverflow.com/questions/52453426/grouping-by-multiple-dimensions/52612073#52612073
"""

import itertools
from collections import defaultdict

import numpy as np
import xarray as xr
from xarray import DataArray

class DataAssembly(DataArray):
    def multi_dim_groupby(self, groups, apply):
        # align
        groups = sorted(groups, key=lambda group: self.dims.index(self[group].dims[0]))
        # build indices
        groups = {group: np.unique(self[group]) for group in groups}
        group_dims = {self[group].dims: group for group in groups}
        indices = defaultdict(lambda: defaultdict(list))
        result_indices = defaultdict(dict)
        for group in groups:
            for index, value in enumerate(self[group].values):
                indices[group][value].append(index)
                if value not in result_indices[group]:  # if captured once, it will be "grouped away"
                    index = max(result_indices[group].values()) + 1 if len(result_indices[group]) > 0 else 0
                    result_indices[group][value] = index

        coords = {coord: (dims, value) for coord, dims, value in walk_coords(self)}

        def simplify(value):
            return value.item() if value.size == 1 else value

        def indexify(dict_indices):
            return [(i,) if isinstance(i, int) else tuple(i) for i in dict_indices.values()]

        # group and apply
        # making this a DataArray right away and then inserting through .loc would slow things down
        result = np.zeros([len(indices) for indices in result_indices.values()])
        result_coords = {coord: (dims, [None] * len(result_indices[group_dims[dims]]))
                         for coord, (dims, value) in coords.items()}
        for values in itertools.product(*groups.values()):
            group_values = dict(zip(groups.keys(), values))
            self_indices = {group: indices[group][value] for group, value in group_values.items()}
            values_indices = indexify(self_indices)
            cells = self.values[values_indices]  # using DataArray would slow things down. thus we pass coords as kwargs
            cells = simplify(cells)
            cell_coords = {coord: (dims, value[self_indices[group_dims[dims]]])
                           for coord, (dims, value) in coords.items()}
            cell_coords = {coord: (dims, simplify(np.unique(value))) for coord, (dims, value) in cell_coords.items()}

            # ignore dims when passing to function
            passed_coords = {coord: value for coord, (dims, value) in cell_coords.items()}
            merge = apply(cells, **passed_coords)
            result_idx = {group: result_indices[group][value] for group, value in group_values.items()}
            result[indexify(result_idx)] = merge
            for coord, (dims, value) in cell_coords.items():
                if isinstance(value, np.ndarray):  # multiple values for coord -> ignore
                    if coord in result_coords:  # delete from result coords if not yet deleted
                        del result_coords[coord]
                    continue
                assert dims == result_coords[coord][0]
                coord_index = result_idx[group_dims[dims]]
                result_coords[coord][1][coord_index] = value

        # re-package
        result = type(self)(result, coords=result_coords, dims=list(itertools.chain(*group_dims.keys())))
        return result