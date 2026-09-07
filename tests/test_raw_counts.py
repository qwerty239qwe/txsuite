import importlib.util
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch


@unittest.skipUnless(
    importlib.util.find_spec("scanpy"), "Requires the single-cell image"
)
class RawCountsTest(unittest.TestCase):
    def test_pseudobulk_rejects_invalid_cells_before_sum(self):
        import numpy as np
        import pandas as pd
        import scanpy as sc
        from scipy import sparse

        from txsuite.resources.single_cell_python import single_cell as script

        with tempfile.TemporaryDirectory() as temporary:
            args = SimpleNamespace(
                input="unused.h5ad",
                output=temporary,
                sample_column="sample",
                design="condition",
                covariate=[],
                group_column=None,
                group_value=None,
                reference=None,
                test=None,
                counts_layer="counts",
            )
            for values in ([-1, 3], [0.5, 1.5], [np.nan, 1], [np.inf, 1], [2**53, 1]):
                for convert in (np.asarray, sparse.csr_matrix):
                    with self.subTest(values=values, storage=convert):
                        data = sc.AnnData(
                            convert(np.array(values, dtype=float).reshape(2, 1)),
                            obs=pd.DataFrame(
                                {"sample": ["s", "s"], "condition": ["c", "c"]},
                                index=["a", "b"],
                            ),
                        )
                        with (
                            patch.object(script.sc, "read_h5ad", return_value=data),
                            self.assertRaises(ValueError),
                        ):
                            script.pseudobulk(args)
            data.X = sparse.csr_matrix([[1], [2]])
            data.layers["raw"] = sparse.csr_matrix([[3], [4]])
            args.counts_layer = "raw"
            with patch.object(script.sc, "read_h5ad", return_value=data):
                script.pseudobulk(args)
                self.assertIn(
                    "\t7", (Path(temporary) / "pseudobulk-counts.tsv").read_text()
                )
                args.counts_layer = "missing"
                with self.assertRaisesRegex(ValueError, "does not exist"):
                    script.pseudobulk(args)
