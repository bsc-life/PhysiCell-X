"""Small synthetic output tests; never launches a simulation."""
import importlib.util
from pathlib import Path
import tempfile
import unittest
import xml.etree.ElementTree as ET

import numpy as np
from scipy.io import savemat

SPEC = importlib.util.spec_from_file_location('comparison', Path(__file__).resolve().parents[2] / 'tools/compare_simulations.py')
comparison = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(comparison)


def fixture(folder, name, time, source='PhysiCell', reorder=False, stamp='2026-10-07T12:00:00Z'):
    root = ET.Element('MultiCellDS')
    metadata = ET.SubElement(root, 'metadata')
    ET.SubElement(metadata, 'current_time', units='min').text = str(time)
    ET.SubElement(metadata, 'last_modified').text = stamp
    domain = ET.SubElement(ET.SubElement(root, 'microenvironment'), 'domain')
    mesh = ET.SubElement(domain, 'mesh')
    ET.SubElement(mesh, 'bounding_box').text = '-20 -20 -10 20 20 10'
    variables = ET.SubElement(domain, 'variables')
    ET.SubElement(variables, 'variable', ID='0', name='oxygen', units='dimensionless')
    ET.SubElement(ET.SubElement(domain, 'data'), 'filename').text = name + '_microenvironment0.mat'
    node = ET.SubElement(root, 'simplified_data', source=source)
    ET.SubElement(node, 'filename').text = name + '_cells.mat'
    labels = ET.SubElement(node, 'labels')
    layout = [('ID', 1), ('position', 3), ('dead', 1), ('radius', 1), ('pressure', 1), ('total_volume', 1), ('current_phase', 1), ('current_cycle_phase_exit_rate', 1)]
    values = [[1, 2], [-5, 5], [0, 0], [0, 0], [0, 1], [6, 6], [.1, .2], [100, 100], [14, 100], [.001, 0]]
    if reorder:
        # Change both row layout and cell order, not just one array's columns.
        layout = layout[2:] + layout[:2]
        values = values[4:] + values[:4]
    index = 0
    for key, size in layout:
        ET.SubElement(labels, 'label', index=str(index), size=str(size), units='none').text = key
        index += size
    cells = np.array(values, dtype=float)
    micro = np.array([[-10, 10, -10, 10], [-10, -10, 10, 10], [0, 0, 0, 0], [8000] * 4, [5, 6, 7, 8]], dtype=float)
    if reorder:
        cells = cells[:, ::-1]; micro = micro[:, [2, 0, 3, 1]]
    savemat(folder / (name + '_cells.mat'), {'cells': cells})
    savemat(folder / (name + '_microenvironment0.mat'), {'multiscale_microenvironment': micro})
    path = folder / (name + '.xml')
    ET.ElementTree(root).write(path)
    return path


class ComparisonTests(unittest.TestCase):
    def test_biofvm_schema_and_reordered_rows_voxels_ids(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp)
            a = comparison.read_snapshot(fixture(p, 'a', 0))
            b = comparison.read_snapshot(fixture(p, 'b', 0, source='BioFVM', reorder=True))
            self.assertTrue(comparison.same_state(a, b))
            np.testing.assert_array_equal(a['fields']['oxygen'], b['fields']['oxygen'])

    def test_initial_final_stale_and_alias_discovery(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp)
            fixture(p, 'initial', 0)
            fixture(p, 'output00000000', 0, reorder=True)
            fixture(p, 'final', 12.01, stamp='2026-10-07T12:01:00Z')
            fixture(p, 'output00000999', 999, stamp='2026-10-06T12:00:00Z')
            warnings = []
            results = comparison.discover(p, 1e-4, warnings)
            self.assertEqual([r['time'] for r in results], [0, 12.01])
            self.assertEqual(len(warnings), 1)

    def test_conflicting_same_time_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp)
            fixture(p, 'initial', 0)
            fixture(p, 'output00000000', 0)
            savemat(p / 'output00000000_microenvironment0.mat', {'multiscale_microenvironment': np.array([
                [-10, 10, -10, 10], [-10, -10, 10, 10], [0, 0, 0, 0], [8000] * 4, [6, 6, 7, 8]])})
            with self.assertRaisesRegex(ValueError, 'Conflicting'):
                comparison.discover(p, 1e-4, [])

    def test_time_matching_and_empty_radial_bins(self):
        pairs, a, b = comparison.match_times([dict(time=0), dict(time=60), dict(time=100)],
                                            [dict(time=0.00001), dict(time=90)], .0001)
        self.assertEqual(len(pairs), 1)
        self.assertEqual(a, [60, 100]); self.assertEqual(b, [90])
        values = comparison.bin_mean(np.array([1., 3.]), np.array([2., 4.]), np.array([0, 2, 4, 6]))
        np.testing.assert_allclose(values[:2], [2, 4])
        self.assertTrue(np.isnan(values[2]))

    def test_planar_gradients_and_nearest_voxel_exposure(self):
        from types import SimpleNamespace
        with tempfile.TemporaryDirectory() as tmp:
            a = comparison.read_snapshot(fixture(Path(tmp), 'a', 0))
            args = SimpleNamespace(center=np.zeros(3), core_radius=100)
            values = comparison.field_comparison(a, a, 'oxygen', args)
            self.assertEqual(values['relative_l2_percent'], 0)
            self.assertEqual(values['gradient_relative_l2_percent'], 0)
            self.assertEqual(values['simulation_negative_voxels'], 0)
            self.assertEqual(len(comparison.nearest_field(a, np.empty((0, 3)))), 0)

    def test_empty_and_single_cell_populations(self):
        from types import SimpleNamespace
        args = SimpleNamespace(center=np.zeros(3), growth_field='current_cycle_phase_exit_rate', growth_threshold=1e-12)
        with tempfile.TemporaryDirectory() as tmp:
            s = comparison.read_snapshot(fixture(Path(tmp), 'a', 0))
            for n in (1, 0):
                subset = dict(s, n=n, labels={k: v[:, :n] for k, v in s['labels'].items()})
                result = comparison.cell_metrics(subset, args)
                self.assertEqual(result['n_cells'], n)
                self.assertTrue(np.isnan(result['nearest_neighbor_mean']))
                self.assertEqual(result['dead_count'], 0)


if __name__ == '__main__':
    unittest.main()
