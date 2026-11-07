"""Paratyphi B marker regressions: both SNP states and unresolved evidence."""
import importlib.util
import pickle
from pathlib import Path
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'bin'))
spec = importlib.util.spec_from_file_location('seqsero_paratyphi', ROOT / 'bin/SeqSero2_package.py')
ss = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ss)

# Public ST28 sequence recovered independently from the reported NCBI genomes.
ST28 = (
    'AATATGCTGACCCGCTATGTTGGTGACCCTTCTGCCTGATCTTTTCGCCTGCGTGAGCATCGCAGGCGCTTCAACCAATT'
    'CTTCAATTCTATTAATCACGTCCAGTAACGATCCAACTACGATCTTGATGTAGATCCCAATCCTGCATCTTTCCTGAGAA'
    'ACGATTTGACTTAGTTTTTGTATGGGAATAATTTCTCCCTGTCAACATTGGATACAATATCATATATCATTGATTGGATA'
    'ATTATAACTATGGAACCCATTACCCTTACCCTATGTTTATTGGTGTTCGCT'
)

class ParatyphiBTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        refs = ss.multifasta_dict(str(ROOT / 'seqsero2_db/H_and_O_and_specific_genes.fasta'))
        cls.refs = {k: v for k, v in refs.items() if k.startswith('gntR-family-regulatory-protein_dt-')}
        with (ROOT / 'seqsero2_db/antigens.pickle').open('rb') as f:
            library = pickle.load(f)
        cls.library = {k: library[k] for k in cls.refs}

    def prediction(self, sequence, scheme='ss2s'):
        kmers = set(ss.createKmerDict_reads([sequence], 27))
        _, _, scores = ss.get_kmer_dict(self.library, kmers)
        return ss.seqsero_from_formula_to_serotypes('4', 'b', '1,2', scores, 'I')

    def test_new_reference_pair_has_only_the_discriminating_snp(self):
        positive = self.refs['gntR-family-regulatory-protein_dt-positive_ST28-context__291']
        negative = self.refs['gntR-family-regulatory-protein_dt-negative_ST28-context__291']
        self.assertEqual(positive, ST28)
        self.assertEqual(len(positive), 291)
        self.assertEqual([(i+1, a, b) for i, (a, b) in enumerate(zip(positive, negative)) if a != b], [(252, 'G', 'A')])

    def test_packaged_kmers_match_all_four_fasta_probes(self):
        self.assertEqual(len(self.refs), 4)
        for name, seq in self.refs.items():
            with self.subTest(probe=name):
                self.assertEqual(set(ss.createKmerDict_reads([seq], 27)), self.library[name])

    def test_both_marker_states_in_both_sequence_backgrounds(self):
        for name, seq in self.refs.items():
            for scheme in ['ss2s', 'ss2']:
                with self.subTest(probe=name, scheme=scheme):
                    _, label, _, note, _ = self.prediction(seq, scheme)
                    if 'dt-positive' in name:
                        self.assertEqual(label, 'Paratyphi B var. L(+) tartrate+')
                        self.assertIn('was not detected', note)
                    else:
                        self.assertEqual(label, 'Paratyphi B')
                        self.assertTrue(note.startswith('Detected the SNP'))

    def test_reverse_complement_st28_has_the_same_variant(self):
        self.assertEqual(self.prediction(ss.reverse_complement(ST28))[1], 'Paratyphi B var. L(+) tartrate+')

    def test_terminal_marker_chunk_is_not_dropped(self):
        target = set().union(*self.library.values())
        for sequence in [ST28, ST28[:251]+'A'+ST28[252:]]:
            with self.subTest(state=sequence[251]), tempfile.TemporaryDirectory() as temp:
                path = Path(temp) / 'short-marker.fasta'
                path.write_text('>marker\n'+sequence+'\n')
                observed = ss.target_multifasta_kmerizer(str(path), 27, target)
                self.assertEqual(observed, set(ss.createKmerDict_reads([sequence], 27)))

    def test_weak_later_probe_cannot_overwrite_stronger_variant_evidence(self):
        prefix = 'gntR-family-regulatory-protein_dt-'
        cases = [
            ([(prefix+'positive_old', 100), (prefix+'negative_old', 90), (prefix+'positive_new', 60)], 'Paratyphi B var. L(+) tartrate+'),
            ([(prefix+'negative_old', 100), (prefix+'positive_old', 90), (prefix+'negative_new', 60)], 'Paratyphi B')]
        for pairs, expected in cases:
            for items in [pairs, list(reversed(pairs))]:
                result = ss.seqsero_from_formula_to_serotypes('4', 'b', '1,2', dict(items), 'I')
                self.assertEqual(result[1], expected)

    def test_mixed_alleles_do_not_become_a_positive_call(self):
        pos = ST28
        neg = ST28[:251] + 'A' + ST28[252:]
        kmers = set(ss.createKmerDict_reads([pos, neg], 27))
        _, _, scores = ss.get_kmer_dict(self.library, kmers)
        result = ss.seqsero_from_formula_to_serotypes('4', 'b', '1,2', scores, 'I')
        self.assertEqual(result[1], 'Paratyphi B (tartrate variant unresolved)')
        self.assertIn('missing or tied', result[3])

    def test_absent_and_shared_only_evidence_are_explicitly_unresolved(self):
        positive = self.library['gntR-family-regulatory-protein_dt-positive_ST28-context__291']
        negative = self.library['gntR-family-regulatory-protein_dt-negative_ST28-context__291']
        for observed in [set(), positive & negative]:
            _, _, scores = ss.get_kmer_dict(self.library, observed)
            result = ss.seqsero_from_formula_to_serotypes('4', 'b', '1,2', scores, 'I')
            self.assertEqual(result[1], 'Paratyphi B (tartrate variant unresolved)')
            self.assertIn('Unable to distinguish', result[3])

if __name__ == '__main__':
    unittest.main()
