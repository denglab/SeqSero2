"""Software regressions; external scientific tools are stubbed, not validated."""
import contextlib
import gzip
import importlib.util
import io
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest import mock

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / 'bin' / 'SeqSero2_package.py'
spec = importlib.util.spec_from_file_location('seqsero2s', SCRIPT)
seqsero = importlib.util.module_from_spec(spec)
# Importing helpers must not inspect the caller's arguments or exit.
with mock.patch.object(sys, 'argv', ['caller', '--unrelated-option']):
    spec.loader.exec_module(seqsero)


class SoftwareTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)

    def cli(self, *args, env=None):
        return subprocess.run([sys.executable, str(SCRIPT), *map(str, args)],
                              cwd=self.root, env=env, capture_output=True, text=True)

    def input(self, name='reads.fastq', content='@r\nACGT\n+\nIIII\n'):
        path = self.root / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(content)
        return path

    def test_help_and_version_without_inputs(self):
        for args in [(), ('-h',), ('--version',)]:
            with self.subTest(args=args):
                result = self.cli(*args)
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertNotIn('Traceback', result.stderr)

    def test_invalid_inputs_fail_before_creating_output(self):
        source = self.input()
        cases = [[], ['-t', '2', '-i', source], ['-t', '3', '-i', source, source],
                 ['-t', '3', '-i', self.root / 'missing.fastq'],
                 ['-t', '4', '-i', source],
                 *[['-t', '3', '-i', source, '-p', p] for p in ['0', '-1', 'abc', '1.5']]]
        for args in cases:
            with self.subTest(args=args):
                output = self.root / 'output'
                result = self.cli(*args, '-d', output)
                self.assertEqual(result.returncode, 2, result.stderr)
                self.assertNotIn('Traceback', result.stderr)
                self.assertFalse(output.exists())

    def test_colliding_input_names_are_rejected(self):
        first = self.input('a/reads.fq')
        second = self.input('b/reads.fq')
        result = self.cli('-t', '2', '-i', first, second)
        self.assertEqual(result.returncode, 2)
        reserved = self.input('data_log.txt')
        result = self.cli('-t', '3', '-i', reserved)
        self.assertEqual(result.returncode, 2)
        self.assertEqual(reserved.read_text(), '@r\nACGT\n+\nIIII\n')

    def test_existing_output_is_never_reused_or_cleaned(self):
        source = self.input()
        output = self.root / 'existing'
        output.mkdir()
        keep = output / 'keep.fastq'
        keep.write_text('keep me')
        result = self.cli('-m', 'k', '-t', '3', '-i', source, '-d', output, '-c')
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('already exists', result.stderr)
        self.assertEqual(keep.read_text(), 'keep me')
        self.assertEqual(list(output.iterdir()), [keep])

    def test_dependency_check_exit_status(self):
        for found, expected in [(None, 1), ('/stub/tool', 0)]:
            with mock.patch.object(seqsero.shutil, 'which', return_value=found), contextlib.redirect_stdout(io.StringIO()):
                self.assertEqual(seqsero.check_dependencies(), expected)

    def test_paths_are_literal_and_diagnostics_are_logged(self):
        output = self.root / 'output (1).txt'
        log = self.root / 'log.txt'
        literal = 'reads (1); touch SHOULD_NOT_EXIST; $(echo bad) [x].fq'
        with seqsero.working_directory(self.root):
            seqsero.run_command([sys.executable, '-c',
                                 'import sys; print(sys.argv[1]); print("diagnostic", file=sys.stderr)', literal],
                                stdout=output, log=log)
            self.assertFalse(Path('SHOULD_NOT_EXIST').exists())
        self.assertEqual(output.read_text().strip(), literal)
        self.assertIn('diagnostic', log.read_text())

    def test_deinterleave_plain_and_gzip_and_restore_directory(self):
        forward = '@r/1\nACGT\n+\nIIII\n'
        reverse = '@r/2\nTGCA\n+\nIIII\n'
        previous = os.getcwd()
        for suffix in ['.fastq', '.fastq.gz']:
            with self.subTest(suffix=suffix):
                output = self.root / ('out' + suffix)
                output.mkdir()
                source = output / ('reads (1); literal' + suffix)
                if suffix.endswith('.gz'):
                    with gzip.open(source, 'wt') as handle:
                        handle.write(forward + reverse)
                else:
                    source.write_text(forward + reverse)
                f1, f2 = seqsero.get_input_files(str(output), [str(source)], '1', str(ROOT / 'bin'))
                self.assertEqual((output / f1).read_text(), forward)
                self.assertEqual((output / f2).read_text(), reverse)
                self.assertEqual(os.getcwd(), previous)

    def test_truncated_pairs_raise_and_restore_directory(self):
        source = self.input()
        previous = os.getcwd()
        with self.assertRaisesRegex(ValueError, 'complete pairs'):
            seqsero.get_input_files(str(self.root), [str(source)], '1', '')
        self.assertEqual(os.getcwd(), previous)

    def test_uncompressed_reads_are_streamed_and_match_gzip(self):
        source = self.input(content='@r\nACGTACGT\n+\nIIIIIIII\n')
        zipped = self.root / 'reads.fastq.gz'
        with gzip.open(zipped, 'wt') as handle:
            handle.write(source.read_text())
        mers = {'ACG', 'CGT', 'GTA', 'TAC'}
        self.assertEqual(seqsero.target_read_kmerizer(str(source), 3, mers),
                         seqsero.target_read_kmerizer(str(zipped), 3, mers))
        real_open = open
        def no_readlines(*args, **kwargs):
            handle = real_open(*args, **kwargs)
            handle.readlines = mock.Mock(side_effect=AssertionError('must stream'))
            return handle
        with mock.patch('builtins.open', side_effect=no_readlines):
            self.assertTrue(seqsero.target_read_kmerizer(str(source), 3, mers))

    def test_single_end_samse_and_samtools_version_branches(self):
        for version, legacy in [('1.2', True), ('1.2.1', False), ('1.20', False), ('1.22.1-7-gabc', False)]:
            commands = []
            def record(args, **kwargs):
                commands.append((args, kwargs))
            with mock.patch.object(seqsero, 'run_command', side_effect=record), \
                 mock.patch.object(seqsero.subprocess, 'run', return_value=subprocess.CompletedProcess([], 1, b'', ('Version: ' + version + '\n').encode())):
                seqsero.map_and_sort('2', 'db.fa', './reads (1).fq', '', 'a.sam', 'a.bam', 'a.sai', '', 'sorted.bam', 'sam')
            samse = next(args for args, _ in commands if args[:2] == ['bwa', 'samse'])
            self.assertEqual(samse[-1], './reads (1).fq')
            sort, options = commands[-1]
            self.assertEqual(sort[:2], ['samtools', 'sort'])
            self.assertEqual('stdout' not in options, legacy)

    def test_no_extracted_reads_returns_without_unbound_variables(self):
        with seqsero.working_directory(self.root), mock.patch.object(seqsero, 'run_command', side_effect=lambda args, **kw: Path(kw['stdout']).write_bytes(b'')):
            result = seqsero.extract_mapped_reads_and_do_assembly_and_blast('now', 'sorted.bam', 'combined.fq', '', '', '1', 'reads.fq', '', 'db', 'mem', 'auto')
        self.assertEqual(result, ('NA', None))

    def test_assembly_commands_preserve_paths_threads_and_quality_option(self):
        commands = []
        def tool(args, **kwargs):
            commands.append(args)
            if 'stdout' in kwargs:
                Path(kwargs['stdout']).write_text('A' * 120)
            if args[:2] == ['samtools', 'bam2fq'] and '-1' in args:
                Path(args[args.index('-1') + 1]).write_text('A' * 120)
                Path(args[args.index('-2') + 1]).write_text('A' * 120)
            if args[0] == 'spades.py':
                for option in ['--pe1-s', '--pe1-1', '--pe1-2']:
                    path = Path(args[args.index(option) + 1])
                    self.assertNotIn(' ', str(path))
                    self.assertEqual(path.read_text(), 'A' * 120)
                output = Path(args[args.index('-o') + 1])
                output.mkdir()
                (output / 'contigs.fasta').write_text('>contig\nACGT\n')
        with seqsero.working_directory(self.root), mock.patch.object(seqsero, 'run_command', side_effect=tool):
            xml, fasta = seqsero.extract_mapped_reads_and_do_assembly_and_blast(
                'now (1)', 'sorted (1).bam', 'combined (1).fq', 'forward (1).fq', 'reverse (1).fq',
                '10', 'reads (1).fq', 'reverse (1).fq', 'db.fa', 'mem', '64')
            self.assertEqual(Path(fasta).read_text(), '>contig\nACGT\n')
            self.assertFalse(Path('now (1)_temp').exists())
        spades = next(args for args in commands if args[0] == 'spades.py')
        self.assertEqual(spades[spades.index('-t') + 1], '4')
        self.assertEqual(spades[spades.index('--phred-offset') + 1], '64')
        self.assertEqual(xml, 'blasted_output.xml')
        blast = next(args for args in commands if args[0] == 'blastn')
        self.assertNotIn(' ', blast[blast.index('-db') + 1])
        self.assertFalse(Path(blast[blast.index('-db') + 1]).parent.exists())

    def test_assembly_failure_is_not_reported_as_a_prediction(self):
        source = self.input()
        output = self.root / 'failed-run'
        failure = subprocess.CalledProcessError(7, ['spades.py'])
        previous = os.getcwd()
        with mock.patch.object(sys, 'argv', ['SeqSero2_package.py', '-t', '3', '-i', str(source), '-d', str(output)]), \
             mock.patch.object(seqsero, 'map_and_sort'), \
             mock.patch.object(seqsero, 'extract_mapped_reads_and_do_assembly_and_blast', side_effect=failure):
            with self.assertRaises(subprocess.CalledProcessError):
                seqsero.main()
        self.assertFalse((output / 'SeqSero_result.tsv').exists())
        self.assertEqual(os.getcwd(), previous)


    def test_gzip_nanopore_fasta_is_decompressed(self):
        source = self.root / 'reads.fasta.gz'
        with gzip.open(source, 'wt') as output:
            output.write('>read\nACGT\n')
        converted = Path(seqsero.format_check(str(source)))
        self.assertEqual(converted.read_text(), '>read\nACGT\n')

    def test_failed_tools_keep_error_status_and_log(self):
        log = self.root / 'error.log'
        with self.assertRaises(subprocess.CalledProcessError) as error:
            seqsero.run_command([sys.executable, '-c', 'import sys; print("tool failed", file=sys.stderr); sys.exit(7)'], log=log)
        self.assertEqual(error.exception.returncode, 7)
        self.assertIn('tool failed', log.read_text())

    def test_bad_nanopore_format_has_clear_error(self):
        source = self.input(content='invalid\n')
        with self.assertRaisesRegex(ValueError, 'FASTA or FASTQ'):
            seqsero.format_check(str(source))

    def test_assembly_cli_with_special_paths_and_clean_mode(self):
        tool_dir = self.root / 'tools'
        tool_dir.mkdir()
        tool = tool_dir / 'mlst'
        tool.write_text('#!' + sys.executable + '\n' +
                        'import json, sys\n' +
                        'json.dump(sys.argv[1:], open("mlst_args.json", "w"))\n' +
                        'json.dump([{"sequence_type":"-", "alleles":{}}], open("mlst.json", "w"))\n')
        tool.chmod(0o755)
        env = dict(os.environ, PATH=str(tool_dir) + os.pathsep + os.environ.get('PATH', ''))
        content = '>contig\n' + 'A' * 300 + '\n'
        source = self.input('Centroidal_assembly_(SRR9043663);literal.fasta', content)
        for clean in [False, True]:
            with self.subTest(clean=clean):
                output = self.root / ('output (1);literal ' + str(clean))
                args = ['-m', 'k', '-t', '4', '-i', source, '-d', output]
                if clean:
                    args += ['-c']
                result = self.cli(*args, env=env)
                self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
                self.assertEqual(source.read_text(), content)
                if clean:
                    self.assertFalse(output.exists())
                else:
                    self.assertTrue((output / 'SeqSero_result.tsv').is_file())
                    self.assertFalse((output / source.name).exists())

    def test_database_updater_uses_packaged_data_directory(self):
        update_spec = importlib.util.spec_from_file_location('updater', ROOT / 'bin' / 'SeqSero2_update_kmer_database.py')
        updater = importlib.util.module_from_spec(update_spec)
        update_spec.loader.exec_module(updater)
        with mock.patch.object(updater, 'parse_args'), \
             mock.patch.object(updater, 'multifasta_to_kmers_dict', return_value={}) as load, \
             mock.patch.object(updater, 'get_salmid_invA_database', return_value={}), \
             mock.patch('builtins.open', mock.mock_open()), mock.patch.object(updater.pickle, 'dump'):
            updater.main()
        self.assertEqual(load.call_args.args[0], str(ROOT / 'seqsero2_db' / 'H_and_O_and_specific_genes.fasta'))

    def test_compressed_assembly_cli_prepares_fasta_for_both_stages(self):
        tools = self.root / 'tools'
        tools.mkdir()
        mlst = tools / 'mlst'
        mlst.write_text('#!' + sys.executable + '\nimport json, sys\n'
                        'assert open(sys.argv[-1]).read().startswith(">contig\\n")\n'
                        'json.dump([{"sequence_type":"-", "alleles":{}}], open("mlst.json", "w"))\n')
        mlst.chmod(0o755)
        source = self.root / 'assembly (1).fna.gz'
        content = '>contig\n' + 'A' * 300 + '\n'
        with gzip.open(source, 'wt') as handle:
            handle.write(content)
        output = self.root / 'compressed-output'
        env = dict(os.environ, PATH=str(tools) + os.pathsep + os.environ.get('PATH', ''))
        result = self.cli('-m', 'k', '-t', '4', '-i', source, '-d', output, env=env)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertTrue((output / 'SeqSero_result.tsv').exists())
        with gzip.open(source, 'rt') as handle:
            self.assertEqual(handle.read(), content)


if __name__ == '__main__':
    unittest.main()
