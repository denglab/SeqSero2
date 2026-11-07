"""Decision regressions for reported issues; fixtures are synthetic, not isolate validation."""
import contextlib,importlib.util,io,itertools,os,subprocess,sys,tempfile,unittest
from pathlib import Path
from unittest import mock
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'bin'))
import SeqSero2_package as ss
import SeqSero2_summarize as summary


def hit(gene,node=1,identity=1.,length=1002):
    return (f'{gene}__{length}___NODE_{node}_length_{length}_cov_30',length*2.,identity,set(range(1,length)))


class IssueTests(unittest.TestCase):
    def test_7_current_default_agrees_with_documented_allele_mode(self):
        with tempfile.TemporaryDirectory() as t:
            p=Path(t)/'reads.fq';p.write_text('@r\nACGT\n+\nIIII\n')
            self.assertEqual(ss.parse_args(['-t','3','-i',str(p)]).m,'a')

    def test_15_exact_reported_wbaV_scores_without_wzy_resolve_O9(self):
        with tempfile.TemporaryDirectory() as t,ss.working_directory(t):
            for score in [54.09836065573771,80.84016393442623]:
                result=ss.call_O_and_H_type({'O-1,3,19_not_in_3,10__130':73.07692307692307,'O-9,46_wbaV__1002':score},{},{},t)
                self.assertEqual(result[0],'O-9')

    def test_43_54_multiple_auxiliary_O_markers_do_not_crash(self):
        rows=[hit('O-1,3,19_not_in_3,10',i,length=130) for i in [1,2,3]]
        result=ss.decide_O_type_and_get_special_genes(rows,rows)
        self.assertEqual(result[0],'-')
        self.assertEqual(result[4],'')

    def test_45_partial_wzy_is_same_family_without_hiding_other_O(self):
        for gene in ['O-9,46_wzy','O-9,46_wzy_partial']:
            rows=[hit('O-9,46_wbaV'),hit(gene,2,length=216)]
            result=ss.decide_O_type_and_get_special_genes(rows,rows)
            self.assertEqual(result[0],'O-9,46')
            self.assertEqual(result[4],'')
            mixed=rows+[hit('O-7_wzx',3)]
            self.assertTrue(ss.decide_O_type_and_get_special_genes(mixed,mixed)[4])

    def test_46_O2_O9_maximum_evidence_is_order_independent(self):
        base=hit('O-9,46_wbaV')
        markers=[hit('tyr-O-9',2,.99),hit('tyr-O-9',3,.4),hit('tyr-O-2',4,.7),hit('tyr-O-2',5,.5)]
        for ordered in itertools.permutations(markers):
            rows=[base,*ordered]
            self.assertEqual(ss.decide_O_type_and_get_special_genes(rows,rows)[0],'O-9')

    def test_57_missing_type_has_actionable_error(self):
        with self.assertRaisesRegex(ValueError,'-t 4'):
            ss.get_input_K('genome.fna',{},None,27)

    def test_13_invalid_salmid_table_is_not_a_prediction(self):
        with tempfile.TemporaryDirectory() as t,ss.working_directory(t):
            for data in [b'',b'warning\n',b'a\tb\n1\t2\n',b'a\tb\tc\td\te\tf\tg\th\n0\t0\t0\t0\tbad\t1\t2\t3\n']:
                with mock.patch.object(ss.subprocess,'run',return_value=subprocess.CompletedProcess([],0,data,b'')):
                    with self.assertRaisesRegex(RuntimeError,'SalmID returned'):
                        ss.judge_subspecies('reads.fq')

    def test_24_O13_database_resolves_from_an_unrelated_directory(self):
        with tempfile.TemporaryDirectory() as t,ss.working_directory(t):
            result=ss.seqsero_from_formula_to_serotypes('13','b','e,n,z',{},'I')
            self.assertEqual(result[0],'13:b:e,n,z')

    def test_salmid_uses_literal_filename_and_propagates_failure(self):
        with tempfile.TemporaryDirectory() as t,ss.working_directory(t):
            with mock.patch.object(ss.subprocess,'run',side_effect=subprocess.CalledProcessError(7,['SalmID.py'])) as run:
                with self.assertRaises(subprocess.CalledProcessError):
                    ss.judge_subspecies('reads (1);literal.fq')
                self.assertEqual(run.call_args.args[0],['SalmID.py','-i','reads (1);literal.fq'])
                self.assertTrue(run.call_args.kwargs['check'])

    def test_49_wrong_input_format_is_rejected_before_staging(self):
        with tempfile.TemporaryDirectory() as t:
            source=Path(t)/'reads.fastq';source.write_text('@r\nACGT\n+\nIIII\n')
            output=Path(t)/'output'
            result=subprocess.run([sys.executable,str(ROOT/'bin/SeqSero2_package.py'),'-m','k','-t','4','-i',str(source),'-d',str(output)],capture_output=True,text=True)
            self.assertEqual(result.returncode,2)
            self.assertIn('input format does not match',result.stderr)
            self.assertFalse(output.exists())

    def test_40_summary_combines_mixed_workflow_columns_and_quotes(self):
        with tempfile.TemporaryDirectory() as t:
            a=Path(t)/'a.tsv';b=Path(t)/'b.tsv'
            a.write_text('Sample name\tPredicted serotype\tNote\nA\tCerro\tclear\n')
            b.write_text('Sample name\tPredicted serotype\tPotential inter-serotype contamination\tNote\nB\tBraenderup\tno\tclear\n')
            out=io.StringIO();self.assertEqual(summary.summarize([a,b],out),2)
            import csv
            rows=list(csv.DictReader(io.StringIO(out.getvalue()),delimiter='\t'))
            self.assertEqual(rows[0]['Sample name'],'A');self.assertEqual(rows[1]['Sample name'],'B')
            self.assertEqual(rows[0]['Potential inter-serotype contamination'],'')
            a.write_text('A\tCerro\tclear\n')
            with self.assertRaisesRegex(ValueError,'missing TSV header'):
                summary.summarize([a],io.StringIO())


if __name__=='__main__':unittest.main()
