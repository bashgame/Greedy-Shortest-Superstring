"""
Test Cases for greedy shortest super-string algorithm
"""
from unittest import TestCase
from modules.FastxFiles import reads_from_fastq as fastq
from modules.overlaps import overlap, build_dictionary, overlap_all_pairs
from modules.greedySuperString import greedySuperString

class TestGreedySS(TestCase):
    """ Testing greedy shortest super-string algorithm """

    @classmethod
    def setUpClass(cls):
        """ Set up variables and data structures """
        cls.fastq_filename = 'data/ads1_week4_reads.fq'
        _, cls.reads, _ = fastq(cls.fastq_filename)

    @classmethod
    def tearDownClass(cls):
        """ Do something to clean up tests, please update """

    def setUp(self):
        """ Ensure each test runs clean, please update """

    def tearDown(self):
        """ Ensure each test runs clean, please update """

    ###########################################################################
    #   T E S T  C A S E S
    ###########################################################################


    def test_fastq(self):
        """ It should return name, reads, and quals """
        name, reads, quals = fastq(self.fastq_filename)
        self.assertIsNotNone(name)
        self.assertIsNotNone(reads)
        self.assertIsNotNone(quals)

    def test_fastq_filename(self):
        """ It should return an error message """
        name, reads, quals = fastq('bob')
        self.assertIn('Error', name)
        self.assertIsNone(reads)
        self.assertIsNone(quals)

    def test_overlap(self):
        """ It should return the length of overlap between 2 strings """
        test_str1 = 'GATTACA'
        test_str2 = 'ACATTAG'
        expected = 3
        min_len = 1
        result = overlap(test_str1, test_str2, min_len)
        self.assertEqual(expected, result)

    def test_bad_overlap(self):
        """ It should return a zero """
        test_str1 = 'GATTACA'
        test_str2 = 'ACATTAG'
        expected = 0
        min_len = 4
        result = overlap(test_str1, test_str2, min_len)
        self.assertEqual(expected, result)
    
    def test_create_dictionary(self):
        """ It should create a dictionary of substrings and sets of reads """
        test_dict = {
            'foo': {'foobarbaz', 'foobaz'}, 'bar': {'foobarbaz'},
            'baz': {'foobarbaz', 'foobaz'}, 'rba': {'foobarbaz'},
            'oba': {'foobarbaz', 'foobaz'}, 'arb': {'foobarbaz'},
            'oob': {'foobarbaz', 'foobaz'}
        }
        min_length = 3
        res_dict = build_dictionary(['foobarbaz', 'foobaz'], min_length)
        self.assertEqual(test_dict, res_dict)

    def test_overlap_all(self):
        """ It should return a set of ordered tuples """
        test_reads = ['CGTACG', 'TACGTA', 'GTACGT', 'ACGTAC', 'GTACGA', 'TACGAT']
        min_len = 4
        results = overlap_all_pairs(test_reads, min_len)
        expected = {
            ('CGTACG', 'TACGTA', 4),
            ('CGTACG', 'GTACGT', 5),
            ('CGTACG', 'GTACGA', 5),
            ('CGTACG', 'TACGAT', 4),
            ('TACGTA', 'ACGTAC', 5),
            ('TACGTA', 'CGTACG', 4),
            ('GTACGT', 'TACGTA', 5),
            ('GTACGT', 'ACGTAC', 4),
            ('ACGTAC', 'GTACGA', 4),
            ('ACGTAC', 'GTACGT', 4),
            ('ACGTAC', 'CGTACG', 5),
            ('GTACGA', 'TACGAT', 5)
        }
        self.assertEqual(expected, results)

    def test_super_string(self):
        """ It should return a super-string """
        strings = ['CCT', 'CTT', 'TGC', 'TGG', 'GAT', 'ATT']
        result = greedySuperString(strings)
        is_super = True
        for s in strings:
            is_super = s in result
        self.assertTrue(is_super)
