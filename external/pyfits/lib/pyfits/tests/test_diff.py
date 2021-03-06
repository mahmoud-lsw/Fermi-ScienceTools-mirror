import textwrap

import numpy as np

import pyfits
from pyfits.column import Column
from pyfits.diff import *
from pyfits.hdu import HDUList, PrimaryHDU, ImageHDU
from pyfits.hdu.table import new_table
from pyfits.header import Header
from pyfits.tests import PyfitsTestCase

from nose.tools import (assert_true, assert_false, assert_equal,
                        assert_not_equal)


class TestDiff(PyfitsTestCase):
    def test_identical_headers(self):
        ha = Header([('A', 1), ('B', 2), ('C', 3)])
        hb = ha.copy()
        assert_true(HeaderDiff(ha, hb).identical)

    def test_slightly_different_headers(self):
        ha = Header([('A', 1), ('B', 2), ('C', 3)])
        hb = ha.copy()
        hb['C'] = 4
        assert_false(HeaderDiff(ha, hb).identical)

    def test_common_keywords(self):
        ha = Header([('A', 1), ('B', 2), ('C', 3)])
        hb = ha.copy()
        hb['C'] = 4
        hb['D'] = (5, 'Comment')
        assert_equal(HeaderDiff(ha, hb).common_keywords, ['A', 'B', 'C'])

    def test_different_keyword_count(self):
        ha = Header([('A', 1), ('B', 2), ('C', 3)])
        hb = ha.copy()
        del hb['B']
        diff = HeaderDiff(ha, hb)
        assert_false(diff.identical)
        assert_equal(diff.diff_keyword_count, (3, 2))

        # But make sure the common keywords are at least correct
        assert_equal(diff.common_keywords, ['A', 'C'])

    def test_different_keywords(self):
        ha = Header([('A', 1), ('B', 2), ('C', 3)])
        hb = ha.copy()
        hb['C'] = 4
        hb['D'] = (5, 'Comment')
        ha['E'] = (6, 'Comment')
        ha['F'] = (7, 'Comment')
        diff = HeaderDiff(ha, hb)
        assert_false(diff.identical)
        assert_equal(diff.diff_keywords, (['E', 'F'], ['D']))

    def test_different_keyword_values(self):
        ha = Header([('A', 1), ('B', 2), ('C', 3)])
        hb = ha.copy()
        hb['C'] = 4
        diff = HeaderDiff(ha, hb)
        assert_false(diff.identical)
        assert_equal(diff.diff_keyword_values, {'C': [(3, 4)]})

    def test_different_keyword_comments(self):
        ha = Header([('A', 1), ('B', 2), ('C', 3, 'comment 1')])
        hb = ha.copy()
        hb.comments['C'] = 'comment 2'
        diff = HeaderDiff(ha, hb)
        assert_false(diff.identical)
        assert_equal(diff.diff_keyword_comments,
                     {'C': [('comment 1', 'comment 2')]})

    def test_different_keyword_values_with_duplicate(self):
        ha = Header([('A', 1), ('B', 2), ('C', 3)])
        hb = ha.copy()
        ha.append(('C', 4))
        hb.append(('C', 5))
        diff = HeaderDiff(ha, hb)
        assert_false(diff.identical)
        assert_equal(diff.diff_keyword_values, {'C': [None, (4, 5)]})

    def test_asymmetric_duplicate_keywords(self):
        ha = Header([('A', 1), ('B', 2), ('C', 3)])
        hb = ha.copy()
        ha.append(('A', 2, 'comment 1'))
        ha.append(('A', 3, 'comment 2'))
        hb.append(('B', 4, 'comment 3'))
        hb.append(('C', 5, 'comment 4'))
        diff = HeaderDiff(ha, hb)
        assert_false(diff.identical)
        assert_equal(diff.diff_keyword_values, {})
        assert_equal(diff.diff_duplicate_keywords,
                     {'A': (3, 1), 'B': (1, 2), 'C': (1, 2)})

    def test_floating_point_tolerance(self):
        ha = Header([('A', 1), ('B', 2.00001), ('C', 3.000001)])
        hb = ha.copy()
        hb['B'] = 2.00002
        hb['C'] = 3.000002
        diff = HeaderDiff(ha, hb)
        assert_false(diff.identical)
        assert_equal(diff.diff_keyword_values,
                     {'B': [(2.00001, 2.00002)], 'C': [(3.000001, 3.000002)]})
        diff = HeaderDiff(ha, hb, tolerance=1e-6)
        assert_false(diff.identical)
        assert_equal(diff.diff_keyword_values, {'B': [(2.00001, 2.00002)]})

    def test_ignore_blanks(self):
        pyfits.STRIP_HEADER_WHITESPACE = False
        try:
            ha = Header([('A', 1), ('B', 2), ('C', 'A       ')])
            hb = ha.copy()
            hb['C'] = 'A'
            assert_not_equal(ha['C'], hb['C'])

            diff = HeaderDiff(ha, hb)
            # Trailing blanks are ignored by default
            assert_true(diff.identical)
            assert_equal(diff.diff_keyword_values, {})

            # Don't ignore blanks
            diff = HeaderDiff(ha, hb, ignore_blanks=False)
            assert_false(diff.identical)
            assert_equal(diff.diff_keyword_values, {'C': [('A       ', 'A')]})
        finally:
            pyfits.STRIP_HEADER_WHITESPACE = True

    def test_ignore_blank_cards(self):
        """Test for #152--ignore blank cards."""

        ha = Header([('A', 1), ('B', 2), ('C', 3)])
        hb = Header([('A', 1), ('', ''), ('B', 2), ('', ''), ('C', 3)])
        hc = ha.copy()
        hc.append()
        hc.append()

        # We now have a header with interleaved blanks, and a header with end
        # blanks, both of which should ignore the blanks
        assert_true(HeaderDiff(ha, hb).identical)
        assert_true(HeaderDiff(ha, hc).identical)
        assert_true(HeaderDiff(hb, hc).identical)

        assert_false(HeaderDiff(ha, hb, ignore_blank_cards=False).identical)
        assert_false(HeaderDiff(ha, hc, ignore_blank_cards=False).identical)

        # Both hb and hc have the same number of blank cards; since order is
        # currently ignored, these should still be identical even if blank
        # cards are not ignored
        assert_true(HeaderDiff(hb, hc, ignore_blank_cards=False).identical)

        hc.append()
        # But now there are different numbers of blanks, so they should not be
        # ignored:
        assert_false(HeaderDiff(hb, hc, ignore_blank_cards=False).identical)

    def test_ignore_keyword_values(self):
        ha = Header([('A', 1), ('B', 2), ('C', 3)])
        hb = ha.copy()
        hb['B'] = 4
        hb['C'] = 5
        diff = HeaderDiff(ha, hb, ignore_keywords=['*'])
        assert_true(diff.identical)
        diff = HeaderDiff(ha, hb, ignore_keywords=['B'])
        assert_false(diff.identical)
        assert_equal(diff.diff_keyword_values, {'C': [(3, 5)]})

        report = diff.report()
        assert_true('Keyword B        has different values' not in report)
        assert_true('Keyword C        has different values' in report)

        # Test case-insensitivity
        diff = HeaderDiff(ha, hb, ignore_keywords=['b'])
        assert_false(diff.identical)
        assert_equal(diff.diff_keyword_values, {'C': [(3, 5)]})

    def test_ignore_keyword_comments(self):
        ha = Header([('A', 1, 'A'), ('B', 2, 'B'), ('C', 3, 'C')])
        hb = ha.copy()
        hb.comments['B'] = 'D'
        hb.comments['C'] = 'E'
        diff = HeaderDiff(ha, hb, ignore_comments=['*'])
        assert_true(diff.identical)
        diff = HeaderDiff(ha, hb, ignore_comments=['B'])
        assert_false(diff.identical)
        assert_equal(diff.diff_keyword_comments, {'C': [('C', 'E')]})

        report = diff.report()
        assert_true('Keyword B        has different comments' not in report)
        assert_true('Keyword C        has different comments' in report)

        # Test case-insensitivity
        diff = HeaderDiff(ha, hb, ignore_comments=['b'])
        assert_false(diff.identical)
        assert_equal(diff.diff_keyword_comments, {'C': [('C', 'E')]})

    def test_trivial_identical_images(self):
        ia = np.arange(100).reshape((10, 10))
        ib = np.arange(100).reshape((10, 10))
        diff = ImageDataDiff(ia, ib)
        assert_true(diff.identical)
        assert_equal(diff.diff_total, 0)

    def test_identical_within_tolerance(self):
        ia = np.ones((10, 10)) - 0.00001
        ib = np.ones((10, 10)) - 0.00002
        diff = ImageDataDiff(ia, ib, tolerance=1.0e-4)
        assert_true(diff.identical)
        assert_equal(diff.diff_total, 0)

    def test_identical_comp_image_hdus(self):
        """Regression test for #189.

        For this test we mostly just care that comparing to compressed images
        does not crash, and returns the correct results.  Two compressed images
        will be considered identical if the decompressed data is the same.
        Obviously we test whether or not the same compression was used by
        looking for (or ignoring) header differences.
        """

        data = np.arange(100.0).reshape((10, 10))
        hdu = pyfits.CompImageHDU(data=data)
        hdu.writeto(self.temp('test.fits'))
        hdula = pyfits.open(self.temp('test.fits'))
        hdulb = pyfits.open(self.temp('test.fits'))
        diff = FITSDiff(hdula, hdulb)
        assert_true(diff.identical)

    def test_different_dimensions(self):
        ia = np.arange(100).reshape((10, 10))
        ib = np.arange(100) - 1

        # Although ib could be reshaped into the same dimensions, for now the
        # data is not compared anyways
        diff = ImageDataDiff(ia, ib)
        assert_false(diff.identical)
        assert_equal(diff.diff_dimensions, ((10, 10), (100,)))
        assert_equal(diff.diff_total, 0)

        report = diff.report()
        assert_true('Data dimensions differ' in report)
        assert_true('a: 10 x 10' in report)
        assert_true('b: 100' in report)
        assert_true('No further data comparison performed.')

    def test_different_pixels(self):
        ia = np.arange(100).reshape((10, 10))
        ib = np.arange(100).reshape((10, 10))
        ib[0,0] = 10
        ib[5,5] = 20
        diff = ImageDataDiff(ia, ib)
        assert_false(diff.identical)
        assert_equal(diff.diff_dimensions, ())
        assert_equal(diff.diff_total, 2)
        assert_equal(diff.diff_ratio, 0.02)
        assert_equal(diff.diff_pixels, [((0, 0), (0, 10)), ((5, 5), (55, 20))])

    def test_identical_tables(self):
        c1 = Column('A', format='L', array=[True, False])
        c2 = Column('B', format='X', array=[[0], [1]])
        c3 = Column('C', format='4I', dim='(2, 2)',
                    array=[[0, 1, 2, 3], [4, 5, 6, 7]])
        c4 = Column('D', format='J', bscale=2.0, array=[0, 1])
        c5 = Column('E', format='A3', array=['abc', 'def'])
        c6 = Column('F', format='E', unit='m', array=[0.0, 1.0])
        c7 = Column('G', format='D', bzero=-0.1, array=[0.0, 1.0])
        c8 = Column('H', format='C', array=[0.0+1.0j, 2.0+3.0j])
        c9 = Column('I', format='M', array=[4.0+5.0j, 6.0+7.0j])
        c10 = Column('J', format='PI(2)', array=[[0, 1], [2, 3]])

        columns = [c1, c2, c3, c4, c5, c6, c7, c8, c9, c10]

        ta = new_table(columns)
        tb = new_table([c.copy() for c in columns])

        diff = TableDataDiff(ta.data, tb.data)
        assert_true(diff.identical)
        assert_equal(len(diff.common_columns), 10)
        assert_equal(diff.common_column_names,
                     set(['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j']))
        assert_equal(diff.diff_ratio, 0)
        assert_equal(diff.diff_total, 0)

    def test_diff_empty_tables(self):
        """
        Regression test for #178.

        Ensure that diffing tables containing empty data doesn't crash.
        """

        c1 = Column('D', format='J')
        c2 = Column('E', format='J')
        thdu = new_table([c1, c2], nrows=0)

        hdula = pyfits.HDUList([thdu])
        hdulb = pyfits.HDUList([thdu])

        diff = FITSDiff(hdula, hdulb)
        assert_true(diff.identical)

    def test_ignore_table_fields(self):
        c1 = Column('A', format='L', array=[True, False])
        c2 = Column('B', format='X', array=[[0], [1]])
        c3 = Column('C', format='4I', dim='(2, 2)',
                    array=[[0, 1, 2, 3], [4, 5, 6, 7]])

        c4 = Column('B', format='X', array=[[1], [0]])
        c5 = Column('C', format='4I', dim='(2, 2)',
                    array=[[1, 2, 3, 4], [5, 6, 7, 8]])

        ta = new_table([c1, c2, c3])
        tb = new_table([c1, c4, c5])

        diff = TableDataDiff(ta.data, tb.data, ignore_fields=['B', 'C'])
        assert_true(diff.identical)

        # The only common column should be c1
        assert_equal(len(diff.common_columns), 1)
        assert_equal(diff.common_column_names, set(['a']))
        assert_equal(diff.diff_ratio, 0)
        assert_equal(diff.diff_total, 0)

    def test_different_table_field_names(self):
        ca = Column('A', format='L', array=[True, False])
        cb = Column('B', format='L', array=[True, False])
        cc = Column('C', format='L', array=[True, False])

        ta = new_table([ca, cb])
        tb = new_table([ca, cc])

        diff = TableDataDiff(ta.data, tb.data)

        assert_false(diff.identical)
        assert_equal(len(diff.common_columns), 1)
        assert_equal(diff.common_column_names, set(['a']))
        assert_equal(diff.diff_column_names, (['B'], ['C']))
        assert_equal(diff.diff_ratio, 0)
        assert_equal(diff.diff_total, 0)

    def test_different_table_field_counts(self):
        """
        Test tables with some common columns, but different number of columns
        overall.
        """

        ca = Column('A', format='L', array=[True, False])
        cb = Column('B', format='L', array=[True, False])
        cc = Column('C', format='L', array=[True, False])

        ta = new_table([cb])
        tb = new_table([ca, cb, cc])

        diff = TableDataDiff(ta.data, tb.data)

        assert_false(diff.identical)
        assert_equal(diff.diff_column_count, (1, 3))
        assert_equal(len(diff.common_columns), 1)
        assert_equal(diff.common_column_names, set(['b']))
        assert_equal(diff.diff_column_names, ([], ['A', 'C']))
        assert_equal(diff.diff_ratio, 0)
        assert_equal(diff.diff_total, 0)

    def test_different_table_rows(self):
        """
        Test tables taht are otherwise identical but one has more rows than the
        other.
        """

        ca1 = Column('A', format='L', array=[True, False])
        cb1 = Column('B', format='L', array=[True, False])
        ca2 = Column('A', format='L', array=[True, False, True])
        cb2 = Column('B', format='L', array=[True, False, True])

        ta = new_table([ca1, cb1])
        tb = new_table([ca2, cb2])

        diff = TableDataDiff(ta.data, tb.data)

        assert_false(diff.identical)
        assert_equal(diff.diff_column_count, ())
        assert_equal(len(diff.common_columns), 2)
        assert_equal(diff.diff_rows, (2, 3))
        assert_equal(diff.diff_values, [])

        report = diff.report()

        assert_true('Table rows differ' in report)
        assert_true('a: 2' in report)
        assert_true('b: 3' in report)
        assert_true('No further data comparison performed.')

    def test_different_table_data(self):
        """
        Test diffing table data on columns of several different data formats
        and dimensions.
        """

        ca1 = Column('A', format='L', array=[True, False])
        ca2 = Column('B', format='X', array=[[0], [1]])
        ca3 = Column('C', format='4I', dim='(2, 2)',
                     array=[[0, 1, 2, 3], [4, 5, 6, 7]])
        ca4 = Column('D', format='J', bscale=2.0, array=[0.0, 2.0])
        ca5 = Column('E', format='A3', array=['abc', 'def'])
        ca6 = Column('F', format='E', unit='m', array=[0.0, 1.0])
        ca7 = Column('G', format='D', bzero=-0.1, array=[0.0, 1.0])
        ca8 = Column('H', format='C', array=[0.0+1.0j, 2.0+3.0j])
        ca9 = Column('I', format='M', array=[4.0+5.0j, 6.0+7.0j])
        ca10 = Column('J', format='PI(2)', array=[[0, 1], [2, 3]])

        cb1 = Column('A', format='L', array=[False, False])
        cb2 = Column('B', format='X', array=[[0], [0]])
        cb3 = Column('C', format='4I', dim='(2, 2)',
                     array=[[0, 1, 2, 3], [5, 6, 7, 8]])
        cb4 = Column('D', format='J', bscale=2.0, array=[2.0, 2.0])
        cb5 = Column('E', format='A3', array=['abc', 'ghi'])
        cb6 = Column('F', format='E', unit='m', array=[1.0, 2.0])
        cb7 = Column('G', format='D', bzero=-0.1, array=[2.0, 3.0])
        cb8 = Column('H', format='C', array=[1.0+1.0j, 2.0+3.0j])
        cb9 = Column('I', format='M', array=[5.0+5.0j, 6.0+7.0j])
        cb10 = Column('J', format='PI(2)', array=[[1, 2], [3, 4]])

        ta = new_table([ca1, ca2, ca3, ca4, ca5, ca6, ca7, ca8, ca9, ca10])
        tb = new_table([cb1, cb2, cb3, cb4, cb5, cb6, cb7, cb8, cb9, cb10])

        diff = TableDataDiff(ta.data, tb.data, numdiffs=20)
        assert_false(diff.identical)
        # The column definitions are the same, but not the column values
        assert_equal(diff.diff_columns, ())
        assert_equal(diff.diff_values[0], (('A', 0), (True, False)))
        assert_equal(diff.diff_values[1], (('B', 1), ([1], [0])))
        assert_equal(diff.diff_values[2][0], ('C', 1))
        assert_true((diff.diff_values[2][1][0] == [[4, 5], [6, 7]]).all())
        assert_true((diff.diff_values[2][1][1] == [[5, 6], [7, 8]]).all())
        assert_equal(diff.diff_values[3], (('D', 0), (0, 2.0)))
        assert_equal(diff.diff_values[4], (('E', 1), ('def', 'ghi')))
        assert_equal(diff.diff_values[5], (('F', 0), (0.0, 1.0)))
        assert_equal(diff.diff_values[6], (('F', 1), (1.0, 2.0)))
        assert_equal(diff.diff_values[7], (('G', 0), (0.0, 2.0)))
        assert_equal(diff.diff_values[8], (('G', 1), (1.0, 3.0)))
        assert_equal(diff.diff_values[9], (('H', 0), (0.0+1.0j, 1.0+1.0j)))
        assert_equal(diff.diff_values[10], (('I', 0), (4.0+5.0j, 5.0+5.0j)))
        assert_equal(diff.diff_values[11][0], ('J', 0))
        assert_true((diff.diff_values[11][1][0] == [0, 1]).all())
        assert_true((diff.diff_values[11][1][1] == [1, 2]).all())
        assert_equal(diff.diff_values[12][0], ('J', 1))
        assert_true((diff.diff_values[12][1][0] == [2, 3]).all())
        assert_true((diff.diff_values[12][1][1] == [3, 4]).all())

        assert_equal(diff.diff_total, 13)
        assert_equal(diff.diff_ratio, 0.65)

    def test_identical_files_basic(self):
        """Test identicality of two simple, extensionless files."""

        a = np.arange(100).reshape((10, 10))
        hdu = PrimaryHDU(data=a)
        hdu.writeto(self.temp('testa.fits'))
        hdu.writeto(self.temp('testb.fits'))
        diff = FITSDiff(self.temp('testa.fits'), self.temp('testb.fits'))
        assert_true(diff.identical)

        report = diff.report()
        # Primary HDUs should contain no differences
        assert_true('Primary HDU' not in report)
        assert_true('Extension HDU' not in report)
        assert_true('No differences found.' in report)

    def test_partially_identical_files1(self):
        """
        Test files that have some identical HDUs but a different extension
        count.
        """

        a = np.arange(100).reshape((10, 10))
        phdu = PrimaryHDU(data=a)
        ehdu = ImageHDU(data=a)
        hdula = HDUList([phdu, ehdu])
        hdulb = HDUList([phdu, ehdu, ehdu])
        diff = FITSDiff(hdula, hdulb)
        assert_false(diff.identical)
        assert_equal(diff.diff_hdu_count, (2, 3))

        # diff_hdus should be empty, since the third extension in hdulb
        # has nothing to compare against
        assert_equal(diff.diff_hdus, [])

        report = diff.report()
        assert_true('Files contain different numbers of HDUs' in report)
        assert_true('a: 2\n b: 3' in report)
        assert_true('No differences found between common HDUs' in report)

    def test_partially_identical_files2(self):
        """
        Test files that have some identical HDUs but one different HDU.
        """

        a = np.arange(100).reshape((10, 10))
        phdu = PrimaryHDU(data=a)
        ehdu = ImageHDU(data=a)
        ehdu2 = ImageHDU(data=(a + 1))
        hdula = HDUList([phdu, ehdu, ehdu])
        hdulb = HDUList([phdu, ehdu2, ehdu])
        diff = FITSDiff(hdula, hdulb)

        assert_false(diff.identical)
        assert_equal(diff.diff_hdu_count, ())
        assert_equal(len(diff.diff_hdus), 1)
        assert_equal(diff.diff_hdus[0][0], 1)

        hdudiff = diff.diff_hdus[0][1]
        assert_false(hdudiff.identical)
        assert_equal(hdudiff.diff_extnames, ())
        assert_equal(hdudiff.diff_extvers, ())
        assert_equal(hdudiff.diff_extension_types, ())
        assert_true(hdudiff.diff_headers.identical)
        assert_false(hdudiff.diff_data is None)

        datadiff = hdudiff.diff_data
        assert_true(isinstance(datadiff, ImageDataDiff))
        assert_false(datadiff.identical)
        assert_equal(datadiff.diff_dimensions, ())
        assert_equal(datadiff.diff_pixels,
                     [((0, y), (y, y + 1)) for y in range(10)])
        assert_equal(datadiff.diff_ratio, 1.0)
        assert_equal(datadiff.diff_total, 100)

        report = diff.report()
        # Primary HDU and 2nd extension HDU should have no differences
        assert_true('Primary HDU' not in report)
        assert_true('Extension HDU 2' not in report)
        assert_true('Extension HDU 1' in report)

        assert_true('Headers contain differences' not in report)
        assert_true('Data contains differences' in report)
        for y in range(10):
            assert_true('Data differs at [%d, 1]' % (y + 1) in report)
        assert_true('100 different pixels found (100.00% different).' in
                    report)

    def test_diff_nans(self):
        """Regression test for #204."""

        # First test some arrays that should be equivalent....
        arr = np.empty((10, 10), dtype=np.float64)
        arr[:5] = 1.0
        arr[5:] = np.nan
        arr2 = arr.copy()

        table = np.rec.array([(1.0, 2.0), (3.0, np.nan), (np.nan, np.nan)],
                             names=['cola', 'colb']).view(pyfits.FITS_rec)
        table2 = table.copy()

        assert_true(ImageDataDiff(arr, arr2).identical)
        assert_true(TableDataDiff(table, table2).identical)

        # Now let's introduce some differences, where there are nans and where
        # there are not nans
        arr2[0][0] = 2.0
        arr2[5][0] = 2.0
        table2[0][0] = 2.0
        table2[1][1] = 2.0

        diff = ImageDataDiff(arr, arr2)
        assert_false(diff.identical)
        assert_equal(diff.diff_pixels[0], ((0, 0), (1.0, 2.0)))
        assert_equal(diff.diff_pixels[1][0], (5, 0))
        assert_true(np.isnan(diff.diff_pixels[1][1][0]))
        assert_equal(diff.diff_pixels[1][1][1], 2.0)

        diff = TableDataDiff(table, table2)
        assert_false(diff.identical)
        assert_equal(diff.diff_values[0], (('cola', 0), (1.0, 2.0)))
        assert_equal(diff.diff_values[1][0], ('colb', 1))
        assert_true(np.isnan(diff.diff_values[1][1][0]))
        assert_equal(diff.diff_values[1][1][1], 2.0)

        # What about in Headers?

        h = pyfits.Header([('A', 1), ('B', np.nan)])
        h2 = h.copy()

        assert_true(HeaderDiff(h, h2).identical)

        h2['B'] = 1.0

        diff = HeaderDiff(h, h2)
        assert_false(diff.identical)
        assert_true(np.isnan(diff.diff_keyword_values['B'][0][0]))
        assert_equal(diff.diff_keyword_values['B'][0][1], 1.0)
