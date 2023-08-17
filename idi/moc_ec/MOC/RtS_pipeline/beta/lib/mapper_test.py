import unittest
from mapper import Mapper

class TestMapper(unittest.TestCase):
    def setUp(self):
        confd = None
        self.mapo = Mapper(confd)

    def test_get_prefix_startswith(self):
        prefix_lst = ['abcde_suf', 'ijklm_suf', 'zzzzz_suff1']
        lsample_id = 'abcde'
        final_prefix_lst = set(['abcde_suf'])
        mapo = self.mapo
        lres = mapo.get_prefix_startswith(prefix_lst, lsample_id)
        self.assertEqual(lres, final_prefix_lst)

    # From https://stackoverflow.com/questions/8672754/how-to-show-the-error-messages-caught-by-assertraises-in-unittest-in-python2-7
    def assertRaisesWithMessage(self, msg, func, *args, **kwargs):
        try:
            func(*args, **kwargs)
            self.assertFail()
        except Exception as inst:
            self.assertEqual(inst.message, msg)

    def test_get_prefix_startswith_2(self):
        prefix_lst = ['abcde_suf', 'abcde_suf2', 'ijklm_suf', 'zzzzz_suff1']
        lsample_id = 'abcde'
        mapo = self.mapo
        lerr_str = 'From get_prefix_startswith, l_int_lst contains multiple prefix'
        self.assertRaisesWithMessage(lerr_str, mapo.get_prefix_startswith, prefix_lst, lsample_id)

    def test_get_prefix_startswith_3(self):
        prefix_lst = ['eabcde_suf', 'babcde_suf2', 'ijklm_suf', 'zzzzz_suff1']
        lsample_id = 'abcde'
        mapo = self.mapo
        lerr_str = 'From get_prefix_startswith, l_int_lst empty'
        self.assertRaisesWithMessage(lerr_str, mapo.get_prefix_startswith, prefix_lst, lsample_id)

    def test_get_prefix_startswith_4(self):
        prefix_lst = ['rabcde_suf', 'ijklm_suf', 'zzzzz_suff1']
        lsample_id = 'abcde'
        mapo = self.mapo
        lerr_str = 'From get_prefix_startswith, l_int_lst empty'
        self.assertRaisesWithMessage(lerr_str, mapo.get_prefix_startswith, prefix_lst, lsample_id)

    def test_get_prefix_startswith_5(self):
        prefix_lst = []
        lsample_id = 'abcde'
        mapo = self.mapo
        lerr_str = 'From get_prefix_startswith, prefix_lst empty'
        self.assertRaisesWithMessage(lerr_str, mapo.get_prefix_startswith, prefix_lst, lsample_id)

    def test_get_prefix_startswith_6(self):
        prefix_lst = ['abcde_suf', 'ijklm_suf', 'zzzzz_suff1']
        lsample_id = ''
        mapo = self.mapo
        lerr_str = 'From get_prefix_startswith, lsample_id empty'
        self.assertRaisesWithMessage(lerr_str, mapo.get_prefix_startswith, prefix_lst, lsample_id)

    def test_get_prefix_startswith_7(self):
        prefix_lst = []
        lsample_id = ''
        mapo = self.mapo
        lerr_str = 'From get_prefix_startswith, both args empty'
        self.assertRaisesWithMessage(lerr_str, mapo.get_prefix_startswith, prefix_lst, lsample_id)



if __name__ == '__main__':
    unittest.main()
    
