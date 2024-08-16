from src import *

import unittest

class TestSrcMethods(unittest.TestCase):

    def setUp(self):
        self.src = src_open(SRC_PROFILE_DEFAULT, MONO,  48000, 48000)

    def tearDown(self):
        src_close(self.src)

    def test_flow(self):
        src_set_scale(self.src, 1.)


if __name__ == '__main__':
    unittest.main()