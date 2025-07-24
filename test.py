import unittest

from functions.masterdata import (
        master_data, 
        read_Surf_Areas, 
        read_config, 
        make_locate_list, 
        enter_locations, 
        read_plate_info,
        get_unique_combos,
        check_plate_info,
        infer_plate_info)


class TestStringMethods(unittest.TestCase):

    def test_upper(self):
        self.assertEqual('foo'.upper(), 'FOO')

    def test_isupper(self):
        self.assertTrue('FOO'.isupper())
        self.assertFalse('Foo'.isupper())

    def test_split(self):
        s = 'hello world'
        self.assertEqual(s.split(), ['hello', 'world'])
        # check that s.split fails when the separator is not a string
        with self.assertRaises(TypeError):
            s.split(2)

class read_plate_info (unittest.TestCase):

    def setUp(self) -> None:
        "Set up an instance before each test"
        pass

    def test_plate(self):

    def tearDown(self) -> None:
        "Clean up after this test"
        pass



class Test_read_config(unittest.TestCase):

    def setUp(self) -> None:
        "Set up an instance before each test"
        pass

    def test_read_config(self):
        result = readConfig()
        print(result)
        self.assertEqual(1,1)

        # Check that minimum configuration parameters are present

        assertIs(type(result['rootDir']), 
            type(str()), 
            msg='RootDir must be a string')
        assertIs(type(result['initialDataPath']), 
            type(str()), 
            msg='initialDataPath must be a string')
        assertIs(type(result['QCDataPath']), 
            type(str()), 
            msg='QCDataPath must be a string')
        assertIs(type(result['labWorksheet01Path']), 
            type(str()), 
            msg='labWorksheet01Path must be a string')
        assertIs(type(result['projectName']), 
            type(str()), 
            msg='projectName must be a string')
        assertIs(type(result['selectedData']), 
            type(list([])), 
            msg='selectedData must be a list')

    def tearDown(self) -> None:
        "Clean up after this test"
        pass


# class TestMasterDataMethods(unittest.TestCase):
#     def setUp(self):
#         self.data = master_data('The widget')

#         configDict = read_config()

# 	def test_read_config(self):
# 		result = readConfig()
# 		print(result)
# 		self.assertEqual(1,1)


if __name__ == '__main__':
    unittest.main()