import unittest


import rdkit
import rdkit.Chem
import rdkit.Chem.AtomPairs.Utils

import extract_fragments


class TestCalc(unittest.TestCase):
    def test_atom_pairs(self):
        molecule = rdkit.Chem.MolFromSmiles("c1ccccn1")
        result = [{"smiles": "CC", "index": 2509202, "type": "AP", "size": 2},
                  {"smiles": "CC", "index": 2509202, "type": "AP", "size": 2},
                  {"smiles": "CC", "index": 2509202, "type": "AP", "size": 2},
                  {"smiles": "CC", "index": 2509202, "type": "AP", "size": 2},
                  {"smiles": "CC", "index": 3557778, "type": "AP", "size": 3},
                  {"smiles": "CC", "index": 3557778, "type": "AP", "size": 3},
                  {"smiles": "CC", "index": 3557778, "type": "AP", "size": 3},
                  {"smiles": "CC", "index": 3557778, "type": "AP", "size": 3},
                  {"smiles": "CC", "index": 4606354, "type": "AP", "size": 4},
                  {"smiles": "CC", "index": 4606354, "type": "AP", "size": 4},
                  {"smiles": "CN", "index": 2574738, "type": "AP", "size": 2},
                  {"smiles": "CN", "index": 2574738, "type": "AP", "size": 2},
                  {"smiles": "CN", "index": 3623314, "type": "AP", "size": 3},
                  {"smiles": "CN", "index": 3623314, "type": "AP", "size": 3},
                  {"smiles": "CN", "index": 4671890, "type": "AP", "size": 4}]

        self.assertEqual(extract_fragments.extract_atompair_fragments(molecule), result)

        molecule = rdkit.Chem.MolFromSmiles("c1nccc2n1ccc2")
        result = [{"smiles": "CC", "index": 2509202, "type": "AP", "size": 2},
                  {"smiles": "CC", "index": 2509202, "type": "AP", "size": 2},
                  {"smiles": "CC", "index": 2509202, "type": "AP", "size": 2},
                  {"smiles": "CC", "index": 3557778, "type": "AP", "size": 3},
                  {"smiles": "CC", "index": 3557778, "type": "AP", "size": 3},
                  {"smiles": "CC", "index": 3557778, "type": "AP", "size": 3},
                  {"smiles": "CC", "index": 3557778, "type": "AP", "size": 3},
                  {"smiles": "CC", "index": 4606354, "type": "AP", "size": 4},
                  {"smiles": "CC", "index": 4606354, "type": "AP", "size": 4},
                  {"smiles": "CC", "index": 4606354, "type": "AP", "size": 4},
                  {"smiles": "CC", "index": 4606354, "type": "AP", "size": 4},
                  {"smiles": "CC", "index": 4606354, "type": "AP", "size": 4},
                  {"smiles": "CC", "index": 4606354, "type": "AP", "size": 4},
                  {"smiles": "CC", "index": 5654930, "type": "AP", "size": 5},
                  {"smiles": "CC", "index": 5654930, "type": "AP", "size": 5},
                  {"smiles": "CC", "index": 2510226, "type": "AP", "size": 2},
                  {"smiles": "CC", "index": 2510226, "type": "AP", "size": 2},
                  {"smiles": "CC", "index": 3558802, "type": "AP", "size": 3},
                  {"smiles": "CC", "index": 3558802, "type": "AP", "size": 3},
                  {"smiles": "CC", "index": 3558802, "type": "AP", "size": 3},
                  {"smiles": "CC", "index": 3558802, "type": "AP", "size": 3},
                  {"smiles": "CN", "index": 2574738, "type": "AP", "size": 2},
                  {"smiles": "CN", "index": 2574738, "type": "AP", "size": 2},
                  {"smiles": "CN", "index": 3623314, "type": "AP", "size": 3},
                  {"smiles": "CN", "index": 4671890, "type": "AP", "size": 4},
                  {"smiles": "CN", "index": 5720466, "type": "AP", "size": 5},
                  {"smiles": "CN", "index": 5720466, "type": "AP", "size": 5},
                  {"smiles": "CN", "index": 4671891, "type": "AP", "size": 4},
                  {"smiles": "CN", "index": 2575762, "type": "AP", "size": 2},
                  {"smiles": "CN", "index": 2575762, "type": "AP", "size": 2},
                  {"smiles": "CN", "index": 3624338, "type": "AP", "size": 3},
                  {"smiles": "CN", "index": 3624338, "type": "AP", "size": 3},
                  {"smiles": "CN", "index": 3624338, "type": "AP", "size": 3},
                  {"smiles": "CN", "index": 4672914, "type": "AP", "size": 4},
                  {"smiles": "CN", "index": 2575763, "type": "AP", "size": 2},
                  {"smiles": "NN", "index": 3624402, "type": "AP", "size": 3}]

        self.assertEqual(extract_fragments.extract_atompair_fragments(molecule), result)

        molecule = rdkit.Chem.MolFromSmiles("CCO")
        result = [{"smiles": "CC", "index": 2492801, "type": "AP", "size": 2},
                  {"smiles": "CO", "index": 3671425, "type": "AP", "size": 3},
                  {"smiles": "CO", "index": 2622850, "type": "AP", "size": 2}]

        self.assertEqual(extract_fragments.extract_atompair_fragments(molecule), result)

    def test_neighbourhood_fragments(self):
        #ECFP
        molecule = rdkit.Chem.MolFromSmiles("c1ccccn1")
        options = {
            "kekule": True,
            "isomeric": True
        }
        size = 3
        result = [{"smiles": "N1:C:C:C:C:C:1", "index": 755035130, "type": "ECFP", "size": 3}]
        self.assertEqual(extract_fragments.extract_neighbourhood_fragments(molecule, size, options, True), result)

        size = 2
        result = [{"smiles": "C(:C:C):C:N", "index": 1207774339, "type": "ECFP", "size": 2},
                  {"smiles": "C(:C:C):C:N", "index": 1207774339, "type": "ECFP", "size": 2},
                  {"smiles": "N(:C:C):C:C", "index": 1343371647, "type": "ECFP", "size": 2},
                  {"smiles": "C(:C:C):N:C", "index": 1821698485, "type": "ECFP", "size": 2},
                  {"smiles": "C(:C:C):N:C", "index": 1821698485, "type": "ECFP", "size": 2},
                  {"smiles": "C(:C:C):C:C", "index": 2763854213, "type": "ECFP", "size": 2}]

        self.assertEqual(extract_fragments.extract_neighbourhood_fragments(molecule, size, options, True), result)

        options = {
            "kekule": False,
            "isomeric": True
        }
        result = [{"smiles": "c(cc)cn", "index": 1207774339, "type": "ECFP", "size": 2},
                  {"smiles": "c(cc)cn", "index": 1207774339, "type": "ECFP", "size": 2},
                  {"smiles": "n(cc)cc", "index": 1343371647, "type": "ECFP", "size": 2},
                  {"smiles": "c(cc)nc", "index": 1821698485, "type": "ECFP", "size": 2},
                  {"smiles": "c(cc)nc", "index": 1821698485, "type": "ECFP", "size": 2},
                  {"smiles": "c(cc)cc", "index": 2763854213, "type": "ECFP", "size": 2}]

        self.assertEqual(extract_fragments.extract_neighbourhood_fragments(molecule, size, options, True), result)

        molecule = rdkit.Chem.MolFromSmiles("c1nccc2n1ccc2")
        options = {
            "kekule": False,
            "isomeric": False
        }
        result = [{"smiles": "c(cn)c(c)n", "index": 201245292, "type": "ECFP", "size": 2},
                  {"smiles": "c(cc)n(c)c", "index": 405194198, "type": "ECFP", "size": 2},
                  {"smiles": "n(cc)(cn)c(c)c", "index": 924977737, "type": "ECFP", "size": 2},
                  {"smiles": "c(cc)nc", "index": 1717044408, "type": "ECFP", "size": 2},
                  {"smiles": "c(cc)(cc)n(c)c", "index": 2345490282, "type": "ECFP", "size": 2},
                  {"smiles": "c(nc)n(c)c", "index": 2558786292, "type": "ECFP", "size": 2},
                  {"smiles": "n(cc)cn", "index": 2910395211, "type": "ECFP", "size": 2},
                  {"smiles": "c(cc)cn", "index": 3428161631, "type": "ECFP", "size": 2},
                  {"smiles": "c(cc)c(c)n", "index": 3896685563, "type": "ECFP", "size": 2}]

        self.assertEqual(extract_fragments.extract_neighbourhood_fragments(molecule, size, options, True), result)

        size = 3
        molecule = rdkit.Chem.MolFromSmiles("C[C@H](O)c1ccccc1")
        options = {
            "kekule": False,
            "isomeric": True
        }
        result = [{"smiles": "c1ccccc1", "index": 742000539, "type": "ECFP", "size": 3},
                  {"smiles": "c1cccc(C)c1", "index": 997097697, "type": "ECFP", "size": 3},
                  {"smiles": "c1ccccc1[C@H](C)O", "index": 1566387358, "type": "ECFP", "size": 3}]

        self.assertEqual(extract_fragments.extract_neighbourhood_fragments(molecule, size, options, True), result)

        options = {
            "kekule": False,
            "isomeric": False
        }
        result = [{"smiles": "c1ccccc1", "index": 742000539, "type": "ECFP", "size": 3},
                  {"smiles": "c1cccc(C)c1", "index": 997097697, "type": "ECFP", "size": 3},
                  {"smiles": "c1ccccc1C(C)O", "index": 1566387358, "type": "ECFP", "size": 3}]

        self.assertEqual(extract_fragments.extract_neighbourhood_fragments(molecule, size, options, True), result)

        options = {
            "kekule": True,
            "isomeric": True
        }
        result = [{"smiles": "C1:C:C:C:C:C:1", "index": 742000539, "type": "ECFP", "size": 3},
                  {"smiles": "C1:C:C:C:C(C):C:1", "index": 997097697, "type": "ECFP", "size": 3},
                  {"smiles": "C1:C:C:C:C:C:1[C@H](C)O", "index": 1566387358, "type": "ECFP", "size": 3}]

        self.assertEqual(extract_fragments.extract_neighbourhood_fragments(molecule, size, options, True), result)

        options = {
            "kekule": True,
            "isomeric": False
        }
        result = [{"smiles": "C1:C:C:C:C:C:1", "index": 742000539, "type": "ECFP", "size": 3},
                  {"smiles": "C1:C:C:C:C(C):C:1", "index": 997097697, "type": "ECFP", "size": 3},
                  {"smiles": "C1:C:C:C:C:C:1C(C)O", "index": 1566387358, "type": "ECFP", "size": 3}]

        self.assertEqual(extract_fragments.extract_neighbourhood_fragments(molecule, size, options, True), result)

        # FCFP
        molecule = rdkit.Chem.MolFromSmiles("c1ccccn1")
        options = {
            "kekule": True,
            "isomeric": True
        }
        size = 3
        result = [{"smiles": "C1:C:C:C:C:N:1", "index": 1067478186, "type": "FCFP", "size": 3}]
        self.assertEqual(extract_fragments.extract_neighbourhood_fragments(molecule, size, options, False), result)

        molecule = rdkit.Chem.MolFromSmiles("c1nccc2n1ccc2")
        options = {
            "kekule": False,
            "isomeric": False
        }
        size = 2
        result = [{"smiles": "c(cc)(cc)n(c)c", "index": 435849959, "type": "FCFP", "size": 2},
                  {"smiles": "n(cc)cn", "index": 1127424909, "type": "FCFP", "size": 2},
                  {"smiles": "c(cc)cn", "index": 1230564256, "type": "FCFP", "size": 2},
                  {"smiles": "c(cc)nc", "index": 1251070542, "type": "FCFP", "size": 2},
                  {"smiles": "n(cc)(cn)c(c)c", "index": 1476508118, "type": "FCFP", "size": 2},
                  {"smiles": "c(nc)n(c)c", "index": 2154510652, "type": "FCFP", "size": 2},
                  {"smiles": "c(cc)n(c)c", "index": 2226952373, "type": "FCFP", "size": 2},
                  {"smiles": "c(cc)c(c)n", "index": 2460461453, "type": "FCFP", "size": 2},
                  {"smiles": "c(cn)c(c)n", "index": 2460461555, "type": "FCFP", "size": 2}]

        self.assertEqual(extract_fragments.extract_neighbourhood_fragments(molecule, size, options, False), result)

        molecule = rdkit.Chem.MolFromSmiles("CCO")
        size = 1
        result = [{"smiles": "CC", "index": 3205495869, "type": "FCFP", "size": 1},
                  {"smiles": "OC", "index": 3205496825, "type": "FCFP", "size": 1},
                  {"smiles": "C(C)O", "index": 3766532901, "type": "FCFP", "size": 1}]

        self.assertEqual(extract_fragments.extract_neighbourhood_fragments(molecule, size, options, False), result)


    def test_path_fragments(self):
        molecule = rdkit.Chem.MolFromSmiles("c1ccccn1")
        options = {
            "kekule": False,
            "isomeric": False
        }
        size = 2
        result = [{"smiles": "cc", "index": 83025, "type": "TT", "size": 2},
                  {"smiles": "cn", "index": 148561, "type": "TT", "size": 2},
                  {"smiles": "cc", "index": 83025, "type": "TT", "size": 2},
                  {"smiles": "cc", "index": 83025, "type": "TT", "size": 2},
                  {"smiles": "cc", "index": 83025, "type": "TT", "size": 2},
                  {"smiles": "cn", "index": 148561, "type": "TT", "size": 2}]

        self.assertEqual(extract_fragments.extract_path_fragments(molecule, size, options), result)

        size = 3
        result = [{"smiles": "ccc", "index": 85016657, "type": "TT", "size": 3},
                  {"smiles": "cnc", "index": 85082193, "type": "TT", "size": 3},
                  {"smiles": "ccn", "index": 152125521, "type": "TT", "size": 3},
                  {"smiles": "ccc", "index": 85016657, "type": "TT", "size": 3},
                  {"smiles": "ccc", "index": 85016657, "type": "TT", "size": 3},
                  {"smiles": "ccn", "index": 152125521, "type": "TT", "size": 3}]

        self.assertEqual(extract_fragments.extract_path_fragments(molecule, size, options), result)

        molecule = rdkit.Chem.MolFromSmiles("c1nccc2n1ccc2")
        options = {
            "kekule": False,
            "isomeric": True
        }
        size = 5;
        result = [{"smiles": "cccnc", "index": 90245936857169, "type": "TT", "size": 5},
                  {"smiles": "cccnc", "index": 89216219430993, "type": "TT", "size": 5},
                  {"smiles": "cccnc", "index": 89216219430993, "type": "TT", "size": 5},
                  {"smiles": "cccnc", "index": 89216218382417, "type": "TT", "size": 5},
                  {"smiles": "ccncn", "index": 159515237499985, "type": "TT", "size": 5},
                  {"smiles": "ccncn", "index": 159515237499985, "type": "TT", "size": 5},
                  {"smiles": "ccncn", "index": 159515237498961, "type": "TT", "size": 5},
                  {"smiles": "ncccn", "index": 160615754711185, "type": "TT", "size": 5},
                  {"smiles": "ccccn", "index": 159515169342545, "type": "TT", "size": 5},
                  {"smiles": "cncnc", "index": 90315730075729, "type": "TT", "size": 5},
                  {"smiles": "cncnc", "index": 89216218447953, "type": "TT", "size": 5},
                  {"smiles": "cccnc", "index": 89216219430993, "type": "TT", "size": 5},
                  {"smiles": "ccccc", "index": 89146426212433, "type": "TT", "size": 5},
                  {"smiles": "ccncn", "index": 160614748078161, "type": "TT", "size": 5},
                  {"smiles": "ccncc", "index": 89147567063121, "type": "TT", "size": 5},
                  {"smiles": "ccccc", "index": 89147498905681, "type": "TT", "size": 5},
                  {"smiles": "c1ccnc1", "index": 90315730010193, "type": "TT", "size": 5}]

        self.assertEqual(extract_fragments.extract_path_fragments(molecule, size, options), result)

        options = {
            "kekule": True,
            "isomeric": False
        }
        size = 3
        result = [{"smiles": "C:N:C", "index": 85082193, "type": "TT", "size": 3},
                  {"smiles": "C:N:C", "index": 86131793, "type": "TT", "size": 3},
                  {"smiles": "C:N:C", "index": 85083217, "type": "TT", "size": 3},
                  {"smiles": "N:C:N", "index": 153174161, "type": "TT", "size": 3},
                  {"smiles": "C:C:N", "index": 152125521, "type": "TT", "size": 3},
                  {"smiles": "C:C:C", "index": 86065233, "type": "TT", "size": 3},
                  {"smiles": "C:C:N", "index": 153175121, "type": "TT", "size": 3},
                  {"smiles": "C:C:C", "index": 85017681, "type": "TT", "size": 3},
                  {"smiles": "C:N:C", "index": 86131793, "type": "TT", "size": 3},
                  {"smiles": "C:C:C", "index": 86065233, "type": "TT", "size": 3},
                  {"smiles": "C:C:N", "index": 153175121, "type": "TT", "size": 3},
                  {"smiles": "C:C:N", "index": 153174097, "type": "TT", "size": 3},
                  {"smiles": "C:C:C", "index": 85016657, "type": "TT", "size": 3}]

        self.assertEqual(extract_fragments.extract_path_fragments(molecule, size, options), result)

        options = {
            "kekule": True,
            "isomeric": True
        }
        size = 8;
        result =[{"smiles": "C:C:N1:C:C:C:N:C:1", "index": 95794032134814304387153, "type": "TT", "size": 8},
                 {"smiles": "C:C:C:C:C:C:N:C", "index": 95794032134814236229713, "type": "TT", "size": 8},
                 {"smiles": "C:C:C1:C:C:C:N:1:C", "index": 95795185056317770383441, "type": "TT", "size": 8},
                 {"smiles": "C:C1:C:C:C:N:1:C:N", "index": 171278182067926592472145, "type": "TT", "size": 8},
                 {"smiles": "N:C:C:C1:C:C:C:N:1", "index": 171278108885601952546897, "type": "TT", "size": 8},
                 {"smiles": "C:N:C:N1:C:C:C:C:1", "index": 95794032206282492035153, "type": "TT", "size": 8},
                 {"smiles": "C:C:C1:C:C:N:C:N:1", "index": 95720317216182156476497, "type": "TT", "size": 8},
                 {"smiles": "C:C:C:N:C:N:C:C", "index": 95720317216182155427921, "type": "TT", "size": 8},
                 {"smiles": "C:C1:C:C:N:C:N:1:C", "index": 95795185126686513513553, "type": "TT", "size": 8}]

        self.assertEqual(extract_fragments.extract_path_fragments(molecule, size, options), result)

        molecule = rdkit.Chem.MolFromSmiles("CCO")
        options = {
            "kekule": False,
            "isomeric": True
        }
        size = 2
        result = [{"smiles": "CC", "index": 66624, "type": "TT", "size": 2},
                  {"smiles": "CO", "index": 196673, "type": "TT", "size": 2}]

        self.assertEqual(extract_fragments.extract_path_fragments(molecule, size, options), result)


if __name__ == "__main__":
    unittest.main()

