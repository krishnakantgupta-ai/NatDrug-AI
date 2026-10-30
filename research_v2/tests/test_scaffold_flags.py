import sys
import unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'code'))
import pandas as pd
from correct_scaffolds import correct_scaffolds

class ScaffoldFlagTests(unittest.TestCase):
    def test_acyclic_not_shared(self):
        train=pd.DataFrame({'Scaffold':['', 'c1ccccc1']})
        screen=pd.DataFrame({'Scaffold':['','c1ccccc1','C1CCCCC1',None],
            'Scaffold_seen_in_train':[True,True,False,True],
            'Predicted_pIC50':[5.0,6.0,7.0,4.0]})
        out=correct_scaffolds(screen,train)
        self.assertEqual(out.Scaffold_status_full_curated_CHEMBL203.tolist(),
            ['NO_MURCKO_FRAMEWORK','SHARED_FRAMEWORK','NOVEL_FRAMEWORK','NO_MURCKO_FRAMEWORK'])
        self.assertEqual(out.Scaffold_status_model_training_subset.iloc[0], 'NO_MURCKO_FRAMEWORK')
        self.assertTrue(pd.isna(out.Scaffold_seen_in_train_CORRECTED.iloc[0]))
        self.assertEqual(out.Predicted_pIC50.tolist(), [5.0,6.0,7.0,4.0])

if __name__=='__main__': unittest.main()
