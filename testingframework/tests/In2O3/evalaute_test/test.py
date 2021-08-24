from testingframework.share.evaluate import evaluate_all_xyz
import os


properties = evaluate_all_xyz(os.path.abspath(os.path.dirname(__file__)), compare=True)
