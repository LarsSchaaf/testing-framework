import os.path
import vacancy_LS as vacancy

properties = vacancy.vacancy_all_xyz(
    os.path.abspath(os.path.dirname(__file__)), nn_cutoff=2.7, dummy=True
)
