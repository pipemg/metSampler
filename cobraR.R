library(reticulate)

py_available()
py_config()

main <- import_main()
builtins <- import_builtins()

cobra.test <- import("cobra.test")
source_python('sample.py')

from cobra.test import create_test_model
from cobra.flux_analysis import sample

model = create_test_model("textbook")
s = sample(model, 100)
s.head()