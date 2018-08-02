from cobra.test import create_test_model
from cobra.flux_analysis import sample

model = create_test_model("textbook")
s = sample(model, 100)
s.head()


