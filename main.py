from file_extractor import FileExtractor

class MainScript:

    def __init__(self, test):
        self.test = test

    def file_extractor(self):
        extractor_object = FileExtractor(self.test)
        extractor_object.extract_reference_genes()

test_object = MainScript("abedyr")

# Extract files and make directories
test_object.file_extractor()
