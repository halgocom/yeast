from interface import AppInterface
from dotenv import load_dotenv
from os import getenv

class CrestInterface(AppInterface):
    app_path = f"{getenv("BIN_PATH")}/crest/crest"

    def __init__(self,save_path):
        self.save_path=save_path

    def gen_xyz(self,smiles:str=None):
        ...

    def opt(self,mol):
        ...
    def protonate(self,mol):
        ...
    def deprotonate(self,mol):
        ...