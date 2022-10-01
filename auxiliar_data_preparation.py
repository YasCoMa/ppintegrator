import os

class Auxiliar_system_preparation:
    def run(self):
        print("Running for predrep")
        os.system("python3 prepare_data_triplification.py -rt 1 -fec params_predrep_5k.json")

        print("Running for hint")
        os.system("python3 prepare_data_triplification.py -rt 3 -fec params_hint.json")

        print("Running for string")
        os.system("python3 prepare_data_triplification.py -rt 2 -fec params_string.json")

a=Auxiliar_system_preparation()
a.run()
